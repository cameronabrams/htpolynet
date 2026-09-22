"""Runs a ring-closing cure: cyclotrimerization rather than pairwise crosslinking.

CURE forms one bond at a time between two reactants, growing a search radius until
enough candidates appear.  A cyclotrimerization forms three bonds among three groups
at once, and the ring has to be pulled shut before any of them exists, so the loop is
its own.  Everything else is shared: the same templates, the same relaxation and
equilibration machinery, the same postcure path.

What one iteration does:

1. find the triangles of unreacted groups whose three new bonds are all short enough
   (:mod:`htpolynet.cure.triplesearch`), and pack them so no group is used twice;
2. restrain all three pairs of every chosen ring and walk the restraints down to a
   bond length, which is what makes the ring closable at all;
3. form the bonds, splice the trimer template over each ring's three residues
   (:mod:`htpolynet.core.productsplice`), and settle the charge per molecule;
4. relax, and count the groups consumed.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

import numpy as np
import yaml

from ..core import projectfilesystem as pfs
from ..core.productsplice import map_product_from_template
from ..cure.triplesearch import (candidate_triples, reactive_sites, residue_map,
                                 select_disjoint, site_name_maps)
from ..external.gromacs import mdp_modify
from ..repair.topology_surgery import neutralize_touched_fragments

logger = logging.getLogger(__name__)


class RingCureState:
    """How far a ring cure has got, and enough to resume it."""

    def __init__(self, iteration=1, rings=0, groups_consumed=0, radius_index=0, total_groups=0):
        self.iter = iteration
        self.rings = rings
        self.groups_consumed = groups_consumed
        self.radius_index = radius_index
        self.total_groups = total_groups

    def to_yaml(self, filename):
        with open(filename, 'w') as f:
            yaml.dump(self.__dict__, f)

    @classmethod
    def from_yaml(cls, filename):
        with open(filename) as f:
            d = yaml.safe_load(f)
        # the attribute reads better as .iter, which is not what the constructor calls it
        d['iteration'] = d.pop('iter', 1)
        return cls(**d)

    @property
    def conversion(self):
        """Fraction of reactive groups consumed, which is what this chemistry reports.

        A cyclotrimerization consumes three groups per ring, so conversion is counted
        in groups rather than bonds -- the quantity experiments measure, and the one a
        bond count has to be mapped onto in the pre-formed-triazine route.
        """
        return self.groups_consumed / self.total_groups if self.total_groups else 0.0


class RingController:
    """Drives a cure whose reaction closes a ring."""

    defaults = {
        'search_radius': 0.5,
        'radial_increment': 0.05,
        'max_radius': 0.0,          # 0 -> half the shortest box vector
        'max_iterations': 100,
        'desired_conversion': 0.95,
        'max_rings_per_iteration': 0,
        'min_rings_per_iteration': 4,
        'same_molecule': False,
        'same_residue': False,
        'relax': [{'ensemble': 'min'},
                  {'ensemble': 'nvt', 'temperature': 300, 'nsteps': 2000}],
        'closure': {'nstages': 8, 'target': 0.15, 'kb': 300000.0,
                    'equilibration': [{'ensemble': 'min'},
                                      {'ensemble': 'nvt', 'temperature': 600, 'nsteps': 1000}]},
    }

    def __init__(self, ringdict=None, state=None):
        d = dict(self.defaults)
        d.update(ringdict or {})
        d['closure'] = {**self.defaults['closure'], **(d.get('closure') or {})}
        self.dicts = d
        self.state = state or RingCureState()
        self.rings_df = None

    # ---- what the loop decides, all of it testable without MD ----

    def setup(self, total_groups, max_radius):
        self.state.total_groups = total_groups
        if not self.dicts['max_radius']:
            self.dicts['max_radius'] = max_radius

    @property
    def radius(self):
        """Current search radius, grown by whole increments as iterations fail."""
        r = self.dicts['search_radius'] + self.state.radius_index * self.dicts['radial_increment']
        return min(r, self.dicts['max_radius']) if self.dicts['max_radius'] else r

    def is_cured(self):
        """True when the target conversion is reached or the iterations run out."""
        if self.state.conversion >= self.dicts['desired_conversion']:
            logger.info(f'Ring cure reached conversion {self.state.conversion:.3f}')
            return True
        if self.state.iter > self.dicts['max_iterations']:
            logger.info(f'Ring cure stopped after {self.dicts["max_iterations"]} iteration(s) '
                        f'at conversion {self.state.conversion:.3f}')
            return True
        return False

    def widen(self):
        """Grows the search radius after an iteration that found nothing.

        Returns:
            bool: False if there is room to grow, True if the radius is at its limit
                and the cure is therefore finished
        """
        if self.radius >= self.dicts['max_radius']:
            logger.info(f'No rings found at the maximum search radius {self.radius:.3f} nm')
            return True
        self.state.radius_index += 1
        logger.info(f'No rings found; widening the search to {self.radius:.3f} nm')
        return False

    def ring_floor(self):
        """How many rings an iteration should find before it stops widening.

        The pairwise cure grows its radius within an iteration until it has at least
        ``min_bonds_per_iteration`` bonds; this is the same rule counted in rings.
        Widening only when an iteration finds *nothing* is not enough: a rigid monomer
        keeps yielding one or two rings at the starting radius, so the radius never
        grows and the cure runs out of iterations far short of its target.  Measured on
        150 bisphenol A dicyanates, that is conversion 0.40 in eight iterations rather
        than the 0.60 asked for.

        The floor is clamped against the rings still needed to reach the target, so the
        last iteration does not widen in search of rings it would decline to form.
        """
        remaining = self.dicts['desired_conversion'] * self.state.total_groups - self.state.groups_consumed
        needed = max(1, int(np.ceil(remaining / 3.0)))
        floor = min(int(self.dicts['min_rings_per_iteration']), needed)
        if self.dicts['max_rings_per_iteration']:
            floor = min(floor, int(self.dicts['max_rings_per_iteration']))
        return max(1, floor)

    def search(self, adf, positions, box, resname, groups):
        """Finds and packs this iteration's rings, widening the radius until it has enough.

        Args:
            adf (pandas.DataFrame): the system's atoms
            positions (dict): global atom index -> position
            box (array-like): box diagonal
            resname (str): residue carrying the reactive groups
            groups (iterable): (donor, acceptor) atom-name pairs

        Returns:
            tuple: (sites, chosen) -- the groups considered and the rings accepted
        """
        sites = reactive_sites(adf, resname, groups=groups)
        floor = self.ring_floor()
        while True:
            cands = candidate_triples(sites, positions, self.radius, box,
                                      same_molecule=self.dicts['same_molecule'],
                                      same_residue=self.dicts['same_residue'])
            chosen = select_disjoint(cands, max_triples=self.dicts['max_rings_per_iteration'])
            logger.info(f'Iteration {self.state.iter}: {len(sites)} unreacted group(s), '
                        f'{len(cands)} candidate ring(s) within {self.radius:.3f} nm, '
                        f'{len(chosen)} accepted')
            if len(chosen) >= floor or self.radius >= self.dicts['max_radius']:
                return sites, chosen
            self.state.radius_index += 1
            logger.info(f'Radius increased to {self.radius:.3f} nm '
                        f'({len(chosen)}/{floor} ring(s) so far)')

    def record(self, chosen):
        """Counts the rings formed and the groups they consumed."""
        self.state.rings += len(chosen)
        self.state.groups_consumed += 3 * len(chosen)
        logger.info(f'Iteration {self.state.iter}: {len(chosen)} ring(s), '
                    f'{self.state.groups_consumed} of {self.state.total_groups} groups consumed '
                    f'(conversion {self.state.conversion:.3f})')

    # ---- the parts that need MD ----

    def close_rings(self, TC, bdf, gromacs_dict=None):
        """Pulls every chosen ring shut before any of its bonds exists.

        All three pairs of each ring are restrained together: with no bonds yet, a
        ladder on one pair lets the other two groups drift while it pulls.

        Args:
            TC (TopoCoord): the system
            bdf (pandas.DataFrame): this iteration's bonds, three rows per ring
            gromacs_dict (dict): gromacs directives

        Returns:
            pandas.DataFrame: the pairs, with their final distances
        """
        c = self.dicts['closure']
        work = bdf[['ai', 'aj']].copy().reset_index(drop=True)
        TC.add_length_attribute(work, attr_name='initial_distance')
        logger.info(f'Closing {work.shape[0] // 3} ring(s) from at most '
                    f'{work["initial_distance"].max():.3f} nm over {c["nstages"]} stages')
        TC.Topology.add_restraints(work, typ=6)
        for stage in range(c['nstages']):
            frac = (stage + 1) / c['nstages']
            lengths = work['initial_distance'] + frac * (c['target'] - work['initial_distance'])
            TC.Topology.set_restraint_parameters(work, lengths, c['kb'])
            # the MD reads the topology from disk, so the restraints have to be
            # written out again every stage or the run pulls on nothing
            TC.write_top(f'ringclose-{self.state.iter}-{stage + 1}.top')
            TC.write_gro(f'ringclose-{self.state.iter}-{stage + 1}.gro')
            self._run_stages(TC, f'ringclose-{self.state.iter}-{stage + 1}', c['equilibration'],
                             gromacs_dict or {})
            TC.add_length_attribute(work, attr_name='current_distance')
            logger.debug(f'ring-close stage {stage + 1}/{c["nstages"]}: '
                         f'max {work["current_distance"].max():.3f} nm')
        TC.Topology.remove_restraints(work)
        TC.add_length_attribute(work, attr_name='final_distance')
        worst = work['final_distance'].max()
        logger.info(f'Rings closed to at most {worst:.3f} nm (target {c["target"]:.3f})')
        if worst > 2.0 * c['target']:
            logger.warning(f'A ring-closing pair is still {worst:.3f} nm apart; its bond will be '
                           f'formed long and relaxed, which may strain the network')
        return work

    def relax(self, TC, gromacs_dict=None):
        """Eases a freshly closed ring's bonds in before anything else runs.

        A ring bond is formed at whatever length the closure ladder reached, a good way
        short of the 1.34 A it wants, so the system is strained the moment the bonds
        exist.  CURE relaxes its new bonds for the same reason; without it here, the
        strain is still there when postcure MD starts, and that run dies.

        Args:
            TC (TopoCoord): the system
            gromacs_dict (dict): gromacs directives
        """
        self._run_stages(TC, f'ringrelax-{self.state.iter}', self.dicts['relax'],
                         gromacs_dict or {})
        logger.info(f'Iteration {self.state.iter}: new ring bonds relaxed')

    def _run_stages(self, TC, deffnm, stages, gromacs_dict):
        """Runs one equilibration sequence, as the cure's ladders do."""
        for stage in stages:
            ensemble = stage['ensemble']
            mdp = f'relax-{ensemble}' if ensemble != 'min' else 'relax-min'
            pfs.checkout(pfs.Dirs.mdp_file(mdp))
            if ensemble != 'min':
                mdp_modify(f'{mdp}.mdp', {'ref_t': stage.get('temperature', 300),
                                          'gen-temp': stage.get('temperature', 300),
                                          'gen-vel': 'yes', 'nsteps': stage.get('nsteps', 1000)})
            TC.grompp_and_mdrun(out=f'{deffnm}-{ensemble}', mdp=mdp, **gromacs_dict)

    def form_rings(self, TC, bdf, sites, chosen, template, template_resids, name_translations=None,
                   atom_names=None):
        """Forms each ring's bonds and gives its residues the template's parameters.

        Args:
            TC (TopoCoord): the system
            bdf (pandas.DataFrame): this iteration's bonds, with a triple column
            sites (pandas.DataFrame): the groups, from the search
            chosen (list): the accepted rings
            template (Molecule): the trimer template
            template_resids (list): its residue numbers, in reactant order
            name_translations (dict): {group label: {template name: that group's name}},
                for residues ringing through a site other than the template's
            atom_names (set): template atom names to splice, one reactive site's share of
                a residue, for a monomer carrying more than one

        Returns:
            dict: ``atoms`` touched and ``rings`` spliced
        """
        pairs = [(int(r.ai), int(r.aj), int(r.order)) for r in bdf.itertuples()]
        # an addition: nothing is lost, so no hydrogen is offered up for any of them
        TC.make_bonds(pairs, explicit_sacH={i: [] for i in range(len(pairs))})
        touched = set()
        maps = residue_map(chosen, sites, template_resids)
        names = site_name_maps(chosen, sites, template_resids, name_translations or {})
        for n, (rmap, nmap) in enumerate(zip(maps, names)):
            ring_bonds = [(int(r.ai), int(r.aj)) for r in bdf[bdf['triple'] == n].itertuples()]
            stats = map_product_from_template(TC, template, rmap, new_bonds=ring_bonds,
                                              name_maps=nmap, atom_names=atom_names)
            touched.update(stats['atoms'])
        # the splice replaces charges wholesale over three residues at a time, so the
        # molecules it touched have to be brought back to neutral
        neutralize_touched_fragments(TC, touched)
        for r in bdf.itertuples():
            for idx in (int(r.ai), int(r.aj)):
                TC.set_gro_attribute_by_attributes('z', 0, {'globalIdx': idx})
                TC.set_gro_attribute_by_attributes('nreactions', 1, {'globalIdx': idx})
        logger.info(f'Formed {len(chosen)} ring(s); {len(touched)} atom(s) took template parameters')
        return {'atoms': sorted(touched), 'rings': len(chosen)}
