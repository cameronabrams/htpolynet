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


def _triple_key(donors):
    """A ring's identity, in terms that survive from one iteration to the next.

    A triple is held as row numbers into the site table, and that table is rebuilt
    every iteration and shrinks as groups react, so those numbers name a different
    ring next time.  The three donor atoms do not move: one per group, unique to it,
    and fixed for the life of the system.

    The two ring orientations deliberately share a key.  :func:`candidate_triples`
    scores both and keeps the better one, so the orientation that was declined is the
    better one, and the other would have farther to go.
    """
    return frozenset(int(d) for d in donors)


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
                  {'ensemble': 'nvt', 'temperature': 600, 'nsteps': 1000},
                  {'ensemble': 'npt', 'temperature': 600, 'pressure': 1, 'nsteps': 2000}],
        'settle': [{'ensemble': 'min'},
                   {'ensemble': 'nvt', 'temperature': 300, 'nsteps': 2500}],
        'closure': {'nstages': 8, 'target': 0.15, 'kb': 300000.0, 'max_accept': 0.30,
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
        self.declined = set()

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
            allowed = self._drop_declined(cands, sites)
            chosen = select_disjoint(allowed, max_triples=self.dicts['max_rings_per_iteration'])
            passed = len(cands) - len(allowed)
            logger.info(f'Iteration {self.state.iter}: {len(sites)} unreacted group(s), '
                        f'{len(cands)} candidate ring(s) within {self.radius:.3f} nm'
                        + (f' ({passed} passed over as already failed)' if passed else '')
                        + f', {len(chosen)} accepted')
            if len(chosen) >= floor or self.radius >= self.dicts['max_radius']:
                if not chosen and cands and self.declined:
                    # the memory is now the only thing keeping this iteration idle, and
                    # the groups have moved since they were declined; try them again
                    logger.info(f'Iteration {self.state.iter}: every ring in reach has '
                                f'failed before; forgetting {len(self.declined)} of them '
                                f'and offering them again')
                    self.declined.clear()
                    chosen = select_disjoint(cands,
                                             max_triples=self.dicts['max_rings_per_iteration'])
                return sites, chosen
            self.state.radius_index += 1
            logger.info(f'Radius increased to {self.radius:.3f} nm '
                        f'({len(chosen)}/{floor} ring(s) so far)')

    def _drop_declined(self, cands, sites):
        """Passes over the triples this cure has already failed to pull shut.

        A declined ring leaves its three groups unreacted by design, and the next search
        offers them again -- so a rigid endgame proposes the same triples, runs the same
        closure ladder, and declines them again.  Measured at conversion 0.97 with 24
        groups left: nine iterations running, 77 declines, about 67 s apiece, before one
        triple finally closed and the cure finished at 0.971.

        A preference, not a prohibition, because the geometry between one decline and the
        next search is not identical -- the inter-iteration relaxation moves the groups,
        which is how that cure escaped.  Skipping a triple saves a ladder that would
        almost certainly fail; refusing it outright would have ended that cure three
        iterations early.  :meth:`search` therefore lifts the whole memory rather than
        come back empty-handed.
        """
        if not self.declined:
            return cands
        return [c for c in cands
                if _triple_key(sites.at[a, 'donor'] for a in c[1]) not in self.declined]

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
        return work

    def accept(self, bdf, work, chosen):
        """Drops any ring the closure ladder could not actually pull shut.

        Late in a cure the matrix is rigid enough that a ring sometimes stops well short
        -- one batch at conversion 0.95 ended 0.327 nm apart against a target of 0.150.
        Forming the bond anyway leaves it long: two ring bonds survived relaxation and a
        500 K anneal still 2.6 A apart, strained but stable, which is a defect the
        network then carries for good.

        Declining costs almost nothing.  The groups stay unreacted and are offered again
        next iteration, by which time the neighbourhood has moved; only if they can never
        close does the cure lose them, which is the right outcome.

        The combination itself is recorded for :meth:`_drop_declined`, which prefers not
        to spend another ladder on it while any untried ring is in reach.  That is what
        a rigid endgame otherwise does: nine iterations running, 77 declines, the same
        few triples every time.

        Args:
            bdf (pandas.DataFrame): this iteration's bonds, three rows per ring
            work (pandas.DataFrame): the same pairs with their final distances
            chosen (list): the rings the search accepted

        Returns:
            tuple: (bdf, chosen) keeping only the rings that closed
        """
        limit = self.dicts['closure']['max_accept']
        if not limit:
            return bdf, chosen
        worst = work['final_distance'].to_numpy().reshape(-1, 3).max(axis=1)
        keep = [n for n in range(len(chosen)) if worst[n] <= limit]
        if len(keep) == len(chosen):
            return bdf, chosen
        for n in range(len(chosen)):
            if n not in keep:
                self.declined.add(_triple_key(bdf[bdf['triple'] == n]['ai']))
                logger.info(f'Iteration {self.state.iter}: declining a ring still '
                            f'{worst[n]:.3f} nm from closing (limit {limit:.3f} nm); its '
                            f'groups stay unreacted, and this combination waits until '
                            f'nothing untried is in reach')
        renumber = {old: new for new, old in enumerate(keep)}
        bdf = bdf[bdf['triple'].isin(keep)].copy()
        bdf['triple'] = [renumber[t] for t in bdf['triple']]
        return bdf.reset_index(drop=True), [chosen[n] for n in keep]

    def relax(self, TC, gromacs_dict=None):
        """Eases a freshly closed ring's bonds in, and lets the box respond.

        A ring bond is formed at whatever length the closure ladder reached, a good way
        short of the 1.34 A it wants, so the system is strained the moment the bonds
        exist.  CURE relaxes its new bonds for the same reason; without it here, the
        strain is still there when postcure MD starts, and that run dies.

        These stages mirror CURE's ``relax.equilibration`` exactly, and the whole
        sequence is above Tg rather than only its last stage.  The NVT ran at 300 K
        until 2.11.4, which cooled the network between the closure ladder at 600 K and
        the next ladder at 600 K, and would have handed the barostat a configuration
        equilibrated cold -- a good way to measure no densification and conclude the
        constant-pressure stage does not work.

        The last stage is constant-pressure, and above Tg, because a cure that never
        runs NPT cannot densify.  Before this, every ring-cure stage was NVT: measured
        over 22 iterations of a 233-triazine build the box stayed at 5.241 nm to four
        decimals and the density at 1155.9 kg/m3 from first ring to last, with the
        entire volume change deferred to a single postcure NPT.  CURE's density climbs
        through its cure, from about 1085 to 1170.  Cure shrinkage -- 0.0324 ml/g for
        BADCy by Snow -- cannot be reproduced at constant volume by construction, and a
        network formed in a box that cannot respond builds up internal stress with
        nowhere to put it.

        Args:
            TC (TopoCoord): the system
            gromacs_dict (dict): gromacs directives
        """
        self._run_stages(TC, f'ringrelax-{self.state.iter}', self.dicts['relax'],
                         gromacs_dict or {})
        logger.info(f'Iteration {self.state.iter}: new ring bonds relaxed')

    def settle(self, TC, gromacs_dict=None):
        """Relaxes the network one last time before the cure hands it on.

        Every batch of rings but the last is relaxed twice: once by :meth:`relax` as
        soon as its bonds exist, and again -- much harder -- by the next iteration's
        closure ladder, which is eight minimizations and 16 ps at 600 K against relax's
        one minimization and 4 ps at 300 K.  The final batch only ever gets the first,
        and whatever the cure ends on goes straight into postcure.

        Measured by htpolynet-study on three 2.11.1 builds, all of which ended on a
        ring-forming iteration, so all of which had been relaxed: bond energy entering
        the anneal tracked how many rings that last iteration formed -- 1, 1 and 4 rings
        for 18849, 26113 and 35959 kJ/mol, against 13619 +/- 226 for a route that enters
        already relaxed.  An inventory of one such structure found five bonds beyond
        2.0 A: one triazine's three ring bonds at 2.24 A, and two monomers whose own
        backbones the ladder had stretched to 2.34 and 3.20 A while pulling their two
        arms toward different rings.  Nothing after the restraints come off relieves
        either, and the strain then has to be absorbed by a 500 K anneal, which is where
        11 of 20 builds died.

        Unconditional on purpose.  The cure has four ways to end -- target reached,
        iterations exhausted, radius exhausted, every ring declined -- and which of them
        leaves strain behind is not worth reasoning about per exit when the remedy is
        this cheap.

        Args:
            TC (TopoCoord): the system
            gromacs_dict (dict): gromacs directives
        """
        stages = self.dicts['settle']
        if not stages:
            return False
        TC.write_top(f'ringsettle-{self.state.iter}.top')
        TC.write_gro(f'ringsettle-{self.state.iter}.gro')
        self._run_stages(TC, f'ringsettle-{self.state.iter}', stages, gromacs_dict or {})
        logger.info(f'Ring cure settled after {self.state.rings} ring(s); the network is '
                    f'relaxed before postcure rather than at 500 K')
        return True

    def _run_stages(self, TC, deffnm, stages, gromacs_dict):
        """Runs one equilibration sequence, as the cure's ladders do."""
        for stage in stages:
            ensemble = stage['ensemble']
            mdp = f'relax-{ensemble}' if ensemble != 'min' else 'relax-min'
            pfs.checkout(pfs.Dirs.mdp_file(mdp))
            if ensemble != 'min':
                mods = {'ref_t': stage.get('temperature', 300),
                        'gen-temp': stage.get('temperature', 300),
                        'gen-vel': 'yes', 'nsteps': stage.get('nsteps', 1000)}
                if ensemble == 'npt':
                    mods['ref_p'] = stage.get('pressure', 1)
                mdp_modify(f'{mdp}.mdp', mods)
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
