"""Manages generation of molecular templates.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging
import os
import shutil

from copy import deepcopy
from itertools import product

import numpy as np
import pandas as pd

from ..core import paramcache
from ..core import projectfilesystem as pfs
from ..core.bondtemplate import BondTemplate, BondTemplateList, ReactionBond, ReactionBondList
from ..core.topocoord import TopoCoord
from ..cure.chain import ChainManager
from ..cure.reaction import Reaction, ReactionList, reaction_stage
from ..external.ambertools import GAFFParameterize
from ..external.command import run
from ..external.gromacs import mdp_modify,gro_from_trr
from ..geometry.matrix4 import Matrix4
from ..io.gro import GRX_ATTRIBUTES

logger=logging.getLogger(__name__)


def _yield_bonds_as_df(R: Reaction, TC: TopoCoord, resid_mapper) -> pd.DataFrame:
    """Returns a DataFrame of reaction bonds by matching template bond records to atom instances.

    Args:
        R (Reaction): a Reaction
        TC (TopoCoord): molecule topology and coordinates
        resid_mapper (list[dict]): per-reactant maps from in-reactant resid to in-product resid

    Returns:
        pd.DataFrame: one row per bond with columns ai, aj, ri, rj, order, reactantName
    """
    rows = []
    for bondrec in R.bonds:
        atom_keys = bondrec['atoms']
        assert len(atom_keys) == 2, f'bond record must have exactly 2 atoms, got {len(atom_keys)}'
        atomrecs = [R.atoms[k] for k in atom_keys]
        atom_names = [rec['atom'] for rec in atomrecs]
        in_product_resids = [resid_mapper[atomrecs[x]['reactant'] - 1][atomrecs[x]['resid']] for x in (0, 1)]
        ai = TC.get_gro_attribute_by_attributes('globalIdx', {'resNum': in_product_resids[0], 'atomName': atom_names[0]})
        aj = TC.get_gro_attribute_by_attributes('globalIdx', {'resNum': in_product_resids[1], 'atomName': atom_names[1]})
        rows.append({'ai': ai, 'aj': aj, 'ri': in_product_resids[0], 'rj': in_product_resids[1],
                     'order': bondrec['order'], 'reactantName': R.product})
    return pd.DataFrame(rows)

_TRANSROT_STEP_DEG_ = 10.0
"""Angular step of the turn scan in Molecule.transrot."""

_TRANSROT_WARN_NM_ = 0.10
"""Closest approach, in nm, below which a template's starting geometry is reported as crowded.

Contacts of 1.2-1.5 A between the new piece and the rest are routine before
minimization; below 1 A they have been seen to stop sqm converging (a TMB ring
methyl H 0.82 A from a triazine N).
"""


def _rotate_about_axis(xyz, point, axis, theta):
    """Rotates points by theta about the line through point along unit vector axis (Rodrigues).

    Args:
        xyz (numpy.ndarray): (n, 3) positions
        point (numpy.ndarray): a point on the axis
        axis (numpy.ndarray): unit vector along the axis
        theta (float): angle in radians

    Returns:
        numpy.ndarray: rotated (n, 3) positions
    """
    if theta == 0.0:
        return xyz.copy()
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)
    return (xyz - point) @ R.T + point


def _min_distance(a, b):
    """Smallest distance between any point of a and any point of b; inf if either is empty."""
    if len(a) == 0 or len(b) == 0:
        return np.inf
    diff = a[:, np.newaxis, :] - b[np.newaxis, :, :]
    return float(np.sqrt((diff * diff).sum(axis=2)).min())


class Molecule:
    def __init__(self, name='', generator:Reaction=None, origin:str=None):
        self.name = name
        self.parentname = name # stereoisomer parent
        self.TopoCoord = TopoCoord()
        self.generator:Reaction = generator
        self.sequence = []
        self.origin = origin
        self.reaction_bonds:ReactionBondList = []
        self.bond_templates:BondTemplateList = []
        self.symmetry_relateds = []
        self.stereocenters = [] # list of atomnames
        self.stereoisomers = {}
        self.nconformers = 0
        self.conformers_dict = {}
        self.conformers = [] # just a list of gro file basenames
        self.zrecs = []
        self.is_reactant = False
        self.chain_manager = ChainManager()

    @classmethod
    def New(cls, mol_name, generator: Reaction, molrec: dict | None = None):
        """Generates a new, partially populated Molecule based on directives in the configuration input.

        Args:
            mol_name (str): name of molecule
            generator (Reaction): reaction that generates this molecule, if applicable
            molrec (dict, optional): dictionary of directives for this molecule, defaults to None

        Returns:
            Molecule: a new Molecule object
        """
        M = cls(name=mol_name, generator=generator)
        if not molrec: return M
        M.symmetry_relateds = molrec.get('symmetry_equivalent_atoms', [])
        M.stereocenters = molrec.get('stereocenters', [])
        # expand stereocenters: any atom in a symmetry-equivalent set with a declared
        # stereocenter pulls the rest of that set in as stereocenters too
        extra = set()
        for stc in M.stereocenters:
            for sc in M.symmetry_relateds:
                if stc in sc:
                    extra.update(sc)
        extra -= set(M.stereocenters)
        M.stereocenters.extend(extra)
        logger.debug(f'{M.name} stereocenters: {M.stereocenters}')
        # generate shells for new stereoisomers
        M.create_new_stereoisomers()
        logger.debug(f'{M.name} stereoisomers: {[s.name for s in M.stereoisomers.values()]}')
        M.conformers_dict = molrec.get('conformers', {})
        logger.debug(f'{M.name} conformers_dict {M.conformers_dict}')
        return M

    def update_zrecs(self, zrecs, moldict):
        """Updates the "z-records" based on z's declared in the input configuration file.

        Args:
            zrecs (dict): zrecs extracted from configuration file for this molecule
            moldict (dict): dictionary of available molecules
        """
        def replace_if_greater(rec, D, matchattr=[], maxattr=[]):
            if not matchattr or not maxattr: return
            for r in D:
                matched = all([r[a] == rec[a] for a in matchattr])
                if matched:
                    replace = all([r[a] < rec[a] for a in maxattr])
                    if replace:
                        D.remove(r)
                        D.append(rec)
                    return True
            return False
        logger.debug(f'Update zrecs in {self.name} from {zrecs}')
        seq = self.sequence
        for zr in zrecs:
            resid = zr['resid'] - 1
            rname = seq[resid]
            target = moldict[rname]
            logger.debug(f'{target.name} {target.zrecs} ->')
            found = replace_if_greater(zr, target.zrecs, matchattr=['resid', 'atom'], maxattr=['z'])
            if not found: target.zrecs.append(zr)
            logger.debug(f'-> {target.name} {target.zrecs}')

    def determine_sequence(self, moldict):
        """Recursively determines the sequence of a molecule using the network of reactions that must be executed to generate it from primitives.

        Args:
            moldict (dict): dictionary of available molecules

        Returns:
            list: list of resnames in order of sequence
        """
        if not self.generator: return [self.parentname]
        R:Reaction = self.generator
        thisseq = []
        for rid, rname in R.reactants.items():
            parentname = moldict[rname].parentname
            thisseq.extend(moldict[parentname].determine_sequence(moldict))
            # logger.debug(thisseq)
        return thisseq
    
    def set_sequence_from_moldict(self, moldict):
        """Sets the sequence of this molecule using the recursive determine_sequence method.

        Args:
            moldict (dict): dictionary of available molecules

        Returns:
            Molecule: self
        """
        self.sequence = self.determine_sequence(moldict)
        return self

    def set_sequence_from_coordinates(self):
        """set_sequence Establish the sequence-list (residue names in order) based on resNum attributes in atom list
        """
        adf = self.TopoCoord.gro_DataFrame('atoms')
        trial_sequence = []
        current_resid = 0
        for i, r in adf.iterrows():
            ri = r['resNum']
            rn = r['resName']
            if ri != current_resid:
                current_resid = ri
                trial_sequence.append(rn)
        assert trial_sequence == self.sequence, f'Error: trial {trial_sequence} does not equal self.seq {self.sequence}'
        return self

    def create_new_stereoisomers(self):
        """Generates new molecules to hold stereoisomers of self.

        Returns:
            None: returns early if no action taken
        """
        if self.generator: return  # we only consider stereoisomers on monomers
        if not self.stereocenters: return
        basename = self.name + '-S'
        b = [[0, 1] for _ in range(len(self.stereocenters))]
        sseq = product(*b)
        next(sseq) # skip the unmodified; it's the base molecule
        for x in sseq:
            s = ''.join([str(_) for _ in x])
            mname = f'{basename}{s}'
            self.stereoisomers[mname] = Molecule.New(mname, None)
            self.stereoisomers[mname].parentname = self.name
            
    def initialize_molecule_rings(self):
        """initialize_molecule_rings generates the dictionary of rings

        the dictionary of rings is keyed on ring size
        """
        TC = self.TopoCoord
        TC.Topology.detect_rings()
        logger.debug(f'Detected {len(TC.Topology.rings)} unique rings.')

    def initialize_monomer_grx_attributes(self):
        """initialize_monomer_grx_attributes initializes all GRX attributes of atoms in molecule
        """
        logger.debug(f'{self.name}')
        TC = self.TopoCoord
        TC.set_gro_attribute('z', 0)
        TC.set_gro_attribute('nreactions', 0)
        TC.set_gro_attribute('molecule', 1)
        TC.set_gro_attribute('molecule_name', self.name)
        for att in ['sea_idx','bondchain','bondchain_idx']:
            TC.set_gro_attribute(att,-1)
        # set symmetry class indices
        sea_idx=1
        logger.debug(f'{self.name}: symmetry_relateds {self.symmetry_relateds}')
        for s in self.symmetry_relateds:
            # logger.debug(f'sea_idx {sea_idx} set for set {s}')
            for atomName in s:
                # logger.debug(f'{atomName} {sea_idx}')
                TC.set_gro_attribute_by_attributes('sea_idx', sea_idx, {'atomName': atomName})
            sea_idx += 1
        # set z and nreactions
        idx = []
        for zr in self.zrecs:
            an = zr['atom']
            rnum = zr['resid']
            if rnum != 1: continue
            z = zr['z']
            logger.debug(f'{self.name} setting z for {an} {rnum} {z}')
            TC.set_gro_attribute_by_attributes('z', z, {'atomName': an, 'resNum': rnum})
            idx.append(TC.get_gro_attribute_by_attributes('globalIdx', {'atomName': an, 'resNum': rnum}))
            for sr in self.symmetry_relateds:
                # logger.debug(f'{self.name}: setting z for {an}, considering sr {sr}')
                if an in sr:
                    for bn in sr:
                        if bn == an: continue
                        # logger.debug(f'{self.name}: setting z for {bn}')
                        idx.append(TC.get_gro_attribute_by_attributes('globalIdx', {'atomName': bn, 'resNum': rnum}))
                        TC.set_gro_attribute_by_attributes('z', z, {'atomName': bn, 'resNum': rnum})

        # TC.idx_lists['bondchain']=[]
        self.chain_manager = ChainManager(create_if_missing=True)
        pairs = product(idx, idx)
        for i, j in pairs:
            if i < j:
                iname = TC.get_gro_attribute_by_attributes('atomName', {'globalIdx': i})
                jname = TC.get_gro_attribute_by_attributes('atomName', {'globalIdx': j})
                if TC.Topology.bondlist.are_bonded(i, j) and iname.startswith('C') and jname.startswith('C'):
                    # this monomer has two carbon atoms capable of reacting
                    # that are bound to each other -- this means that
                    # the two originated as a double-bond.
                    # *If* there is one with three hydrogens
                    # (remember this is an activated monomer)
                    # then it is the "tail"; the other is the "head".
                    i_nH = TC.count_H(i)
                    j_nH = TC.count_H(j)
                    if i_nH == 3 and j_nH != 3:
                        # i is the tail
                        entry = [j, i]
                    elif i_nH != 3 and j_nH == 3:
                        # j is the tail
                        entry = [i, j]
                    else:
                        logger.warning(f'In molecule {self.name}, cannot identify bonded reactive head and tail atoms\nAssuming {j} is head and {i} is tail')
                        entry = [j, i]
                    # logger.debug(f'Adding {entry} to chainlist of {self.name}')
                    self.chain_manager.injest_bond(entry[0], entry[1])
                    # TC.idx_lists['bondchain'].append(entry)
        # TC.reset_grx_attributes_from_idx_list('bondchain')
        self.chain_manager.to_dataframe(TC.Coordinates.A)
        self.initialize_molecule_rings()

    def previously_parameterized(self):
        """Returns True if a gro file exists in the project molecule/parameterized directory for this molecule.

        Returns:
            bool: True if gro file found
        """
        rval = True
        for ext in ['gro']:
            rval = rval and pfs.exists(os.path.join(pfs.Dirs.molecules_parameterized, f'{self.name}.{ext}'))
        return rval

    def parameterize(self, outname='', input_structure_format='mol2', ambertools={}):
        """Manages GAFF parameterization of this molecule.

        Args:
            outname (str, optional): output file basename, defaults to ''
            input_structure_format (str, optional): input structure format, defaults to 'mol2' ('pdb' is other possibility)
            ambertools (dict): ambertools configuration directives, defaults to {}
        """
        assert os.path.exists(f'{self.name}.{input_structure_format}'),f'Cannot parameterize molecule {self.name} without {self.name}.{input_structure_format} as input'
        if outname == '':
            outname = f'{self.name}'
        if input_structure_format == 'pdb':
            # Convert PDB → mol2 before antechamber so that repeated CONECT entries
            # (which encode double bonds in PDB format) are written as proper bond
            # orders in the mol2 bond table rather than duplicated single bonds.
            from .coordinates import Coordinates
            coords = Coordinates.read_pdb(f'{self.name}.pdb')
            coords.write_mol2(f'{self.name}.mol2', molname=self.name)
            logger.info(f'Converted {self.name}.pdb → {self.name}.mol2 with bond orders from CONECT records')
            input_structure_format = 'mol2'
        GAFFParameterize(self.name, outname, input_structure_format=input_structure_format, ambertools=ambertools)
        self.load_top_gro(f'{outname}.top', f'{outname}.gro', mol2filename=f'{outname}.mol2', wrap_coords=False)
        self.initialize_molecule_rings()
        self.TopoCoord.write_tpx(f'{outname}.tpx')

    def minimize(self, outname=''):
        """Manages invocation of vacuum minimization.

        Args:
            outname (str, optional): output file basename, defaults to ''
        """
        if outname == '':
            outname = f'{self.name}'
        self.TopoCoord.vacuum_minimize(outname)

    def relax(self, relax_dict):
        """Manages invocation of MD relaxations.

        Args:
            relax_dict (dict): dictionary of simulation directives
        """
        deffnm = relax_dict.get('deffnm', f'{self.name}-relax')
        nsteps = relax_dict.get('nsteps', 10000)
        temperature = relax_dict.get('temperature', 10000)
        n = self.name
        boxsize = np.array(self.TopoCoord.maxspan()) + 2 * np.ones(3)
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix = 'single-molecule-nvt'
        pfs.checkout(pfs.Dirs.mdp_file(mdp_prefix))
        mdp_modify(f'{mdp_prefix}.mdp', {'nsteps': nsteps, 'gen-vel': 'yes', 'ref_t': temperature, 'gen-temp': temperature})
        logger.info(f'In vacuo equilibration of {self.name}.gro for {nsteps} steps at {temperature} K')
        self.TopoCoord.grompp_and_mdrun(out=deffnm, mdp=mdp_prefix, boxSize=boxsize)

    def center_coords(self, new_boxsize: np.ndarray = None):
        """Wrapper for the TopoCoord.center_coords method.

        Args:
            new_boxsize (np.ndarray, optional): new box size, defaults to None
        """
        self.TopoCoord.center_coords(new_boxsize)

    def generate(self, outname='', available_molecules={}, gaff={}, ambertools={}):
        """Manages generating topology and coordinates for self.

        Args:
            outname (str, optional): output file basename, defaults to ''
            available_molecules (dict, optional): dictionary of available molecules, defaults to {}
            gaff (dict): GAFF configuration directives, defaults to {}
            ambertools (dict): ambertools configuration directives, defaults to {}
        """
        logger.debug(f'{self.name}.generate() begins')
        if outname == '': outname = f'{self.name}'
        do_minimization = gaff.get('minimize_molecules', True)
        do_parameterization = False
        if self.generator:
            R = self.generator
            if R.stage in [reaction_stage.cure, reaction_stage.param, reaction_stage.cap, reaction_stage.repair]: do_parameterization = True
            self.TopoCoord = TopoCoord()
            logger.debug(f'Using reaction {R.name} ({str(R.stage)}) to generate {self.name} parent {self.parentname}')
            isf = 'mol2'
            resid_mapper = []
            for ri in R.reactants.values():
                logger.debug(f'Adding {ri}')
                new_reactant = deepcopy(available_molecules[ri])
                new_reactant.TopoCoord.write_mol2(filename=f'{self.name}-reactant{ri}-prebonding.mol2', molname=self.name)
                rnr = len(new_reactant.sequence)
                shifts = self.TopoCoord.merge(new_reactant.TopoCoord, self_cm=self.chain_manager, other_cm=new_reactant.chain_manager)
                # for ln in self.TopoCoord.Coordinates.A.head().to_string().split('\n'): logger.debug(ln)
                resid_mapper.append({k:v for k,v in zip(range(1,rnr+1),range(1+shifts[2],1+rnr+shifts[2]))})
            # logger.debug(f'{self.name}: resid_mapper {resid_mapper}')
            # logger.debug(f'{self.TopoCoord.idx_lists}')
            # logger.debug(f'\n{self.TopoCoord.Coordinates.A.to_string()}')
            # logger.debug(f'composite prebonded molecule in box {self.TopoCoord.Coordinates.box}')
            self.TopoCoord.write_mol2(filename=f'{self.name}-prebonding.mol2', molname=self.name)
            self.set_sequence_from_coordinates()
            bdf = _yield_bonds_as_df(R, self.TopoCoord, resid_mapper)
            # logger.debug(f'Generation of {self.name}: composite molecule has {len(self.sequence)} resids')
            # logger.debug(f'generation of {self.name}: composite molecule:\n{composite_mol.TopoCoord.Coordinates.A.to_string()}')
            self.make_bonds(bdf, available_molecules, R.stage)
            # self.TopoCoord.set_gro_attribute('reactantName',R.product)
            self.TopoCoord.set_gro_attribute('sea_idx', -1) # turn off symmetry-equivalence for multimers
            self.TopoCoord.set_gro_attribute('molecule', 1)
            self.TopoCoord.set_gro_attribute('molecule_name', self.name)
            self.write_gro_attributes(GRX_ATTRIBUTES, f'{R.product}.grx')
            # if pfs.exists(f'molecules/inputs/{self.name}.mol2'): # an override structure is present
            #     logger.debug(f'Using override input molecules/inputs/{self.name}.{isf} as a generator')
            #     pfs.checkout(f'molecules/inputs/{self.name}.{isf}')
            # else:
            self.TopoCoord.write_mol2(filename=f'{self.name}.mol2', molname=self.name)
            if not do_parameterization:
                self.TopoCoord.write_gro(f'{self.name}.gro', grotitle=self.name)
                self.TopoCoord.write_top(f'{self.name}.top')
            # if pfs.exists(f'molecules/inputs/{self.name}.pdb'): # an override structure is present
            #     isf='pdb'
            #     logger.debug(f'Using override input molecules/inputs/{self.name}.{isf} as a generator')
            #     pfs.checkout(f'molecules/inputs/{self.name}.{isf}')
        else:
            # this molecule was not assigned a generator: implies it is a monomer and must have a parameterization
            input_structure_formats = ['mol2','pdb']
            isf = None
            for isf in input_structure_formats:
                if pfs.exists(f'{pfs.Dirs.molecules_inputs}/{self.name}.{isf}'):
                    logger.debug(f'Using input {pfs.Dirs.molecules_inputs}/{self.name}.{isf} as a generator')
                    pfs.checkout(f'{pfs.Dirs.molecules_inputs}/{self.name}.{isf}')
                    break
            assert isf, 'Error: no valid input structure file found'
            do_parameterization = True

        reactantName = self.name
        if do_parameterization:
            self.parameterize(outname, input_structure_format=isf, ambertools=ambertools)
        else:
            if self.name != self.parentname:
                logger.info(f'Built {self.name} using topology of {self.parentname}; copying {self.parentname}.top to {self.name}.top')
                self.load_top_gro(f'{self.parentname}.top', f'{self.name}.gro', tpxfilename=f'{self.parentname}.tpx', wrap_coords=False)
                shutil.copy(f'{self.parentname}.top', f'{self.name}.top')
                shutil.copy(f'{self.parentname}.grx', f'{self.name}.grx')
                shutil.copy(f'{self.parentname}.tpx', f'{self.name}.tpx')
                # The isomer carries the parent's parameters, so it carries the
                # parent's provenance too; without this it would look like a
                # library entry of unknown origin on every subsequent run.
                parent_key = paramcache.key_filename(self.parentname)
                if os.path.exists(parent_key):
                    shutil.copy(parent_key, paramcache.key_filename(self.name))

        if do_minimization:
            self.minimize(outname)
        self.set_sequence_from_coordinates()
        if not self.generator:
            self.TopoCoord.set_gro_attribute('reactantName', reactantName)
            self.initialize_monomer_grx_attributes()
            self.write_gro_attributes(GRX_ATTRIBUTES, f'{reactantName}.grx')
        else:
            if do_parameterization or self.name != self.parentname:
                grx=f'{reactantName}.grx'
                if (os.path.exists(grx)):
                    self.TopoCoord.read_gro_attributes(grx)
                #self.reset_chains_from_attributes()
        # logger.debug(f'{self.name} gro\n{self.TopoCoord.Coordinates.A.to_string()}')
        self.prepare_new_bonds(available_molecules=available_molecules)

        # for ln in self.TopoCoord.Coordinates.A.head().to_string().split('\n'): logger.debug(ln)
        logger.info(f'{self.name}: {self.get_molecular_weight():.2f} g/mol')
        logger.debug('Done.')

    def get_molecular_weight(self):
        """Returns the molecular weight of self.

        Returns:
            float: molecular weight in g/mol
        """
        mass = self.TopoCoord.Topology.total_mass(units='gromacs') # g
        return mass

    def prepare_new_bonds(self, available_molecules={}):
        """Populates the bond templates and reaction bonds for self.

        Args:
            available_molecules (dict, optional): dictionary of available molecules, defaults to {}
        """
        # logger.debug(f'set_reaction_bonds: molecules {list(available_molecules.keys())}')
        R = self.generator
        if not R:
            return
        self.reaction_bonds = []
        self.bond_templates = []
        TC = self.TopoCoord
        # logger.debug(f'prepare_new_bonds {self.name}: chainlists {TC.idx_lists["bondchain"]}')
        for bondrec in R.bonds:
            atom_keys = bondrec['atoms']
            order = bondrec['order']
            assert len(atom_keys) == 2
            atomrecs = [R.atoms[x] for x in atom_keys]
            atom_names = [x['atom'] for x in atomrecs]
            reactant_keys = [x['reactant'] for x in atomrecs]
            in_reactant_resids = [x['resid'] for x in atomrecs]
            if reactant_keys[0] == reactant_keys[1]:  # this is an intraresidue bond
                reactant_names = [R.reactants[reactant_keys[0]]]
            else:
                reactant_names = [R.reactants[x] for x in reactant_keys]
            reactant_sequences = [available_molecules[x].sequence for x in reactant_names]
            product_sequence = []
            for seq in reactant_sequences:
                product_sequence.extend(seq)
            # logger.debug(f'product_sequence {product_sequence}')
            sequence_residue_idx_origins = [0, 0]
            if len(reactant_sequences) == 2:
                sequence_residue_idx_origins[1] = len(reactant_sequences[0])
            in_product_resids = [in_reactant_resids[x] + sequence_residue_idx_origins[x] for x in [0, 1]]
            # logger.debug(f'in_product_resids {in_product_resids}')
            in_product_resnames = [product_sequence[in_product_resids[x] - 1] for x in [0, 1]]
            atom_idx = [TC.get_gro_attribute_by_attributes('globalIdx', {'resNum': in_product_resids[x], 'atomName': atom_names[x]}) for x in [0, 1]]
            # logger.debug(f'{R.name} names {atom_names} in_product_resids {in_product_resids} idx {atom_idx}')
            bystander_resids, bystander_resnames, bystander_atomidx, bystander_atomnames = TC.get_bystanders(atom_idx)
            oneaway_resids, oneaway_resnames, oneaway_atomidx, oneaway_atomnames = TC.get_oneaways(atom_idx, chain_manager=self.chain_manager)
            # logger.debug(f'{self.name} bystanders {bystander_resids} {bystander_resnames} {bystander_atomidx} {bystander_atomnames}')
            # logger.debug(f'{self.name} oneaways {oneaway_resids} {oneaway_resnames} {oneaway_atomidx} {oneaway_atomnames}')
            self.reaction_bonds.append(ReactionBond(atom_idx, in_product_resids, order, bystander_resids, bystander_atomidx, oneaway_resids, oneaway_atomidx))
            intraresidue = in_product_resids[0] == in_product_resids[1]
            self.bond_templates.append(BondTemplate(atom_names, in_product_resnames, intraresidue, order, bystander_resnames, bystander_atomnames, oneaway_resnames, oneaway_atomnames, siblings=TC.get_siblings(atom_idx)))

    def idx_mappers(self, otherTC: TopoCoord, other_bond, bystanders, oneaways, uniq_atom_idx: set):
        """Computes the mapping dictionary from molecule template index to instance index in the other TopoCoord.

        Args:
            otherTC (TopoCoord): the other TopoCoord
            other_bond: 2 atom indices of the bond in the other TopoCoord
            bystanders: bystander lists, one for each reacting atom
            oneaways: oneaways, one for each atom in the bond
            uniq_atom_idx (set): set of unique atoms in template that must be mapped to instance

        Raises:
            Exception: if there is a buggy double-counting of one or more indexes

        Returns:
            tuple: two-way dictionaries of index mappers instance<->template
        """
        assert len(other_bond) == 2
        assert len(bystanders) == 2
        assert len(oneaways) == 2
        ut = uniq_atom_idx.copy()
        logger.debug(f'Template name {self.name}')
        i_idx, j_idx = other_bond
        i_resName, i_resNum, i_atomName = otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i_idx})
        j_resName, j_resNum, j_atomName = otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j_idx})
        logger.debug(f'i_idx {i_idx} i_resName {i_resName} i_resNum {i_resNum} i_atomName {i_atomName}')
        logger.debug(f'j_idx {j_idx} j_resName {j_resName} j_resNum {j_resNum} j_atomName {j_atomName}')
        # identify the template bond represented by the other_bond parameter
        ij = []
        for RB, BT in zip(self.reaction_bonds, self.bond_templates):
            temp_resids = RB.resids
            temp_iname, temp_jname = BT.names
            temp_iresname, temp_jresname = BT.resnames
            temp_bystander_resids = RB.bystander_resids
            temp_oneaway_resids = RB.oneaway_resids
            # logger.debug(f'temp_iresname {temp_iresname} temp_iname {temp_iname}')
            # logger.debug(f'temp_jresname {temp_jresname} temp_jname {temp_jname}')
            if (i_atomName, i_resName) == (temp_iname, temp_iresname):
                ij = [0, 1]
                break # found it -- stop looking
            elif (i_atomName, i_resName) == (temp_jname, temp_jresname):
                ij = [1, 0]
                break
        assert len(ij) == 2, f'Mappers using template {self.name} unable to map from instance bond {i_resName}-{i_resNum}-{i_atomName}---{j_resName}-{j_resNum}-{j_atomName}'
        inst_resids = [i_resNum, j_resNum]
        inst_resids = [inst_resids[ij[x]] for x in [0, 1]]
        inst_bystander_resids = [bystanders[ij[x]] for x in [0, 1]]
        inst_oneaway_resids = [oneaways[ij[x]] for x in [0, 1]]
        # Bystander counts on each side may legitimately differ when the
        # template (a small parameterization fragment) carries fewer
        # bystanders than the in-system instance — e.g. cyanate-ester cure,
        # where the CY+CY dimer template has no BPA bystander but every
        # BCY-embedded instance does (find_template's subset semantics
        # let the small-fragment template match in this case).  Allow
        # the template to be the shorter list on each side; the extra
        # instance bystanders are residues outside the template and don't
        # need to be mapped here.  When both lists carry chain context
        # the counts agree as before.
        # Pair region-by-region so that a length difference on one side
        # doesn't shift alignment of the other side's bystanders or the
        # oneaways.  Within each region zip naturally truncates the longer
        # (always the instance) list.
        instdf = otherTC.Coordinates.A
        tempdf = self.TopoCoord.Coordinates.A
        inst2temp = {}
        temp2inst = {}
        resid_pairs = list(zip(inst_resids, temp_resids))
        resid_pairs.extend(zip(inst_bystander_resids[0], temp_bystander_resids[0]))
        resid_pairs.extend(zip(inst_bystander_resids[1], temp_bystander_resids[1]))
        resid_pairs.extend(zip(inst_oneaway_resids, temp_oneaway_resids))
        for inst, temp in resid_pairs:
            if inst and temp:  # None's in the bystander lists and oneaways lists should be ignored
                logger.debug(f'map inst resid {inst} to template resid {temp}')
                idf = instdf[instdf['resNum'] == inst][['globalIdx','atomName']].copy()
                # logger.debug(f'idf res {inst}:\n{idf.to_string()}')
                tdf = tempdf[tempdf['resNum'] == temp][['globalIdx','atomName']].copy()
                # logger.debug(f'tdf res {temp}:\n{tdf.to_string()}')
                tdf = tdf.merge(idf, on='atomName', how='inner', suffixes=('_template','_instance'))
                # logger.debug(f'merged\n{tdf.to_string()}')
                for i, r in tdf.iterrows():
                    temp_idx = r['globalIdx_template']
                    inst_idx = r['globalIdx_instance']
                    # logger.debug(f't {temp_idx} <-> i {inst_idx}')
                    # only map template atoms that are identified in the passed in set
                    if temp_idx in ut:
                        ut.remove(temp_idx)
                        temp2inst[temp_idx] = inst_idx
                        inst2temp[inst_idx] = temp_idx
                    if temp_idx in temp2inst and temp2inst[temp_idx] != inst_idx:
                        raise Exception(f'Error: temp_idx {temp_idx} already claimed in temp2inst; bug')
        assert len(inst2temp) == len(temp2inst), f'Error: could not establish two-way dict of atom globalIdx'
        return (inst2temp, temp2inst)

    def get_angles_dihedrals(self, bond):
        """Returns copies of selections from the Topology interaction-type dataframes that contain the two atoms indicated in the bond.

        Args:
            bond: 2-element list-like container of ints

        Raises:
            Exception: if a NaN is found in any selection

        Returns:
            tuple: copies of the dataframe selections for angles, dihedrals, and 1-4 pairs
        """
        ai, aj = bond
        d = self.TopoCoord.Topology.D['angles']
        ad = d[((d.ai == ai) & (d.aj == aj)) |
               ((d.ai == aj) & (d.aj == ai)) |
               ((d.aj == ai) & (d.ak == aj)) |
               ((d.aj == aj) & (d.ak == ai))].copy()
        d = self.TopoCoord.Topology.D['dihedrals']
        td = d[((d.ai == ai) & (d.aj == aj)) |
               ((d.ai == aj) & (d.aj == ai)) |
               ((d.aj == ai) & (d.ak == aj)) |
               ((d.aj == aj) & (d.ak == ai)) |
               ((d.ak == ai) & (d.al == aj)) |
               ((d.ak == aj) & (d.al == ai))].copy()
        check = True
        for a in ['ai', 'aj', 'ak', 'al']:
            check = check and td[a].isnull().values.any()
        if check:
            logger.error('NAN in molecule/dihedrals')
            raise ValueError('NAN in molecule/dihedrals')

        d = self.TopoCoord.Topology.D['pairs']
        paird = pd.DataFrame()
        for ai, al in zip(td.ai, td.al):
            tpair = d[((d.ai == ai) & (d.aj == al)) |
                      ((d.ai == al) & (d.aj == ai))].copy()
            paird = pd.concat((paird, tpair), ignore_index=True)
        check = True
        for a in ['ai', 'aj']:
            check = check and paird[a].isnull().values.any()
        if check:
            logger.error('NAN in molecule/pairs')
            raise ValueError('NAN in molecule/pairs')
        return ad, td, paird

    def get_resname(self, internal_resid):
        """Returns the residue name at position internal_resid in the molecule's sequence.

        Args:
            internal_resid (int): molecule-internal residue index

        Returns:
            str: residue name
        """
        return self.sequence[internal_resid - 1]

    # def inherit_attribute_from_reactants(self,attribute,available_molecules,increment=True,no_increment_if_negative=True):
    #     """inherit_attribute_from_reactants populate certain atom attributes in molecule from its constituent reactants

    #     :param attribute: _description_
    #     :type attribute: _type_
    #     :param available_molecules: _description_
    #     :type available_molecules: _type_
    #     :param increment: _description_, defaults to True
    #     :type increment: bool, optional
    #     :param no_increment_if_negative: _description_, defaults to True
    #     :type no_increment_if_negative: bool, optional
    #     """
    #     adf=self.TopoCoord.Coordinates.A
    #     ordered_attribute_idx=[]
    #     curr_max=0
    #     # logger.debug(f'{self.name}({adf.shape[0]}) inheriting {attribute} from {self.sequence}')
    #     # logger.debug(f'available molecules {list(available_molecules.keys())}')
    #     for i,r in enumerate(self.sequence):
    #         '''
    #         for this residue number, read the list of unique atom names
    #         '''
    #         namesinres=list(adf[adf['resNum']==(i+1)]['atomName'])
    #         '''
    #         access coordinates of standalone residue template with this name 'r' on the list of available molecules
    #         '''
    #         rdf=available_molecules[r].TopoCoord.Coordinates.A
    #         '''
    #         get the attribute values from residue template
    #         '''
    #         x=list(rdf[rdf['atomName'].isin(namesinres)][attribute])
    #         # logger.debug(f'{r}->{len(x)}')
    #         '''
    #         increment these attribute value based on residue number in this molecule
    #         '''
    #         if increment:
    #             i_x=[]
    #             for y in x:
    #                 if y>0 or (y<0 and not no_increment_if_negative):
    #                     i_x.append(y+curr_max)
    #                 else:
    #                     i_x.append(y)
    #             curr_max=max(i_x)
    #         ordered_attribute_idx.extend(i_x)
    #     assert len(ordered_attribute_idx)==adf.shape[0]
    #     adf[attribute]=ordered_attribute_idx

    def merge(self, other):
        """Merges TopoCoord from other into self's TopoCoord.

        Args:
            other (Molecule): another Molecule

        Returns:
            tuple: a shift tuple (returned by Coordinates.merge())
        """
        shifts = self.TopoCoord.merge(other.TopoCoord, self_cm=self.chain_manager, other_cm=other.chain_manager)
        return shifts

    def load_top_gro(self, topfilename, grofilename, tpxfilename='', mol2filename='', wrap_coords=False, ignore_bonds=False, overwrite_coordinates=False):
        """Generates a new TopoCoord member object for this molecule by reading in a Gromacs topology file and a Gromacs gro file.

        Args:
            topfilename (str): Gromacs topology file
            grofilename (str): Gromacs gro file
            tpxfilename (str, optional): extended topology file, defaults to ''
            mol2filename (str, optional): alternative coordinate mol2 file, defaults to ''
            wrap_coords (bool): if True, wrap coordinates into the box after reading gro, defaults to False
            ignore_bonds (bool): if True, skip reading bonds from mol2, defaults to False
            overwrite_coordinates (bool): if True, overwrite existing coordinates when reading mol2, defaults to False
        """
        self.TopoCoord = TopoCoord(topfilename=topfilename, grofilename=grofilename, tpxfilename=tpxfilename, mol2filename=mol2filename, wrap_coords=wrap_coords, ignore_bonds=ignore_bonds, overwrite_coordinates=overwrite_coordinates)

    def set_gro_attribute(self, attribute, srs):
        """Sets attribute of atoms to srs (drills through to Coordinates.set_atomset_attributes()).

        Args:
            attribute (str): name of attribute
            srs: scalar or list-like attribute values in same ordering as self.A
        """
        self.TopoCoord.set_gro_attribute(attribute, srs)

    def read_gro_attributes(self, grxfilename, attribute_list=[]):
        """Reads attributes from file into self.TopoCoord.Coordinates.A.

        Args:
            grxfilename (str): name of input file
            attribute_list (list, optional): list of attributes to take, defaults to [] (take all)
        """
        self.TopoCoord.read_gro_attributes(grxfilename, attribute_list=attribute_list, chain_manager=self.chain_manager)

    def write_gro_attributes(self, attribute_list, grxfilename):
        """Writes atomic attributes to a file.

        Args:
            attribute_list (list): list of attributes to write
            grxfilename (str): name of output file
        """
        self.TopoCoord.write_gro_attributes(attribute_list, grxfilename)

    def make_bonds(self, bdf: pd.DataFrame, moldict, stage):
        """Adds new bonds to the molecule's topology and deletes any sacrificial hydrogens.

        Args:
            bdf (pd.DataFrame): pandas dataframe identifying new bonds
            moldict (dict): dictionary of available molecular templates
            stage (reaction_stage): enumerated parameter indicating reaction_stage
        """
        TC = self.TopoCoord
        explicit_sacrificial_Hs = {}
        for i, r in bdf.iterrows():
            aname, bname = [TC.get_gro_attribute_by_attributes('atomName', {'globalIdx': x}) for x in [r.ai, r.rj]]
            logger.debug(f'generating {self.name} bond {r.ri}:{aname}:{r.ai}-{r.rj}:{bname}:{r.aj} order {r.order}')
            if r.ri != r.rj:
                resid_sets = TC.get_resid_sets([r.ai, r.aj])
                hxi, hxj = self.transrot(r.ai, r.ri, r.aj, r.rj, connected_resids=resid_sets[1])
                explicit_sacrificial_Hs[i] = [hxi, hxj]
        if stage in [reaction_stage.cure, reaction_stage.param, reaction_stage.cap, reaction_stage.repair]:
            template_source = 'ambertools'
        else:
            template_source = 'internal'  # signals that a template molecule should be identified to parameterize this bond
        TC.update_topology_and_coordinates(bdf, moldict, explicit_sacH=explicit_sacrificial_Hs, template_source=template_source, chain_manager=self.chain_manager)
        self.initialize_molecule_rings()

    def transrot(self, at_idx, at_resid, from_idx, from_resid, connected_resids=[]):
        """Given a composite molecule, translates and rotates the piece downstream of the yet-to-be-created bond to minimize steric overlaps and identify the best two sacrificial hydrogens.

        Args:
            at_idx (int): global index of left-hand atom in new bond
            at_resid (int): resid of left-hand residue
            from_idx (int): global index of right-hand atom in new bond
            from_resid (int): resid of right-hand residue
            connected_resids (list, optional): list of all other resids attached to right-hand residue that should move with it, defaults to []

        Returns:
            tuple: 2-tuple containing global indices of the sacrificial hydrogens
        """
        # Rotate and translate
        if at_resid == from_resid:
            return # should never happen but JIC
        logger.debug(f'{self.name} connected resids {connected_resids}')
        TC = self.TopoCoord
        ATC = TopoCoord()
        BTC = TopoCoord()
        C = TC.gro_DataFrame('atoms')
        ATC.Coordinates.A = C[C['resNum'] == at_resid].copy()
        bresids = connected_resids.copy()
        bresids.append(from_resid)
        BTC.Coordinates.A = C[C['resNum'].isin(bresids)].copy()
        for ln in BTC.Coordinates.A.to_string().split('\n'): logger.debug(ln)
        NONROT = C[~C['resNum'].isin(bresids)].shape[0]
        logger.debug(f'{self.TopoCoord.Coordinates.A.shape[0]} atoms')
        logger.debug(f'holding {at_resid} ({NONROT})')
        logger.debug(f'rotating/translating {bresids} ({BTC.Coordinates.A.shape[0]})')
        assert self.TopoCoord.Coordinates.A.shape[0] == (NONROT + BTC.Coordinates.A.shape[0])
        mypartners = TC.Topology.bondlist.partners_of(at_idx)
        otpartners = TC.Topology.bondlist.partners_of(from_idx)
        logger.debug(f'Partners of {at_idx} {mypartners}')
        logger.debug(f'Partners of {from_idx} {otpartners}')
        myHpartners = {k: v for k, v in zip(mypartners, [C[C['globalIdx'] == i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')}
        otHpartners = {k: v for k, v in zip(otpartners, [C[C['globalIdx'] == i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')}
        myHighestH = {k: v for k, v in myHpartners.items() if v == max([k for k in myHpartners.values()], key=lambda x: int(x.split('H')[1] if x.split('H')[1] != '' else '0'))}
        otHighestH = {k: v for k, v in otHpartners.items() if v == max([k for k in otHpartners.values()], key=lambda x: int(x.split('H')[1] if x.split('H')[1] != '' else '0'))}
        assert len(myHighestH) == 1
        assert len(otHighestH) == 1
        logger.debug(f'Highest-named H partner of {at_idx} is {myHighestH}')
        logger.debug(f'Highest-named H partner of {from_idx} is {otHighestH}')
        assert len(myHpartners) > 0, f'Error: atom {at_idx} does not have a deletable H atom!'
        assert len(otHpartners) > 0, f'Error: atom {from_idx} does not have a deletable H atom!'

        Ri = TC.get_R(at_idx)
        Rj = TC.get_R(from_idx)
        logger.debug(f'Ri {at_idx} {Ri} {type(Ri)} {Ri.dtype}')
        logger.debug(f'Rj {from_idx} {Rj} {type(Rj)} {Rj.dtype}')
        # The two sacrificial H's are made to coincide, which puts at_idx, the H's
        # and from_idx on one line.  The moving piece can still turn freely about
        # that line, and that freedom is what clears a crowded site: a third
        # tetramethyl bisphenol on a triazine lands its ortho methyl inside the
        # ring at the unrotated placement.  So each H pairing is tried at a
        # range of turns, and the one farthest from the residues that stay put
        # wins.  (Measuring against all of TC would include the moving piece's
        # own pre-move positions.)
        fixed = C[~C['resNum'].isin(bresids)]
        overall_maximum = (-1.e9, -1, -1, 0.0)
        coord_trials = {}
        for myH, myHnm in myHpartners.items():  # keys are globalIdx's, values are names
            coord_trials[myH] = {}
            Rh = TC.get_R(myH)
            logger.debug(f'  Rh {myH} {Rh} {Rh.dtype}')
            Rih = Ri - Rh
            Rih *= 1.0 / np.linalg.norm(Rih)
            fixed_xyz = fixed.loc[fixed['globalIdx'] != myH, ['posX', 'posY', 'posZ']].to_numpy(dtype=float)
            for otH, otHnm in otHpartners.items():
                logger.debug(f'{self.name}: Considering {myH} {otH}')
                aligned = deepcopy(BTC)
                Rk = aligned.get_R(otH)
                Rkj = Rk - Rj
                Rkj *= 1.0 / np.linalg.norm(Rkj)
                cp = np.cross(Rkj,Rih)
                c = np.dot(Rkj,Rih)
                v = np.array([[0,-cp[2],cp[1]],[cp[2],0,-cp[0]],[-cp[1],cp[0],0]])
                v2 = np.dot(v,v)
                I = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                # R is the rotation matrix that will rotate donb to align with accb
                R = I + v + v2 / (1. + c)
                aligned.rotate(R)
                Rk = aligned.get_R(otH)
                # overlap the two H atoms by translation
                aligned.translate(Rh - Rk)
                A = aligned.Coordinates.A
                moving_mask = (A['globalIdx'] != otH).to_numpy()
                xyz = A[['posX', 'posY', 'posZ']].to_numpy(dtype=float)
                best = None
                for theta in np.deg2rad(np.arange(0.0, 360.0, _TRANSROT_STEP_DEG_)):
                    trial = _rotate_about_axis(xyz, Rh, Rih, theta)
                    minD = _min_distance(trial[moving_mask], fixed_xyz)
                    if best is None or minD > best[0]:
                        best = (minD, theta, trial)
                minD, theta, trial = best
                A[['posX', 'posY', 'posZ']] = trial
                coord_trials[myH][otH] = aligned
                logger.debug(f'{self.name}: {myH}/{otH} best turn {np.rad2deg(theta):.0f} deg, minD {minD}')
                if minD > overall_maximum[0]:
                    overall_maximum = (minD, myH, otH, theta)
        logger.debug(f'{self.name}: overall_maximum {overall_maximum}')
        minD, myH, otH, _ = overall_maximum
        if minD < _TRANSROT_WARN_NM_:
            logger.warning(f'{self.name}: the best placement found for the new bond still leaves two atoms '
                           f'{minD*10:.2f} A apart; a semi-empirical charge method (bcc, abcg2) may fail to converge '
                           f'on this starting geometry')
        BTC = coord_trials[myH][otH]
        TC.overwrite_coords(BTC)
        TC.swap_atom_names(myH, list(myHighestH.keys())[0])
        TC.swap_atom_names(otH, list(otHighestH.keys())[0])
        return myH, otH

    def atoms_w_same_attribute_as(self, find_dict={}, same_attribute='', return_attribute=''):
        """Returns a list of atom attributes named in return_attribute from atoms that share an attribute named in same_attribute with the atom identified by find_dict.

        Args:
            find_dict (dict, optional): dictionary of attribute:value pairs that should uniquely identify an atom, defaults to {}
            same_attribute (str, optional): name of attribute used to screen atoms, defaults to ''
            return_attribute (str, optional): attribute value to return a list of from the screened atoms, defaults to ''

        Returns:
            list: list of attribute values
        """
        att_val = self.TopoCoord.get_gro_attribute_by_attributes(same_attribute, find_dict)
        return self.TopoCoord.get_gro_attributelist_by_attributes(return_attribute, {same_attribute: att_val})

    def flip_stereocenters(self, idxlist):
        """Flips stereochemistry of atoms in idxlist.

        Args:
            idxlist (list): global indices of chiral atoms
        """
        self.TopoCoord.flip_stereocenters(idxlist)

    def rotate_bond(self, a, b, deg):
        """Rotates all atoms in molecule on b-side of a-b bond by deg degrees.

        Args:
            a (int): index of a
            b (int): index of b
            deg (float): angle of rotation in degrees
        """
        TC = self.TopoCoord
        A = TC.Coordinates.A
        branch = TopoCoord()
        bl = deepcopy(self.TopoCoord.Topology.bondlist)
        branchidx = bl.half_as_list((a, b), 99)
        ra = TC.get_R(a)
        rb = TC.get_R(b)
        Rab = ra - rb
        rab = Rab / np.linalg.norm(Rab)
        O = rb
        TC.translate(-1 * O)
        branch.Coordinates.A = A[A['globalIdx'].isin(branchidx)].copy()
        # do stuff
        rx, ry, rz = rab
        az = np.arctan2(ry, rx)   # azimuthal angle of bond axis (radians)
        ay = np.arctan2(np.sqrt(rx**2 + ry**2), rz)  # polar angle of bond axis (radians)
        M = (Matrix4()
           .rot(np.degrees(az),  'z')
           .rot(np.degrees(ay),  'y')
           .rot(deg,             'z')
           .rot(np.degrees(-ay), 'y')
           .rot(np.degrees(-az), 'z'))
        branch.rotate(M.m[:3,:3])
        bdf = branch.Coordinates.A
        A.loc[A['globalIdx'].isin(branchidx), ['posX', 'posY', 'posZ']] = bdf[['posX', 'posY', 'posZ']]
        TC.translate(O)

    # def sea_of(self,idx):
    #     clu=self.atoms_w_same_attribute_as(find_dict={'globalIdx':idx},
    #                                             same_attribute='sea_idx',
    #                                             return_attribute='globalIdx')
    #     return list(clu)

    def generate_stereoisomers(self):
        """Generates list of Molecule shells, one for each stereoisomer.

        Returns:
            None: returns early if no stereoisomers need to be generated
        """
        if self.TopoCoord.Topology.D['atoms'].shape[0] == 0: return  # self has not yet acquired topology/coordinates
        if len(self.stereoisomers) == 0: return
        flip = [[0, 1] for _ in range(len(self.stereocenters))]
        st_idx = [self.TopoCoord.get_gro_attribute_by_attributes('globalIdx', {'atomName': n}) for n in self.stereocenters]
        P = product(*flip)
        next(P) # one with no flips is the original molecule, so skip it
        for p in P:
            si_name = self.name + '-S' + ''.join([str(_) for _ in p])
            if not si_name in self.stereoisomers:
                logger.debug(f'{si_name} not found in dict of stereoisomers of {self.name}')
            logger.debug(f'Stereocenter sequence {p} generates stereoisomer {si_name}')
            M = self.stereoisomers[si_name]
            M.origin = self.origin
            M.TopoCoord = deepcopy(self.TopoCoord)
            fsc = [st_idx[i] for i in range(len(self.stereocenters)) if p[i]]
            M.flip_stereocenters(fsc)
            M.TopoCoord.write_gro(f'{si_name}.gro')

    def generate_conformers(self):
        """generate_conformers generates this molecule's conformer instances using either gromacs or obabel
        """
        # only generates gro files
        default_gromacs_params = {'ensemble': 'nvt', 'temperature': 600, 'ps': 100, 'begin_at': 50}
        # if self.nconformers==0: return
        cd = self.conformers_dict
        if not cd: return
        logger.debug(f'{self.name} conformer_dict {cd}')
        self.nconformers = cd['count']
        minimize = cd.get('minimize', False)
        generator = cd.get('generator', {})
        if generator == 'obabel' and not minimize:
            logger.debug(f'Confomers generated by obabel should be energy-minimized.  Indicate "mimimize: True" in the confomers directive for {self.name}')
            minimize = True
        if not generator: return
        logger.info(f'Generating {self.nconformers*(1+len(self.stereoisomers))} conformers for {self.name}')
        gronames = [f'{self.name}']
        nd = cd.get('nzeros', 2)
        for si in self.stereoisomers:
            gronames.append(f'{si}')
        for gro in gronames:
            pfx = f'{gro}-C'
            fmt = r'{A}{B:0'+str(nd) + r'd}'  # the trjconv command in gro_from_trr must generate these files
            cfnl = [fmt.format(A=pfx,B=x) for x in range(self.nconformers)]
            if all(os.path.exists(f'{mname}.gro') for mname in cfnl):
                logger.info(f'Reusing {self.nconformers} existing conformer gro files for {gro}')
                self.conformers.extend(cfnl)
                continue
            if generator['name'] == 'obabel':
                compfile = f'{gro}-obabel-confs.gro'
                run(f'obabel -igro {gro}.gro -O {compfile} --conformer --nconf {self.nconformers} --writeconformers')
                out,err = run(f'wc -l {compfile}')
                tok = out.split()
                lpf = int(tok[0]) // self.nconformers
                run(f'split -d -n {nd} -l {lpf} {compfile} {pfx} --additional-suffix=".gro"')
            elif generator['name'] == 'gromacs':
                params = generator.get('params', default_gromacs_params)
                begin_at = params.get('begin_at', 0.0)
                compfile = f'{gro}-gromacs-confs'
                TC = self.TopoCoord if gro == self.name else self.stereoisomers[gro].TopoCoord
                TC.vacuum_simulate(outname=f'{compfile}', nsamples=cd['count'], params=params)
                gro_from_trr(compfile, nzero=nd, outpfx=pfx, b=begin_at)
            # os.remove(f'{gro}-confs.gro')
            for mname in cfnl:
                assert os.path.exists(f'{mname}.gro'), f'Error: Conformer coordinates file {mname}.gro not found'
            logger.debug(f'Conformer coordinate filenames {cfnl}')
            self.conformers.extend(cfnl)
        if minimize:
            saveTC = deepcopy(self.TopoCoord)
            for mname in self.conformers:
                self.TopoCoord.read_gro(f'{mname}.gro')
                logger.info(f'Minimizing conformer {mname}')
                self.TopoCoord.vacuum_minimize(outname=f'{mname}')
            self.TopoCoord = saveTC

MoleculeDict = dict[str,Molecule]
MoleculeList = list[Molecule]