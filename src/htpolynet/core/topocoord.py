"""Class for jointly handling Topology and Coordinate objects.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging
import os
import shutil

from copy import deepcopy
from enum import Enum
from itertools import product

import networkx as nx
import numpy as np
import pandas as pd

from ..core import projectfilesystem as pfs
from ..core.bondtemplate import BondTemplate, ReactionBond
from ..core.coordinates import Coordinates
from ..core.topology import Topology
from ..utils.convergence import series_converged
from ..external.gromacs import grompp_and_mdrun,mdp_get, mdp_modify, gmx_energy_trace
from ..geometry.matrix4 import Matrix4
from ..io.gro import GRX_ATTRIBUTES, GRX_GLOBALLY_UNIQUE, GRX_UNSET_DEFAULTS

logger=logging.getLogger(__name__)

class BTRC(Enum):
    """Bond test return codes: bond tests are applied to those bond-candidates that are within search radius of each other.
    """
    passed = 0
    failed_pierced_ring = 1       # candidate bond pierces a ring
    failed_shortcircuit = 2       # candidate bond creates a short-circuit
    unset = 99

class TopoCoord:
    """Container for Topology and Coordinates, along with methods that
        use either or both of them
    """
    def __init__(self, topfilename='', tpxfilename='', grofilename='', grxfilename='', mol2filename='', system_name='htpolynet', wrap_coords=False, ignore_bonds=False, overwrite_coordinates=False):
        """Constructor method for TopoCoord.

        Args:
            topfilename (str): name of Gromacs-format topology file (top), defaults to ''
            tpxfilename (str): name of custom-format extended topology file (tpx), defaults to ''
            grofilename (str): name of Gromacs-format coordinate file (gro), defaults to ''
            grxfilename (str): name of custom-format extended coordinate file (grx), defaults to ''
            mol2filename (str): name of SYBYL MOL2-format coordinate/bonds file, defaults to ''
            system_name (str): name of the system, defaults to 'htpolynet'
            wrap_coords (bool): if True, wrap coordinates into the box after reading gro, defaults to False
            ignore_bonds (bool): if True, skip reading bonds from mol2, defaults to False
            overwrite_coordinates (bool): if True, overwrite existing coordinates when reading mol2, defaults to False
        """
        self.files = {}
        self.files['top'] = os.path.abspath(topfilename)
        self.files['tpx'] = os.path.abspath(tpxfilename)
        self.files['gro'] = os.path.abspath(grofilename)
        self.files['grx'] = os.path.abspath(grxfilename)
        self.files['mol2'] = os.path.abspath(mol2filename)
        self.grxattr = []

        # self.idx_lists={}
        # self.idx_lists['bondchain']=[]
        if grofilename != '':
            self.read_gro(grofilename, wrap_coords=wrap_coords)
        else:
            self.Coordinates = Coordinates()  # empty
        if topfilename != '':
            self.read_top(topfilename)
        else:
            self.Topology = Topology(system_name=system_name) # empty
        if tpxfilename != '':
            self.read_tpx(tpxfilename)
        if mol2filename != '':
            self.read_mol2(mol2filename, ignore_bonds=ignore_bonds, overwrite_coordinates=overwrite_coordinates)
        if grxfilename != '':
            self.read_gro_attributes(grxfilename)
        self.Coordinates.claim_parent(self)

    @classmethod
    def from_top_gro(cls, top, gro):
        """Constructs a TopoCoord object from a top / gro pair.

        Args:
            top (str): Gromacs topology file
            gro (str): Gromacs gro coordinate file

        Returns:
            TopoCoord: TopoCoord instance with the top / gro read in
        """
        return cls(topfilename=top, grofilename=gro)

    def _sacH(self, ai, aj, rename=False):
        """Finds the two H atoms (one bound to ai, one to aj) that are closest to each other.

        Args:
            ai (int): index of one bond atom
            aj (int): index of the other bond atom
            rename (bool): if True, renames remaining H atoms so highest-sorted names appear sacrificed, defaults to False

        Returns:
            list: global indexes of the two sacrificial H atoms
        """
        bondlist = self.Topology.bondlist
        i_partners = bondlist.partners_of(ai)
        j_partners = bondlist.partners_of(aj)
        A = self.Coordinates.A
        i_Hpartners = {k: v for k, v in zip(i_partners, [A[A['globalIdx'] == i]['atomName'].values[0] for i in i_partners]) if v.startswith('H')}
        j_Hpartners = {k: v for k, v in zip(j_partners, [A[A['globalIdx'] == i]['atomName'].values[0] for i in j_partners]) if v.startswith('H')}
        assert len(i_Hpartners) > 0, f'Error: atom {ai} does not have a deletable H atom ligand!'
        assert len(j_Hpartners) > 0, f'Error: atom {aj} does not have a deletable H atom ligand!'
        minHH = (1.e9, -1, -1)
        for ih in i_Hpartners:
            RiH = self.Coordinates.get_R(ih)
            for jh in j_Hpartners:
                RjH = self.Coordinates.get_R(jh)
                RijH = RiH - RjH
                rijh = np.sqrt(RijH.dot(RijH))
                if rijh < minHH[0]:
                    minHH = (rijh, ih, jh)
        if rename:
            i_avails = list(sorted(i_Hpartners.values(), key=lambda x: int(x.split('H')[1] if x.split('H')[1] != '' else '0')))[:-1]
            j_avails = list(sorted(j_Hpartners.values(), key=lambda x: int(x.split('H')[1] if x.split('H')[1] != '' else '0')))[:-1]
            logger.debug(f'i_avails {i_avails}')
            logger.debug(f'j_avails {j_avails}')
            del i_Hpartners[ih]
            del j_Hpartners[jh]
            Top = self.Topology.D['atoms']
            Cor = self.Coordinates.A
            for h in i_Hpartners:
                i_Hpartners[h] = i_avails.pop(0)
                Top.iloc[h-1, Top.columns == 'atom'] = i_Hpartners[h]
                Cor.iloc[h-1, Cor.columns == 'atomName'] = i_Hpartners[h]
                logger.debug(f'i: changed name of {h} to {i_Hpartners[h]}')
            for h in j_Hpartners:
                j_Hpartners[h] = j_avails.pop(0)  
                Top.iloc[h-1, Top.columns == 'atom'] = j_Hpartners[h]
                Cor.iloc[h-1, Cor.columns == 'atomName'] = j_Hpartners[h]
                logger.debug(f'j: changed name of {h} to {j_Hpartners[h]}')
        return [ih, jh]

    def _find_sacrificial_H(self, pairs, rename=False, explicit_sacH={}):
        """Identifies all appropriate sacrificial hydrogen atoms for the given bond pairs.

        Args:
            pairs (list): list of (ai, aj, order) tuples indicating bonds
            rename (bool): if True, renames remaining H atoms so highest-order names appear sacrificed, defaults to False
            explicit_sacH (dict): pre-chosen sacrificial H atoms keyed by pair index, defaults to {}

        Returns:
            list: global atom indices to delete
        """
        idx_to_delete = []
        for i,b in enumerate(pairs):
            if not i in explicit_sacH:
                ai, aj, o = b
                idx_to_delete.extend(self._sacH(ai, aj, rename=rename))
            else:
                idx_to_delete.extend(explicit_sacH[i])
        return idx_to_delete

    def make_bonds(self, pairs, explicit_sacH={}, chain_manager=None):
        """Adds new bonds to the global topology.

        Args:
            pairs (list): list of pairs of atom global indices indicating each new bond
            explicit_sacH (dict): explicit mapping of sacrificial H atoms, defaults to {}
            chain_manager: optional ChainManager owned by the caller; updated in place if provided

        Returns:
            list: list of indexes of atoms that must now be deleted (sacrificial H's)
        """
        idx_to_ignore = self._find_sacrificial_H(pairs, explicit_sacH=explicit_sacH)
        logger.debug(f'idx_to_ignore {idx_to_ignore}')
        self.Topology.add_bonds(pairs)
        if chain_manager is not None:
            logger.debug(f'Prior to injesting bonds, chainmanager reports {len(chain_manager.chains)} chains')
            chain_manager.injest_bonds(pairs)
            logger.debug(f'After injesting bonds, chainmanager reports {len(chain_manager.chains)} chains')
            chain_manager.to_dataframe(self.Coordinates.A)
        # self.bondchainlist_update(pairs,msg='TopoCoord.make_bonds')
        self.Topology.null_check(msg='add_bonds')
        rename = False if len(explicit_sacH) > 0 else False
        idx_to_delete = self._find_sacrificial_H(pairs, rename=rename, explicit_sacH=explicit_sacH)
        assert type(idx_to_delete) == list
        return idx_to_delete

    def add_restraints(self, pairdf, typ=6):
        """Adds bonds of type typ to the topology from the dataframe pairdf.

        Args:
            pairdf (pandas.DataFrame): dataframe with ai, aj, initial-distance columns
            typ (int): bond type to add, defaults to 6 (non-topological restraint)
        """
        self.Topology.add_restraints(pairdf, typ=typ)

    def delete_atoms(self, atomlist):
        """Deletes atoms from both the Topology and Coordinates instances.

        Args:
            atomlist (list): list of global indexes of atoms to delete

        Returns:
            dict: old-to-new atom global index mapper dictionary resulting from reindexing remaining atoms to make sure global indexes are sequential
        """
        # logger.debug(f'delete_atoms: {atomlist}')
        self.Coordinates.delete_atoms(atomlist)
        idx_mapper = self.Topology.delete_atoms(atomlist)
        assert type(idx_mapper) == dict
        # logger.debug(f'idx_mapper: {idx_mapper}')
        # for list_name in ['bondchain']:
        #     # logger.debug(f'remapping idxs in stale {list_name} lists: {self.idx_lists[list_name]}')
        #     self.reset_idx_list_from_grx_attributes(list_name)
            # self.remap_idx_list(list_name,idx_mapper)
        # logger.debug(f'finished')
        return idx_mapper

    def count_H(self,idx):
        """Counts the number of hydrogens bound to atom with index idx; any atom whose name begins with the character 'H' or 'h' is assumed to be a hydrogen.

        Args:
            idx (int): atom global index

        Returns:
            int: number of H atoms bound to this atom
        """
        aneigh = self.Topology.bondlist.partners_of(idx)
        aneighnames = [self.get_gro_attribute_by_attributes('atomName', {'globalIdx': y}) for y in aneigh]
        anH = sum([int(x.upper().startswith('H')) for x in aneighnames])
        return anH

    def map_from_templates(self, bdf, moldict, overcharge_threshhold=0.1, chain_manager=None, adjust_charges=True):
        """Updates angles, pairs, dihedrals, atom types, and charges, based on product templates associated with each bond in 'bdf'.

        Args:
            bdf (pandas.DataFrame): dataframe with columns 'ai', 'aj', 'reactantName'
            moldict (dict): dictionary of template Molecules keyed by name
            overcharge_threshhold (float): threshold for charge adjustment, defaults to 0.1
            chain_manager: optional ChainManager owned by the caller; passed to get_oneaways
            adjust_charges (bool): if True (the default), spread any resulting system overcharge
                uniformly over every mapped atom; a caller that deletes atoms afterwards, and so
                can only settle the charge once those are gone, passes False

        Returns:
            list: global indices of every system atom that received template attributes

        Raises:
            Exception: if nan found in any attribute of any new system angle
            Exception: if nan found in any attribute of any new system dihedral
            Exception: if nan found in any attribute of any new system pair that came along with a dihedral
            Exception: if nan found in ai attribute of any template pair
            Exception: if nan found in aj attribute of any template pair
            Exception: if nan found in any system pair
        """
        atdf = self.Topology.D['atoms']
        # ij=self.Topology.D['bondtypes'].set_index(['i','j'])
        #mb=self.D['mol2_bonds']
        # bmi=self.Topology.D['bonds'].set_index(['ai','aj']).sort_index().index
        grodf = self.Coordinates.A
        grodf['old_reactantName'] = grodf['reactantName'].copy()
        logger.debug(f'Mapping {bdf.shape[0]} bonds.')
        premapping_total_charge = self.Topology.total_charge()
        logger.debug(f'Must compensate for an overcharge of {premapping_total_charge:.4f}')
        mapped_inst_atoms = []
        for b in bdf.itertuples():
            logger.debug(f'Mapping bond {b}')
            # for ln in b.to_string().split('\n'):
            #     logger.debug(f'  -> {ln}')
            bb = [b.ai,b.aj]
            order = b.order
            names = [self.get_gro_attribute_by_attributes('atomName', {'globalIdx': x}) for x in bb]
            resnames = [self.get_gro_attribute_by_attributes('resName', {'globalIdx': x}) for x in bb]
            logger.debug(f'{resnames}')
            resids = [self.get_gro_attribute_by_attributes('resNum', {'globalIdx': x}) for x in bb]
            # this is the product name of the reaction used to identify this bond
            product_name = b.reactantName
            if product_name in moldict:
                P=moldict[product_name]
                if P.is_reactant:
                    logger.debug(f'{P.name} is a reactant; updating \"reactantName\" attributes of {bb}')
                    # if the product of this reaction is also a reactant in a cure reaction, we should
                    # change the reactantName attribute of the two atoms to match this so either of
                    # these atoms can be found in a later bond search
                    for j in bb:
                        self.set_gro_attribute_by_attributes('reactantName',product_name,{'globalIdx':j})
                else:
                    logger.debug(f'{P.name} is not a reactant; no update of \"reactantName\" attributes')
            bystander_resids, bystander_resnames, bystander_atomidx, bystander_atomnames = self.get_bystanders(bb)
            oneaway_resids, oneaway_resnames, oneaway_atomidx, oneaway_atomnames = self.get_oneaways(bb, chain_manager=chain_manager)
            intraresidue = resids[0] == resids[1]
            BT = BondTemplate(names, resnames, intraresidue, order, bystander_resnames, bystander_atomnames, oneaway_resnames, oneaway_atomnames, siblings=self.get_siblings(bb))
            RB = ReactionBond(bb, resids, order, bystander_resids, bystander_atomidx, oneaway_resids, oneaway_atomidx)
            logger.debug(f'apparent bond template {str(BT)}')
            logger.debug(f'apparent bond instance {str(RB)}')
            T, rb, reverse_bond = find_template(BT, moldict)
            if reverse_bond: RB.reverse()
            temp_i_idx, temp_j_idx = rb.idx
            d = T.TopoCoord.Topology.D['bonds']
            # copy all bond records matching these two bonds; should be only one!!
            d = d[((d.ai == temp_i_idx) & (d.aj == temp_j_idx)) |
                ((d.ai == temp_j_idx) & (d.aj == temp_i_idx))].copy()
            assert d.shape[0] == 1, f'Using {T.name} is sent inst-bond {RB.idx[0]}-{RB.idx[1]} which is claimed to map to {temp_i_idx}-{temp_j_idx}, but no such unique bond is found:\n{T.TopoCoord.Topology.D["bonds"].to_string()}'
            ''' check passed '''
            temp_angles, temp_dihedrals, temp_pairs = T.get_angles_dihedrals((temp_i_idx, temp_j_idx))
            logger.debug(f'Mapping {temp_angles.shape[0]} angles, {temp_dihedrals.shape[0]} dihedrals, and {temp_pairs.shape[0]} pairs from template {T.name}')
            # determine the set of unique atoms in template that must be mapped to instance -- it is all
            # atoms involved in these bonded interactions
            uniq_temp_idx = set()
            for df in [temp_angles, temp_dihedrals, temp_pairs]:
                for label in ['ai','aj','ak','al']:
                    if label in df:
                        uniq_temp_idx = uniq_temp_idx.union(set(df[label].to_list()))
            logger.debug(f'Template atom indexes that must be mapped: {uniq_temp_idx}')
            # get the bidirectional instance<->template mapping dictionaries
            inst2temp, temp2inst = T.idx_mappers(self, RB.idx, RB.bystander_resids, RB.oneaway_resids, uniq_temp_idx)
            # some hard checks on compatibility of the dicts
            assert len(inst2temp) == len(temp2inst)
            check = True
            for k, v in inst2temp.items():
                check = check and (k == temp2inst[v])
            for k, v in temp2inst.items():
                check = check and (k == inst2temp[v])
            assert check,f'Error: bidirectional dicts are incompatible; bug\n{inst2temp}\b{temp2inst}'
            # logger.debug(f'inst2temp {inst2temp}')
            # logger.debug(f'temp2inst {temp2inst}')
            # RB, not bb: a template that matched the reversed bond has reversed RB to line up with it
            i_idx, j_idx = RB.idx
            _temp_i_idx, _temp_j_idx = inst2temp[i_idx], inst2temp[j_idx]
            assert temp_i_idx == _temp_i_idx, f'mapping mismatch -- bug'
            assert temp_j_idx == _temp_j_idx, f'mapping mismatch -- bug'

            need_new_bond_parameters = False
            total_dcharge = 0.0
            temp_atdf = T.TopoCoord.Topology.D['atoms']
            for temp_atom, inst_atom in temp2inst.items():
                assert inst_atom in atdf['nr'].values, f'Error: mapped atom {inst_atom} not found in [ atoms ]'
                inst_type, inst_charge, inst_name, inst_resn, inst_rnam = atdf[atdf['nr'] == inst_atom][['type','charge','atom','resnr','residue']].values[0]
                temp_type, temp_charge, temp_name, temp_resn, temp_rnam = temp_atdf[temp_atdf['nr'] == temp_atom][['type','charge','atom','resnr','residue']].values[0]
                # logger.debug(f'temp {temp_atom} {temp_name} {temp_rnam} {temp_resn} {temp_type} {temp_charge}')
                # logger.debug(f'inst {inst_atom} {inst_name} {inst_rnam} {inst_resn} {inst_type} {inst_charge}')
                if inst_type != temp_type:
                    logger.debug(f'changing type of inst atom {inst_atom} ({inst_resn} {inst_rnam} {inst_name}) from {inst_type} to {temp_type}')
                    atdf.loc[atdf['nr'] == inst_atom, 'type'] = temp_type
                    if inst_atom == b.ai or inst_atom == b.aj:
                        # changed type of one or both of the bond atoms
                        need_new_bond_parameters = True
                if inst_charge != temp_charge:
                    logger.debug(f'charge {inst_atom} ({inst_resn} {inst_rnam} {inst_name}) from {inst_charge} to {temp_charge}')
                    atdf.loc[atdf['nr'] == inst_atom, 'charge'] = temp_charge
                    dcharge = temp_charge - inst_charge
                    total_dcharge += dcharge
            mapped_inst_atoms.extend(list(temp2inst.values()))
            if need_new_bond_parameters:
                self.Topology.reset_override_from_type('bonds', 'bondtypes', inst_idx=(b.ai, b.aj))

            # get all angle, dihedrals, and pairs from template that result from the existence of the specified bond
            # temp_angles,temp_dihedrals,temp_pairs=T.get_angles_dihedrals((temp_i_idx,temp_j_idx))
            # logger.debug(f'Mapping {temp_angles.shape[0]} angles, {temp_dihedrals.shape[0]} dihedrals, and {temp_pairs.shape[0]} pairs from template {T.name}')
            # map from template atom indicies to system atom indicies in angles
            # logger.debug(f'Template angles:')
            # for ln in temp_angles.to_string().split('\n'):
            #     logger.debug(ln)
            inst_angles = temp_angles.copy()
            inst_angles.ai = temp_angles.ai.map(temp2inst)
            inst_angles.aj = temp_angles.aj.map(temp2inst)
            inst_angles.ak = temp_angles.ak.map(temp2inst)
            # Drop rows that reference template atoms with no system counterpart.
            # This happens when the cure-reactive atom's H count differs between
            # the fresh-from-SMILES cure template and the post-build system —
            # e.g. cyanate-ester CY where the build consumes one H from C before
            # cure, so a system CY has one fewer H on the cure-side than the
            # template's fresh CY does.  Template force-field entries for those
            # phantom H's describe parameters for atoms that don't exist in the
            # cured system and are correctly skipped.
            angle_mask = inst_angles[['ai','aj','ak']].notna().all(axis=1)
            n_skipped_ang = int((~angle_mask).sum())
            if n_skipped_ang:
                logger.debug(f'Skipping {n_skipped_ang} template angle(s) referencing unmapped atoms')
            inst_angles = inst_angles[angle_mask]
            # add new angles to the system topology
            d = self.Topology.D['angles']
            self.Topology.D['angles'] = pd.concat((d, inst_angles), ignore_index=True)
            # hard check for any nan's in any atom index attribute in any angle
            d = self.Topology.D['angles']
            check = True
            for a in ['ai', 'aj', 'ak']:
                check = check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in angles')
                raise ValueError('NAN in angles')

            # logger.debug(f'Template dihedrals:')
            # for ln in temp_dihedrals.to_string().split('\n'):
            #     logger.debug(ln)
            # map from template atom indicies to system atom indicies in dihedrals
            inst_dihedrals = temp_dihedrals.copy()
            inst_dihedrals.ai = temp_dihedrals.ai.map(temp2inst)
            inst_dihedrals.aj = temp_dihedrals.aj.map(temp2inst)
            inst_dihedrals.ak = temp_dihedrals.ak.map(temp2inst)
            inst_dihedrals.al = temp_dihedrals.al.map(temp2inst)
            # Same filter as angles — drop rows where any leg of the dihedral
            # references a template atom with no system counterpart.
            dh_mask = inst_dihedrals[['ai','aj','ak','al']].notna().all(axis=1)
            n_skipped_dh = int((~dh_mask).sum())
            if n_skipped_dh:
                logger.debug(f'Skipping {n_skipped_dh} template dihedral(s) referencing unmapped atoms')
            inst_dihedrals = inst_dihedrals[dh_mask]
            # add new dihedrals to global topology
            d = self.Topology.D['dihedrals']
            self.Topology.D['dihedrals'] = pd.concat((d, inst_dihedrals), ignore_index=True)
            # hard check for no nans
            d = self.Topology.D['dihedrals']
            check = True
            for a in ['ai','aj','ak','al']:
                check = check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in dihedrals')
                raise ValueError('NAN in dihedrals')
            d = self.Topology.D['pairs']
            check = True
            for a in ['ai','aj']:
                check = check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in pairs premapping')
                raise ValueError('NAN in pairs premapping')

            # Map pair atoms.  As with angles/dihedrals, a pair entry that
            # references a template atom without a system counterpart (an
            # H present in the fresh-CY template but absent in the post-build
            # system CY) is silently skipped — its parameters describe an
            # interaction the cured system doesn't have.
            temp_pairs.ai = temp_pairs.ai.map(temp2inst)
            temp_pairs.aj = temp_pairs.aj.map(temp2inst)
            pair_mask = temp_pairs[['ai','aj']].notna().all(axis=1)
            n_skipped_pr = int((~pair_mask).sum())
            if n_skipped_pr:
                logger.debug(f'Skipping {n_skipped_pr} template pair(s) referencing unmapped atoms')
            temp_pairs = temp_pairs[pair_mask]
            # add these pairs to the topology
            self.Topology.D['pairs'] = pd.concat((d,temp_pairs),ignore_index=True)
            d = self.Topology.D['pairs']
            # check AGAIN for nans
            check = True
            for a in ['ai','aj']:
                check = check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in pairs post mapping')
                raise ValueError('NAN in pairs post mapping')
            # return temp_pairs
        mapped_inst_atoms = list(set(mapped_inst_atoms))
        logger.debug(f'System overcharge after mapping: {self.Topology.total_charge():.4f}')
        if adjust_charges:
            self.Topology.adjust_charges(atoms=mapped_inst_atoms, overcharge_threshhold=overcharge_threshhold, msg=f'overcharge magnitude exceeds {overcharge_threshhold}')
        return mapped_inst_atoms

    def enumerate_1_4_pairs(self, at_idx):
        """Enumerates all 1-4 pair interactions resulting from new bonds in at_idx.

        Args:
            at_idx (list): list of 2-tuples, each containing global indices of atoms that bond to each other

        Returns:
            pandas.DataFrame: dataframe of all pairs
        """
        pdf = self.Topology.D['pairs']
        # each of these bonds results in 1-4 pair interactions
        bl = self.Topology.bondlist
        pai = []
        paj = []
        pri = []
        prj = []
        prsource = []
        for p in at_idx:
            # determine all nj--nk pairs from j-k bond, where nj are neighbors of j not including k, and nk are neighbors of k not including j
            # and all nnj--k pairs from all nj-j bonds, where nnj are next-nearest neighbors of j excluding j
            # and all j--nnk pairs from all k-nk bonds, where nnk are next-nearest neighbors of k excluding k
            j, k = p
            nj = bl.partners_of(j)
            nj.remove(k)
            nnj = []
            for nn in nj:
                t = bl.partners_of(nn)
                t.remove(j)
                nnj.append(t)
            nk = bl.partners_of(k)
            nk.remove(j)
            nnk = []
            for nn in nk:
                t = bl.partners_of(nn)
                t.remove(k)
                nnk.append(t)
            logger.debug(f'{j} {nj} {k} {nk}')
            logger.debug(f'{nnj} {nnk}')
            for pp in product(nj, nk):
                jj, kk = pp
                pai.append(jj)
                paj.append(kk)
                pri.append(j)
                prj.append(k)
                prsource.append('c')
                for nn in nnj:
                    for ii in nn:
                        pai.append(ii)
                        paj.append(k)
                        pri.append(j)
                        prj.append(k)
                        prsource.append('l')
                for nn in nnk:
                    for ll in nn:
                        pai.append(j)
                        paj.append(ll)
                        pri.append(j)
                        prj.append(k)
                        prsource.append('r')

        pi_df = pd.DataFrame({'ai': pai, 'aj': paj, 'source': prsource, 'bi': pri, 'bj': prj})
        pi_df.drop_duplicates(inplace=True, ignore_index=True)
        return pi_df

    def update_topology_and_coordinates(self,bdf,template_dict={},write_mapper_to=None,chain_manager=None,**kwargs):
        """Updates global topology and necessary atom attributes in the configuration to reflect formation of all bonds listed in bdf.

        Args:
            bdf (pandas.DataFrame): bonds dataframe with columns 'ai', 'aj', 'reactantName'
            template_dict (dict): dictionary of molecule templates keyed on molecule name, defaults to {}
            write_mapper_to (str): filename to write index mapper to, defaults to None
            chain_manager: optional ChainManager owned by the caller; updated in place if provided

        Returns:
            tuple: bonds dataframe and pairs dataframe with atom indices updated to reflect any atom deletions
        """
        explicit_sacH=kwargs.get('explicit_sacH',{})
        template_source=kwargs.get('template_source','internal')
        overcharge_threshhold=kwargs.get('overcharge_threshhold',0.1)
        logger.debug(f'begins.')
        if bdf.shape[0]>0:
            assert bdf['ai'].dtype==int
            assert bdf['aj'].dtype==int
            # pull out just the atom index pairs (first element of each tuple)
            at_idx=[(int(x.ai),int(x.aj),x.order) for x in bdf.itertuples()]
            logger.debug(f'Making {len(at_idx)} bonds.')
            idx_to_delete=self.make_bonds(at_idx,explicit_sacH=explicit_sacH,chain_manager=chain_manager)
            logger.debug(f'Deleting {len(idx_to_delete)} atoms.')
            idx_mapper=self.delete_atoms(idx_to_delete) # will result in full reindexing
            if chain_manager is not None:
                chain_manager.remap(idx_mapper)
            self.Topology.null_check(msg='delete_atoms')
            # reindex all atoms in the list of bonds sent in, and write it out
            logger.debug(f'z-decrement, nreactions increment')
            ri_bdf=bdf.copy()
            ri_bdf.ai=ri_bdf.ai.map(idx_mapper)
            ri_bdf.aj=ri_bdf.aj.map(idx_mapper)
            at_idx=[(x.ai,x.aj) for x in ri_bdf.itertuples()]
            for idx_pair in at_idx:
                for idx in idx_pair:
                    self.decrement_gro_attribute_by_attributes('z',{'globalIdx':idx})
                    self.increment_gro_attribute_by_attributes('nreactions',{'globalIdx':idx})
            if template_source=='internal':
                logger.debug(f'calling map_from_templates')
                self.map_from_templates(ri_bdf,template_dict,overcharge_threshhold=overcharge_threshhold,chain_manager=chain_manager)
            logger.debug(f'1-4 pair update')
            pi_df=self.enumerate_1_4_pairs(at_idx)
            self.Topology.null_check(msg='update_topology_and_coordinates')
            if write_mapper_to:
                tdf=pd.DataFrame({'old':list(idx_mapper.keys()),'new':list(idx_mapper.values())})
                tdf.to_csv(write_mapper_to,sep=' ',index=False)
            logger.debug('finished')
            self.Topology.rings.remap(idx_mapper)
            # self.bondchainlist_remap(idx_mapper)
            return ri_bdf,pi_df

    def read_top(self,topfilename):
        """Creates a new Topology member by reading from a Gromacs-style top file. Just a wrapper for the read_top method of the Topology class.

        Args:
            topfilename (str): name of topology file
        """
        self.files['top']=os.path.abspath(topfilename)
        self.Topology=Topology.read_top(topfilename)

    def read_tpx(self,tpxfilename):
        """Reads in extended topology information.

        Args:
            tpxfilename (str): name of tpx file
        """
        self.files['tpx']=os.path.abspath(tpxfilename)
        self.Topology.read_tpx(tpxfilename)

    def read_gro(self,grofilename,preserve_box=False,wrap_coords=False):
        """Creates a new Coordinates member by reading from a Gromacs-style coordinates file. Just a wrapper for the read_gro method of Coordinates.

        Args:
            grofilename (str): name of gro file
            preserve_box (bool): if True, preserve the existing box, defaults to False
            wrap_coords (bool): if True, wrap coordinates into the box, defaults to False
        """
        self.files['gro']=os.path.abspath(grofilename)
        if preserve_box:
            savebox=self.Coordinates.box.copy()
        self.Coordinates=Coordinates.read_gro(grofilename,wrap_coords=wrap_coords)
        self.Coordinates.claim_parent(self)
        if preserve_box:
            self.Coordinates.box=savebox
        # logger.debug(f'box: {self.Coordinates.box}')

    def read_mol2(self,mol2filename,ignore_bonds=False,overwrite_coordinates=False):
        """Creates a new Coordinates member by reading from a SYBYL-style MOL2 file. A wrapper for read_mol2 from Coordinates, but also sets the 'mol2_bonds' dataframe in the Topology if the parameter ignore_bonds is False. If the mol2_bonds dataframe is created, and the Topology already has a 'bonds' dataframe, a consistency check is performed.

        Args:
            mol2filename (str): name of mol2 file
            ignore_bonds (bool): if True, skip reading bonds from mol2, defaults to False
            overwrite_coordinates (bool): if True, overwrite existing coordinates, defaults to False
        """
        self.files['mol2']=os.path.abspath(mol2filename)
        temp_coords=Coordinates.read_mol2(mol2filename)
        if self.Coordinates.empty or overwrite_coordinates:
            save_box=self.Coordinates.box.copy()
            self.Coordinates=temp_coords
            self.Coordinates.box=save_box
            temp_coords=self.Coordinates
        if not ignore_bonds:
            self.Topology.D['mol2_bonds']=temp_coords.mol2_bonds.copy()
            if 'bonds' in self.Topology.D:
                self.Topology.bond_source_check()

    def swap_atom_names(self,ai,aj):
        """Swaps the names of the two atoms with global indices ai and aj. This is used when automatically selecting which of several possible sacrificial H's will actually be selected. Surviving H's are renamed so it always appears that the H with the "least important" name (lowest order if sorted) is the sacrificial H, giving perfect control of the names of the atoms that survive a bond.

        Args:
            ai (int): global index of first atom
            aj (int): global index of second atom
        """
        T=self.Topology.D['atoms']
        C=self.Coordinates.A
        l1=T.columns=='atom'
        l2=C.columns=='atomName'
        iname=T.iloc[ai-1,l1].values[0]
        jname=T.iloc[aj-1,l1].values[0]
        # logger.debug(f'Swapping names of atoms {ai}({iname}) and {aj}({jname})')
        tmpNm=T.iloc[ai-1,l1].values[0]
        T.iloc[ai-1,l1]=T.iloc[aj-1,l1]
        T.iloc[aj-1,l1]=tmpNm
        C.iloc[ai-1,l2]=C.iloc[aj-1,l2]
        C.iloc[aj-1,l2]=tmpNm

    def read_top_gro(self,topfilename,grofilename):
        """Wrapper for read_top and read_gro; generates new Topology and Coordinates members.

        Args:
            topfilename (str): name of topology file
            grofilename (str): name of coordinates file (Gromacs format)
        """
        self.read_top(topfilename)
        self.read_gro(grofilename)

    def write_top(self,topfilename):
        """Writes a Gromacs-format topology file; this will only write an in-line version, no itp; wrapper for Topology.write_top().

        Args:
            topfilename (str): name of file to write
        """
        self.Topology.write_top(topfilename)
        self.files['top']=os.path.abspath(topfilename)

    def write_tpx(self,tpxfilename):
        self.Topology.write_tpx(tpxfilename)
        self.files['tpx']=os.path.abspath(tpxfilename)

    def write_gro(self,grofilename,grotitle=''):
        """Writes a Gromacs-format coordinate file; wrapper for Coordinates.write_gro().

        Args:
            grofilename (str): name of file to write
            grotitle (str): title to put in gro file, defaults to ''
        """
        self.Coordinates.write_gro(grofilename,grotitle=grotitle)
        self.files['gro']=os.path.abspath(grofilename)

    # def write_top_gro(self,topfilename,grofilename):
    #     """write_top_gro Writes both a Gromacs top file and Gromacs coordinate file

    #     :param topfilename: name of topology file to write
    #     :type topfilename: str
    #     :param grofilename: name of coordinate file to write
    #     :type grofilename: str
    #     """
    #     self.write_top(topfilename)
    #     self.write_gro(grofilename)

    def return_bond_lengths(self,bdf):
        """Returns the length of all bonds in bdf.

        Args:
            bdf (pandas.DataFrame): bonds dataframe with columns 'ai', 'aj', 'reactantName'

        Returns:
            list: list of lengths parallel to bonds
        """
        return self.Coordinates.return_bond_lengths(bdf)

    def add_length_attribute(self,bdf:pd.DataFrame,attr_name='length'):
        """Computes bond lengths based on bonds indicated by the parallel 'ai' and 'aj' columns of the parameter dataframe bdf and stores result in a new column called attr_name.

        Args:
            bdf (pd.DataFrame): a pandas dataframe with 'ai' and 'aj' columns of atom indices indicating bonds
            attr_name (str): name of length attribute column, defaults to 'length'
        """
        bdf[attr_name]=self.Coordinates.return_bond_lengths(bdf)

    def copy_coords(self,other):
        """Copies coordinates and box size from other to self.

        Args:
            other (TopoCoord): a TopoCoord instance
        """
        self.Coordinates.copy_coords(other.Coordinates)
        self.Coordinates.box=other.Coordinates.box.copy()
        self.files['gro']=os.path.abspath(other.files['gro'])

    def set_grx_attributes(self,attributes=[]):
        """Overrides the global GRX_ATTRIBUTES.

        Args:
            attributes (list): new GRX attributes to use, defaults to []
        """
        if len(attributes)==0:
            self.grxattr=GRX_ATTRIBUTES
        else:
            self.grxattr=attributes
        logger.debug(f'grxattr set to {self.grxattr}')
    
    def write_gro_attributes(self,attributes_list,grxfilename):
        """Writes atomic attributes to a file.

        Args:
            attributes_list (list): list of attributes to write
            grxfilename (str): name of output file
        """
        self.Coordinates.write_atomset_attributes(attributes_list,grxfilename)
        self.files['grx']=os.path.abspath(grxfilename)

    def write_grx_attributes(self,grxfilename):
        """Writes GRX attributes to a file.

        Args:
            grxfilename (str): name of output file
        """
        self.write_gro_attributes(self.grxattr,grxfilename)

    def read_gro_attributes(self,grxfilename,attribute_list=[],chain_manager=None):
        """Reads attributes from file into self.Coordinates.A.

        Args:
            grxfilename (str): name of input file
            attribute_list (list): list of attributes to take, defaults to [] (take all)
            chain_manager: optional ChainManager owned by the caller; reconstructed from the data if provided
        """
        self.files['grx']=os.path.abspath(grxfilename)
        logger.debug(f'Reading {grxfilename}')
        attributes_read=self.Coordinates.read_atomset_attributes(grxfilename,attributes=attribute_list)
        if chain_manager is not None:
            chain_manager.from_dataframe(self.Coordinates.A)
        if attributes_read!=self.grxattr:
            self.grxattr=attributes_read

    def set_gro_attribute(self,attribute,srs):
        """Sets attribute of atoms to srs (drills through to Coordinates.set_atomset_attributes()).

        Args:
            attribute (str): name of attribute
            srs: scalar or list-like attribute values in same ordering as self.A
        """
        self.Coordinates.set_atomset_attribute(attribute,srs)

    def set_gro_attribute_by_attributes(self,att_name,att_value,attribute_dict):
        """Sets the attribute named att_name to att_value for the set of atoms specified by attribute_dict.

        Args:
            att_name (str): name of attribute to set
            att_value: value to set the attribute to
            attribute_dict (dict): dictionary of attribute:value pairs that specify the set of atoms to be considered
        """
        self.Coordinates.set_atom_attribute(att_name,att_value,attribute_dict)

    def get_gro_attribute_by_attributes(self,att_name,attribute_dict):
        """Returns values of attributes listed in att_name from atoms specified by attribute:value pairs in attribute_dict.

        Args:
            att_name: name or list of names of attributes whose values are to be returned
            attribute_dict (dict): dictionary of attribute:value pairs that specify the set of atoms to be considered

        Returns:
            scalar or list: one or more return attribute values (list if att_name is a list; scalar otherwise)
        """
        return self.Coordinates.get_atom_attribute(att_name,attribute_dict)

    def increment_gro_attribute_by_attributes(self,att_name,attribute_dict):
        """Adds one to attribute att_name of all atoms identified by attribute:value pairs in attribute_dict.

        Args:
            att_name (str): name of attribute to increment
            attribute_dict (dict): attribute:value pairs that specify atoms to which to apply this incrementation
        """
        val=self.get_gro_attribute_by_attributes(att_name,attribute_dict)
        val+=1
        self.set_gro_attribute_by_attributes(att_name,val,attribute_dict)

    def decrement_gro_attribute_by_attributes(self,att_name,attribute_dict):
        """Subtracts one from attribute att_name of all atoms identified by attribute:value pairs in attribute_dict.

        Args:
            att_name (str): name of attribute to decrement
            attribute_dict (dict): attribute:value pairs that specify atoms to which to apply this decrementation
        """
        val=self.get_gro_attribute_by_attributes(att_name,attribute_dict)
        val-=1
        self.set_gro_attribute_by_attributes(att_name,val,attribute_dict)

    def get_gro_attributelist_by_attributes(self,attribute_list,attribute_dict):
        """Returns all rows of atoms dataframe and columns named in attribute_list for atoms identified by the attribute_dict.

        Args:
            attribute_list (list): list of attributes to be in the rows that are returned
            attribute_dict (dict): dictionary of attribute:value pairs that specify the set of atoms to be considered

        Returns:
            pd.DataFrame: a dataframe segment
        """
        return self.Coordinates.get_atoms_w_attribute(attribute_list,attribute_dict)

    def get_R(self,idx):
        """Returns the cartesian position of atom with global index idx.

        Args:
            idx (int): atom global index

        Returns:
            numpy.ndarray: position of atom as array of shape (3,)
        """
        return self.Coordinates.get_R(idx)

    def rotate(self,R):
        """Applies rotation matrix R to all atom positions.

        Args:
            R (numpy.ndarray): rotation matrix of shape (3, 3)
        """
        self.Coordinates.rotate(R)

    def translate(self,L):
        """Applies translation vector L to all atom positions.

        Args:
            L (numpy.ndarray): translation vector of shape (3,)
        """
        self.Coordinates.translate(L)

    def resid_partners_of(self,ri):
        """Returns list of resid partners of resid ri.

        Args:
            ri (int): residue index

        Returns:
            list: list of partner residue indices (two residues are partners if there is at least one interatomic bond joining them)
        """
        result=[]
        adf=self.Coordinates.A
        radf=adf[adf['resNum']==ri]
        for at in radf['globalIdx'].to_list():
            bl=self.Topology.bondlist.partners_of(at)
            for j in bl:
                theirresid=adf.iloc[j-1]['resNum']
                if theirresid!=ri and not theirresid in result:
                    result.append(theirresid)
        return result

    def interresidue_partners_of(self,i):
        """Returns list of atom indices that are bonded partners of atom i and are not in the residue of atom i.

        Args:
            i (int): atom global index

        Returns:
            list: list of atom indices of atoms that are partners of i not in i's residue
        """
        result=[]
        bl=self.Topology.bondlist.partners_of(i)
        # logger.debug(f'{i} partners {bl}')
        myresid=self.Coordinates.A.iloc[i-1]['resNum']
        for j in bl:
            theirresid=self.Coordinates.A.iloc[j-1]['resNum']
            if theirresid!=myresid:
                result.append(j)
        return result

    def get_siblings(self, atom_idx):
        """Identifies already-reacted atoms two bonds from each atom of a bond, inside its own residue.

        A template splice overwrites the charges of every template atom that
        shares an angle, dihedral or 1-4 pair with the new bond, which reaches
        two bonds into each residue.  If one of those atoms has already bonded
        to another residue, a template that shows it unreacted gives it the
        wrong charge: in a triazine, the second and third ring carbons to react
        each reset the carbons that reacted before them.  Recording these
        atoms lets find_template prefer a template that shows them reacted.

        Only the second shell counts.  A reacted atom directly bonded to the
        bonding atom is the C=C chain case, which the one-away context already
        describes; counting it too would change which template every chain
        polymerization uses.

        Args:
            atom_idx: container of two global atom indices specifying a bond

        Returns:
            list: two sorted lists of (atomName, partner resName) tuples, one per atom in
                atom_idx, excluding the bond itself
        """
        A = self.Coordinates.A
        bl = self.Topology.bondlist
        result = [[], []]
        for x in (0, 1):
            a, b = atom_idx[x], atom_idx[1 - x]
            resnum = A.iloc[a - 1]['resNum']
            seen = {a}
            frontier = [a]
            for _ in range(2):
                reached = []
                for i in frontier:
                    for j in bl.partners_of(i):
                        if j not in seen and A.iloc[j - 1]['resNum'] == resnum:
                            seen.add(j)
                            reached.append(j)
                frontier = reached
            for sib in frontier:
                for j in bl.partners_of(sib):
                    if j != b and A.iloc[j - 1]['resNum'] != resnum:
                        result[x].append((A.iloc[sib - 1]['atomName'], A.iloc[j - 1]['resName']))
            result[x].sort()
        return result

    def minimum_distance(self, other, self_excludes=None, other_excludes=None):
        """Computes the distance of closest approach between self's Coordinates and other's Coordinates.

        Args:
            other (TopoCoord): another TopoCoord object
            self_excludes (list, optional): atom globalIdx values to ignore in self, defaults to None
            other_excludes (list, optional): atom globalIdx values to ignore in other, defaults to None

        Returns:
            float: distance of closest approach (nm)
        """
        return self.Coordinates.minimum_distance(other.Coordinates, self_excludes=self_excludes, other_excludes=other_excludes)

    def gro_DataFrame(self,name):
        """Returns the appropriate Coordinates dataframe based on the directive in 'name'.

        Args:
            name (str): either 'atoms' or 'mol2_bonds'

        Returns:
            pandas.DataFrame: the desired dataframe
        """
        if name=='atoms':
            return self.Coordinates.A
        elif name=='mol2_bonds':
            return self.Coordinates.mol2_bonds
        else:
            return None

    def overwrite_coords(self,other):
        """Overwrites coordinates in self by those in other.

        Args:
            other (TopoCoord): another TopoCoord object
        """
        # logger.debug(f'Overwriting {other.Coordinates.A.shape[0]} coordinates')
        C=self.Coordinates.A
        C=C.set_index('globalIdx')
        # logger.debug(f'before update:\n{C.to_string()}')
        B=other.Coordinates.A
        B=B.set_index('globalIdx')
        B=B[['posX','posY','posZ']].copy()
        # logger.debug(f'new coordinates:\n{B.to_string()}')
        C.update(B)
        self.Coordinates.A=C.reset_index()
        # logger.debug(f'after update:\n{self.Coordinates.A.to_string()}')

    def linkcell_initialize(self,cutoff,force_repopulate=True):
        """Initializes the linkcell structure; a wrapper for Coordinates.

        Args:
            cutoff (float): cell side length (nm)
            force_repopulate (bool): if True, ignore any cached linkcell_idx file, defaults to True
        """
        self.Coordinates.linkcell_initialize(cutoff,populate=True,force_repopulate=force_repopulate)

    def linkcell_cleanup(self):
        """Removes linkcell_idx attribute from Coordinate.A.
        """
        self.Coordinates.A.drop(columns=['linkcell_idx'],inplace=True)

    def atom_count(self):
        """Checks to be sure the Coordinate and Topology members contain the same number of atoms.

        Returns:
            int: the number of atoms
        """
        assert self.Coordinates.A.shape[0]==self.Topology.D['atoms'].shape[0]
        return self.Coordinates.A.shape[0]

    def total_volume(self,units='SI'):
        """Returns total volume represented by the system's Coordinates.

        Args:
            units (str): unit designation, defaults to 'SI'

        Returns:
            float: box volume
        """
        return self.Coordinates.box_volume(units=units)

    def density(self,units='SI'):
        """Returns system density.

        Args:
            units (str): unit designation, defaults to 'SI'

        Returns:
            float: density
        """
        return self.Topology.total_mass(units)/self.total_volume(units)

    def wrap_coords(self):
        """Wraps all coordinates into center periodic box.
        """
        self.Coordinates.wrap_coords()

    def inherit_grx_attributes_from_molecules(self,molecule_dict,initial_composition,globally_unique=GRX_GLOBALLY_UNIQUE,unset_defaults=GRX_UNSET_DEFAULTS,overall_default=0):
        """Copies non-Gromacs-standard atom attributes from molecule templates in molecule_dict according to molecule counts in dict initial_composition.

        Args:
            molecule_dict (dict): dictionary of available molecules (name:Molecule)
            initial_composition (dict): dictionary of initial composition (name:count)
            globally_unique (list): boolean list indicating attributes that must be globally unique, defaults to GRX_GLOBALLY_UNIQUE
            unset_defaults (list): list of values, one per attribute, that signify UNSET, defaults to GRX_UNSET_DEFAULTS
            overall_default (int): default UNSET value for all attributes if unset_defaults is empty, defaults to 0
        """
        # logger.debug(f'inherit grx {attributes} {unset_defaults} {globally_unique}')
        ''' set up the globally_unique and unset_defaults list if necessary '''
        if len(globally_unique)!=len(self.grxattr):
            globally_unique=[False for _ in range(len(self.grxattr))]
        if len(unset_defaults)!=len(self.grxattr):
            unset_defaults=[overall_default for _ in range(len(self.grxattr))]

        ''' drop all attribute values from current global atom dataframe '''
        adf=self.Coordinates.A
        self.Coordinates.A=adf.drop(columns=[a for a in self.grxattr if a in adf])

        logger.debug(f'{adf.shape[0]} atoms inheriting these attributes from molecular templates:')
        logger.debug(f'    Attribute name     Default value   Globally unique?')
        for attname,defval,gu in zip(self.grxattr,unset_defaults,globally_unique):
            logger.debug(f'    {attname:<15s}    {str(defval):<13s}   {gu}')
        if len(unset_defaults)==len(self.grxattr):
            att_running_maxval={}
            for k,v in zip(self.grxattr,unset_defaults):
                if type(v)==int or type(v)==float:
                    att_running_maxval[k]=0
                else:
                    att_running_maxval[k]='0'
        else:
            att_running_maxval={k:overall_default for k in self.grxattr}

        attribute_lists={k:[] for k in self.grxattr}
        value_counts={k:0 for k in self.grxattr}
        mol_idx_counter=0
        for icdict in [cc for cc in initial_composition if 'count' in cc]:
            molecule=icdict['molecule']
            count=icdict['count']
            mol_adf=molecule_dict[molecule].TopoCoord.Coordinates.A
            logger.debug(f'Inheriting from {molecule}')
            for ln in mol_adf.head().to_string().split('\n'):
                 logger.debug(ln)
            mol_attr_df=mol_adf[self.grxattr]
            for molecule_number in range(count):
                for i,k in enumerate(self.grxattr):
                    tra=mol_attr_df[k].to_list()
                    if k=='molecule': 
                        utra=[mol_idx_counter for _ in tra]
                        mol_idx_counter+=1
                    else:
                    # nv=len(tra)-tra.count(unset_defaults[i])
                        nuv=len(list(set([x for x in tra if x != unset_defaults[i]])))
                        utra=[]
                        for x in tra:
                            if globally_unique[i] and (type(x)==int or type(x)==float):
                                xx=x+value_counts[k] if x!=unset_defaults[i] else unset_defaults[i]
                            else:
                                xx=x
                            utra.append(xx)
                        value_counts[k]+=nuv
                    attribute_lists[k].extend(utra)

        for k,L in attribute_lists.items():
            self.Coordinates.A[k]=L
        logger.debug(f'postinherit adf columns {self.Coordinates.A.columns}')

    def maxspan(self):
        """Returns the maxspan of the Coordinates (dimensions of orthorhombic convex hull enclosing Coordinates). Just a wrapper.

        Returns:
            numpy.ndarray: array of x-span, y-span, z-span
        """
        return self.Coordinates.maxspan()

    def minmax(self):
        """Returns the coordinates of the atoms at the lower-leftmost and upper-rightmost positions in the constellation of points in the atoms dataframe.

        Returns:
            tuple: two numpy arrays, lower-leftmost and upper-rightmost positions, respectively
        """
        return self.Coordinates.minmax()

    def checkbox(self):
        """Checks that the entire constellation of points in the atoms dataframe fits within the designated box for this Configuration object.

        Returns:
            tuple: (bool, bool), True for each of lower-leftmost and upper-rightmost points if they are within the box
        """
        return self.Coordinates.checkbox()

    def write_mol2(self,filename,molname='',element_names_as_types=False):
        """Writes a SYBYL MOL2-format file using Coordinates, with certain atom attributes borrowed from the Topology.

        Args:
            filename (str): name of file to write
            molname (str): name of molecule to put in mol2 file, defaults to ''
            element_names_as_types (bool): if True, use element names as atom types, defaults to False
        """
        if molname=='':
            molname='This Molecule has no name'
        other_attributes=pd.DataFrame()
        if element_names_as_types:
            element_names=[x[0] for x in self.Topology.D['atoms']['atom']]
            other_attributes['type']=element_names
        else:
            other_attributes['type']=self.Topology.D['atoms']['type']
        other_attributes['charge']=self.Topology.D['atoms']['charge']
        self.files['mol2']=os.path.abspath(filename)
        # logger.debug(f'write_mol2, other_attributes:\n{other_attributes.to_string()}')
        if 'mol2_bonds' in self.Topology.D:
            self.Coordinates.write_mol2(filename,molname=molname,bondsDF=self.Topology.D['mol2_bonds'],other_attributes=other_attributes)
        elif 'bonds' in self.Topology.D and not self.Topology.D['bonds'].empty:
            # A molecule built by copying a parent's topology -- a symmetry
            # sibling such as GMAS-4 -- is loaded from top+gro and never gets a
            # mol2 bond table, so this used to write a bondless MOL2 and warn.
            # Its topology does carry every bond; write those instead.  Bond
            # orders are not in a GROMACS topology, so they are written as single,
            # the same assumption Topology.add_bonds makes for a new bond.
            b=self.Topology.D['bonds'][['ai','aj']].reset_index(drop=True)
            topo_bonds=pd.DataFrame({'bondIdx':range(1,len(b)+1),'ai':b['ai'],'aj':b['aj'],'order':'1'})
            logger.debug(f'{filename}: no mol2 bond table; writing {len(b)} bonds from the topology as single bonds')
            self.Coordinates.write_mol2(filename,molname=molname,bondsDF=topo_bonds,other_attributes=other_attributes)
        else:
            self.Coordinates.write_mol2(filename,molname=molname,other_attributes=other_attributes)

    def merge(self,other,self_cm=None,other_cm=None):
        """Merges the TopoCoord instance "other" to self.

        Args:
            other (TopoCoord): another TopoCoord instance
            self_cm: optional ChainManager for self, updated in place if provided
            other_cm: optional ChainManager for other, shifted in place if provided

        Returns:
            tuple: a shift tuple (returned by Coordinates.merge())
        """
        self.Topology.merge(other.Topology)
        shifts=self.Coordinates.merge(other.Coordinates)
        if self_cm is not None and other_cm is not None:
            other_cm.shift(shifts[0]) # updates atom idx only
            self_cm.injest_other(other_cm)
            self_cm.to_dataframe(self.Coordinates.A)
        # for name,idx_lists in other.idx_lists.items():
        #     # logger.debug(f'TopoCoord merge: list_name {name} lists {idx_lists}')
        #     for a_list in idx_lists:
        #         self.idx_lists[name].append([x+shifts[0] for x in a_list])
        #     self.reset_grx_attributes_from_idx_list(name)
        return shifts

    def bondtest_df(self,df:pd.DataFrame,pbc=[1,1,1],show_piercings=True):
        """Applies bond filters to all bonds in the dataframe.

        Args:
            df (pd.DataFrame): dataframe of possible bonds with columns 'ai', 'aj' and 'r'; this method adds the column 'results'
            pbc (list): flags indicating dimensions in which pbc are used, defaults to [1,1,1]
            show_piercings (bool): toggles diagnostic output for pierced rings, defaults to True

        Returns:
            pandas.DataFrame: input data frame with new 'results' column
        """
        if df.empty:
            return df
        self.Topology.rings.injest_coordinates(self.Coordinates.A)
        results=[]
        for r in df.itertuples():
            result,dummy=self.bondtest((r.ai,r.aj,r.r),pbc=pbc,show_piercings=show_piercings)
            results.append(result)
        df['result']=results
        for ln in df[df['result']==BTRC.passed].to_string().split('\n'):
            logger.debug(ln)
        return df

    def bondtest(self,b,pbc=[1,1,1],show_piercings=True):
        """Determines if bond b is to be allowed to form based on geometric and topological criteria.

        Args:
            b (tuple): bond as tuple (ai, aj, rij)
            pbc (list): periodic boundary condition flags in each direction, defaults to [1,1,1]
            show_piercings (bool): flag indicating you want to write gro files showing pierced rings, defaults to True

        Returns:
            BTRC: BTRC enum instance
        """
        i,j,rij=b
        # i=int(i)
        # j=int(j)
        if self.makes_shortcircuit(i,j):
            return BTRC.failed_shortcircuit,0
        # if self.makes_bondcycle(i,j):
        #     return BTRC.failed_bondcycle,0
        if self.pierces_ring(i,j,pbc=pbc,show_piercings=show_piercings):
            return BTRC.failed_pierced_ring,0
        # logger.debug(f'passed: {i:>7d} {j:>7d} {rij:>6.3f} nm')
        return BTRC.passed,rij

    def pierces_ring(self,i,j,pbc=[1,1,1],show_piercings=True):
        """Checks to see if bond i-j would pierce any covalent ring structure.

        Args:
            i (int): global index of an atom
            j (int): global index of another atom
            pbc (list): flags indicating which dimensions have pbc applied, defaults to [1,1,1]
            show_piercings (bool): toggles diagnostic output of ring piercings, defaults to True

        Returns:
            bool: True if a ring is pierced; False otherwise
        """
        adf=self.Coordinates.A
        LC=self.Coordinates.linkcell
        assert 'linkcell_idx' in adf.columns,f'Error: atoms have no linkcell_idx attribute - bug!'
        i_lcidx=int(self.get_gro_attribute_by_attributes('linkcell_idx',{'globalIdx':i}))
        j_lcidx=int(self.get_gro_attribute_by_attributes('linkcell_idx',{'globalIdx':j}))
        nearby_cells=LC.neighbor_cell_set(i_lcidx) | LC.neighbor_cell_set(j_lcidx)
        nearby_atom_ids=LC.nearby_atom_ids(adf,nearby_cells)
        nearby_rings=self.Topology.rings.filter(nearby_atom_ids)
        logger.debug(f'Ring-pierce check for bond {i}-{j} will consider {len(nearby_rings)} rings')
        B=adf.iloc[[i-1,j-1]].copy()
        seg=np.array(B[['posX','posY','posZ']].values)
        # saveB=B.copy()
        # logger.debug(f'bond {i} {j} subframe')
        # for ln in B.to_string().split('\n'):
        #     logger.debug(ln)
        for r in nearby_rings:
            ru=r.unwrap(seg[0],pbc=pbc,unwrapf=self.Coordinates.unwrap)
            seg[1]=self.Coordinates.unwrap(seg[1],seg[0],pbc=pbc)
            is_pierced,pierce_point=ru.pierced_by(seg)
            if is_pierced:
                if show_piercings:
                    C=adf[adf['globalIdx'].isin(r.idx)].copy()
                    J=pd.concat([B,C])
                    sub=self.Coordinates.subcoords(J)
                    sub.write_gro(f'ring-{str(r)}=bond-{i}-{j}'+'.gro')
                #     sub=self.Coordinates.subcoords(pd.concat([saveB,saveC]))
                #     sub.write_gro(f'ring-orig-{i}-{j}'+'.gro')
                #     logger.debug(f'Cycle {str} pierced by bond candidate {i}-{j}')
                #     for ln in sub.A.to_string().split('\n'):
                #         logger.debug(ln)
                return True
        return False

    def makes_shortcircuit(self,i,j):
        """Determines whether atoms i and j, if bonded, would produce a short circuit, defined as an instance in which i and j belong to residues that are already bonded to each other.

        Args:
            i (int): global index of first atom
            j (int): global index of second atom

        Returns:
            bool: True if a short circuit would happen, False otherwise
        """
        i_resName,i_resNum,i_atomName,i_molNum=self.get_gro_attribute_by_attributes(['resName','resNum','atomName','molecule'],{'globalIdx':i})
        j_resName,j_resNum,j_atomName,j_molNum=self.get_gro_attribute_by_attributes(['resName','resNum','atomName','molecule'],{'globalIdx':j})
        '''
        In a cure reaction, atoms that react should be in different residues
        '''
        assert i_resNum!=j_resNum,f'shortcircuit test error {i}-{j} both in residue {i_resNum}?'
        i_neighbors=self.Topology.bondlist.partners_of(i)
        j_neighbors=self.Topology.bondlist.partners_of(j)
        # if i in j_neighbors or j in i_neighbors:
        #     # logger.debug(f'atoms {i} and {j} already on each other\'s list of bonded partners')
        #     return True
        '''
        Should never test a proposed bond that already exists
        '''
        assert not j in i_neighbors
        assert not i in j_neighbors
        '''
        Set up a DataFrame for heavy atoms in each molecule
        '''
        ADF=self.Coordinates.A
        R1DF=ADF[ADF['molecule']==i_molNum]
        R1DF=R1DF[[(not (b.startswith('H') or b.startswith('h'))) for b in R1DF['atomName']]]
        R2DF=ADF[ADF['molecule']==j_molNum]
        R2DF=R2DF[[(not (b.startswith('H') or b.startswith('h'))) for b in R2DF['atomName']]]
        assert R1DF.shape[0]>0
        assert R2DF.shape[0]>0
        '''
        Test to see if there exists any bond between these two molecules
        '''
        for i,j in product(R1DF['globalIdx'].to_list(),R2DF['globalIdx'].to_list()):
            if self.Topology.bondlist.are_bonded(i,j):
                return True
        return False

    def reset_grx_attributes_from_idx_list(self,list_name):
        """Uses information in the "index lists" to repopulate appropriate GRX attributes. Currently, there is only one index list, for bondchains. Each index list is a list of lists; each element corresponds to a unique structure (bondchain) and is a list of global atom indices for atoms that make up that structural instance.

        Args:
            list_name (str): only allowed value currently is 'bondchain'
        """
        self.set_gro_attribute(list_name,-1)
        self.set_gro_attribute(f'{list_name}_idx',-1)
        for i,c in enumerate(self.idx_lists[list_name]):
            for j,x in enumerate(c):
                self.set_gro_attribute_by_attributes(list_name,i,{'globalIdx':x})
                self.set_gro_attribute_by_attributes(f'{list_name}_idx',j,{'globalIdx':x})

    def reset_idx_list_from_grx_attributes(self,list_name):
        """Is the inverse of reset_grx_attributes_from_idx_list: it uses the GRX attributes to rebuild index lists.

        Args:
            list_name (str): only allowed value currently is 'bondchain'
        """
        adf=self.Coordinates.A
        logger.debug(f'reset: columns {adf.columns}')
        tmp_dict={}
        for i,r in adf.iterrows():
            gix=r['globalIdx']
            cid=r[list_name]
            cix=r[f'{list_name}_idx']
            if cid!=-1:
                if not cid in tmp_dict:
                    tmp_dict[cid]={}
                tmp_dict[cid][cix]=gix
        # logger.debug(f'{list_name} tmp_dict item count: {len(tmp_dict)}')
        # logger.debug(f'{tmp_dict}')
        if tmp_dict:
            consec_test=[a in tmp_dict for a in range(len(tmp_dict))]
            # logger.debug(f'{consec_test}')
            assert all(consec_test),f'{list_name} reset_idx_list for group attribute {list_name} has non-consecutive integer keys -- bug\n{[a in tmp_dict for a in range(len(tmp_dict))]}'
            ngroups=len(tmp_dict)
            self.idx_lists[list_name]=[[] for _ in range(ngroups)]
            for i in range(ngroups):
                for j in range(len(tmp_dict[i])):
                    self.idx_lists[list_name][i].append(tmp_dict[i][j])
        # logger.debug(f'-> idx_lists[{list_name}]: {self.idx_lists[list_name]}')

    # def bondchainlist_remap(self,idx_mapper):
    #     for c in self.idx_lists['bondchain']:
    #         newlist=[idx_mapper[x] for x in c]
            
    #         for i in range(len(c)):
    #             c[i]=idx_mapper.get(c[i],c[i])

    # def bondchainlist_update(self,new_bond_recs,msg=''):
    #     """bondchainlist_update updates the bondchain index lists due to generation of new bonds

    #     :param new_bond_recs: list of bond records, each a tuple of two ints corresponding to atom indices
    #     :type new_bond_recs: list
    #     :param msg: a nice message (unused), defaults to ''
    #     :type msg: str, optional
    #     """
    #     chainlists=self.idx_lists['bondchain']
    #     if len(chainlists)==0: return
    #     # logger.debug(f'pre {msg} chainlists')
    #     # for i,c in enumerate(chainlists):
    #     #     logger.debug(f'  {i} {c}')
    #     for b in new_bond_recs:
    #         aidx,bidx=b[0],b[1]
    #         ar=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':aidx})
    #         br=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':bidx})
    #         if ar==br: continue # ignore intramolecular bonds
    #         # logger.debug(f'bondchainlist_update pair {aidx} {bidx}')
    #         ac=self.get_gro_attribute_by_attributes('bondchain',{'globalIdx':aidx})
    #         bc=self.get_gro_attribute_by_attributes('bondchain',{'globalIdx':bidx})
    #         logger.debug(f'ac {ac} bc {bc}')
    #         if ac==-1 or bc==-1:
    #             # neither of these newly bonded atoms is already in a bondchain, so
    #             # there is no possibility that this new bond can join two chains.
    #             continue
    #         # logger.debug(f'bondchain of bidx {bidx}: {chainlists[bc]}')
    #         # logger.debug(f'bondchain of aidx {aidx}: {chainlists[ac]}')
    #         aci=self.get_gro_attribute_by_attributes('bondchain_idx',{'globalIdx':aidx})
    #         bci=self.get_gro_attribute_by_attributes('bondchain_idx',{'globalIdx':bidx})
    #         logger.debug(f' -> {aidx}-{bidx}: ac {ac} #{len(chainlists[ac])} aci {aci} bc {bc} #{len(chainlists[bc])} bci {bci}')
    #         # one must be a head and the other a tail
    #         if aci==0: # a is a head
    #             if len(chainlists[bc])-1!=bci:
    #                 logger.warning(f'Atom {bidx} has index {bci} in bondchain {bc} but that bondchain has {len(chainlists[bc])} elements, so it cannot be the tail of the bondchain. The tail appears to be atom {chainlists[bc][-1]}. So something is wrong, and I am not merging these chains.')
    #                 logger.debug(f'bondchain of {aidx} {chainlists[ac]}')
    #                 logger.debug(f'bondchain of {bidx} {chainlists[bc]}')
    #                 continue
    #             c1=ac
    #             c2=bc
    #         elif bci==0: # b is a head
    #             if len(chainlists[ac])-1!=aci:
    #                 logger.warning(f'Atom {aidx} has index {aci} in bondchain {ac} but that bondchain has {len(chainlists[ac])} elements, so it cannot be the tail of the bondchain. The tail appears to be atom {chainlists[ac][-1]}. So something is wrong, and I am not merging these chains.')
    #                 logger.debug(f'bondchain of {aidx} {chainlists[ac]}')
    #                 logger.debug(f'bondchain of {bidx} {chainlists[bc]}')
    #                 continue
    #             c1=bc
    #             c2=ac
    #         else:
    #             logger.debug(f'Neither of {aidx} or {bidx} are heads of their respective chains.')
    #             logger.debug(f'This is most likely happening because of branching, which is permitted but')
    #             logger.debug(f'does not allow for merging of 1-D chains.')
    #             continue
    #         chainlists[c2].extend(chainlists[c1])
    #         for aidx in chainlists[c1]:
    #             self.set_gro_attribute_by_attributes('bondchain',c2,{'globalIdx':aidx})
    #             self.set_gro_attribute_by_attributes('bondchain_idx',chainlists[c2].index(aidx),{'globalIdx':aidx})
    #         # logger.debug(f'removing bondchain {c1}')
    #         chainlists.remove(chainlists[c1])
    #         # since we remove c1, all indices greater than c1 must decrement
    #         dec_us=np.array(self.Coordinates.A['bondchain'])
    #         bad_chain_idx=np.where(dec_us>c1)
    #         # logger.debug(f'bad_chain_idx: {bad_chain_idx}')
    #         # logger.debug(f'{dec_us[bad_chain_idx]}')
    #         dec_us[bad_chain_idx]-=1
    #         self.Coordinates.A['bondchain']=dec_us
    #     # cnms=[]
    #     # for c in self.idx_lists['bondchain']:
    #     #     cnms.append([self.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in c])
    #     # logger.debug(f'post {msg} chains {self.idx_lists["bondchain"]} {cnms}')

    # def makes_bondcycle(self,aidx,bidx):
    #     """makes_bondcycle checks the current 'bondchain' index lists to see if a bond between aidx and bidx (global atom indices) would generate a C-C' bondcycle, which they do if they have the same 'bondchain' attribute, shorter than the minimum allowable bondcycle length, if this is not 0

    #     :param aidx: global index of an atom
    #     :type aidx: int
    #     :param bidx: global index of another atom
    #     :type bidx: int
    #     :return: True if a bond betwen aidx and bidx would form a C-C' bondcycle
    #     :rtype: bool
    #     """
    #     # is there a bondchain with aidx as head and bidx as tail, or vice versa?
    #     for c in self.idx_lists['bondchain']:
    #         if (aidx==c[0] and bidx==c[-1]) or (aidx==c[-1] and bidx==c[0]):
    #             if self.min_bondcycle_length<=0 or len(c)<self.min_bondcycle_length:
    #                 return True
    #     return False

    def bondcycle_collective(self,bdf:pd.DataFrame,chain_manager=None):
        """Checks to see if, when considered as a collective, this set of bond records leads to one or more cyclic bondchains; if so, longest bonds that break bondcycles are removed from the bond records list and the resulting list is returned.

        Args:
            bdf (pandas.DataFrame): bonds data frame, sorted in ascending order by bondlength
            chain_manager: optional ChainManager owned by the caller; used read-only (deep-copied for simulation)

        Returns:
            pandas.DataFrame: bond records that results in no new bondcycles
        """
        # class working_bondlist:
        #     def __init__(self,my_idx,idx_list):
        #         self.idx=my_idx
        #         self.atidx_list=idx_list
        #         self.added_bonds=[]
        #         self.is_new_cycle=False
        #     def idx_of(self,ai):
        #         return self.atidx_list.index(ai)
        #     def is_head(self,ai):
        #         return self.idx_of(ai)==0
        #     def len(self):
        #         return len(self.atidx_list)
        #     def is_tail(self,ai):
        #         return self.idx_of(ai)==self.len()-1
        #     def addbond(self,bidx):
        #         self.added_bonds.append(bidx)
        #     def cycletag(self):
        #         self.is_new_cycle=True
        #     def merge(self,other):
        #         self.atidx_list.extend(other.atidx_list)
        #         other.atidx_list=[]
        #         self.added_bonds.extend(other.added_bonds)
        #         other.added_bonds=[]
        #         other.absorber=self

        # def addlink(bondchains,b_idx,ai,aj):
        #     ibl=None
        #     jbl=None
        #     for bl in bondchains:
        #         if ai in bl.atidx_list:
        #             ibl=bl
        #         if aj in bl.atidx_list:
        #             jbl=bl
        #     if not ibl or not jbl:
        #         # evidently this is a system that cannot form bondchains
        #         return
        #     logger.debug(f'i-atom {ai} is in bondchain {ibl.idx} (length {ibl.len()}) at index {ibl.idx_of(ai)}')
        #     logger.debug(f'j-atom {aj} is in bondchain {jbl.idx} (length {jbl.len()}) at index {jbl.idx_of(aj)}')
        #     if ibl==jbl:
        #         # same bondchain
        #         assert (ibl.is_head(ai) and ibl.is_tail(aj)) or (ibl.is_tail(ai) and ibl.is_head(aj))
        #         ibl.cycletag()
        #         ibl.addbond(b_idx)
        #         logger.debug(f'bondchain {ibl.idx} is cyclized!')
        #         # return dict(makescycle=True,cyclizedchain=i_cidx["bondchainidx"])
        #     # otherwise, we unify the two bondchains
        #     else:
        #         if ibl.is_head(ai):
        #             assert jbl.is_tail(aj)
        #             jbl.addbond(b_idx)
        #             jbl.merge(ibl)
        #         else:
        #             assert ibl.is_tail(ai) and jbl.is_head(aj)
        #             ibl.addbond(b_idx)
        #             ibl.merge(jbl)

        logger.debug(f'Checking set of {bdf.shape[0]} bonds for cyclic C-C bondchains')
        new_bdf=bdf.copy()
        new_bdf['remove-to-uncyclize']=[False for _ in range(new_bdf.shape[0])]
        if chain_manager is None or len(chain_manager.chains)==0:
            logger.debug(f'System has no bondchains; cyclic C-C bondchain checking is skipped')
            return new_bdf
        # make a temporary working copy of the bondchains structure from its snapshot
        tmp_local_cm=deepcopy(chain_manager)
        existing_cyclic_chains=[x for x in tmp_local_cm.chains if x.is_cyclic]
        logger.debug(f'Existing cyclic C-C chains: {[c.idx for c in existing_cyclic_chains]}')
        # bondchains=[working_bondlist(i,self.idx_lists['bondchain'][i][:]) for i in range(len(self.idx_lists['bondchain']))]
        for i,r in bdf.iterrows():
            tmp_local_cm.injest_bond(r.ai,r.aj)
            c=tmp_local_cm.chain_of(r.ai)
            if c is None:
                # Neither endpoint participates in a vinyl C-C bondchain
                # (e.g. urethane O-C cure bonds), so this bond cannot
                # contribute to a C-C bondcycle — skip it.
                continue
            if not hasattr(c,'added_bonds'):
                c.added_bonds=[]
            c.added_bonds.append(i)
            # addlink(bondchains,i,r.ai,r.aj)
        for c in tmp_local_cm.chains:
            if c.is_cyclic and not c in existing_cyclic_chains:
                logger.debug(f'Bondchain {c.idx} ({"-".join([str(x) for x in c.idx_list])}) is a new cycle')
                logger.debug(f'-> New bonds added: {",".join([str(x) for x in c.added_bonds])}')
                logger.debug(f'-> Cycle length is {len(c.idx_list)} atoms')
                if self.min_bondcycle_length<=0 or len(c.idx_list)<self.min_bondcycle_length:
                    longest_bond_index=max(c.added_bonds) # remember, bdf is sorted by bondlength
                    logger.debug(f'-> limit is {self.min_bondcycle_length}, so we will break the longest added bond')
                    longest_bond_length=new_bdf.loc[longest_bond_index,'r']
                    logger.debug(f'-> bond index {longest_bond_index} (length {longest_bond_length:.3f} nm) is the longest added bond in the cycle')
                    new_bdf.loc[longest_bond_index,'remove-to-uncyclize']=True

        #     if rdict.get('makescycle',False):
        #         chains_of_bonds[i]=rdict["cyclizedchain"]
        #         cyclized_chains.append(rdict["cyclizedchain"])
        #     else:
        #         if rdict.get('retainedchain',-1)>-1:
        #             # a chain merger happened
        #             chains_of_bonds[i]=rdict["retainedchain"]
        #             for k in chains_of_bonds.keys():
        #                 if chains_of_bonds[k]==rdict["removedchain"]:
        #                     chains_of_bonds[k]=rdict["retainedchain"]
        #         else:
        #             # no chain merger happened; do nothing
        #             pass
        # # for each bondcycle, identify the new bonds that collectively cause the cyclization
        # cyclizing_bonds={}
        # for k,v in chains_of_bonds.items():
        #     if v in cyclized_chains:
        #         if not v in cyclizing_bonds:
        #             cyclizing_bonds[v]=[]
        #         cyclizing_bonds[v].append(k)
        
        # for cycleidx,newbonds in cyclizing_bonds.items():
        #     cyclelength=len(chains[cycleidx])
        # #     # logger.debug(f'bond {i} {r.ai}-{r.aj} ({r.reactantName}) adds to bondchain {chain_idx}')
        # #     chains_of_bonds.append(chain_idx)
        # #     if cyclized and not chain_idx in cyclized_chains:
        # #         cyclized_chains[chain_idx]=[]
        # # assert len(chains_of_bonds)==bdf.shape[0]
        # for i,r in bdf.iterrows():
        #     # logger.debug(f'compiling bonds of cyclized chains: bond {i}')
        #     cidx=chains_of_bonds[i]
        #     # logger.debug(f'compiling bonds of cyclized chains: cidx {chains_of_bonds[i]}')
        #     if cidx in cyclized_chains:
        #         cyclized_chains[cidx].append(i)
        # for k,v in cyclized_chains.items():
        #     if self.min_bondcycle_length<=0 or len(v)<self.min_bondcycle_length:
        #         lb=v[-1]  # last bondrec in list
        #         new_bdf.loc[lb,'remove-to-uncyclize']=True
        return new_bdf

    def get_bystanders(self,atom_idx):
        """Identifies and returns bystanders at a particular proposed bond specified by atom_idx.

        Args:
            atom_idx: container of two global atom indices specifying a proposed bond

        Returns:
            tuple: 4-tuple of special bystander info containers
        """
        bystander_atomidx=[[int],[int]] # elements correspond to atoms in atom_idx
        bystander_atomnames=[[str],[str]]
        bystander_resids=[[int],[int]]
        bystander_resnames=[[str],[str]]
        for x in [0,1]:
            bystander_atomidx[x]=self.interresidue_partners_of(atom_idx[x])
            if atom_idx[1-x] in bystander_atomidx[x]: # if new bond has not yet formed
                bystander_atomidx[x].remove(atom_idx[1-x])
            bystander_atomnames[x]=[self.get_gro_attribute_by_attributes('atomName',{'globalIdx':y}) for y in bystander_atomidx[x]]
            bystander_resids[x]=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':y}) for y in bystander_atomidx[x]]
            bystander_resnames[x]=[self.get_gro_attribute_by_attributes('resName',{'globalIdx':y}) for y in bystander_atomidx[x]]
        return bystander_resids,bystander_resnames,bystander_atomidx,bystander_atomnames
    
    def get_oneaways(self,atom_idx,chain_manager=None):
        """Identifies and returns one-aways for a particular proposed bond specified by atom_idx.

        Args:
            atom_idx (list): two-element container of atom indices specifying proposed bond
            chain_manager: optional ChainManager owned by the caller; returns empty lists if None

        Returns:
            tuple: tuple of oneaways-lists
        """
        oneaway_resids=[None,None]
        oneaway_resnames=[None,None]
        oneaway_atomidx=[None,None]
        oneaway_atomnames=[None,None]
        if chain_manager is None:
            return oneaway_resids,oneaway_resnames,oneaway_atomidx,oneaway_atomnames
        resids=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in atom_idx]
        chains=[self.get_gro_attribute_by_attributes('bondchain',{'globalIdx':x}) for x in atom_idx]
        chain_idx=[self.get_gro_attribute_by_attributes('bondchain_idx',{'globalIdx':x}) for x in atom_idx]
        logger.debug(f'atoms {atom_idx} chains {chains} chain_idx {chain_idx}')
        logger.debug(f'chain manager has {len(chain_manager.chains)} chains')
        # assert chain_idx[0]==0 or chain_idx[1]==0 # one must be a head
        if chains!=[-1,-1]:
            chainlists_idx=[chain_manager.chains[chains[x]].idx_list for x in [0,1]]
            logger.debug(f'chainlists_idx {chainlists_idx}')
            # chainlists_idx=[self.idx_lists['bondchain'][chains[x]] for x in [0,1]]
            chainlists_atomnames=[]
            chainlists_resids=[]
            chainlists_resnames=[]
            for x in [0,1]:
                atomnamelist=[]
                residlist=[]
                resnamelist=[]
                for y in chainlists_idx[x]:
                    rnum,rname,aname=self.get_gro_attribute_by_attributes(['resNum','resName','atomName'],{'globalIdx':y})
                    residlist.append(rnum)
                    resnamelist.append(rname)
                    atomnamelist.append(aname)
                chainlists_atomnames.append(atomnamelist)
                chainlists_resids.append(residlist)
                chainlists_resnames.append(resnamelist)
            # logger.debug(f'chainlists_idx: {chainlists_idx}')
            # logger.debug(f'chainlists_atomnames: {chainlists_atomnames}')
            # logger.debug(f'chainlists_resids: {chainlists_resids}')
            # logger.debug(f'chainlists_resnames: {chainlists_resnames}')
            oa_chain_idx=[]
            if chain_idx[0]==0 or chain_idx[0]>chain_idx[1]:
                oa_chain_idx=[chain_idx[0]+2,chain_idx[1]-2]
            elif chain_idx[1]==0 or chain_idx[1]>chain_idx[0]:
                oa_chain_idx=[chain_idx[0]-2,chain_idx[1]+2]
            # logger.debug(f'oa_chain_idx {oa_chain_idx}')
            for i in range(len(oa_chain_idx)):
                if oa_chain_idx[i]<0:
                    oa_chain_idx[i]=None
            oa_idx=[None,None]
            for i in [0,1]:
                try:
                    oa_idx[i]=chainlists_idx[i][oa_chain_idx[i]]
                except:
                    oa_idx[i]=None
                if oa_idx[i]:
                    oneaway_atomidx[i]=chainlists_idx[i][oa_chain_idx[i]]
                    oneaway_atomnames[i]=chainlists_atomnames[i][oa_chain_idx[i]]
                    oneaway_resids[i]=chainlists_resids[i][oa_chain_idx[i]]
                    oneaway_resnames[i]=chainlists_resnames[i][oa_chain_idx[i]]
        # logger.debug(f'result oneaway_resids {oneaway_resids} oneaway_resnames {oneaway_resnames}')
        # logger.debug(f'oneaway_atomidx {oneaway_atomidx} oneaway_atomnames {oneaway_atomnames}')
        return oneaway_resids,oneaway_resnames,oneaway_atomidx,oneaway_atomnames

    def grab_files(self):
        """Using absolute pathname information, grabs the most up-to-date gromacs files for this system and deposits them into the cwd.
        """
        cwd=os.getcwd()
        for ext in ['top','tpx','gro','grx']:
            filename=self.files[ext]
            logger.debug(f'{ext} grabbing {pfs.proj_abspath(filename)} into {pfs.cwd()}')  
            if os.path.commonprefix([cwd,filename])!=cwd:
                shutil.copy(filename,cwd)
            self.files[ext]=os.path.abspath(os.path.basename(filename))

    def grompp_and_mdrun(self,out,mdp,mylogger=logger.debug,**kwargs):
        """Manages invoking a single Gromacs run using the current TopoCoord.

        Args:
            out (str): output filename basename
            mdp (str): name of mdp file
            mylogger: a logger, defaults to logger.debug

        Returns:
            str: any message generated by gromacs.grompp_and_mdrun
        """
        wrap_coords=kwargs.get('wrap_coords',False)
        self.grab_files()
        # make sure required files exist in this directory
        top=os.path.basename(self.files['top']).replace('.top','')
        gro=os.path.basename(self.files['gro']).replace('.gro','')
        assert os.path.exists(f'{top}.top')
        assert os.path.exists(f'{gro}.gro')
        assert os.path.exists(f'{mdp}.mdp')
        mylogger(f'TopoCoord: running Gromacs {top}.top, {gro}.gro, {mdp}.mdp')
        msg=grompp_and_mdrun(gro=gro,top=top,out=out,mdp=mdp,**kwargs)
        self.copy_coords(TopoCoord(grofilename=f'{out}.gro',wrap_coords=wrap_coords))
        mylogger(f'after grompp_and_run: gro {self.files["gro"]}')
        return msg

    def load_files(self,filenames:dict):
        """Loads all gromacs files into TopoCoord.

        Args:
            filenames (dict): dictionary of extension:filename
        """
        for e,n in filenames.items():
            if not e in ['gro','top','grx','tpx','mol2']: continue
            bn,ext=os.path.splitext(n)
            logger.debug(f'bn {bn} ext {ext}')
            if   ext=='.gro':  self.read_gro(n)
            elif ext=='.top':  self.read_top(n)
            elif ext=='.tpx':  self.read_tpx(n)
            elif ext=='.grx':  self.read_gro_attributes(n)
            elif ext=='.mol2': self.read_mol2(n)
            else:
                logger.debug(f'Warning: file {n} has unknown file extension.  Skipped.')

    def center_coords(self,new_boxsize:np.ndarray=None):
        """Centers all coordinates in box.

        Args:
            new_boxsize (numpy.ndarray): new boxsize if desired, defaults to None
        """
        if type(new_boxsize)==np.ndarray:
            if new_boxsize.shape==(3,):
                box_vectors=new_boxsize*np.identity(3,dtype=float)
            elif new_boxsize.shape==(3,3):
                box_vectors=new_boxsize
            self.Coordinates.box=box_vectors
        center=self.Coordinates.box.diagonal()/2.0
        gc=self.Coordinates.geometric_center()
        addme=center-gc
        self.Coordinates.translate(addme)

    def vacuum_minimize(self,outname='minimized',**kwargs):
        """The minimize analog to grompp_and_mdrun; performs an energy minimization using mdrun.

        Args:
            outname (str): output file basename, defaults to 'minimized'
        """
        pad=kwargs.get('pad',5)
        boxsize=np.array(self.maxspan())+pad*np.ones(3)
        logger.debug(f'{self.maxspan()} -> {boxsize}')
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix='single-molecule-min'
        pfs.checkout(pfs.Dirs.mdp_file(mdp_prefix))
        # gromacs_dict={'nt':1,'nb':'cpu','pme':'cpu','pmefft':'cpu','bonded':'cpu','update':'cpu'}
        self.grompp_and_mdrun(out=f'{outname}',
            mdp=mdp_prefix,boxSize=boxsize,single_molecule=True,wrap_coords=False) #,**gromacs_dict)

    def vacuum_simulate(self,outname='simulated',**kwargs):
        """Performs a vacuum MD simulation using mdrun.

        Args:
            outname (str): output file basename, defaults to 'simulated'
        """
        pad=kwargs.get('pad',5)
        boxsize=np.array(self.maxspan())+pad*np.ones(3)
        logger.debug(f'{self.maxspan()} -> {boxsize}')
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix='single-molecule-nvt'
        pfs.checkout(pfs.Dirs.mdp_file(mdp_prefix))
        nsamples=kwargs.get('nsamples',10)
        params=kwargs.get('params',{})
        T=params.get('temperature',300.0)
        ps=params.get('ps',0.0)
        nsteps=params.get('nsteps',-1)
        begin_at=params.get('begin_at',0.0)
        assert ps!=0.0 or nsteps!=-1
        dt=float(mdp_get(f'{mdp_prefix}.mdp','dt'))
        if ps!=0.0:
            nsteps=int(float(ps)/dt)
        else:
            ps=nsteps*dt
        sampling_duration_ps=ps-begin_at
        sampling_duration_nsteps=int(sampling_duration_ps/dt)
        sample_interval=sampling_duration_nsteps//nsamples
        logger.debug(f'nsteps {nsteps} begin_at {begin_at} ps {ps} sampling_duration_ps {sampling_duration_ps} ({sampling_duration_nsteps} steps) sample_interval {sample_interval}')
        # nsteps=nsamples*(sample_interval+1)
        # nsteps=mdp_get(f'{mdp_prefix}.mdp','nsteps')
        # nstxout=mdp_get(f'{mdp_prefix}.mdp','nstxout')
        mdp_modify(f'{mdp_prefix}.mdp',{'nsteps':nsteps,'nstxout':sample_interval,'ref_t':T})
        # gromacs_dict={'nt':1,'nb':'cpu','pme':'cpu','pmefft':'cpu','bonded':'cpu','update':'cpu'}
        self.grompp_and_mdrun(out=f'{outname}',
            mdp=mdp_prefix,boxSize=boxsize,single_molecule=True,wrap_coords=False) #,**gromacs_dict)

    def equilibrate(self,deffnm='equilibrate',edict={},gromacs_dict={}):
        """Performs an MD simulation using mdrun.

        Args:
            deffnm (str): output file basename, defaults to 'equilibrate'
            edict (dict): dictionary of simulation directives, defaults to {}
            gromacs_dict (dict): dictionary of gromacs directives, defaults to {}

        Returns:
            list: list of edr files this equilibration generates
        """
        mod_dict={}
        edr_list=[]
        ens=edict['ensemble']
        repeat=edict.get('repeat',0)
        assert ens in ['min','npt','nvt'],f'Bad ensemble: {ens}'
        pfs.checkout(pfs.Dirs.mdp_file(ens)) # plain jane?
        if ens=='min':
            msg='minimization'
        else:
            msg=f'{ens} ensemble'
        if ens in ['nvt','npt']:
            dt=float(mdp_get(f'{ens}.mdp','dt'))
            mod_dict['ref_t']=edict['temperature']
            mod_dict['gen-temp']=edict['temperature']
            mod_dict['gen-vel']='yes'
            nsteps=edict.get('nsteps',0)
            if not nsteps:
                ps=edict.get('ps',0.0)
                if not ps: return
                nsteps=int(float(ps)/dt)
            else:
                ps=nsteps*dt
            mod_dict['nsteps']=nsteps
            msg+=f'; {ps:6.2f} ps, {float(edict["temperature"]):7.2f} K'
            if ens=='npt': 
                mod_dict['ref_p']=edict['pressure']
                msg+=f', {float(edict["pressure"]):6.2f} bar'
        new_mdp=f'{deffnm}-{ens}'
        mdp_modify(f'{ens}.mdp',mod_dict,new_filename=f'{new_mdp}.mdp')
        logger.info(f'Running Gromacs: {msg}')
        self.grompp_and_mdrun(out=f'{deffnm}-{ens}',mdp=new_mdp,quiet=False,**gromacs_dict)
        edr_list=[f'{deffnm}-{ens}']
        density=None
        if ens=='npt':
            self._log_box()
            density=self._density_series(f'{deffnm}-{ens}',gromacs_dict)
        for rep in range(repeat):
            logger.info(f'Repeat {rep+1} out of {repeat}')
            this_deffnm=f'{deffnm}-repeat-{rep+1}-{ens}'
            self.grompp_and_mdrun(out=this_deffnm,mdp=new_mdp,quiet=False,**gromacs_dict)
            edr_list.append(this_deffnm)
            if ens=='npt':
                self._log_box()
            density=self._density_series(this_deffnm,gromacs_dict)
        converge=edict.get('converge',{}) or {}
        if converge and ens=='npt':
            edr_list+=self._converge_density(deffnm,new_mdp,converge,density,gromacs_dict)
        elif converge:
            logger.warning(f'ignoring "converge" on a {ens} stage; density convergence needs npt')
        return edr_list

    def _log_box(self):
        """Logs the current box side lengths."""
        box=self.Coordinates.box.diagonal()
        logger.info(f'Current box side lengths: {box[0]:.3f} nm x {box[1]:.3f} nm x {box[2]:.3f} nm')

    def _density_series(self,edr_pfx,gromacs_dict={}):
        """Density trace from one edr, as a plain array, or None.

        Args:
            edr_pfx (str): edr basename, without the '.edr' extension
            gromacs_dict (dict): gromacs directives passed through to gmx energy

        Returns:
            numpy.ndarray: the Density column in kg/m^3, or None if unavailable
        """
        try:
            df=gmx_energy_trace(edr_pfx,['Density'],report_averages=True,**gromacs_dict)
        except Exception as e:
            logger.debug(f'no density available from {edr_pfx}: {e}')
            return None
        if df.empty or 'Density' not in df.columns: return None
        return df['Density'].to_numpy(dtype=float)

    def _converge_density(self,deffnm,mdp,converge,density,gromacs_dict={}):
        """Repeats an NPT stage until its density settles, or until a ceiling.

        This is the fixed-``repeat`` mechanism with a measured stopping rule
        instead of a counted one.  Densification is the right place for it:
        the 200 -> ~1100 kg/m^3 compaction is one-shot, involves a large volume
        change, and otherwise runs for however many steps someone guessed.

        Each extension is a fresh ``mdrun``, so velocities are regenerated
        between segments exactly as ``repeat`` already does; the criterion is
        therefore applied to the *last* segment rather than to a concatenation
        of segments that do not share a velocity history.

        A build that never settles says so, loudly, rather than silently
        reporting whatever density it reached.

        Args:
            deffnm (str): output file basename of the stage being extended
            mdp (str): name of the mdp file to rerun
            converge (dict): 'tolerance' (kg/m^3), 'max_repeats', 'min_samples', 'drift_sems'
            density (numpy.ndarray): density trace of the stage just run, or None
            gromacs_dict (dict): gromacs directives

        Returns:
            list: edr basenames of any extension segments run
        """
        tol=float(converge.get('tolerance',1.0))
        ceiling=int(converge.get('max_repeats',10))
        min_samples=int(converge.get('min_samples',50))
        drift_sems=float(converge.get('drift_sems',2.0))
        extra=[]
        for attempt in range(ceiling+1):
            if density is None:
                logger.warning('density convergence requested but no Density trace is available; not gating')
                return extra
            r=series_converged(density,tol,min_samples=min_samples,drift_sems=drift_sems)
            if r['converged']:
                logger.info(
                    f'Density converged after {attempt} extension(s): '
                    f'{r["mean"]:.1f} +/- {r["sem"]:.2f} kg/m^3 '
                    f'(tau_int {r["tau_int"]:.0f} samples, drift {r["drift"]:+.2f})'
                )
                return extra
            if attempt==ceiling:
                logger.warning(
                    f'Density did NOT converge after {ceiling} extension(s) -- {r["reason"]}. '
                    f'The last window averaged {r["mean"]:.1f} kg/m^3; treat that as an '
                    f'unsettled box, not as this system\'s density'
                )
                return extra
            logger.info(f'Density not settled ({r["reason"]}); extending, {attempt+1} of {ceiling}')
            this_deffnm=f'{deffnm}-converge-{attempt+1}-npt'
            self.grompp_and_mdrun(out=this_deffnm,mdp=mdp,quiet=False,**gromacs_dict)
            extra.append(this_deffnm)
            self._log_box()
            density=self._density_series(this_deffnm,gromacs_dict)
        return extra

    def get_resid_sets(self,atom_pair):
        """Identifies individual sets of separate resids owned by unbonded atoms i and j.

        Args:
            atom_pair (tuple): tuple of i, j atom indices

        Returns:
            tuple: tuple containing the two apposing residue sets
        """
        # assertion: i and j are not bonded and they represent two separate sets of residues in this topocoord.
        i,j=atom_pair
        ri=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':i})
        rj=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':j})
        ci=self.Topology.local_resid_cluster(ri)
        cj=self.Topology.local_resid_cluster(rj)
        if set(ci)==set(cj):
            logger.debug(f'resid sets overlap: {i} {ri} : {ci} and {j} {rj} : {cj}')
            return []
        return [ci,cj]

    def check_your_topology(self):
        """Checks topology for duplicate 1-4 pair interactions and deletes them.
        """
        T=self.Topology
        C=self.Coordinates
        aT={}
        pdf=T.D['pairs']
        aT['dx']=[]
        aT['dy']=[]
        aT['dz']=[]
        checked=[]
        drops=[]
        for i,r in pdf.iterrows():
            ai,aj=min([r.ai,r.aj]),max([r.ai,r.aj])
            assert ai!=aj
            if not (ai,aj) in checked: 
                checked.append((ai,aj))
            else:
                # print(f'repeated pair {ai} {aj}')
                drops.append(i)
            ri=self.get_R(ai)
            rj=self.get_R(aj)
            D=C.mic(ri-rj,pbc=[1,1,1])
            aT['dx'],aT['dy'],aT['dz']=D
        print(f'{len(drops)} duplicate pairs detected')
        T.D['pairs']=pdf.drop(drops)
        pdf['dx']=aT['dx']
        pdf['dy']=aT['dy']
        pdf['dz']=aT['dz']
        print(pdf.sort_values(by='dx').head(3).to_string())
        print(pdf.sort_values(by='dy').head(3).to_string())
        print(pdf.sort_values(by='dz').head(3).to_string())
        self.write_top('checked.top')

    def flip_stereocenters(self,idxlist):
        """Flips stereochemistry of atoms in idxlist.

        Args:
            idxlist (list): global indices of chiral atoms
        """
        A=self.Coordinates.A
        g=self.Topology.bondlist.graph()
        for idx in idxlist:
            ligand_idx=self.Topology.bondlist.partners_of(idx)
            O=self.get_R(idx)
            name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':idx})
            if len(ligand_idx)!=4:
                logger.debug(f'Atom {idx} cannot be a stereocenter; it only has {len(ligand_idx)} ligands.')
            else:
                logger.debug(f'Flipping stereochemistry at atom {idx}')
                # do stuff
                cg=g.copy()
                # remove the stereocenter and see what unconnected components are generated
                cg.remove_node(idx)
                cc=[cg.subgraph(c).copy() for c in nx.connected_components(cg)]
                ncc=len(cc)
                logger.debug(f'At stereocenter {idx} there are {ncc} unconnected ligands')
                # for each ligand-group, assign one or more stereocenter ligand atoms
                associates={k:[] for k in range(ncc)}
                for i in ligand_idx:
                    for icc in range(ncc):
                        tcc=cc[icc]
                        if i in tcc:
                            associates[icc].append(i)
                sizes=[]
                for icc in range(ncc):
                    sizes.append(len(cc[icc]))
                sizes=np.array(sizes)
                sarg=sizes.argsort()
                logger.debug(f'Atom {idx} grouped associate ligand subgraphs {associates}')
                logger.debug(f'Subgroups sorted by size:')
                for sa in sarg:
                    logger.debug(f'  {sa}: {" ".join([self.get_gro_attribute_by_attributes("atomName",{"globalIdx":x}) for x in cc[sa]])}')

                # The flip is a 180-degree rotation of "half" of the atoms connected to the stereocenter
                # about an axis defined by (a) the location of the stereocenter and (b) the midpoint on
                # the vector connecting two specially-selected ligand atoms.
                #
                # How the two ligand atoms are selected:
                # 1. If there are four unconnected ligand groups, then the two groups with the fewest atoms
                #    are selected.
                # 2. If there are three unconnected groups, then one of them has two ligand atoms.
                #    These two atoms are selected if that group is smaller than the union of
                #    the other two groups.  If that group is not smaller than the union of the
                #    other two groups, then the ligand atoms of those other two groups are selected.
                # 3. If there are two unconnected groups, AND both have two ligand atoms, the group with
                #    the smaller number of atoms has its ligand atoms selected.  If one group has three and
                #    the other has one, we cannot flip this stereocenter using this method.  However,
                #    we could perform a reflection of all atoms in this group through the plane defined by
                #    (a) the stereocenter, (b) one of three ligand atoms, and (c) the midpoint between the
                #    the other two.  However, this would also flip all stereocenters in that group.
                #    So for now, we throw an error.
                # 4. If there is only one unconnected group, there is no way we can do this.
                if ncc==4:
                    l1idx=associates[sarg[0]][0]
                    l2idx=associates[sarg[1]][0]
                    l1name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':l1idx})
                    l2name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':l1idx})
                    logger.debug(f'ncc {ncc}: Atoms {l1name} and {l2name} are the selected ligands of stereocenter {name}')
                    p1=self.get_R(l1idx)
                    p2=self.get_R(l2idx)
                    mp=0.5*(p1+p2)
                    ax=O-mp
                    g1idx=list(cc[sarg[0]])
                    g2idx=list(cc[sarg[1]])
                    gidx=g1idx+g2idx # set of indices of rotating group
                    M=Matrix4()
                    M.translate(-O).rotate_axis(180.0,ax).translate(O)
                    self.Coordinates.homog_trans(M,indices=gidx)
                elif ncc==3:
                    # largest group must have the two connected atoms?  not necessarily.
                    the_idx_of_the_double=[x for x in sarg if len(associates[x])==2][0]
                    other_idx=[x for x in sarg if len(associates[x])==1]
                    assert len(other_idx)==2
                    double_size=len(cc[the_idx_of_the_double])
                    joint_size=sum([len(cc[x]) for x in other_idx])
                    logger.debug(f'ncc {ncc}: Size of UG with two ligands: {double_size}; side of union of UGs with one ligand each: {joint_size}')
                    if joint_size<double_size:
                        l1idx=associates[other_idx[0]][0]
                        l2idx=associates[other_idx[1]][0]
                        gidx=list(cc[other_idx[0]])+list(cc[other_idx[1]])
                    else:
                        l1idx=associates[the_idx_of_the_double][0]
                        l1idx=associates[the_idx_of_the_double][1]
                        gidx=list(cc[the_idx_of_the_double])
                    l1name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':l1idx})
                    l2name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':l2idx})
                    logger.debug(f'ncc {ncc}: Atoms {l1name} and {l2name} are the selected ligands of stereocenter {name}')
                    p1=self.get_R(l1idx)
                    p2=self.get_R(l2idx)
                    mp=0.5*(p1+p2)
                    ax=O-mp
                    M=Matrix4()
                    M.translate(-O).rotate_axis(180.0,ax).translate(O)
                    self.Coordinates.homog_trans(M,indices=gidx)
                elif ncc==2:
                    if all([len(x)==2 for x in associates.values()]):
                        l1idx=associates[sarg[0]][0]
                        l2idx=associates[sarg[0]][1]
                        l1name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':l1idx})
                        l2name=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':l2idx})
                        logger.debug(f'ncc {ncc}: Atoms {l1name} and {l2name} are the selected ligands of stereocenter {name}')
                        p1=self.get_R(l1idx)
                        p2=self.get_R(l2idx)
                        mp=0.5*(p1+p2)
                        ax=mp-O
                        logger.debug(f'p1 {p1} p2 {p2} mp {mp} ax {ax} O {O}')
                        gidx=list(cc[sarg[0]])
                        logger.debug(f'The rotating group is {" ".join([self.get_gro_attribute_by_attributes("atomName",{"globalIdx":x}) for x in gidx])}')
                        M=Matrix4()
                        M.translate(-O).rotate_axis(180.0,ax).translate(O)
                        self.Coordinates.homog_trans(M,indices=gidx)
                    else:
                        logger.debug(f'Arrangment of unconnected ligand groups at stereocenter {idx} does not permit rotation.')
                    pass
                elif ncc==4:
                    logger.debug(f'Cannot flip at stereocenter {idx}.')


def find_template(BT:BondTemplate,moldict):
    """Searches the dictionary of available molecules to identify a bond template that matches the passed-in template, returning the corresponding template molecule and reaction-bond.

    Args:
        BT (BondTemplate): bond template to search for
        moldict (MoleculeDict): dictionary of available molecules

    Raises:
        Exception: if no matching template is found

    Returns:
        tuple: template Molecule object, corresponding ReactionBond object from that template, and a boolean flag indicating whether or not the match required a symmetric-reversal of the template
    """
    # Collect every (template, bond-index, reverse) tuple whose bond_template
    # matches the instance under subset-bystander semantics, then pick the
    # most-specific one (most bystanders declared).  This keeps small-fragment
    # cure templates (no bystanders) matching instances embedded in larger
    # assembled monomers while still letting chain-extension templates win
    # against the bare-dimer template when their additional bystander/oneaway
    # context is exactly the in-chain instance's.
    best = None  # (specificity, template, b_idx, reverse_bond)
    for template_name, T in moldict.items():
        for b_idx in range(len(T.bond_templates)):
            bt = T.bond_templates[b_idx]
            if bt.matches(BT):
                spec = bt.bystander_count()
                if best is None or spec > best[0]:
                    best = (spec, T, b_idx, False)
            elif bt.matches_reverse_of(BT):
                spec = bt.bystander_count()
                if best is None or spec > best[0]:
                    best = (spec, T, b_idx, True)
    if best is None:
        logger.error(f'No template is found for {str(BT)}')
        raise Exception('you have a bond for which I cannot find a template')
    _, use_T, b_idx, reverse_bond = best
    logger.debug(f'Using template {use_T.name} and bond index {b_idx}'
                 f' (specificity {best[0]}, reversed={reverse_bond})')
    rb = use_T.reaction_bonds[b_idx]
    return use_T, rb, reverse_bond

