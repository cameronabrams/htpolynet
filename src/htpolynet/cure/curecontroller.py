"""Manages execution of the CURE algorithm.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging
import os

from enum import Enum
from functools import partial
from itertools import product
from multiprocessing import Pool

import yaml

import numpy as np
import pandas as pd

from ..analysis.plot import trace
from ..core import projectfilesystem as pfs
from ..core.molecule import MoleculeDict
from ..core.topocoord import TopoCoord, BTRC
from ..cure.reaction import ReactionList
from ..cure.reaction import reaction_stage
from ..external.gromacs import gmx_energy_trace, gromacs_distance, mdp_modify
from ..utils import checkpoint as cp
from ..utils.stringthings import my_logger

logger = logging.getLogger(__name__)

def residue_reaction_counts(adf):
    """Counts, per residue, how many bonds that residue's atoms have already formed.

    Args:
        adf (pandas.DataFrame): the global atom dataframe; must carry the 'resNum' and 'nreactions' attributes

    Returns:
        pandas.Series: number of bonds already formed, indexed by residue number
    """
    return adf.groupby('resNum')['nreactions'].sum()

def cure_site_mask(adf,RL,MD):
    """Marks the atoms a cure-stage reaction can ever bond.

    A residue's ``z + nreactions`` counts every reactive site it was built
    with, and that includes bonds made at the param and build stages and
    valences no cure reaction can consume -- in example 2 each HIE carries a
    param-stage bond and a leftover valence on C4, so its four "sites" include
    two the cure never touches.  Functionality for asking whether the cure made
    a network has to be counted on cure sites only.

    The site set is read from the runtime reaction list, which has already been
    symmetry- and stereo-expanded: a reaction written against ``TAZ:C1`` has
    siblings for C2 and C3 there.  Reading the configuration's reactions instead
    would name one atom per symmetry class and undercount -- that version was
    tried and silenced example 6's genuine check.  Residue names are resolved the
    way the bond search resolves them, through the reactant's sequence.

    Args:
        adf (pandas.DataFrame): the global atom dataframe; must carry 'resName' and 'atomName'
        RL (ReactionList): the runtime, symmetry-expanded reaction list
        MD (MoleculeDict): molecule templates, for resolving a reactant's residue names

    Returns:
        pandas.Series: boolean per atom, aligned to adf; None if no cure sites can be determined
    """
    if RL is None or MD is None or not all(c in adf.columns for c in ('resName','atomName')):
        return None
    sites=set()
    for R in RL:
        if R.stage!=reaction_stage.cure: continue
        for rec in R.atoms.values():
            try:
                resname=MD[R.reactants[rec['reactant']]].get_resname(rec['resid'])
            except (KeyError,IndexError,AttributeError,TypeError):
                continue
            sites.add((resname,rec['atom']))
    if not sites:
        return None
    return pd.Series([k in sites for k in zip(adf['resName'],adf['atomName'])],index=adf.index)

def residue_functionality(adf):
    """Counts, per residue, how many reactive sites it started with.

    Forming a bond decrements an atom's 'z' and increments its 'nreactions' in
    the same operation, so their sum is conserved per atom for the whole cure.
    Summed over a residue it is that residue's initial functionality, and it
    reads the same at iteration 0 as at the end.

    Args:
        adf (pandas.DataFrame): the global atom dataframe; must carry the 'resNum', 'z' and 'nreactions' attributes

    Returns:
        pandas.Series: initial reactive sites, indexed by residue number
    """
    return adf.groupby('resNum')[['z','nreactions']].sum().sum(axis=1)

def minimum_image_displacements(before,after,box):
    """Per-atom displacement magnitudes across a window, under minimum image.

    Args:
        before (numpy.ndarray): (N,3) positions before the window, nm
        after (numpy.ndarray): (N,3) positions after the window, nm
        box (numpy.ndarray): (3,3) box matrix; only its diagonal is used

    Returns:
        numpy.ndarray: (N,) displacement magnitudes, nm; empty if N is 0
    """
    b=np.asarray(before,dtype=float).reshape(-1,3)
    a=np.asarray(after,dtype=float).reshape(-1,3)
    assert a.shape==b.shape,f'shape mismatch {b.shape} vs {a.shape}'
    d=a-b
    L=np.asarray(box,dtype=float)
    L=L.diagonal() if L.ndim==2 else L.reshape(-1)
    # an unset box would make every displacement nan; leave the raw
    # displacement alone instead, which is right for a non-periodic snapshot
    if d.size and np.all(L>0):
        d=d-L*np.round(d/L)
    return np.sqrt((d**2).sum(axis=1))

def varshney_criterion(displacements,search_radius):
    """Summarizes reactive-species mobility against the bond-search capture radius.

    Varshney sized the original 40 ps relaxation window on the requirement that
    unreacted reactive species diffuse far enough during it to find new
    partners.  That requirement decays over a cure: as those species are bonded
    into the growing network they become topologically constrained rather than
    merely slow, and the fixed window silently stops satisfying the criterion it
    was designed against.  Nothing in the build reports this, so a window that
    has stopped working looks exactly like one that still works.

    'crossing_fraction' is the share of still-reactive atoms that moved at least
    one capture radius during the window; it is the quantity that decays.

    Args:
        displacements (numpy.ndarray): per-atom displacement magnitudes, nm, from :func:`minimum_image_displacements`
        search_radius (float): bond-search capture radius, nm

    Returns:
        dict: 'n', 'rmsd', 'mean', 'max' (nm) and 'crossing_fraction'; all but 'n' are None when there are no reactive atoms left
    """
    r=np.asarray(displacements,dtype=float).reshape(-1)
    if r.size==0:
        return {'n':0,'rmsd':None,'mean':None,'max':None,'crossing_fraction':None}
    return {
        'n':int(r.size),
        'rmsd':float(np.sqrt((r**2).mean())),
        'mean':float(r.mean()),
        'max':float(r.max()),
        'crossing_fraction':float((r>=search_radius).mean())
    }

def _relax_stage_density(edr_pfx):
    """Mean density over one relax stage's NPT segment, or None if unavailable.

    Instrumentation must never be able to fail a build, so every failure path --
    a missing edr, a gmx energy that will not run, an ensemble whose menu has no
    Density -- returns None, and the stage reports '--' instead.

    Args:
        edr_pfx (str): edr basename, without the '.edr' extension

    Returns:
        float: mean density in kg/m^3, or None
    """
    try:
        df=gmx_energy_trace(edr_pfx,['Density'])
    except Exception as e:
        logger.debug(f'no density available from {edr_pfx}: {e}')
        return None
    if df.empty or 'Density' not in df.columns: return None
    return float(df['Density'].mean())

def rank_bond_candidates(bdf,reaction_counts=None):
    """Orders bond candidates ahead of the greedy one-bond-per-residue downselection.

    The default ordering is by pair separation alone.  Separation is
    uncorrelated with how many bonds a crosslinker already carries, so bonds
    spread evenly over all crosslinkers rather than completing the
    partly-reacted ones.  Supplying ``reaction_counts`` promotes candidates
    whose B-side residue already carries the most bonds, with each such group
    still ordered by separation; because the downselection that follows is
    greedy and admits at most one bond per residue per iteration, this ordering
    is what decides which residues react.

    Args:
        bdf (pandas.DataFrame): bond candidates, with column 'r' and, when biasing, column 'rj'
        reaction_counts (pandas.Series): bonds already formed per residue, as returned by :func:`residue_reaction_counts`; None (the default) selects the distance-only ordering

    Returns:
        pandas.DataFrame: bdf reordered; its columns are unchanged
    """
    if reaction_counts is None:
        return bdf.sort_values('r',axis=0,ignore_index=True).reset_index(drop=True)
    ranked=bdf.copy()
    ranked['_nrx']=ranked['rj'].map(reaction_counts).fillna(0).astype(int)
    ranked=ranked.sort_values(['_nrx','r'],ascending=[False,True],ignore_index=True)
    return ranked.drop(columns=['_nrx']).reset_index(drop=True)

class cure_step(Enum):
    """Enumerated CURE step
    """
    cure_bondsearch=0
    cure_drag=1
    cure_update=2
    cure_relax=3
    cure_equilibrate=4
    cap_bondsearch=5
    cap_update=6  # dragging is not needed/allowed in capping since these are likely intramolecular bonds...
    cap_relax=7
    cap_equilibrate=8
    finished=9
    unknown=99
    def __str__(self):
        return self.name
    def basename(self):
        name=str(self)
        if 'cap_' in name:
            return name[len('cap_'):]
        if 'cure_' in name:
            return name[len('cure_'):]

class CureState:
    def __init__(self):
        """Generates a new, initialized CureState object.
        """
        self.iter=0
        self.max_nxlinkbonds=0
        self.cum_nxlinkbonds=0
        self.desired_nxlinkbonds=0
        self.max_search_radius=0.0
        self.max_radidx=0
        self.step=cure_step.unknown
        self.curr_nxlinkbonds=0

        self.current_stage={}
        self.current_stage['drag']=0
        self.current_stage['relax']=0

        self.current_radidx=0
        self.current_radius=0.0

        self.bonds_file=''

    def _to_yaml(self,filename='../cure_state.yaml'):
        """Writes CureState object to YAML file.

        Args:
            filename (str): name of output YAML file, defaults to '../cure_state.yaml'
        """
        with open(filename,'w') as f:
            f.write(yaml.dump(self))

    @classmethod
    def from_yaml(cls,filename='cure_state.yaml'):
        """Returns a new CureState object generated by reading in the YAML file.

        Args:
            filename (str): input file name, defaults to 'cure_state.yaml'

        Returns:
            CureState: new CureState object returned by yaml.load
        """
        with open(filename,'r') as f:
            yaml_string=f.read()
        # FullLoader was tightened in recent PyYAML versions to reject Python
        # object tags; use yaml.Loader to round-trip yaml.dump(self) output.
        return yaml.load(yaml_string,Loader=yaml.Loader)

    def reset(self):
        """Resets this CureState object to begin a new CURE iteration.
        """
        self.iter=1
        self.step=cure_step.cure_bondsearch
        self.curr_nxlinkbonds=0
        self.current_stage['drag']=0
        self.current_stage['relax']=0
        self.current_radidx=0
        self.current_radius=0.0
        self.bonds_file=''

class CureController:
    default_equilibration_sequence = [ { 'ensemble': 'min' }, 
        { 'ensemble': 'nvt', 'temperature': 600, 'nsteps': 1000 },
        { 'ensemble': 'npt', 'temperature': 600, 'pressure': 1, 'nsteps': 2000 }
    ]
    curedict_defaults={
        'output': {
            'bonds_file': 'bonds.csv'
        },
        'controls': {
            'max_conversion_per_iteration': 1.0,
            'search_radius': 0.5,
            'radial_increment': 0.05,
            'late_threshold': 0.85,
            'max_iterations': 100,
            'desired_conversion': 0.5,
            'min_allowable_bondcycle_length':-1, # not set
            'min_bonds_per_iteration': 10, # grow radius until >= this many bonds (or max radius)
            'completion_bias': False, # rank bond candidates by B-side completion before distance
            'ncpu' : os.cpu_count()
        },
        'drag': {
            'limit': 0.0,
            'kb': 300000.0, # kJ/mol/nm
            'trigger_distance': 0.0,
            'nstages': 0,
            'increment': 0.0,
            'cutoff_pad': 0.2,
            'equilibration': default_equilibration_sequence
        },
        'relax': {
            'nstages': 6,
            # not 0.0: _distance_attenuation derives its stage count as
            # int(maxL/increment) with no guard, so a config that omits this
            # key used to raise ZeroDivisionError at the first relax.  Drag
            # can default to 0.0 because its own 0.0 limit disables it.
            'increment': 0.08,
            'cutoff_pad': 0.2,
            'equilibration': default_equilibration_sequence
        },
        'equilibrate': {
            'temperature': 300,
            'pressure': 1,
            'nsteps': 50000,
            'ensemble': 'npt'
        },
        'gromacs': {
            'rdefault': 0.9
        }
    }

    def __init__(self,curedict={}):
        """Generates a new CureController object populated from directives in curedict.

        Args:
            curedict (dict): dictionary of cure directives, defaults to {}
        """
        self.state=CureState()
        self.bonds_df:pd.DataFrame=None
        self.bonds_are='nonexistent!'
        self.search_failed=False
        self.dicts={}
        for k,v in self.curedict_defaults.items():
            self.dicts[k]=curedict.get(k,v)  # loads entire default if dict is just missing
            if type(v)==dict: # assign subitem defaults if their keys are missing
                for kk,vv in v.items():
                    if not kk in self.dicts[k]:
                        self.dicts[k][kk]=vv

        self.chain_manager=None
        # once-per-run: has the completion bias been checked against which side
        # of the reaction actually carries the more functional residue?
        self.bias_side_checked=False
        self.dragging_enabled=False
        d=self.dicts['drag']
        if (d['nstages']>0 or d['increment']>0.0) and d['limit']>0.0:
            self.dragging_enabled=True

    def setup(self,max_nxlinkbonds=0,desired_nxlinkbonds=0,max_search_radius=0.0):
        """Sets up this CureController object.

        Args:
            max_nxlinkbonds (int): maximum allowable number of bonds, defaults to 0
            desired_nxlinkbonds (int): desired number of bonds, defaults to 0
            max_search_radius (float): maximum search radius (nm), defaults to 0.0
        """
        st=self.state
        st.max_nxlinkbonds=max_nxlinkbonds
        st.desired_nxlinkbonds=desired_nxlinkbonds
        st.max_search_radius=max_search_radius
        d=self.dicts['controls']
        st.max_radidx=int((st.max_search_radius-d['search_radius'])/d['radial_increment'])

    @cp.enableCheckpoint
    def do_iter(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,gromacs_dict={}):
        """Performs one CURE iteration.

        Args:
            TC (TopoCoord): TopoCoord object; all topological and coordinate information is in here
            RL (ReactionList): list of all reaction types
            MD (MoleculeDict): dictionary of all molecules, keyed by name
            gromacs_dict (dict): dictionary of custom gromacs parameters, defaults to {}

        Returns:
            dict: dictionary of resulting output file names
        """
        my_logger(f'Iteration {self.state.iter} begins',logger.info,fill='~')
        TC.grab_files() # copy files locally
        self.state.step=cure_step.cure_bondsearch
        self._do_bondsearch(TC,RL,MD)
        self._do_preupdate_dragging(TC,gromacs_dict)
        self._do_topology_update(TC,MD)
        self._do_relax(TC,gromacs_dict)
        self._do_equilibrate(TC,gromacs_dict)
        self.state.cum_nxlinkbonds+=self.bonds_df.shape[0]
        logger.info(f'Iteration {self.state.iter} current conversion {self._curr_conversion():.3f} or {self.state.cum_nxlinkbonds} bonds')
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    @cp.enableCheckpoint
    def do_capping(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,gromacs_dict={}):
        """Manages generation of all capping bonds.

        Args:
            TC (TopoCoord): TopoCoord object containing all topology and coordinate information
            RL (ReactionList): list of all reaction types
            MD (MoleculeDict): dictionary of all molecule information
            gromacs_dict (dict): dictionary of custom gromacs parameters, defaults to {}

        Returns:
            dict: dictionary of resulting output files, keyed by extension ('gro', 'top', 'grx')
        """
        self._do_cap_bondsearch(TC,RL,MD)
        self._do_topology_update(TC,MD)
        self._do_relax(TC)
        self._do_equilibrate(TC,gromacs_dict)
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    def is_cured(self):
        """Returns True if system is cured.

        Returns:
            bool: True if system is cured; False otherwise
        """
        st=self.state
        logger.debug(f'search_failed? {self.search_failed}')
        logger.debug(f'cumxlinks {st.cum_nxlinkbonds} maxxlinks {st.max_nxlinkbonds}: {(st.cum_nxlinkbonds>=st.max_nxlinkbonds)}')
        finished=self.search_failed or (st.cum_nxlinkbonds>=st.desired_nxlinkbonds)
        if finished:
            st.step=cure_step.cap_bondsearch
        return finished

    def _curr_conversion(self):
        """Returns current conversion expressed as a fraction of the maximum number of possible bonds.

        Returns:
            float: conversion
        """
        if not self.state.max_nxlinkbonds: return 0
        return float(self.state.cum_nxlinkbonds)/self.state.max_nxlinkbonds

    def conversion(self):
        """Returns the current bond conversion: the fraction of all possible crosslink bonds that have been formed.

        This is the number the cure iterates against and reports.  For a
        chemistry whose crosslinkers only count as reacted once every one of
        their sites is filled, it is an upper bound on the conversion an
        experiment would measure; see :mod:`htpolynet.repair.cyanate_cap`.

        Returns:
            float: bond conversion
        """
        return self._curr_conversion()

    def reset(self):
        """Resets this CureController object by resetting its CureState, its bonds dataframe, and its search_failed member.
        """
        self.state.reset()
        self.bonds_df=pd.DataFrame()
        self.bonds_are='nonexistent!'
        self.search_failed=False

    def next_iter(self):
        """Resets this CureController for the next CURE iteration.

        Returns:
            bool: True if new iteration index exceeds the maximum number of iterations allowed by runtime
        """
        i=self.state.iter+1
        self.reset()
        self.state.iter=i
        return self.state.iter>=self.dicts['controls']['max_iterations']

    def _do_bondsearch(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,reentry=False):
        """Manages the search for new bonds.

        Args:
            TC (TopoCoord): TopoCoord object containing all topology and coordinate information
            RL (ReactionList): list of all reaction types
            MD (MoleculeDict): dictionary of all molecular templates, keyed by name
            reentry (bool): True if we are restarting a bondsearch in this directory, defaults to False
        """
        if self.state.step!=cure_step.cure_bondsearch: return
        opfx=self._pfx()
        d=self.dicts['controls']
        TC.min_bondcycle_length=d['min_allowable_bondcycle_length']
        self.state.current_radius=d['search_radius']+self.state.current_radidx*d['radial_increment']
        logger.info(f'Bond search using radius {self.state.current_radius} nm initiated')
        curr_conversion=self._curr_conversion()
        apply_probabilities=curr_conversion<d['late_threshold']
        bond_limit=int(d['max_conversion_per_iteration']*self.state.max_nxlinkbonds)
        bond_target=int((d['desired_conversion']-curr_conversion)*self.state.max_nxlinkbonds)
        bond_limit=min([bond_limit,bond_target])
        logger.debug(f'Iteration limited to at most {bond_limit} new bonds')
        # Minimum bonds to find before we stop growing the search radius.
        # Clamped against bond_target (don't demand more than physically
        # possible) and bond_limit (don't demand more than we'd accept).
        # If we hit max radius with fewer than this many but at least one,
        # we still proceed -- the post-loop branch handles that.
        min_floor=max(1,min(int(d.get('min_bonds_per_iteration',1)),bond_target,bond_limit))
        logger.debug(f'Iteration will grow radius until >= {min_floor} bonds (or max radius)')
        nbdf=pd.DataFrame()
        nbonds=0
        while nbonds<min_floor and self.state.current_radidx<self.state.max_radidx:
            nbdf=self._searchbonds(TC,RL,MD,stage=reaction_stage.cure,abs_max=bond_limit,apply_probabilities=apply_probabilities,reentry=reentry)
            nbonds=nbdf.shape[0]
            logger.debug(f'Search at radius {self.state.current_radius:.3f} nm found {nbonds} bond(s) (floor {min_floor})')
            if nbonds<min_floor and self.state.current_radidx<self.state.max_radidx:
                self.state.current_radidx+=1
                self.state.current_radius+=d['radial_increment']
                logger.info(f'Radius increased to {self.state.current_radius} nm ({nbonds}/{min_floor} eligible bonds so far)')
        if nbonds>0:
            ess='' if nbonds==1 else 's'
            logger.info(f'Iteration {self.state.iter} will generate {nbdf.shape[0]} new bond{ess}')
            pairs=pd.DataFrame() # empty placeholder
            ''' compute all lengths '''
            TC.add_length_attribute(nbdf,attr_name='initial_distance')
            ''' register all new bonds '''
            self._register_bonds(nbdf,pairs,f'{opfx}-bonds.csv',bonds_are='identified')
            ''' declare intention to go to dragging or topology-update '''
            self.state.step=cure_step.cure_drag if self.dragging_enabled else cure_step.cure_update
        else:
            self.search_failed=True
            ''' declare intention to proceed to bond capping '''
            self.state.step=cure_step.cap_bondsearch
        self.state._to_yaml()
        logger.debug(f'next: {self.state.step}')

    def _do_preupdate_dragging(self,TC:TopoCoord,gromacs_dict={}):
        """Manages the preupdate dragging CURE step.

        Args:
            TC (TopoCoord): the global system topology and coordinates
            gromacs_dict (dict): dictionary of optional gromacs directives, defaults to {}
        """
        if self.state.step!=cure_step.cure_drag: return
        nbdf=self.bonds_df
        assert nbdf.shape[0]>0
        d=self.dicts['drag']
        nogos=[not self.dragging_enabled,nbdf['initial_distance'].max()<d['trigger_distance']]
        logger.debug(f'{nogos} {any(nogos)}')
        if any(nogos):
            logger.debug(f'No dragging is necessary')
        else:
            self._distance_attenuation(TC,mode='drag',gromacs_dict=gromacs_dict)
        self.state.step=cure_step.cure_update
        self.state._to_yaml()
        logger.debug(f'next: {self.state.step}')

    def _do_relax(self,TC:TopoCoord,gromacs_dict={}):
        """Manages relaxation steps within the CURE algorithm.

        Args:
            TC (TopoCoord): global topology and coordinates
            gromacs_dict (dict): dictionary of optional gromacs directives, defaults to {}
        """
        if self.state.step!=cure_step.cure_relax and self.state.step!=cure_step.cap_relax: return
        nbdf=self.bonds_df
        if nbdf.shape[0]==0: # no bonds identified
            if self.state.step==cure_step.cure_relax:
                self.state.step=cure_step.cap_bondsearch # drop out of CURE loop
            else:
                self.state.step=cure_step.finished
            return
        self._distance_attenuation(TC,mode='relax',gromacs_dict=gromacs_dict)
        if self.state.step==cure_step.cure_relax:
            self.state.step=cure_step.cure_equilibrate
        else:
            self.state.step=cure_step.cap_equilibrate
        self.state._to_yaml()
        logger.debug(f'next: {self.state.step}')

    def _do_topology_update(self,TC:TopoCoord,MD:MoleculeDict):
        """Manages the topology update in CURE.

        Args:
            TC (TopoCoord): global system topology and coordinates
            MD (MoleculeDict): dictionary of all molecular templates
        """
        if self.state.step!=cure_step.cure_update and self.state.step!=cure_step.cap_update: return
        opfx=self._pfx()
        logger.debug(f'Topology update')
        bonds_df,pairs_df=TC.update_topology_and_coordinates(self.bonds_df,template_dict=MD,write_mapper_to=f'{opfx}-idx-mapper.csv',chain_manager=self.chain_manager)
        TC.add_length_attribute(bonds_df,attr_name='initial_distance')
        TC.add_length_attribute(pairs_df,attr_name='initial_distance')
        self._register_bonds(bonds_df,pairs_df,f'{opfx}-bonds.csv',bonds_are='unrelaxed')
        TC.write_gro(f'{opfx}.gro')
        TC.write_top(f'{opfx}.top')
        TC.write_tpx(f'{opfx}.tpx')
        TC.write_grx_attributes(f'{opfx}.grx')
        if self.state.step==cure_step.cure_update:
            self.state.step=cure_step.cure_relax
        else:
            self.state.step=cure_step.cap_relax
        self.state._to_yaml()

    def _do_equilibrate(self,TC:TopoCoord,gromacs_dict={}):
        """Manages pre- and post-cure equilibration runs.

        Args:
            TC (TopoCoord): global system topology and coordinates
            gromacs_dict (dict): dictionary of optional gromacs directives, defaults to {}
        """
        if self.state.step!=cure_step.cure_equilibrate and self.state.step!=cure_step.cap_equilibrate: return
        d=self.dicts['equilibrate']
        opfx=self._pfx()
        edr_list=TC.equilibrate(deffnm=f'{opfx}',edict=d,gromacs_dict=gromacs_dict) #,plot_pfx=f'iter-{self.state.iter}-{str(self.state.step)}')
        ens=d['ensemble']
        if ens=='npt':
            plot_pfx=f'iter-{self.state.iter}-{str(self.state.step)}'
            trace('Density',edr_list,outfile=os.path.join(pfs.proj(),f'plots/{plot_pfx}-density.png'),yunits=r'kg/m$^3$')
        if self.state.step==cure_step.cure_equilibrate:
            # go to next iteration -- this whole method is skipped if nbonds==0 in relax
            self.state.step=cure_step.cure_bondsearch # if not self.search_failed else state.cap_bondsearch
        elif self.state.step==cure_step.cap_equilibrate:
            self.state.step=cure_step.finished
        self.state._to_yaml()

    def _do_cap_bondsearch(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict):
        """Manages the post-cure capping reaction bond search.

        Args:
            TC (TopoCoord): global system topology and coordinates
            RL (ReactionList): list of all reactions
            MD (MoleculeDict): dictionary of all molecular templates
        """
        if self.state.step!=cure_step.cap_bondsearch: return
        opfx=self._pfx()
        # multi: nbdf=self.make_cadidates(...,stage='cap')
        nbdf=self._searchbonds(TC,RL,MD,stage=reaction_stage.cap)
        nbonds=nbdf.shape[0]
        ess='' if nbonds==1 else 's'
        logger.info(f'Capping will generate {nbdf.shape[0]} new bond{ess}')
        if nbonds>0:
            cwd=pfs.go_to(pfs.Dirs.systems_capping)
            pairs=pd.DataFrame() # empty placeholder
            TC.add_length_attribute(nbdf,attr_name='initial_distance')
            self._register_bonds(nbdf,pairs,f'{opfx}-bonds.csv',bonds_are='identified')
            self.state.step=cure_step.cap_update
        else:
            self.state.step=cure_step.finished
        self.state._to_yaml()

    def _write_bonds_df(self,bondsfile='bonds.csv'):
        """Writes all bonds in the bonds dataframe to a CSV file.

        Args:
            bondsfile (str): name of csv file to write, defaults to 'bonds.csv'
        """
        self.bonds_df.to_csv(bondsfile,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def _read_bonds_df(self,bonds_file_override=''):
        """Reads bonds into bonds dataframe from a CSV file.

        Args:
            bonds_file_override (str): name of file to read from (overrides default), defaults to ''
        """
        infile=self.dicts['output']['bonds_file'] if not bonds_file_override else bonds_file_override
        assert os.path.exists(infile),f'Error: {infile} not found'
        self.bonds_df=pd.read_csv(infile,sep=r'\s+',header=0)
        self.dicts['output']['bonds_file']=os.path.abspath(infile)

    def _register_bonds(self,bonds,pairs,bonds_file,bonds_are='unrelaxed'):
        """Registers the bonds and 1-4 pairs from the input parameters.

        Args:
            bonds (pd.DataFrame): dataframe of bonds
            pairs (pd.DataFrame): dataframe of 1-4 pairs
            bonds_file (str): file to which bonds are written
            bonds_are (str): a little flag to remind us the state of these bonds, defaults to 'unrelaxed'
        """
        self.bonds_df=bonds
        self.pairs_df=pairs
        self.bonds_are=bonds_are
        self.state.bonds_file=bonds_file
        self._write_bonds_df(bonds_file)
        self.dicts['output']['bonds_file']=os.path.abspath(bonds_file)

    def _pfx(self):
        """Generates a little string that can be used as a state-specific filename prefix.

        Returns:
            str: the little string
        """
        return f'{self.state.step.value}-{self.state.step}'

    def _reactive_positions(self,TC:TopoCoord):
        """Snapshots the positions of every atom that still has a reactive site.

        'z' is the count of a reactive atom's unused sites, so ``z > 0`` is the
        set of species still looking for a partner -- the same selector the bond
        search itself uses.  Bonds only form during the topology update, so this
        set is constant across one relax ladder and the snapshot can be compared
        atom for atom with a later one.

        Args:
            TC (TopoCoord): global topology and coordinates

        Returns:
            tuple: (globalIdx array, (N,3) position array), or None if the atom dataframe carries no 'z' attribute or nothing is still reactive
        """
        adf=TC.Coordinates.A
        if adf is None or not all(c in adf.columns for c in ('z','globalIdx','posX','posY','posZ')): return None
        sel=adf[adf['z']>0]
        if sel.shape[0]==0: return None
        return sel['globalIdx'].to_numpy(),sel[['posX','posY','posZ']].to_numpy(dtype=float)

    def _report_reactive_mobility(self,TC:TopoCoord,pre_reactive):
        """Logs Varshney's criterion across the relax window just completed.

        See :func:`varshney_criterion` for why this is worth reporting.  It is
        pure observation: any failure is logged and swallowed, because a
        diagnostic must not be able to fail a build.

        Args:
            TC (TopoCoord): global topology and coordinates, as left by the relax ladder
            pre_reactive (tuple): the snapshot returned by :meth:`_reactive_positions` before the ladder ran
        """
        try:
            idx,before=pre_reactive
            adf=TC.Coordinates.A.set_index('globalIdx')
            after=adf.loc[idx,['posX','posY','posZ']].to_numpy(dtype=float)
            r=minimum_image_displacements(before,after,TC.Coordinates.box)
            stats=varshney_criterion(r,self.dicts['controls']['search_radius'])
        except Exception as e:
            logger.debug(f'could not report reactive mobility: {e}')
            return
        if not stats['n']: return
        rcap=self.dicts['controls']['search_radius']
        logger.info(
            f'Reactive-species mobility over relax: rmsd {stats["rmsd"]:.3f} nm, '
            f'max {stats["max"]:.3f} nm, {stats["crossing_fraction"]*100:.0f}% of '
            f'{stats["n"]} crossed the {rcap:.2f} nm search radius'
        )
        if stats['crossing_fraction']<0.25:
            logger.warning(
                f'reactive-species mobility is low: only {stats["crossing_fraction"]*100:.0f}% of '
                f'still-reactive atoms moved a full search radius during relax. The '
                f'relax window is no longer satisfying the criterion it was sized '
                f'against, so later bonds are being chosen from a nearly frozen '
                f'neighborhood'
            )

    def _distance_attenuation(self,TC:TopoCoord,mode='drag',gromacs_dict={}):
        """Manages the progressive attenuation of bond parameters for new bonds.

        Args:
            TC (TopoCoord): global system topology and coordinates
            mode (str): string indicating whether this is prebonding or postbonding attenuation; defaults to 'drag' (prebonding); other option is 'relax' (postbonding)
            gromacs_dict (dict): dictionary of optional gromacs directives, defaults to {}
        """
        assert mode in ['drag','relax']
        statename=self.state.step.basename()
        opfx=self._pfx()
        nbdf=self.bonds_df.copy()
        pdf=self.pairs_df.copy()
        d=self.dicts[mode]
        nbdf['current_lengths']=nbdf['initial_distance'].copy()
        maxL,minL,meanL=nbdf['current_lengths'].max(),nbdf['current_lengths'].min(),nbdf['current_lengths'].mean()
        logger.debug(f'Lengths avg/min/max: {meanL:.3f}/{minL:.3f}/{maxL:.3f}')
        ess='' if nbdf.shape[0]==1 else 's'
        logger.info(f'Step "{str(self.state.step)}" initiated on {nbdf.shape[0]} distance{ess} (max {maxL:.3f} nm)')
        roptions=[self.dicts['gromacs']['rdefault'],maxL]
        if mode=='drag':
            TC.add_restraints(nbdf,typ=6)
            logger.info('     Stage  Max-distance (nm)')
        else:
            pdf['current_lengths']=pdf['initial_distance'].copy()
            pmaxL,pminL,pmeanL=pdf['current_lengths'].max(),pdf['current_lengths'].min(),pdf['current_lengths'].mean()
            logger.debug(f'1-4 distances lengths avg/min/max: {pmeanL:.3f}/{pminL:.3f}/{pmaxL:.3f}')
            roptions.append(pmaxL)
            logger.info('     Stage  Max-distance (nm)  Max-1-4-distance (nm)  Density (kg/m3)')
        rcommon=max(roptions)
        for stg_dict in d['equilibration']:
            ensemble=stg_dict['ensemble']
            impfx=f'{statename}-{ensemble}' # e.g., drag-min, drag-nvt, drag-npt
            pfs.checkout(pfs.Dirs.mdp_file(impfx))
            mdp_mods_dict={'rvdw':rcommon,'rcoulomb':rcommon,'rlist':rcommon}
            if ensemble!='min':
                mdp_mods_dict.update({'gen-temp':stg_dict['temperature'],'ref_t':stg_dict['temperature'],'gen-vel':'yes','nsteps':stg_dict['nsteps']})
            if ensemble=='npt':
                mdp_mods_dict['ref_p']=stg_dict['pressure']
            mdp_modify(f'{impfx}.mdp',mdp_mods_dict,add_if_missing=(ensemble!='min'))
        this_nstages=int(maxL/d['increment'])
        this_firststage=self.state.current_stage[mode]
        # Varshney's criterion is measured across the whole relax window, so a
        # build resuming mid-ladder cannot report it -- the window it would
        # measure is only the tail.  Skip rather than report a partial number.
        pre_reactive=None
        if mode=='relax':
            if this_firststage==0:
                pre_reactive=self._reactive_positions(TC)
            else:
                logger.debug(f'resuming at stage {this_firststage}; skipping the mobility report')
        logger.debug(f'{self.state.step} {this_nstages} {this_firststage}')
        saveT=TC.Topology.copy_bond_parameters(self.bonds_df)
        for i in range(this_firststage,this_nstages):
            self.state.current_stage[mode]=i
            if mode=='drag':
                TC.Topology.attenuate_bond_parameters(self.bonds_df,i,this_nstages,minimum_distance=d['limit'],init_colname='initial_distance')
            else:
                TC.Topology.attenuate_bond_parameters(self.bonds_df,i,this_nstages,init_colname='initial_distance')
            stagepfx=f'{opfx}-stage-{i+1}'
            TC.write_top(f'{stagepfx}.top')
            stage_density=None
            for stg_dict in d['equilibration']:
                ensemble=stg_dict['ensemble']
                impfx=f'{statename}-{ensemble}' # e.g., drag-min, drag-nvt, drag-npt
                TC.grompp_and_mdrun(out=f'{stagepfx}-{ensemble}',mdp=f'{impfx}',**gromacs_dict)
                # logger.debug(f'{TC.files["gro"]}')
                # the relax stages are the only above-Tg constant-pressure time in
                # a cure, and until now nothing looked at the density they produce.
                # Drag runs under restraints, so its density is not comparable.
                if mode=='relax' and ensemble=='npt':
                    stage_density=_relax_stage_density(f'{stagepfx}-npt')
            TC.Topology.restore_bond_parameters(saveT)
            TC.add_length_attribute(nbdf,attr_name='current_lengths')
            maxL,minL,meanL=nbdf['current_lengths'].max(),nbdf['current_lengths'].min(),nbdf['current_lengths'].mean()
            logger.debug(f'Distances avg/min/max: {meanL:.3f}/{minL:.3f}/{maxL:.3f}')
            if mode=='drag':
                logger.info(f'{i+1:>10d}  {maxL:>17.3f}')
            else:
                TC.add_length_attribute(pdf,attr_name='current_lengths')
                pmaxL,pminL,pmeanL=pdf['current_lengths'].max(),pdf['current_lengths'].min(),pdf['current_lengths'].mean()
                logger.debug(f'1-4 distances lengths avg/min/max: {pmeanL:.3f}/{pminL:.3f}/{pmaxL:.3f}')
                dstr='--' if stage_density is None else f'{stage_density:.1f}'
                logger.info(f'{i+1:>10d}  {maxL:>17.3f}  {pmaxL:>21.3f}  {dstr:>15}')
            self.state._to_yaml()
        if mode=='drag':
            TC.Topology.remove_restraints(self.bonds_df)
        if pre_reactive is not None:
            self._report_reactive_mobility(TC,pre_reactive)
        TC.write_top(f'{opfx}-complete.top')
        self._register_bonds(nbdf,pdf,f'{opfx}-{mode}-bonds.csv',bonds_are=('relaxed' if mode=='relax' else 'dragged'))
        self.state.current_stage[mode]=0
        self.state._to_yaml()

    def check_iterations_vs_functionality(self,TC:TopoCoord,RL:ReactionList=None,MD:MoleculeDict=None):
        """Warns if the cure produced few or no complete crosslinkers.

        A system whose crosslinkers did not complete has no junctions: it is a
        monomer melt that equilibrates, reports a sensible density, and looks
        in every other respect like a cured network.  Nothing else surfaces
        this.  The conversion the cure reports is a *bond* conversion and is
        perfectly happy with it; only a postcure repair stage would notice, and
        only if one is configured.

        There are two ways to arrive there and this checks for both.

        The first is a counting constraint.  The bond downselection admits at
        most one bond per residue per iteration, so a residue with ``f``
        reactive sites cannot be fully reacted in fewer than ``f`` iterations —
        whatever the bond conversion says — and a cure that reaches its target
        in two iterations leaves *every* trifunctional crosslinker incomplete.

        The second is ordinary shortfall, and it is not bounded away by any
        iteration count.  The bonds a cure forms are not spread evenly across
        its iterations: the per-iteration count decays steeply, and it is the
        *last* iteration that has to supply a crosslinker's final bond.  A cure
        that runs exactly ``f`` iterations and spends its last one on almost
        nothing ends up in the same place as one that ran too few, and the
        counting test above passes it in silence.  This is not hypothetical —
        a build has reached exactly one complete crosslinker in 240 that way.
        So the completed count itself is checked, which is exact by the time
        this runs and needs no threshold on a proxy for it.

        ``completion_bias`` is not the remedy for either — it changes which
        residues react, not the one-per-residue-per-iteration rule — which
        matters because it is the first thing a user would reach for on seeing
        a low crosslinker conversion.

        Args:
            TC (TopoCoord): global system topology and coordinates
            RL (ReactionList): runtime reaction list; when given with MD, functionality and completion are counted on cure-reactive sites only (see :func:`cure_site_mask`)
            MD (MoleculeDict): molecule templates, used with RL
        """
        adf=TC.gro_DataFrame('atoms')
        if adf is None or not all(c in adf.columns for c in ('resNum','resName','z','nreactions')): return
        mask=cure_site_mask(adf,RL,MD)
        sites=adf if mask is None else adf[mask]
        if sites.empty: return
        func=residue_functionality(sites)
        if func.empty: return
        fmax=int(func.max())
        # A crosslink junction needs a residue with at least three reactive
        # sites.  When the most functional residue has two, every chain is linear
        # (chain-growth polystyrene, say): a fully reacted difunctional unit is a
        # chain interior, not a junction, and nothing can percolate however far
        # the cure goes.  Both checks below would then misfire -- the gel-point
        # threshold (1/(f-1))**f is exactly 100% at f=2, so it fired unless every
        # monomer was mid-chain.
        if fmax<3:
            if mask is not None:
                logger.info('No residue carries three or more cure-reactive sites, so crosslink percolation is not assessed. That is correct for a linear polymer; a system that crosslinks through a multi-residue molecule -- a dimethacrylate, a multi-arm prepolymer -- is outside what this residue-level check can see.')
            return
        names=adf.loc[adf['resNum']==func.idxmax(),'resName']
        resname=names.iloc[0] if len(names) else f'residue {func.idxmax()}'
        iterations=self.state.iter
        if iterations<fmax:
            logger.warning(f'The cure ran {iterations} iteration(s), but {resname} carries {fmax} reactive sites and the bond search admits at most one bond per residue per iteration. No {resname} can have reacted at all {fmax} of its sites, whatever bond conversion was reached, so this system contains no fully-reacted {resname} and therefore no crosslink junctions. completion_bias does not change this. To spend more iterations reaching the same conversion, lower CURE.controls.max_conversion_per_iteration or min_bonds_per_iteration.')
            return
        crosslinkers=func[func==fmax].index
        if not len(crosslinkers): return
        # completion counts cure-site bonds only, or a bond made at the param
        # stage (HIE's C4) would make a half-reacted residue look complete
        counts=residue_reaction_counts(sites).reindex(crosslinkers).fillna(0)
        complete=int((counts>=fmax).sum())
        fraction=complete/len(crosslinkers)
        logger.info(f'{complete} of {len(crosslinkers)} {resname} are fully reacted ({fraction:.1%})')
        # The crosslinker conversion an ideal randomly-branching network of this
        # functionality would show at its gel point, f_c**f with f_c=1/(f-1).
        # An idealization, used only to decide when to raise the volume: below
        # it a run has too few complete crosslinkers for a percolating network
        # to be assumed, whatever the bond conversion says.
        gel_point=(1.0/(fmax-1))**fmax
        if fraction<gel_point:
            logger.warning(f'The cure ran {iterations} iteration(s) and formed enough bonds to report its bond conversion, but only {complete} of {len(crosslinkers)} {resname} ({fraction:.1%}) ended up reacted at all {fmax} of their sites, and only those are crosslink junctions. That is below the {gel_point:.1%} an ideal randomly-branching network of functionality {fmax} would reach at its gel point, so do not assume this system percolates. The usual cause is a cure whose last iteration formed very few bonds, since that is the iteration that has to supply each crosslinker its final bond; the per-iteration bond count decays steeply, so reaching the target conversion is not enough on its own. completion_bias does not change this. To spread the same conversion over more iterations, lower CURE.controls.max_conversion_per_iteration or min_bonds_per_iteration.')

    def _completion_bias_counts(self,adf):
        """Returns per-residue reaction counts when the completion bias is enabled, else None.

        Args:
            adf (pandas.DataFrame): the global atom dataframe

        Returns:
            pandas.Series: per-residue bond counts, or None to rank by distance alone
        """
        if not self.dicts['controls'].get('completion_bias',False):
            return None
        if 'nreactions' not in adf.columns:
            logger.warning('completion_bias is enabled but the atom dataframe carries no "nreactions" attribute; ranking bond candidates by distance alone')
            return None
        return residue_reaction_counts(adf)

    def _check_bias_side(self,adf,bdf,RL=None,MD=None):
        """Warns, once per run, if the completion bias is ranking the less functional side.

        The bias ranks on the ``B`` reactant because htpolynet's A2+B3 idiom
        puts the crosslinker there.  Declare the crosslinker as ``A`` instead
        and the bias silently starts completing the *bridges*, which is a
        different claim about which species is a reactive intermediate and
        almost certainly not the one that was wanted.  It cannot be caught by
        looking for an all-zero bias key -- a difunctional bridge accumulates
        reactions perfectly well, so the key is non-zero and the ranking looks
        healthy.  What distinguishes the two cases is which side is the more
        functional one, and that is fixed at setup.

        Args:
            adf (pandas.DataFrame): the global atom dataframe
            bdf (pandas.DataFrame): the bond candidates for this iteration
        """
        if self.bias_side_checked: return
        if not all(c in adf.columns for c in ('z','nreactions')): return
        self.bias_side_checked=True
        mask=cure_site_mask(adf,RL,MD)
        func=residue_functionality(adf if mask is None else adf[mask])
        a=bdf['ri'].map(func).dropna()
        b=bdf['rj'].map(func).dropna()
        if a.empty or b.empty: return
        if a.max()>b.max():
            logger.warning(f'completion_bias ranks candidates on the B-side residue, but the A side is the more functional one here ({int(a.max())} reactive sites vs {int(b.max())}). The bias will complete B-side residues, which is probably not what was intended; declare the crosslinker as the second reactant to bias toward it.')

    # TODO: move to searchbonds.by; split into make_candidates() and apply_filters()
    def _searchbonds(self,TC:TopoCoord,RL:ReactionList,MD:MoleculeDict,stage=reaction_stage.cure,abs_max=0,apply_probabilities=False,reentry=False):
        """Manages the search for bonds.

        Args:
            TC (TopoCoord): global system topology and coordinates
            RL (ReactionList): list of all reactions
            MD (MoleculeDict): dictionary of all molecular templates
            stage (reaction_stage): which reaction stage, defaults to reaction_stage.cure, other choice is reaction_stage.cap
            abs_max (int): upper limit to number of new bonds allowed, defaults to 0, signifying no limit
            apply_probabilities (bool): if true, bond probabilities are applied, defaults to False
            reentry (bool): if true, existing linkcell data stored on disk is used, defaults to False

        Returns:
            pd.DataFrame: the dataframe of proposed bonds
        """
        adf=TC.gro_DataFrame('atoms')
        gro=TC.files['gro']
        ncpu=self.dicts['controls']['ncpu']
        if stage==reaction_stage.cure:
            TC.linkcell_initialize(self.state.current_radius,force_repopulate=reentry)
        raset=adf[adf['z']>0]  # this view will be used for downselecting to potential A-B partners
        bdf=pd.DataFrame()
        Rlist=[x for x in RL if (x.stage==stage and x.probability>0.0)]
        logger.debug(f'reactioncount {len(Rlist)} atomscount {raset.shape[0]}')
        for R in Rlist:
            logger.debug(f'Reaction {R.name} with {len(R.bonds)} bond(s)')
            prob=R.probability
            for bond in R.bonds:
                A=R.atoms[bond['atoms'][0]]
                B=R.atoms[bond['atoms'][1]]
                order=bond['order']
                aname=A['atom']
                areactantname_template=R.reactants[A['reactant']]
                aresid_template=A['resid']
                aresname=MD[areactantname_template].get_resname(aresid_template)
                az=A['z']
                bname=B['atom']
                breactantname_template=R.reactants[B['reactant']]
                bresid_template=B['resid']
                bresname=MD[breactantname_template].get_resname(bresid_template)
                bz=B['z']
                if stage==reaction_stage.cap:
                    assert areactantname_template==breactantname_template,f'Error: capping reaction {R.name} lists a bond whose atoms are in different reactants'
                    assert aresname==bresname,f'Error: capping reaction {R.name} lists a bond whose atoms are in different residues'

                Aset=raset[(raset['atomName']==aname)&(raset['resName']==aresname)&(raset['z']==az)&(raset['reactantName']==areactantname_template)]
                Bset=raset[(raset['atomName']==bname)&(raset['resName']==bresname)&(raset['z']==bz)&(raset['reactantName']==breactantname_template)]
                logger.debug(f'Aset {Aset.shape[0]} atoms')
                logger.debug(f'Bset {Bset.shape[0]} atoms')
                alist=list(zip(Aset['globalIdx'].to_list(),Aset['resNum'].to_list(),Aset['molecule'].to_list()))
                blist=list(zip(Bset['globalIdx'].to_list(),Bset['resNum'].to_list(),Bset['molecule'].to_list()))
                all_possible_pairs=list(product(alist,blist))

                idf=pd.DataFrame({'ai':           [int(x[0][0]) for x in all_possible_pairs],
                                  'ri':           [int(x[0][1]) for x in all_possible_pairs],
                                  'mi':           [int(x[0][2]) for x in all_possible_pairs],
                                  'aj':           [int(x[1][0]) for x in all_possible_pairs],
                                  'rj':           [int(x[1][1]) for x in all_possible_pairs],
                                  'mj':           [int(x[1][2]) for x in all_possible_pairs],
                                  'prob':         [prob for _ in all_possible_pairs],
                                  'reactantName': [R.product for _ in  all_possible_pairs],
                                  'order':        [order for _ in all_possible_pairs]})
                if stage==reaction_stage.cure:
                    # exclude atom pairs that have same resid or molid
                    idf=idf[(idf['ri']!=idf['rj'])&(idf['mi']!=idf['mj'])].copy()
                    logger.debug(f'Examining {idf.shape[0]} bond-candidates of order {order}')
                    if idf.shape[0]>0:
                        ess='' if idf.shape[0]!=1 else 's'
                        bondtestoutcomes={k:0 for k in BTRC}
                        p=Pool(processes=ncpu)
                        idf_split=[idf.iloc[s] for s in np.array_split(np.arange(len(idf)),ncpu)]
                        packets=[(i,idf_split[i]) for i in range(ncpu)]
                        logger.debug(f'Decomposed dataframe lengths: {", ".join([str(x.shape[0]) for x in idf_split])}')
                        results=p.map(partial(gromacs_distance,gro=gro,new_column_name='r'),packets)
                        p.close()
                        p.join()
                        # reassemble final dataframe:
                        # logger.debug(f'Checking dataframe lengths: {", ".join([str(x.shape[0]) for x in results])}')
                        idf=pd.concat(results,ignore_index=True)
                        # gromacs_distance(idf,gro,new_column_name='r') # use "gmx distance" to very quickly get all lengths
                        logger.debug(f'{idf.shape[0]} bond-candidate length{ess} avg/min/max: {idf["r"].mean():0.3f}/{idf["r"].min():0.3f}/{idf["r"].max():0.3f}')
                        idf=idf[idf['r']<self.state.current_radius].copy().reset_index(drop=True)
                        ess='' if idf.shape[0]!=1 else 's'
                        logger.debug(f'{idf.shape[0]} bond-candidate{ess} with lengths below {self.state.current_radius} nm')
                        p=Pool(processes=ncpu)
                        idf_split=[idf.iloc[s] for s in np.array_split(np.arange(len(idf)),ncpu)]
                        # logger.debug(f'Decomposed dataframe lengths: {", ".join([str(x.shape[0]) for x in idf_split])}')
                        results=p.map(partial(TC.bondtest_df),idf_split)
                        p.close()
                        p.join()
                        # reassemble final dataframe:
                        # logger.debug(f'Checking dataframe lengths: {", ".join([str(x.shape[0]) for x in results])}')
                        idf=pd.concat(results,ignore_index=True)
                        if not idf.empty:
                            logger.debug(f'Bond-candidate test outcomes:')
                            for k in bondtestoutcomes:
                                bondtestoutcomes[k]=idf[idf['result']==k].shape[0]
                                logger.debug(f'   {str(k)}: {bondtestoutcomes[k]}')
                            idf=idf[idf['result']==BTRC.passed].copy()
                elif stage==reaction_stage.cap:
                    idf=idf[idf['ri']==idf['rj']].copy()
                    logger.debug(f'Examining {idf.shape[0]} bond-candidates of order {order}')
                if not idf.empty:
                    bdf=pd.concat((bdf,idf),ignore_index=True)
        # logger.debug('Filtered (by single-bond tests) bonds:')
        # for ln in bdf.to_string().split('\n'):
        #     logger.debug(ln)

        if stage==reaction_stage.cure and bdf.shape[0]>0:
            reaction_counts=self._completion_bias_counts(adf)
            bdf=rank_bond_candidates(bdf,reaction_counts)
            if reaction_counts is not None:
                self._check_bias_side(adf,bdf,RL=RL,MD=MD)
                nrx=bdf['rj'].map(reaction_counts).fillna(0).astype(int)
                logger.debug(f'Completion bias: {int((nrx>0).sum())} of {bdf.shape[0]} candidates extend an already-reacted B-side residue (most bonds already carried: {int(nrx.max())})')
            bdf['allowed']=[True for x in range(bdf.shape[0])]
            unique_atomidx=set(bdf.ai.to_list()).union(set(bdf.aj.to_list()))
            unique_resids=set(bdf.ri.to_list()).union(set(bdf.rj.to_list()))
            for i,r in bdf.iterrows():
                if r.ai in unique_atomidx:
                    unique_atomidx.remove(r.ai)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated atom index {r.ai}')
                    bdf.loc[i,'allowed']=False
                if r.aj in unique_atomidx:
                    unique_atomidx.remove(r.aj)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated atom index {r.aj}')
                    bdf.loc[i,'allowed']=False
                if r.ri in unique_resids:
                    unique_resids.remove(r.ri)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated residue index {r.ri}')
                    bdf.loc[i,'allowed']=False
                if r.rj in unique_resids:
                    unique_resids.remove(r.rj)
                else:
                    # logger.debug(f'Disallowing bond {i} due to repeated residue index {r.rj}')
                    bdf.loc[i,'allowed']=False

            logger.debug(f'{bdf[bdf["allowed"]==False].shape[0]} out of {bdf.shape[0]} bonds disallowed due to repeated atom indexes or residue indexes')

            bdf=bdf[bdf['allowed']==True].copy().reset_index(drop=True)
            # logger.debug('Allowed bonds:')
            # for ln in bdf.to_string().split('\n'):
            #     logger.debug(ln)

            bdf=TC.bondcycle_collective(bdf,chain_manager=self.chain_manager)
            logger.debug(f'{bdf[bdf["remove-to-uncyclize"]==True].shape[0]} out of {bdf.shape[0]} bonds removed to break nascent cycles')
            bdf=bdf[bdf['remove-to-uncyclize']==False].copy().reset_index(drop=True)
            # logger.debug('Non-cyclizing bonds:')
            # for ln in bdf.to_string().split('\n'):
            #     logger.debug(ln)

            ''' roll the dice '''
            bdf['lucky']=[True for _ in range(bdf.shape[0])]
            if apply_probabilities:
                for i,r in bdf.iterrows():
                    x=np.random.random()
                    if x>r.prob:
                        bdf.loc[i,'lucky']=False

            logger.debug(f'{bdf[bdf["lucky"]==True].shape[0]} bonds survive probability application')
            bdf=bdf[bdf['lucky']==True].copy().reset_index(drop=True)
            # logger.debug('Lucky bonds:')
            # for ln in bdf.to_string().split('\n'):
            #     logger.debug(ln)

            ''' apply the stated limit '''
            if abs_max>-1:
                if abs_max<bdf.shape[0]:
                    bdf=bdf.loc[:abs_max].copy().reset_index(drop=True)
                    logger.debug(f'Limiting to {bdf.shape[0]} allowed bonds')
            logger.debug('Final bonds:')
            for ln in bdf.to_string().split('\n'):
                logger.debug(ln)

        if stage==reaction_stage.cure:
            TC.linkcell_cleanup()

        return bdf