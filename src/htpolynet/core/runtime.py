"""Handles the HTPolyNet runtime workflow.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging
import os
import random
import shutil

import yaml

from collections import namedtuple
from copy import deepcopy

import networkx as nx
import numpy as np

from networkx.algorithms import isomorphism

from ..analysis.plot import trace
from ..core import paramcache
from ..core import projectfilesystem as pfs
from ..core.configuration import Configuration
from ..core.molecule import Molecule, MoleculeDict
from ..core.topocoord import TopoCoord
from ..core.topology import conflicting_types, select_topology_type_option
from ..cure.chain import ChainManager
from ..cure.curecontroller import CureController, CureState
from ..cure.expandreactions import bondchain_expand_reactions, generate_stereo_reactions, generate_symmetry_reactions, sibling_expand_reactions
from ..cure.reaction import (Reaction, ReactionList, parse_reaction_list, extract_molecule_reactions,
                             is_reactant, is_ring_closing, reaction_stage)
from ..core.productsplice import nearer_site_atoms
from ..cure.ringcontroller import RingController, RingCureState
from ..cure.triplesearch import reactive_sites, triple_bonds_dataframe
from ..external import software as software
from ..external.gromacs import insert_molecules, mdp_modify, mdp_get
from ..external.smiles_input import materialize_smiles_inputs
from ..utils import profiling
from ..utils import checkpoint as cp
from ..utils.stringthings import my_logger

logger = logging.getLogger(__name__)

def logrotate(filename):
    """Renames any existing log file with name 'filename' using a pattern that prevents overwriting any old logs.

    Args:
        filename (str): name of log file
    """
    if os.path.exists(filename):
        n = 1
        while os.path.exists(f'#{n}#{filename}'):
            n += 1
        shutil.copyfile(filename, f'#{n}#{filename}')

_directives=['ps','nsteps','ncycles']
def _nonempty_directives(dirlist):
    """Checks the list of directives for existence of at least one of the strings in the _directives list.

    Args:
        dirlist (list): list of dictionaries

    Returns:
        bool: True if one of the strings in _directives is present in one of the dictionaries
    """
    amts=[]
    for d in dirlist:
        ta=[d.get(x,0) for x in _directives]
        amts.append(any(ta))
    return any(amts)

class Runtime:
    """ Runtime class manages all aspects of a system build using HTPolyNet.
    """
    default_edict={ 'ensemble': 'npt', 'temperature': 300, 'pressure': 10, 'ps': 200, 'nsteps': -2, 'repeat': 0 }
    runtime_defaults={
        'gromacs': {
            'gmx': 'gmx',
            'gmx_options': '-quiet -nobackup',
            'mdrun': 'gmx mdrun'
        },
        'ambertools': {
            'charge_method': 'gas'
        },
        'densification': {
            'initial_density': 200.0,  # kg/m3
            'aspect_ratio': np.array([1.0,1.0,1.0]),
            'equilibration': [
                { 'ensemble': 'min' },
                { 'ensemble': 'nvt', 'temperature': 300, 'ps': 10 },
                { 'ensemble': 'npt', 'temperature': 300, 'pressure': 10, 'ps': 200 }
            ]
        },
        'precure': {  
            'preequilibration': {
                'ensemble': 'npt',
                'temperature': 300,        # K
                'pressure': 1,             # bar
                'ps': 100  # we should probably bring to a desired pressure after densification
            },
            'anneal': {
                'ncycles': 0,  # default none
                'initial_temperature': 300,
                'cycle_segments': [
                    { 'T': 300, 'ps': 0 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 300, 'ps': 20 },
                    { 'T': 300, 'ps': 20 }
                ]
            },
            'postequilibration': {
                'ensemble': 'npt',
                'temperature': 300,        # K
                'pressure': 1,             # bar
                'ps': 0 # default none
            }
        },
        'CURE': {},
        'postcure': {
            'anneal': {
                'ncycles': 0, # default none
                'initial_temperature': 300,
                'cycle_segments': [
                    { 'T': 300, 'ps': 0 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 300, 'ps': 20 },
                    { 'T': 300, 'ps': 20 }
                ]
            },
            'postequilibration': {
                'ensemble': 'npt',
                'temperature': 300,       # K
                'pressure': 1,            # bar
                'ps':  0 # default none
            }
        }
    }

    def __init__(self,cfgfile='',restart=False):
        """Generates a new Runtime object based on the cfg file.

        Args:
            cfgfile (str): name of YAML cfg file, defaults to ''
            restart (bool): flag indicating if this is a restart
        """
        self.cfgfile=cfgfile
        if cfgfile=='':
            logger.error('HTPolyNet requires a configuration file.')
            raise RuntimeError('HTPolyNet requires a configuration file')
        logger.info(f'Configuration: {cfgfile}')
        self.cfg=Configuration.read(os.path.join(pfs.root(),cfgfile))
        self._apply_runtime_defaults()
        software.set_gmx_preferences(self.cfg.gromacs)
        self.TopoCoord=TopoCoord(system_name='htpolynet')
        self.restart=restart
        if restart:
            logger.info(f'Restarting in {pfs.proj()}')
        self.frcmods=self._read_frcmods()
        self.molecules:MoleculeDict={}
        # molecules reused from a library entry that carries no provenance
        # record, collected for a stage-end summary; see generate_molecules
        self.unverified_parameterizations=[]
        self.reactions:ReactionList=[]
        self.molecule_report={}
        self.maxconv=0.0
        self.ncpu=self.cfg.ncpu
        if self.cfg.cure:
            logger.debug('Setting up cure controller')
            self.cc=CureController(self.cfg.cure)
        # Write generated mol2 into the user library (rootPath/lib/...) so that
        # pfs.checkout() finds them later when molecule.generate() runs from
        # inside the project directory.  Defaulting to a relative path here
        # would resolve against projPath instead, leaving the files unfindable.
        if pfs._PFS_ and pfs._PFS_.userlibrary:
            inputs_dir = str(pfs._PFS_.userlibrary.root / pfs.Dirs.molecules_inputs)
        else:
            inputs_dir = os.path.join(pfs.root(), 'lib', pfs.Dirs.molecules_inputs)
        generated=materialize_smiles_inputs(self.cfg.constituents, inputs_dir=inputs_dir)
        if generated:
            logger.info(f'Generated mol2 inputs from SMILES for: {", ".join(generated)}')
        self._build_molecules_and_reactions()

    def _apply_runtime_defaults(self):
        """Merges runtime_defaults into the cfg sub-objects for any missing keys."""
        mapping=[
            ('gromacs',     'gromacs'),
            ('ambertools',  'ambertools'),
            ('densification','densification'),
            ('precure',     'precure'),
            ('CURE',        'cure'),
            ('postcure',    'postcure'),
        ]
        for cfg_key, attr in mapping:
            defaults=self.runtime_defaults.get(cfg_key,{})
            current=getattr(self.cfg,attr)
            if not current:
                setattr(self.cfg,attr,defaults.copy() if isinstance(defaults,dict) else defaults)
            elif isinstance(defaults,dict):
                for k,v in defaults.items():
                    if k not in current:
                        current[k]=v

    def _report_parameterization_summary(self):
        """Logs how many molecule templates were reused and how many were parameterized.

        Counted after every derived template (symmetry, C=C chain, and
        second-shell) exists, so the numbers cover everything this run needs;
        a job that must not parameterize anything, e.g. one reading a
        pre-filled read-only cache, can check for "0 freshly parameterized".
        """
        n_cached=sum(1 for M in self.molecules.values() if M.origin=='previously parameterized')
        n_new=sum(1 for M in self.molecules.values() if (M.origin or '').startswith('newly parameterized'))
        if not (n_cached or n_new): return
        where=[]
        if pfs._PFS_ and pfs._PFS_.userlibrary:
            where.append(f'{pfs._PFS_.userlibrary.root}')
        if pfs._USER_CACHE_ is not None:
            where.append(f'{pfs._USER_CACHE_.root}')
        logger.info(f'Parameterization summary: {n_cached} of {n_cached+n_new} molecule templates reused from the library '
                    f'({" or ".join(where) or "unknown location"}), {n_new} freshly parameterized.')
        if n_cached:
            logger.info('  to force re-parameterization on the next run, use '
                        '`htpolynet run --force-parameterization --force-checkin ...`')

    def _read_frcmods(self):
        """Reads the user frcmod file each constituent names, if any.

        Returns:
            dict: {constituent name: frcmod text}

        Raises:
            FileNotFoundError: if a named frcmod file does not exist
        """
        frcmods={}
        for name,rec in (self.cfg.constituents or {}).items():
            fn=(rec or {}).get('frcmod')
            if not fn: continue
            path=os.path.join(pfs.root(),fn)
            if not os.path.isfile(path):
                raise FileNotFoundError(f'frcmod {fn!r} for constituent {name} not found (looked for {path})')
            with open(path) as f:
                frcmods[name]=f.read()
            logger.info(f'{name}: user frcmod {fn} will be loaded for {name} and every template containing it')
        return frcmods

    def _ambertools_for(self,M:Molecule):
        """Returns the AmberTools directives for parameterizing M.

        These are the configuration's, plus the user frcmod files of every
        constituent residue M contains, so that a template built around a
        monomer carries the same override as the monomer.

        Args:
            M (Molecule): molecule about to be parameterized, or checked against the cache

        Returns:
            dict: ambertools directives
        """
        amb=dict(self.cfg.ambertools or {})
        residues=set(getattr(M,'sequence',None) or [M.name])
        frc=[(n,t) for n,t in sorted(getattr(self,'frcmods',{}).items()) if n in residues]
        if frc:
            amb['frcmod']=frc
        return amb

    def _check_frcmod_type_conflicts(self):
        """Refuses a set of molecules that a user frcmod leaves with conflicting bonded types.

        Raises:
            RuntimeError: if a type is parameterized differently in molecules built with and
                without (or with different) user frcmod files
        """
        if not getattr(self,'frcmods',None): return
        tops={n:(M.TopoCoord.Topology,paramcache.frcmod_digest(self._ambertools_for(M).get('frcmod')))
              for n,M in self.molecules.items() if M.TopoCoord.Topology.D.get('atoms') is not None}
        conflicts=conflicting_types(tops)
        if not conflicts: return
        for typename,key,a,b in conflicts[:10]:
            logger.error(f'{typename} {"-".join(key[:-1])} (func {key[-1]}) differs between {a} and {b}')
        raise RuntimeError(f'{len(conflicts)} bonded type(s) are parameterized differently in molecules built with and '
                           f'without a user frcmod.  A system has one table per type, so the override cannot apply to '
                           f'some molecules and not others.  Give the overridden atoms their own atom type, or apply '
                           f'the frcmod to every constituent that contains the type.')

    def _build_molecules_and_reactions(self):
        """Builds Molecule and Reaction objects from the parsed configuration.

        Populates self.molecules, self.reactions, and self.molecule_report.
        """
        base_reactions=[Reaction(r) for r in self.cfg.reaction_specs]
        self.reactions=parse_reaction_list(base_reactions)

        mol_reac_detected=extract_molecule_reactions(self.reactions)
        self.molecule_report['explicit']=len(mol_reac_detected)
        self.molecule_report['implied by stereochemistry']=0
        self.molecule_report['implied by symmetry']=0

        for mname,gen in mol_reac_detected:
            self.molecules[mname]=Molecule.New(mname,gen,self.cfg.constituents.get(mname,{}))
            for si,S in self.molecules[mname].stereoisomers.items():
                self.molecules[si]=S
                self.molecule_report['implied by stereochemistry']+=1

        for mname,M in self.molecules.items():
            M.set_sequence_from_moldict(self.molecules)

        self.molecule_report['implied by stereochemistry']+=generate_stereo_reactions(self.reactions,self.molecules)
        self.molecule_report['implied by symmetry']+=generate_symmetry_reactions(self.reactions,self.molecules)

        for R in self.reactions:
            for rnum,rname in R.reactants.items():
                zrecs=[]
                for atnum,atrec in R.atoms.items():
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        zrecs.append(cprec)
                self.molecules[rname].update_zrecs(zrecs,self.molecules)

        for molec in self.cfg.initial_composition:
            mname=molec['molecule']
            if mname not in self.molecules:
                self.molecules[mname]=Molecule.New(mname,None,self.cfg.constituents.get(mname,{}))
                self.molecules[mname].set_sequence_from_moldict(self.molecules)

    def _calculate_maximum_conversion(self):
        """Calculates and stores the maximum number of polymerization bonds that can form."""
        Atom=namedtuple('Atom',['name','resid','reactantKey','reactantName','z'])
        Bond=namedtuple('Bond',['ai','aj'])
        N={}
        for item in self.cfg.initial_composition:
            molecule_name=item['molecule']
            molecule_count=item.get('count',0)
            if molecule_count:
                for res in self.molecules[molecule_name].sequence:
                    if res not in N:
                        N[res]=0
                    N[res]+=molecule_count
        Bonds=[]
        Atoms=[]
        for R in [x for x in self.reactions if x.stage==reaction_stage.cure]:
            for b in R.bonds:
                A,B=b['atoms']
                a,b=R.atoms[A],R.atoms[B]
                aan,ban=a['atom'],b['atom']
                ari,bri=a['resid'],b['resid']
                arnum,brnum=a['reactant'],b['reactant']
                arn,brn=R.reactants[arnum],R.reactants[brnum]
                if arnum==brnum: continue
                az,bz=a['z'],b['z']
                ia=Atom(aan,ari,arnum,arn,az)
                ib=Atom(ban,bri,brnum,brn,bz)
                b=Bond(ia,ib)
                if ia not in Atoms and arn in N: Atoms.append(ia)
                if ib not in Atoms and brn in N: Atoms.append(ib)
                if b not in Bonds and arn in N and brn in N: Bonds.append(b)
        Z=[a.z*N[a.reactantName] for a in Atoms]
        MaxB=[]
        for B in Bonds:
            az=Z[Atoms.index(B.ai)]
            bz=Z[Atoms.index(B.aj)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(B.ai)]-=MaxB[-1]
            Z[Atoms.index(B.aj)]-=MaxB[-1]
        self.maxconv=sum(MaxB)

    def generate_molecules(self,force_parameterization=False,force_checkin=False):
        """Manages creation and parameterization of all monomers and oligomer templates.

        Args:
            force_parameterization (bool): forces AmberTools to run parameterizations, defaults to False
            force_checkin (bool): forces HTPolyNet to overwrite molecule library, defaults to False
        """
        for mname,M in self.molecules.items():
            M.origin='unparameterized'
        self.unverified_parameterizations=[]
        pfs.go_to(pfs.Dirs.molecules_parameterized)
        my_logger(f'Templates in {pfs.cwd()}',logger.info)
        ess='' if len(self.molecules)==1 else 's'
        logger.info(f'{len(self.molecules)} molecule{ess} detected in {self.cfgfile}')
        replns=[f'{msg:>30s}: {count:<5d}' for msg,count in self.molecule_report.items()]
        for ln in replns:
            logger.info(ln)
        logger.debug(f'Generating: {list(self.molecules.keys())}')
        for mname,M in self.molecules.items():
            self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
        new_reactions,new_molecules=bondchain_expand_reactions(self.molecules)
        if len(new_molecules)>0:
            ess='' if len(new_molecules)==1 else 's'
            logger.info(f'{len(new_molecules)} molecule{ess} implied by C-C bondchains')
            self.reactions.extend(new_reactions)
            make_molecules={k:v for k,v in new_molecules.items() if k not in self.molecules}
            for mname,M in make_molecules.items():
                self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
                assert M.origin!='unparameterized'
                self.molecules[mname]=M
                logger.debug(f'Generated {mname}')

        new_reactions,new_molecules=sibling_expand_reactions(self.molecules,self.reactions)
        if len(new_molecules)>0:
            ess='' if len(new_molecules)==1 else 's'
            logger.info(f'{len(new_molecules)} molecule{ess} implied by cure bonds next to already-reacted atoms')
            self.reactions.extend(new_reactions)
            for mname,M in new_molecules.items():
                self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
                assert M.origin!='unparameterized'
                self.molecules[mname]=M
                logger.debug(f'Generated {mname}')

        self._report_parameterization_summary()
        self._check_frcmod_type_conflicts()

        for M in self.molecules:
            self.molecules[M].is_reactant=is_reactant(M,self.reactions,stage=reaction_stage.cure)

        resolve_type_discrepancies=self.cfg.gaff.get('resolve_type_discrepancies',[]) or self.cfg.resolve_type_discrepancies
        if resolve_type_discrepancies:
            for resolve_directive in resolve_type_discrepancies:
                logger.info(f'Resolving type discrepancies for directive {resolve_directive}')
                self._type_consistency_check(typename=resolve_directive['typename'],funcidx=resolve_directive.get('funcidx',4),selection_rule=resolve_directive['rule'])

        self._report_unverified_parameterizations()

        ess='' if len(self.molecules)==1 else 's'
        logger.info(f'Generated {len(self.molecules)} molecule template{ess}')
        if self.cfg.initial_composition:
            cmp_msg=", ".join([(x["molecule"]+" "+str(x["count"])) for x in self.cfg.initial_composition if x["count"]>0])
            logger.info(f'Initial composition is {cmp_msg}')
            self._calculate_maximum_conversion()
            logger.info(f'100% conversion is {self.maxconv} bonds')

        logger.debug(f'Reaction bond(s) in each molecular template:')
        for M in self.molecules.values():
            if len(M.reaction_bonds)>0:
                logger.debug(f'Bond {M.name}:')
                for b in M.reaction_bonds:
                    logger.debug(f'   {str(b)}')

        logger.debug(f'Bond template(s) in each molecular template:')
        for M in self.molecules.values():
            if len(M.bond_templates)>0:
                logger.debug(f'Template {M.name}:')
                for b in M.bond_templates:
                    logger.debug(f'   {str(b)}')

        for k,v in self.cfg.constituents.items():
            relaxdict=v.get('relax',{})
            if relaxdict:
                self.molecules[k].relax(relaxdict)

    @cp.enableCheckpoint
    def do_initialization(self,inpfnm='init'):
        """Manages creation of the global system topology file and the execution of 'gmx insert-molecules' to create the initial simulation box.

        Args:
            inpfnm (str): basename for output files, defaults to 'init'

        Returns:
            dict: dictionary of Gromacs files
        """
        my_logger(f'Initialization in {pfs.cwd()}',logger.info)
        top,tpx=self._initialize_topology(inpfnm)
        gro,grx=self._initialize_coordinates(inpfnm)
        return {'top':top,'tpx':tpx,'gro':gro,'grx':grx}

    @cp.enableCheckpoint
    def do_densification(self,deffnm='densified'):
        """Manages execution of gmx mdrun to perform minimization and NPT MD simulation of the initial liquid system. Final coordinates are loaded into the global TopoCoord.

        Args:
            deffnm (str): deffnm prefix fed to gmx mdrun, defaults to 'densified'

        Returns:
            dict: dictionary of Gromacs file names
        """
        my_logger(f'Densification in {pfs.cwd()}',logger.info)
        densification_dict=self.cfg.densification
        assert len(densification_dict)>0,'"densification" directives missing'
        equilibration=densification_dict.get('equilibration',[])
        assert len(equilibration)>0,'equilibration directives missing'
        TC=self.TopoCoord
        infiles=[TC.files[x] for x in ['gro','top','tpx','grx']]
        assert all([os.path.exists(x) for x in infiles]),f'One or more of {infiles} not found'
        self._do_equilibration_series(equilibration,deffnm=f'{deffnm}',plot_pfx='densification')
        logger.info(f'Densified coordinates in {pfs.cwd()}/{os.path.basename(TC.files["gro"])}')
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    @cp.enableCheckpoint
    def do_precure(self):
        """Manages execution of mdrun to perform any pre-cure MD simulation(s).

        Returns:
            dict: dictionary of Gromacs file names
        """
        self._do_pap(pfx='precure')
        return {c:os.path.basename(x) for c,x in self.TopoCoord.files.items() if c!='mol2'}

    def do_cure(self):
        """Manages CURE algorithm execution.

        Returns:
            None: only returns None if CURE cannot execute
        """
        if not hasattr(self,'cc'): 
            logger.debug(f'no cure controller')
            return  # no cure controller
        # A ring-closing reaction is the ring cure's business, and the pairwise search
        # skips it.  With nothing else to form, CURE would iterate to its limit forming
        # nothing, so say so and leave.  (The CURE block is present even when the
        # configuration omits it: the schema supplies its defaults.)
        pairwise=[R for R in self.reactions
                  if R.stage==reaction_stage.cure and not is_ring_closing(R)]
        if not pairwise:
            if self._ring_reactions():
                my_logger('No pairwise cure reactions; crosslinking is the ring cure\'s',logger.info)
            return
        cc=self.cc
        TC=self.TopoCoord
        RL=self.reactions
        MD=self.molecules
        gromacs_dict=self.cfg.gromacs
        pfs.go_proj()
        ''' read cure state or initialize new cure '''
        cc.chain_manager=self.chain_manager
        if os.path.exists(f'{pfs.Dirs.systems}/cure_state.yaml'):
            cc.state=CureState.from_yaml(f'{pfs.Dirs.systems}/cure_state.yaml')
            my_logger(f'Connect-Update-Relax-Equilibrate (CURE) resumes',logger.info)
            my_logger(f'at iteration {cc.state.iter} and {cc.state.cum_nxlinkbonds} bonds',logger.info)
        else:
            cc.setup(max_nxlinkbonds=self.maxconv,desired_nxlinkbonds=int(self.maxconv*cc.dicts['controls']['desired_conversion']),max_search_radius=float(min(TC.Coordinates.box.diagonal()/2)))
            cc.state.iter=1
            my_logger('Connect-Update-Relax-Equilibrate (CURE) begins',logger.info)
        cure_finished=cc.is_cured()
        if cure_finished: 
            logger.debug('cure finished even before loop')
            return
        ''' perform CURE iterations '''
        logger.info(f'Attempting to form {cc.state.desired_nxlinkbonds} bonds')
        while not cure_finished:
            pfs.go_to(pfs.Dirs.systems_iter(cc.state.iter))
            with profiling.stage(f'iter-{cc.state.iter}'):
                cc.do_iter(TC,RL,MD,gromacs_dict=gromacs_dict)
            cure_finished=cc.is_cured()
            if not cure_finished:
                cure_finished=cc.next_iter()
        cc.check_iterations_vs_functionality(TC,RL,MD)
        ''' perform capping if necessary '''
        my_logger(f'Capping begins',logger.info)
        pfs.go_to(pfs.Dirs.systems_capping)
        with profiling.stage('capping'):
            cc.do_capping(TC,RL,MD,gromacs_dict=gromacs_dict)
        self._report_molecule_charges('after cure')
        my_logger('Connect-Update-Relax-Equilibrate (CURE) ends',logger.info)

    def _ring_reactions(self):
        """The cure reactions that close a ring, which CURE's pairwise search skips."""
        return [R for R in self.reactions
                if R.stage == reaction_stage.cure and is_ring_closing(R)]

    def _ring_recipe(self, R):
        """Works out what the ring search needs from one ring-closing reaction.

        Returns:
            tuple: (residue name, ((donor, acceptor), ...), template Molecule,
                [template resid per reactant])
        """
        names = set(R.reactants.values())
        if len(names) != 1:
            raise NotImplementedError(f'reaction "{R.name}" closes a ring among different '
                                      f'reactants ({sorted(names)}); only one kind is supported')
        reactant = next(iter(names))
        if len(self.molecules[reactant].sequence) != 1:
            raise NotImplementedError(f'reaction "{R.name}" closes a ring among multi-residue '
                                      f'reactants; only single-residue ones are supported')
        resname = self.molecules[reactant].get_resname(1)
        donors = {b['atoms'][0] for b in R.bonds}
        pairs = set()
        for key in R.reactants:
            d = [R.atoms[k]['atom'] for k in R.atoms if R.atoms[k]['reactant'] == key and k in donors]
            a = [R.atoms[k]['atom'] for k in R.atoms if R.atoms[k]['reactant'] == key and k not in donors]
            if len(d) != 1 or len(a) != 1:
                raise NotImplementedError(f'reaction "{R.name}" gives reactant {key} '
                                          f'{len(d)} donor(s) and {len(a)} acceptor(s); expected one each')
            pairs.add((d[0], a[0]))
        template = self.molecules[R.product]
        resids = list(range(1, len(R.reactants) + 1))
        if len(pairs) != 1:
            raise NotImplementedError(f'reaction "{R.name}" uses {len(pairs)} different sites; '
                                      f'declare one and let symmetry supply the rest')
        home = next(iter(pairs))
        groups, names = self._symmetry_images(self.molecules[reactant], home)
        # A monomer with more than one reactive site can join more than one ring, and
        # this template shows one site reacted.  Splice only that site's share of each
        # residue, or a second ring would reset the first one's atoms.
        others = [g for g in groups if g != home]
        arm = None
        if others:
            arm = nearer_site_atoms(template, resids[0], home, [n for g in others for n in g])
            logger.info(f'{reactant} carries {len(groups)} reactive sites; splicing {len(arm)} '
                        f'atom(s) per residue, the reacting one\'s share')
        return resname, groups, template, resids, names, arm

    @staticmethod
    def _symmetry_images(M, group):
        """Expands one reactive site over its equivalent sites, with the name translation each needs.

        A dicyanate carries two cyanate arms, and the trimer template is built with one
        of them reacted.  A residue ringing through the other is the same chemistry with
        the two arms exchanged, so one template serves both -- provided the splice is
        told how to read it.

        That exchange is a symmetry of the molecule itself, and it moves *every* atom of
        the arm: the oxygen, the methylenes between them, their hydrogens.  Declaring all
        of that by hand is impractical, and a missed atom maps a bond onto one the
        instance does not have.  So the translation is derived as an automorphism of the
        monomer's own bond graph -- one that carries the declared site onto the other --
        and `symmetry_equivalent_atoms` is consulted only for which other sites exist.

        Args:
            M (Molecule): the monomer carrying the sites
            group (tuple): (donor, acceptor) atom names, as the reaction declares them

        Returns:
            tuple: (groups, names) -- every (donor, acceptor) pair, and
                {group label: {template name: that group's name}}

        Raises:
            NotImplementedError: if a declared equivalent site is not reachable by any
                symmetry of the molecule
        """
        donor, acceptor = group
        sets = [list(s) for s in (M.symmetry_relateds or [])]
        home = next((s.index(donor) for s in sets if donor in s), None)
        if home is None:
            return (group,), {}
        width = min(len(s) for s in sets)
        A = M.TopoCoord.Coordinates.A
        name_of = dict(zip(A['globalIdx'].astype(int), A['atomName']))
        idx_of = {v: k for k, v in name_of.items()}
        G = nx.Graph()
        types = dict(zip(M.TopoCoord.Topology.D['atoms']['nr'].astype(int),
                         M.TopoCoord.Topology.D['atoms']['type']))
        for i, t in types.items():
            G.add_node(i, type=t)
        G.add_edges_from(zip(M.TopoCoord.Topology.D['bonds']['ai'].astype(int),
                             M.TopoCoord.Topology.D['bonds']['aj'].astype(int)))
        groups, names = [], {}
        for k in range(width):
            pair = (next(s[k] for s in sets if donor in s), next(s[k] for s in sets if acceptor in s))
            groups.append(pair)
            if k == home:
                continue
            want = {idx_of[donor]: idx_of[pair[0]], idx_of[acceptor]: idx_of[pair[1]]}
            matcher = isomorphism.GraphMatcher(G, G, node_match=lambda a, b: a['type'] == b['type'])
            found = next((m for m in matcher.isomorphisms_iter()
                          if all(m[a] == b for a, b in want.items())), None)
            if found is None:
                raise NotImplementedError(f'{M.name}: no symmetry of the molecule carries '
                                          f'{group} onto {pair}, so one template cannot serve both; '
                                          f'declare a reaction per site instead')
            names[f'{pair[0]}-{pair[1]}'] = {name_of[a]: name_of[b] for a, b in found.items()}
        return tuple(groups), names

    @cp.enableCheckpoint
    def do_ring_cure(self):
        """Manages a cure whose reaction closes a ring (cyclotrimerization).

        Returns:
            dict: dictionary of Gromacs file basenames, or empty if there is nothing to do
        """
        rxns = self._ring_reactions()
        if not rxns:
            return {}
        if len(rxns) > 1:
            raise NotImplementedError(f'{len(rxns)} ring-closing reactions; only one is supported')
        R = rxns[0]
        resname, groups, template, template_resids, name_translations, arm = self._ring_recipe(R)
        TC = self.TopoCoord
        gromacs_dict = self.cfg.gromacs
        statefile = f'{pfs.Dirs.systems}/ring_state.yaml'
        pfs.go_proj()
        state = RingCureState.from_yaml(statefile) if os.path.exists(statefile) else None
        rc = RingController(self.cfg.ring_cure, state=state)
        adf = TC.Coordinates.A
        total = reactive_sites(adf, resname, groups=groups, require_z=False).shape[0]
        rc.setup(total_groups=total, max_radius=float(min(TC.Coordinates.box.diagonal() / 2)))
        my_logger(f'Ring cure begins: {total} reactive {resname} group(s), '
                  f'{groups}, template {template.name}', logger.info)
        while not rc.is_cured():
            pfs.go_to(pfs.Dirs.systems_iter(rc.state.iter))
            TC.grab_files()
            adf = TC.Coordinates.A
            positions = {int(r.globalIdx): np.array([r.posX, r.posY, r.posZ])
                         for r in adf.itertuples()}
            sites, chosen = rc.search(adf, positions, TC.Coordinates.box.diagonal(), resname, groups)
            if not chosen:
                if rc.widen():
                    break
                rc.state.iter += 1
                continue
            bdf = triple_bonds_dataframe(chosen, sites, R.product, order=R.bonds[0].get('order', 1))
            work = rc.close_rings(TC, bdf, gromacs_dict=gromacs_dict)
            # a ring the ladder could not pull shut would be bonded long and stay that way
            bdf, chosen = rc.accept(bdf, work, chosen)
            if not chosen:
                rc.state.iter += 1
                continue
            rc.form_rings(TC, bdf, sites, chosen, template, template_resids,
                          name_translations=name_translations, atom_names=arm)
            rc.relax(TC, gromacs_dict=gromacs_dict)
            rc.record(chosen)
            TC.write_gro(f'ring-{rc.state.iter}.gro')
            TC.write_top(f'ring-{rc.state.iter}.top')
            TC.write_grx_attributes(f'ring-{rc.state.iter}.grx')
            TC.write_tpx(f'ring-{rc.state.iter}.tpx')
            rc.state.iter += 1
            pfs.go_proj()
            rc.state.to_yaml(statefile)
        self._report_molecule_charges('after ring cure')
        my_logger(f'Ring cure ends: {rc.state.rings} ring(s), conversion '
                  f'{rc.state.conversion:.3f}', logger.info)
        return {c: os.path.basename(x) for c, x in TC.files.items() if c != 'mol2'}

    @cp.enableCheckpoint
    def do_repair(self):
        """Apply postcure topology-repair operations between cure and postcure.

        Each spec in ``self.cfg.postcure_repair`` is dispatched by ``type:``
        to a driver in :mod:`htpolynet.repair`; drivers may sever bonds,
        delete atoms, transfer atoms between residues, and re-template the
        affected linkages — operations the monotonic cure/cap machinery
        cannot perform.  After all drivers run, the modified gro/top is
        written so subsequent stages (postcure MD) start from the repaired
        system.

        Returns:
            dict: dictionary of Gromacs file basenames (or empty if no
                repair specs were configured).
        """
        specs = getattr(self.cfg, 'postcure_repair', None) or []
        if not specs:
            return {}
        from ..repair import run_repair
        TC = self.TopoCoord
        my_logger(f'Postcure repair in {pfs.cwd()}', logger.info)
        n_repaired, repair_stats = run_repair(TC, self.molecules, specs, self.reactions)
        my_logger(f'Postcure repair performed {n_repaired} dismantle operations', logger.info)
        self._report_repair_conversion(repair_stats)
        self._report_molecule_charges('after repair')
        # Write out the repaired system so the relaxation can pick it up
        TC.write_grx_attributes('repaired.grx')
        TC.write_gro('repaired.gro')
        TC.write_top('repaired.top')
        TC.write_tpx('repaired.tpx')
        # Relocated -C#N groups can land close to neighbouring atoms;
        # a steepest-descent minimization plus a short, soft NVT settle
        # is needed before the postcure MD ensemble can take over.
        if n_repaired:
            logger.info('Relaxing repaired geometry')
            self._do_equilibration({'ensemble': 'min'}, deffnm='repair-min')
            self._do_equilibration({'ensemble': 'nvt', 'temperature': 300, 'ps': 5},
                                   deffnm='repair-nvt')
        return {c: os.path.basename(x) for c, x in TC.files.items() if c != 'mol2'}

    def _report_repair_conversion(self, repair_stats):
        """Reports the conversion a repaired system actually carries, next to the one the cure reported.

        The cure iterates on bond conversion: bonds formed over bonds
        possible.  Repair then dismantles every crosslinker that did not fill
        all of its sites, so the structure that leaves this stage contains
        only complete crosslinkers, and the fraction of those is what an
        experiment measures.  The two numbers differ -- for a trifunctional
        crosslinker under random placement the second is roughly the cube of
        the first -- and reporting only the cure's number invites a run at a
        bond conversion of 0.90 to be written up as a 0.90 when it is closer
        to 0.76.  Both numbers, and the counts behind them, are also written
        to ``repair-summary.yaml`` in this stage's directory.

        Args:
            repair_stats (list): statistics dicts returned by :func:`htpolynet.repair.run_repair`
        """
        if not any(st['n_crosslinkers'] for st in repair_stats):
            return
        bond_conversion = None
        cc = getattr(self, 'cc', None)
        if cc is not None:
            bond_conversion = cc.conversion()
        record = {'bond_conversion': bond_conversion, 'repairs': repair_stats}
        for st in repair_stats:
            if not st['n_crosslinkers']:
                continue
            msg = (f'Crosslinker conversion after repair: '
                   f'{st["crosslinker_conversion"]:.3f} '
                   f'({st["n_complete"]}/{st["n_crosslinkers"]} {st["residue"]} complete)')
            if bond_conversion is not None:
                msg += f'; the cure reported a bond conversion of {bond_conversion:.3f}'
            my_logger(msg, logger.info)
            # the transfer count is what predicts whether this stage survives,
            # so it belongs next to the conversion rather than mid-log
            if st.get('n_transferred'):
                tmsg = f'Repair transferred {st["n_transferred"]} cap fragments'
                if st.get('min_clearance_nm') is not None:
                    tmsg += f'; tightest placement {st["min_clearance_nm"]:.3f} nm'
                if st.get('n_below_target'):
                    tmsg += f'; {st["n_below_target"]} below the clearance target'
                my_logger(tmsg, logger.info)
        try:
            with open('repair-summary.yaml', 'w') as f:
                yaml.dump(record, f, default_flow_style=False)
        except Exception as e:
            logger.warning(f'Could not write repair-summary.yaml: {e}')

    @cp.enableCheckpoint
    def do_postcure(self):
        """Manages execution of mdrun to perform any post-cure MD simulation(s).

        Returns:
            dict: dictionary of Gromacs file names
        """
        self._do_pap(pfx='postcure')
        return {c:os.path.basename(x) for c,x in self.TopoCoord.files.items() if c!='mol2'}

    @cp.enableCheckpoint
    def save_data(self,result_name='final'):
        """Writes 'gro', 'top', 'tpx', and 'grx' files for system, plus a PSF
        + TCL pair for VMD visualization.

        Args:
            result_name (str): output file basename, defaults to 'final'

        Returns:
            dict: dictionary of Gromacs file basenames
        """
        TC=self.TopoCoord
        my_logger(f'Final data to {pfs.cwd()}',logger.info)
        dropped = TC.Topology.prune_stale_14_pairs()
        if dropped:
            logger.info(f'Pruned {dropped} stale [ pairs ] entries left over from cure-induced path shortening')
        self._report_molecule_charges(f'in {result_name}')
        TC.write_grx_attributes(f'{result_name}.grx')
        TC.write_gro(f'{result_name}.gro')
        TC.write_top(f'{result_name}.top')
        TC.write_tpx(f'{result_name}.tpx')
        self._write_vmd_viz_files(result_name)
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    def _report_molecule_charges(self, when, tol=0.01):
        """Logs whether every covalently bonded molecule in the system is neutral.

        The system total is zero by construction, so a charged molecule is
        otherwise invisible: a template or repair step that moves charge
        between molecules leaves no other trace.  A percolated network is one
        molecule and neutral by construction, so this cannot see charge
        misplaced inside it.

        Args:
            when (str): stage label for the message, e.g. 'after cure'
            tol (float): largest net charge, in e, treated as neutral, defaults to 0.01

        Returns:
            list: (number of atoms, net charge) of each molecule carrying more than tol
        """
        molecules = self.TopoCoord.Topology.molecule_charges()
        if not molecules:
            return []
        charged = [m for m in molecules if abs(m[1]) > tol]
        ess = '' if len(molecules) == 1 else 's'
        if not charged:
            logger.info(f'Molecule charges {when}: all {len(molecules)} molecule{ess} neutral to within {tol} e')
            return charged
        worst = max(charged, key=lambda m: abs(m[1]))
        logger.warning(f'Molecule charges {when}: {len(charged)} of {len(molecules)} molecule{ess} carry more than '
                       f'{tol} e of net charge; the largest is {worst[1]:+.3f} e on a {worst[0]}-atom molecule.  '
                       f'The system total is {self.TopoCoord.Topology.total_charge():+.4f} e, so this is charge '
                       f'moved between molecules, which parameterization does not intend.')
        return charged

    def _write_vmd_viz_files(self, result_name='final'):
        """Writes a PSF (real bond topology for VMD) and a TCL helper that
        trims bonds spanning periodic boundaries so VMD doesn't render the
        long "wrap-around" bonds inherent to a crosslinked network in PBC.
        Also emits a ``.viz.macros.tcl`` of constituent-keyed atomselect
        macros, sourced automatically from the main ``.viz.tcl``.

        Usage afterward:
            vmd final.viz.psf final.gro -e final.viz.tcl
        """
        from ..utils.vmd_viz import write_viz_files
        try:
            write_viz_files(f'{result_name}.top', f'{result_name}.gro',
                            prefix=result_name, grx=f'{result_name}.grx')
        except ImportError:
            logger.warning('parmed not importable; skipping VMD viz files')
        except Exception as e:
            logger.warning(f'Could not write {result_name}.viz.psf: {e}')

    def do_workflow(self,**kwargs):
        """do_workflow manages runtime for one entire system-build workflow
        """
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        TC=self.TopoCoord
        self.profile = profiling.RunProfile()
        profiling.set_active(self.profile)
        try:
            pfs.go_proj()
            with profiling.stage('setup'):
                self.generate_molecules(
                    force_parameterization=force_parameterization,  # force antechamber/GAFF parameterization
                    force_checkin=force_checkin                     # force check-in to system libraries
                )
            pfs.go_proj()
            last_data=cp.read_checkpoint()
            logger.debug(f'Checkpoint last_data {last_data}')
            TC.load_files(last_data)
            # On a restart, do_initialization is skipped by the checkpoint decorator,
            # so the in-memory chain_manager it would have built isn't there.
            # Rebuild it from the reloaded coordinates before any later stage tries
            # to consume it (do_cure at minimum).
            if last_data and getattr(TC, 'Coordinates', None) is not None and getattr(TC.Coordinates, 'A', None) is not None and not getattr(self, 'chain_manager', None):
                self.chain_manager = ChainManager()
                self.chain_manager.from_dataframe(TC.Coordinates.A)
                logger.info('Rebuilt chain_manager from checkpoint-loaded coordinates')
            pfs.go_to(pfs.Dirs.systems_init)
            with profiling.stage('initialization'):
                self.do_initialization()
            pfs.go_to(pfs.Dirs.systems_densification)
            with profiling.stage('densification'):
                self.do_densification()
            pfs.go_to(pfs.Dirs.systems_precure)
            with profiling.stage('precure'):
                self.do_precure()
            with profiling.stage('cure'):
                self.do_cure()
            if self._ring_reactions():
                with profiling.stage('ring-cure'):
                    self.do_ring_cure()
            if getattr(self.cfg, 'postcure_repair', None):
                pfs.go_to(pfs.Dirs.systems_repair)
                with profiling.stage('repair'):
                    self.do_repair()
            pfs.go_to(pfs.Dirs.systems_postcure)
            with profiling.stage('postcure'):
                self.do_postcure()
            pfs.go_to(pfs.Dirs.systems_final)
            with profiling.stage('final'):
                self.save_data()
            self._write_profile_report()
            pfs.go_proj()
        finally:
            profiling.set_active(None)

    def _write_profile_report(self):
        """Emit the run profile to the log and to ``proj-N/profile.json``."""
        p = getattr(self, 'profile', None)
        if p is None:
            return
        # Log raw (one line per logger.info call) so the table isn't wrapped
        # in my_logger's asterisk fill — that would crowd a multi-line report.
        for line in p.format_report().split('\n'):
            logger.info(line)
        try:
            json_path = os.path.join(pfs.proj(), 'profile.json')
            p.write_json(json_path)
            logger.info(f'Wrote {json_path}')
        except Exception as e:
            logger.warning(f'Could not write profile.json: {e}')

    def _report_unverified_parameterizations(self):
        """Summarizes, at the end of the parameterization stage, every molecule
        reused from a library entry whose provenance could not be checked.

        The per-molecule warning is right for grepping but wrong as the only
        surface: it sits in a log that runs to thousands of lines, alongside
        the ordinary "Using cached parameterization" line that this whole
        mechanism exists because people read straight past.  A block at the end
        of the stage is what someone scanning the tail of a log actually hits,
        and it arrives while the expensive part of the build is still ahead of
        them rather than behind.

        Both messages have to say that the cached parameters are being *used*,
        not merely that they could not be checked.  A reader who upgrades with
        an existing library, asks for one charge method, and gets another one's
        numbers will otherwise conclude the cache guard does not work -- the
        correct behavior here is indistinguishable from the bug it replaced
        unless the message says which is happening.
        """
        if not self.unverified_parameterizations:
            return
        names = sorted(set(self.unverified_parameterizations))
        ess = '' if len(names) == 1 else 's'
        requested = paramcache.build_key(self.cfg.ambertools)
        my_logger(f'{len(names)} parameterization{ess} reused without provenance', logger.warning)
        logger.warning('These were taken from the library, which records nothing about how they were')
        logger.warning('built. This build therefore carries whatever charges those entries hold, and')
        logger.warning(f'they may not be {requested["charge_method"]!r}:')
        my_logger(names, logger.warning)
        logger.warning('This is expected for entries written before htpolynet 2.3, and is not a')
        logger.warning('failure. Rebuild them with --force-parameterization if this build\'s charges')
        logger.warning('need to be certain.')

    @staticmethod
    def _checkin_parameterization(mname, overwrite):
        """Checks a molecule's parameterization files into the user library.

        The provenance record is checked in only when the data files beside it
        are.  pfs.checkin() will not replace a file that already exists unless
        told to, and a library entry can be partially populated, so a record
        checked in unconditionally could attach itself to files it does not
        describe -- leaving the library asserting that a gas parameterization
        was built with bcc, which is worse than the missing record it replaced.

        Args:
            mname (str): molecule name
            overwrite (bool): replace files already in the library
        """
        entry = f'{pfs.Dirs.molecules_parameterized}/{mname}'
        wrote_data = overwrite or not pfs.exists(f'{entry}.gro')
        for ex in ['mol2', 'top', 'tpx', 'itp', 'gro', 'grx']:
            pfs.checkin(f'{entry}.{ex}', overwrite=overwrite)
        if wrote_data:
            pfs.checkin(f'{entry}.{paramcache.CACHE_KEY_EXT}', overwrite=True)

    def _cached_parameterization_mismatch(self, M):
        """Returns descriptions of every way a cached parameterization of M
        disagrees with what this run asks AmberTools for.

        The cache is keyed on molecule name alone, so without this check a
        configuration requesting one charge method silently reuses a
        parameterization built with another, producing a system with mixed
        charge methods and no indication of it anywhere in the output.

        A library entry written before provenance records existed carries no
        record to check.  That is reported and treated as a match, so that an
        existing library keeps working; --force-parameterization remains the
        way to be certain of what a build used.

        Args:
            M (Molecule): the molecule whose cached parameterization is in question

        Returns:
            list: descriptions of the differing directives; empty if the cached
                parameterization agrees with this run or carries no record
        """
        requested = paramcache.build_key(self._ambertools_for(M))
        # Discard any record left in the project directory by an earlier run so
        # that what we check is the record belonging to the library entry that
        # previously_parameterized() just found.
        local = paramcache.key_filename(M.name)
        if os.path.exists(local):
            os.remove(local)
        keyfile = os.path.join(pfs.Dirs.molecules_parameterized, local)
        if pfs.exists(keyfile):
            pfs.checkout(keyfile)
        stored = paramcache.read_key(M.name)
        if stored is None:
            logger.warning(f'Reusing cached parameterization for {M.name}, which carries no record of '
                           f'how it was built: this build takes whatever charges that entry holds, and '
                           f'they may not be {requested["charge_method"]!r}. Expected for an entry '
                           f'written before htpolynet 2.3; rebuild it with --force-parameterization '
                           f'if this build\'s charges need to be certain.')
            self.unverified_parameterizations.append(M.name)
            return []
        diffs = paramcache.describe_mismatch(stored, requested)
        if diffs:
            logger.info(f'Cached parameterization for {M.name} was built differently than this run requests '
                        f'({"; ".join(diffs)}); re-parameterizing')
        return diffs

    def _generate_molecule(self,M:Molecule,**kwargs):
        """Generates and parameterizes a single molecule based on the partially complete instance in parameter M.

        Args:
            M (Molecule): partially-complete molecule instance (populated by config reader)
        """
        if M.origin!='unparameterized' and M.origin!='symmetry_product': return
        mname=M.name
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        # A cached parameterization is keyed on molecule name, which says nothing
        # about the charge method, net charge or atom types that produced its
        # numbers.  Reject one that answers a different question than this run is
        # asking; see _cached_parameterization_mismatch.
        use_cache = (not force_parameterization) and M.previously_parameterized()
        if use_cache and self._cached_parameterization_mismatch(M):
            use_cache = False
        # Either perform a fresh parameterization or fetch parameterization files
        if not use_cache:
            logger.debug(f'Parameterization of {mname} requested -- can we generate {mname}?')
            generatable=(not M.generator) or (all([m in self.molecules for m in M.generator.reactants.values()]))
            if generatable:
                logger.debug(f'Generating {mname}')
                M.generate(available_molecules=self.molecules,gaff=self.cfg.gaff,ambertools=self._ambertools_for(M))
                self._checkin_parameterization(mname, force_checkin)
                M.origin='newly parameterized'
            else:
                logger.debug(f'Error: could not generate {mname}')
                logger.debug(f'not ({mname}.generator) = {bool(not M.generator)}')
                if M.generator:
                    logger.debug(f'reactants {list(M.generator.reactants.values())}')
                return
        else:
            logger.info(f'Using cached parameterization for {mname}')
            exts=pfs.fetch_molecule_files(mname)
            logger.debug(f'fetched {mname} exts {exts}')
            mol2fn = f'{mname}.mol2' if 'mol2' in exts else ''
            M.load_top_gro(f'{mname}.top',f'{mname}.gro',tpxfilename=f'{mname}.tpx',mol2filename=mol2fn,wrap_coords=False)
            M.TopoCoord.read_gro_attributes(f'{mname}.grx')
            M.set_sequence_from_coordinates()
            M.chain_manager.from_dataframe(M.TopoCoord.Coordinates.A)
            stale_cache = False
            if M.generator:
                M.prepare_new_bonds(available_molecules=self.molecules)
                # Detect stale chain data on a cached build product.  When the
                # current run's monomer chain_managers are populated (after the
                # monomer fix above), every chain atom in each reactant should
                # show up in this build product's chain_manager — either as
                # part of an existing chain, or merged via the new bond into a
                # single longer chain.  If the cached product carries fewer
                # chain atoms than the sum across its reactants, the cache
                # predates the current run's reactivity metadata.  Catches
                # both the zero-chains case (homo-dimers from cache written
                # against a no-z monomer) and the short-chain case (hetero-
                # dimers whose chain_manager came out length-3 instead of
                # length-4 because the bond-side chain was never fully merged).
                expected_chain_atoms = sum(
                    sum(len(c.idx_list) for c in self.molecules[r].chain_manager.chains)
                    for r in M.sequence if r in self.molecules
                )
                actual_chain_atoms = sum(len(c.idx_list) for c in M.chain_manager.chains)
                if expected_chain_atoms > 0 and actual_chain_atoms < expected_chain_atoms:
                    logger.info(f'Cached chain data for {mname} appears stale '
                                f'({actual_chain_atoms} chain atom(s); expected {expected_chain_atoms}); re-parameterizing')
                    stale_cache = True
            else:
                # Reactivity-related grx attributes (z, sea_idx, bondchain) are
                # YAML-dependent — they reflect the *current* reaction set, not
                # something intrinsic to the molecule.  A cache entry written by
                # a YAML with no reactions (or different reactions) for this
                # monomer would carry stale z=0 / sea_idx=-1, silently producing
                # 0 reactive atoms in cure.  Re-derive from the current run.
                M.TopoCoord.set_gro_attribute('reactantName', mname)
                M.initialize_monomer_grx_attributes()
            if stale_cache:
                # Reset the molecule state populated by cache-load so the
                # subsequent M.generate() call starts from the same blank
                # slate the fresh-parameterize branch above would have given
                # it.  Otherwise the half-loaded TopoCoord / stale
                # chain_manager interact with merge() and prepare_new_bonds()
                # in ways that produce float-promoted globalIdx columns.
                from .topocoord import TopoCoord as _TopoCoord
                from ..cure.chain import ChainManager as _ChainManager
                M.TopoCoord = _TopoCoord()
                M.chain_manager = _ChainManager()
                M.bond_templates = []
                M.reaction_bonds = []
                # Re-derive sequence from generator (this is what
                # set_sequence_from_moldict normally does at __init__ time;
                # set_sequence_from_coordinates inside generate() will
                # assert it matches what falls out of the merged dataframe).
                M.set_sequence_from_moldict(self.molecules)
                M.generate(available_molecules=self.molecules,
                           gaff=self.cfg.gaff, ambertools=self._ambertools_for(M))
                self._checkin_parameterization(mname, True)
                M.origin='newly parameterized (stale cache refreshed)'
            else:
                M.origin='previously parameterized'

        ''' Generate any stereoisomers and/or conformers '''
        M.generate_stereoisomers()
        M.generate_conformers()

    def _type_consistency_check(self,typename='dihedraltypes',funcidx=4,selection_rule='stiffest'):
        """Checks current topology for any inconsistencies and corrects them.

        Args:
            typename (str): type of Gromacs interaction to check, defaults to 'dihedraltypes'
            funcidx (int): function type of interaction, defaults to 4
            selection_rule (str): one-word designation for correction mechanism, defaults to 'stiffest'
        """
        logger.debug(f'Consistency check of {typename} func {funcidx} on all {len(self.molecules)} molecules requested')
        mnames=list(self.molecules.keys())
        types_duplicated=[]
        for i in range(len(mnames)):
            logger.debug(f'{mnames[i]}...')
            for j in range(i+1,len(mnames)):
                mol1topo=self.molecules[mnames[i]].TopoCoord.Topology
                mol2topo=self.molecules[mnames[j]].TopoCoord.Topology
                this_duptypes=mol1topo.report_duplicate_types(mol2topo,typename=typename,funcidx=funcidx)
                logger.debug(f'...{mnames[j]} {len(this_duptypes)}')
                for d in this_duptypes:
                    if not d in types_duplicated:
                        types_duplicated.append(d)
        logger.debug(f'Duplicated {typename}: {types_duplicated}')
        options={}
        for t in types_duplicated:
            logger.debug(f'Duplicate {t}:')
            options[t]=[]
            for i in range(len(mnames)):
                moltopo=self.molecules[mnames[i]].TopoCoord.Topology
                this_type=moltopo.report_type(t,typename='dihedraltypes',funcidx=4)
                if len(this_type)>0 and not this_type in options[t]:
                    options[t].append(this_type)
                logger.debug(f'{mnames[i]} reports {this_type}')
            logger.debug(f'Conflicting options for this type: {options[t]}')
            selected_type=select_topology_type_option(options[t],typename,rule=selection_rule)
            logger.debug(f'Under selection rule "{selection_rule}", preferred type is {selected_type}')
            for i in range(len(mnames)):
                logger.debug(f'resetting {mnames[i]}')
                TC=self.molecules[mnames[i]].TopoCoord
                moltopo=TC.Topology
                moltopo.reset_type(typename,t,selected_type)

    def _initialize_topology(self,inpfnm='init'):
        """Creates a full gromacs topology that includes all directives necessary for an initial liquid simulation. This will NOT use any #include's; all types will be explicitly in-lined.

        Args:
            inpfnm (str): input file basename, defaults to 'init'

        Returns:
            str: name of top file
        """
        # if cp.passed('initialize_topology'): return
        if os.path.isfile(f'{inpfnm}.top'):
            logger.debug(f'{inpfnm}.top already exists in {pfs.cwd()} but we will rebuild it anyway!')

        ''' for each monomer named in the cfg, either parameterize it or fetch its parameterization '''
        TC=self.TopoCoord
        already_merged=[]
        for item in self.cfg.initial_composition:
            M=self.molecules[item['molecule']]

            N=item['count']
            if not N: continue
            t=deepcopy(M.TopoCoord.Topology)
            logger.debug(f'Merging {N} copies of {M.name}\'s topology into global topology')
            t.adjust_charges(atoms=t.D['atoms']['nr'].to_list(),desired_charge=0.0,overcharge_threshhold=0.1,msg='')
            t.rep_ex(N)
            TC.Topology.merge(t)
            already_merged.append(M.name)
        for othermol,M in self.molecules.items():
            if not othermol in already_merged:
                logger.debug(f'Merging types from {othermol}\'s topology into global topology')
                TC.Topology.merge_types(M.TopoCoord.Topology)
        logger.info(f'Topology "{inpfnm}.top" in {pfs.cwd()}')
        TC.write_top(f'{inpfnm}.top')
        TC.write_tpx(f'{inpfnm}.tpx')
        return f'{inpfnm}.top',f'{inpfnm}.tpx'

    def _initialize_coordinates(self,inpfnm='init'):
        """Builds initial top and gro files for initial liquid simulation.

        Args:
            inpfnm (str): input file basename, defaults to 'init'

        Returns:
            tuple: 'gro' and 'grx' file names
        """
        densification_dict=self.cfg.densification
        TC=self.TopoCoord
        dspec=['initial_density' in densification_dict,'initial_boxsize' in densification_dict]
        # logger.debug(f'{dspec} {any(dspec)} {not all(dspec)}')
        assert any(dspec),'Neither "initial_boxsize" nor "initial_density" are specfied'
        assert not all(dspec),'Cannot specify both "initial_boxsize" and "initial_density"'
        if 'initial_boxsize' in densification_dict:
            boxsize=densification_dict['initial_boxsize']
        else:
            mass_kg=TC.Topology.total_mass(units='SI')
            V0_m3=mass_kg/densification_dict['initial_density']
            ar=densification_dict.get('aspect_ratio',np.array([1.,1.,1.]))
            assert ar[0]==1.,f'Error: parameter aspect_ratio must be a 3-element-list with first element 1'
            ar_yx=ar[1]*ar[2]
            L0_m=(V0_m3/ar_yx)**(1./3.)
            L0_nm=L0_m*1.e9
            boxsize=np.array([L0_nm,L0_nm,L0_nm])*ar
            logger.info(f'Initial density: {densification_dict["initial_density"]} kg/m^3')
            logger.info(f'Total mass: {mass_kg:.3e} kg')
            logger.info(f'Box aspect ratio: {" x ".join([str(x) for x in ar])}')
        logger.info(f'Initial box side lengths: {boxsize[0]:.3f} nm x {boxsize[1]:.3f} nm x {boxsize[2]:.3f} nm')

        clist=[cc for cc in self.cfg.initial_composition if 'count' in cc]
        c_togromacs={}
        for cc in clist:
            M=self.molecules[cc['molecule']]
            tc=cc['count']
            if not tc: continue
            if len(M.conformers)>0:
                nconf=len(M.conformers)
                if tc<nconf:
                    c_togromacs.update({u:1 for u in random.sample(M.conformers,tc)})
                else:
                    q,r=divmod(tc,nconf)
                    logger.debug(f'divmod({nconf},{tc}) gives {q}, {r}')
                    c_togromacs.update({u:q for u in M.conformers})
                    k=0
                    while r:
                        c_togromacs[M.conformers[k]]+=1
                        k+=1
                        r-=1
                for gro in M.conformers:
                    pfs.checkout(f'{pfs.Dirs.molecules_parameterized}/{gro}.gro',altpath=[pfs.subpath(pfs.Dirs.molecules)])        
            else:
                ''' assuming racemic mixture of any stereoisomers '''
                total_isomers=len(M.stereoisomers)+1
                count_per_isomer=tc//total_isomers
                leftovers=tc%total_isomers
                c_togromacs[M.name]=count_per_isomer
                if leftovers>0:
                    c_togromacs[M.name]+=1
                    leftovers-=1
                pfs.checkout(f'{pfs.Dirs.molecules_parameterized}/{M.name}.gro',altpath=[pfs.subpath(pfs.Dirs.molecules)])
                for isomer in M.stereoisomers:
                    c_togromacs[isomer]=count_per_isomer
                    if leftovers>0:
                        c_togromacs[isomer]+=1
                        leftovers-=1
                    pfs.checkout(f'{pfs.Dirs.molecules_parameterized}/{isomer}.gro',altpath=[pfs.subpath(pfs.Dirs.molecules)])
        logger.debug(f'Sending to insert_molecules: {c_togromacs}')
        msg=insert_molecules(c_togromacs,boxsize,inpfnm,inputs_dir='../../molecules/parameterized',scale=self.cfg.densification.get('scale',0.4))
        TC.read_gro(f'{inpfnm}.gro')
        TC.atom_count()
        TC.set_grx_attributes()
        TC.inherit_grx_attributes_from_molecules(self.molecules,self.cfg.initial_composition)
        self.chain_manager=ChainManager()
        self.chain_manager.from_dataframe(TC.Coordinates.A)
        TC.Topology.make_resid_graph()
        TC.write_grx_attributes(f'{inpfnm}.grx')
        logger.info(f'Coordinates "{inpfnm}.gro" in {pfs.cwd()}')
        logger.info(f'Extended attributes "{inpfnm}.grx" in {pfs.cwd()}')
        return f'{inpfnm}.gro',f'{inpfnm}.grx'

    def _do_equilibration(self, edict={}, deffnm='equilibrate', plot_pfx=''):
        """Manages execution of TopoCoord.equilibrate using the directive passed in parameter 'edict'.

        Args:
            edict (dict): equilibration directives, defaults to {}
            deffnm (str): mdrun default file basename, defaults to 'equilibrate'
            plot_pfx (str): basename for any output plots, defaults to ''
        """
        gromacs_dict = self.cfg.gromacs
        TC = self.TopoCoord
        for k in edict.keys():
            if not k in self.default_edict:
                logger.debug(f'ignoring unknown equilibration directive "{k}"')
        edr_list = TC.equilibrate(deffnm, edict, gromacs_dict)
        ens = edict['ensemble']
        if ens == 'npt' and plot_pfx != '':
            trace('Density', edr_list, outfile=os.path.join(pfs.proj(), f'plots/{plot_pfx}-density.png'), yunits=r'kg/m$^3$')

    def _do_equilibration_series(self,eq_stgs=[],deffnm='equilibrate',plot_pfx=''):
        """Manages a series of equilibrations.

        Args:
            eq_stgs (list): list of stage designations, defaults to []
            deffnm (str): mdrun default file basename, defaults to 'equilibrate'
            plot_pfx (str): basename for any output plots, defaults to ''
        """
        if not eq_stgs: return
        for stg in eq_stgs:
            self._do_equilibration(stg,deffnm,plot_pfx)
                
    def _do_pap(self, pfx='precure'):
        """Manages a preanneal-anneal-postanneal series of MD simulations.

        Args:
            pfx (str): either 'precure' or 'postcure', defaults to 'precure'
        """
        if not pfx: return
        assert pfx in ['precure', 'postcure']
        pap_dict = getattr(self.cfg,pfx,{})
        if not pap_dict: return
        preequil = pap_dict.get('preequilibration', {})
        anneal = pap_dict.get('anneal', {})
        postequil = pap_dict.get('postequilibration', {})
        if not any([preequil, anneal, postequil]): return
        if not _nonempty_directives([x for x in [preequil, anneal, postequil] if x]): return
        my_logger(f'{pfx.capitalize()} in {pfs.cwd()}', logger.info)
        # Every branch needs the _nonempty_directives guard, not just anneal:
        # _apply_runtime_defaults injects default preequil/anneal/postequil
        # blocks for any cfg key the user omitted, and the default postequil
        # carries ps=0 (i.e. "no postequilibration").  Without this guard,
        # TC.equilibrate returns None for ps=0 and the trace() that follows
        # crashes on 'NoneType' object is not iterable.
        if preequil and _nonempty_directives([preequil]):
            self._do_equilibration(preequil, deffnm='preequilibration', plot_pfx=f'{pfx}-preequilibration')
        if anneal and _nonempty_directives([anneal]):
            self._do_anneal(anneal, deffnm='annealed')
            trace('Temperature', ['annealed'], outfile=os.path.join(pfs.proj(), f'plots/{pfx}-anneal-T.png'), yunits='K')
        if postequil and _nonempty_directives([postequil]):
            self._do_equilibration(postequil, deffnm='postequilibration', plot_pfx=f'{pfx}-postequilibration')

    def _do_anneal(self, anneal_dict={}, deffnm='anneal'):
        """Manages execution of an annealing MD simulation.

        Args:
            anneal_dict (dict): annealing directives, defaults to {}
            deffnm (str): mdrun default file basename, defaults to 'anneal'

        Raises:
            Exception: if no run duration is indicated in any of the annealing segments
        """
        if not anneal_dict: return
        ncycles = anneal_dict.get('ncycles', 0)
        if not ncycles: return
        mdp_pfx = 'nvt' # assume all annealing run at NVT
        TC = self.TopoCoord
        gromacs_dict = self.cfg.gromacs
        pfs.checkout(pfs.Dirs.mdp_file(mdp_pfx))
        timestep = float(mdp_get(f'{mdp_pfx}.mdp', 'dt'))
        cycle_segments = anneal_dict.get('cycle_segments', [])
        temps = [str(r['T']) for r in cycle_segments]
        durations = [r.get('ps', 0) for r in cycle_segments]
        if not any(durations):
            durations = ['{:.2f}'.format(r.get('nsteps', 0) * timestep) for r in cycle_segments]
            if not any(durations):
                raise Exception(f'No "ps" or "nsteps" found.')
        cycle_duration = sum(durations)
        total_duration = cycle_duration*ncycles
        nsteps = int(total_duration / timestep)
        cum_time = durations.copy()
        for i in range(1, len(cum_time)):
            cum_time[i] += cum_time[i-1]
        logger.info(f'Annealing: {len(durations)} points for {ncycles} cycles over {total_duration} ps')
        mod_dict={
            'ref_t': anneal_dict.get('initial_temperature', 300.0),
            'gen-temp': anneal_dict.get('initial_temperature', 300.0),
            'gen-vel': 'yes',
            'annealing-npoints': len(cycle_segments),
            'annealing-temp': ' '.join(temps),
            'annealing-time': ' '.join([f'{x:.2f}' for x in cum_time]),
            'annealing': 'periodic' if ncycles > 1 else 'single',
            'nsteps': nsteps
            }
        mdp_modify(f'{mdp_pfx}.mdp', mod_dict)
        t_lo = min(float(t) for t in temps)
        t_hi = max(float(t) for t in temps)
        logger.info(f'Running Gromacs: anneal; {total_duration:.2f} ps, {ncycles} cycle(s), {t_lo:.2f}-{t_hi:.2f} K')
        msg = TC.grompp_and_mdrun(out=deffnm, mdp=mdp_pfx, quiet=False, **gromacs_dict)
        logger.info(f'Annealed coordinates in {deffnm}.gro')

