"""Handles the postsim subcommand.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import json
import logging
import os
import shutil

import yaml

import numpy as np

from ..analysis.plot import scatter
from ..core import projectfilesystem as pfs
from ..core.configuration import Configuration
from ..core.topocoord import TopoCoord
from ..external import software as software
from ..external.gromacs import mdp_get, mdp_modify, gmx_energy_trace
from ..utils.logsetup import setup_logging

logger=logging.getLogger(__name__)

class PostSimMD:
    """ Generic class for handling post-cure md simulations; this one just does simple NPT MD equilibration;
    Classes that inherit from this class should define their own default_params and build_npt
    """
    default_params={
        'subdir': f'{pfs.Dirs.postsim}/equilibrate',
        'input_top': f'{pfs.Dirs.systems_final}/final.top',
        'input_gro': f'{pfs.Dirs.systems_final}/final.gro',
        'input_grx': f'{pfs.Dirs.systems_final}/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'ps': 1000,
        'T': 300,
        'P':1,
        'output_deffnm': 'equilibrate',
        'traces': ['Temperature','Density','Volume'],
        'scatter': ('time(ps)',['Density'],'rho_v_ns.png')
    }
    def __init__(self,indict,strict=True):
        self.params={}
        for p,v in self.default_params.items():
            self.params[p]=indict.get(p,v)
        for p,v in indict.items():
            if not p in self.default_params:
                if strict:
                    logger.info(f'Ignoring directive \'{p}\' in yaml input file')
                else:
                    self.params[p]=v
                    
    def do(self,mdp_pfx='npt',**gromacs_dict):
        """Handles executing the postsim MD simulation.

        Args:
            mdp_pfx (str): filename prefix for output files, defaults to 'npt'
        """
        p=self.params
        logger.info(f'do {p}')
        # if a gromacs dict is passed in, assume this overrides the one read in from the file
        if gromacs_dict:
            software.set_gmx_preferences(gromacs_dict)
        else:
            software.set_gmx_preferences(p['gromacs'])
        logger.info(f'going to {p["subdir"]}')
        pfs.go_to(p['subdir'])
        for input_file in ['input_top','input_gro','input_grx']:
            srcnm=os.path.join(pfs.proj(),p[input_file])
            shutil.copy(srcnm,'.')
        local_top=os.path.basename(p['input_top'])
        local_gro=os.path.basename(p['input_gro'])
        local_grx=os.path.basename(p['input_grx'])
        TC=TopoCoord(topfilename=local_top,grofilename=local_gro,grxfilename=local_grx)
        logger.info(f'{TC.Coordinates.A.shape[0]} atoms {TC.Topology.total_mass(units="gromacs"):.2f} amu')
        box=TC.Coordinates.box
        pfs.checkout(pfs.Dirs.mdp_file('npt'))
        os.rename('npt.mdp',f'{mdp_pfx}.mdp')
        self.build_mdp(f'{mdp_pfx}.mdp',box=box)
        msg=TC.grompp_and_mdrun(out=p['output_deffnm'],mdp=mdp_pfx,quiet=False,mylogger=logger.info,**gromacs_dict)
        df=gmx_energy_trace(p['output_deffnm'],p['traces'])
        use_scatter=p['scatter']
        temp_x=None
        temp_y=None
        scat_file=use_scatter[2]
        for sznm in ['Box-X','Box-Y','Box-Z']:
            if sznm in p['traces']:
                L0=box[0][0] if sznm=='Box-X' else box[1][1] if sznm=='Box-Y' else box[2][2]
                temp_x=f'{sznm}-strain'
                df[temp_x]=df[sznm]/L0-1.0
        # shear: the off-diagonal box element starts at zero, so the engineering shear
        # strain is the element itself over the length of the sheared face -- no -1
        for sznm,ref in [('Box-YX',1),('Box-ZX',2),('Box-ZY',2)]:
            if sznm in p['traces']:
                temp_x=f'{sznm}-strain'
                df[temp_x]=df[sznm]/box[ref][ref]
        for sznm in ['Pres-XX','Pres-YY','Pres-ZZ','Pres-XY','Pres-XZ','Pres-ZY']:
            if sznm in p['traces']:
                df[f'{sznm}-stress']=df[sznm]*(-1)
                temp_y=[f'{sznm}-stress']
        if temp_x and temp_y:
            use_scatter=(temp_x,temp_y,scat_file)
        df.to_csv(f'{p["output_deffnm"]}.csv',header=True,index=False)
        scatter(df,*use_scatter)
        logger.info(f'Final coordinates in {p["output_deffnm"]}.gro')
        logger.info(f'Traces saved in {p["output_deffnm"]}.csv')

    def build_mdp(self,mdpname,**kwargs):
        """Builds the GROMACS mdp file required for an NPT equilibration.

        Args:
            mdpname (str): name of mdp file
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        duration=params['ps']
        nsteps=int(duration/timestep)
        mod_dict={
            'ref_t':params['T'],
            'ref_p':params['P'],
            'gen-temp':params['T'],
            'gen-vel':'yes',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(mdpname,mod_dict)
        
class PostSimAnneal(PostSimMD):
    """ a class to handle temperature annealing MD simulation 
    """
    default_params={
        'subdir': f'{pfs.Dirs.postsim}/anneal',
        'input_top': f'{pfs.Dirs.systems_final}/final.top',
        'input_gro': f'{pfs.Dirs.systems_final}/final.gro',
        'input_grx': f'{pfs.Dirs.systems_final}/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun'
        },
        'output_deffnm':'anneal',
        'traces': ['Temperature','Density','Volume'],
        'scatter': ('time(ps)',['Density'],'rho_v_ns.png'),
        'T0': 300,
        'T1': 600,
        'ncycles': 1,
        'T0_to_T1_ps': 1000,
        'T1_ps': 1000,
        'T1_to_T0_ps': 1000,
        'T0_ps': 1000,
        'P':1
    }
    def build_mdp(self,mdpname,**kwargs):
        """Builds the GROMACS mdp file required for an annealing MD simulation.

        Args:
            mdpname (str): name of mdp file
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        timeints=[0.0,params['T0_to_T1_ps'],params['T1_ps'],params['T1_to_T0_ps'],params['T0_ps']]
        timepoints=[0.0]
        temppoints=[params['T0'],params['T1'],params['T1'],params['T0'],params['T0']]
        for i in range(1,len(timeints)):
            timepoints.append(timepoints[-1]+timeints[i])
        duration=timepoints[-1]
        nsteps=int(duration/timestep)*params['ncycles']
        mod_dict={
            'ref_t':params['T0'],
            'ref_p':params['P'],
            'gen-temp':params['T0'],
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'periodic' if params['ncycles']>1 else 'single',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(mdpname,mod_dict)

class PostSimLadder(PostSimMD):
    """ a class to handle a temperature-ladder MD simulation
    """
    default_params={
        'subdir': f'{pfs.Dirs.postsim}/ladder',
        'input_top': f'{pfs.Dirs.systems_final}/final.top',
        'input_gro': f'{pfs.Dirs.systems_final}/final.gro',
        'input_grx': f'{pfs.Dirs.systems_final}/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'output_deffnm':'ladder',
        'traces': ['Temperature','Density','Volume'],
        'scatter': ('time(ps)',['Density'],'rho_v_ns.png'),
        'Tlo':300.0,
        'Thi':600.0,
        'deltaT':5,
        'ps_per_run':1000,
        'ps_per_rise':1000,
        'warmup_ps':5000,
        'P':1
    }

    def build_mdp(self,mdpname,**kwargs):
        """Builds the GROMACS mdp file required for a temperature-ladder MD simulation.

        Args:
            mdpname (str): name of mdp file
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        Tdomain=params['Thi']-params['Tlo']
        nSteps=int(Tdomain/np.abs(params['deltaT']))+1
        T0=params['Tlo'] if params['deltaT']>0.0 else params['Thi']
        T1=params['Thi'] if params['deltaT']>0.0 else params['Tlo']
        Tladder=np.linspace(T0,T1,nSteps)
        timepoints=[0.0,params['warmup_ps']]
        temppoints=[T0,T0]
        for i in range(nSteps):
            cT=Tladder[i]
            timepoints.append(timepoints[-1]+params['ps_per_rise'])
            temppoints.append(cT)
            timepoints.append(timepoints[-1]+params['ps_per_run'])
            temppoints.append(cT)
        duration=timepoints[-1]
        nMDsteps=int(duration/timestep)
        mod_dict={
            'ref_t':T0,
            'ref_p':params['P'],
            'gen-temp':T0,
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'single',
            'nsteps':nMDsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(mdpname,mod_dict)

def _strain_advice(edot,ps,logger=logger,what='strain'):
    """Says what strain a deformation will reach, and what a modulus from it is worth.

    The measured warning is not about the rate, it is about replicates.  Eight
    independent ramps of one protocol on one network gave shear moduli from 1.29 to
    1.74 GPa --- 10% scatter --- while each run's own fit reported a standard error far
    smaller than that.  So a single ramp gives a modulus to roughly 10%, which is often
    good enough, but it cannot give the uncertainty on it, and a fit's own error bar is
    not that uncertainty.  Averaging replicates also does something a single run cannot:
    it separates a real trend in the modulus with strain from one draw wandering.

    Rate does matter past yield --- 1.22 +/- 0.06 GPa at ``edot = 1e-3`` against
    0.96 +/- 0.02 at ``1e-4`` over the same window --- and a fast ramp spends little of
    its time at the low strains a tangent modulus would come from.  Whether that biases
    the modulus is not established: the one experiment that appeared to show it was a
    single trajectory and a replicate refuted it.

    A known bias that is *not* about rate: :class:`PostSimShear` has to switch the
    barostat off for the Cartesian component it shears, because Gromacs refuses the run
    otherwise, and that clamp raises *G* by about 29% --- 1.51 GPa clamped against 1.17
    unclamped on one network, where only the unclamped value agrees with the same
    network's *E* and Poisson ratio.  See ROADMAP.md.
    """
    total=edot*ps
    logger.info(f'This run reaches {what} {total:.3f} at {edot:g} ps^-1 over {ps} ps')
    logger.info('One ramp gives a modulus to about 10% but not its uncertainty: eight '
                'replicates of one protocol spanned 1.29 to 1.74 GPa while each fit '
                'reported a much smaller error. Run several seeds and take the scatter '
                'between them as the error bar.')
    if edot>=1e-3:
        logger.warning(
            f'{what} passes 0.02 after only {0.02/edot:.0f} ps at edot={edot:g}, so '
            f'little of this run samples the low strains a modulus would come from; what '
            f'is fitted is a secant over a large window.  edot=1e-4 reaches 0.03 in the '
            f'same wall time if that is what you want.  Rate dependence past yield is '
            f'real and measured; whether the rate biases the modulus itself is not '
            f'established.')


class PostSimDeform(PostSimMD):
    """ a class to handle a uniaxial deformation MD simulation
    """
    default_params={
        'subdir': f'{pfs.Dirs.postsim}/deform-x',
        'input_top': f'{pfs.Dirs.systems_final}/final.top',
        'input_gro': f'{pfs.Dirs.postsim}/equilibrate/equilibrate.gro',
        'input_grx': f'{pfs.Dirs.systems_final}/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'output_deffnm':'deform-x',
        'traces': ['Box-X','Pres-XX'],
        'scatter': ('Box-X',['Pres-XX'],'tension_v_xlength.png'),
        'direction':'x',
        'T':300.0,
        'P':1.0,
        'ps':1000,
        'edot': 0.001 # strain rate in ps^-1
    }

    def build_mdp(self,mdpname,**kwargs):
        """Builds the GROMACS mdp file required for a uniaxial deformation MD simulation.

        Args:
            mdpname (str): name of mdp file
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        duration=params['ps']
        nsteps=int(duration/timestep)
        box=kwargs.get('box',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]))
        edot=params.get('edot',0.0)
        direction=params.get('direction','x')
        mod_dict={
            'ref_t':params['T'],
            'ref_p':params['P'],
            'gen-temp':params['T'],
            'gen-vel':'yes',
            'tcoupl':'v-rescale',
            'nsteps': nsteps,
            'rlist': 1.2,
            'rcoulomb': 1.2,
            'rvdw': 1.2,
            'tau_t':0.5,
            'tau_p':1.0,
            'refcoord_scaling': 'com',
            'pcoupltype': 'anisotropic',
            # stress is what this stage exists to measure, so sample it far more often
            # than the packaged npt.mdp does.  At 1 ps the frames are already
            # near-independent, so ten times as many cut the standard error on a fitted
            # modulus by sqrt(10) for nothing but a larger edr.
            'nstcalcenergy': 50,
            'nstenergy': 50,
            # required from Gromacs 2025 whenever deform is combined with generated
            # velocities: the initial velocities must carry the flow profile the
            # deformation implies, and grompp refuses the run outright without it
            'deform-init-flow': 'yes'
            }
        if direction=='x':
            strain_vel=box[0][0]*edot
            mod_dict['ref_p']='0.0 1.0 1.0 0 0 0'
            mod_dict['compressibility']='0.0 4.5e-5 4.5e-5 0 0 0'
            mod_dict['deform']=f'{strain_vel:.3e} 0 0 0 0 0'
            params['output_deffnm'] = 'deform-x'
            params['traces']=['Box-X','Pres-XX']
            params['scatter']=('Box-X',['Pres-XX'],'tension_v_xlength.png')
        elif direction=='y':
            strain_vel=box[1][1]*edot
            mod_dict['ref_p']='1.0 0.0 1.0 0 0 0'
            mod_dict['compressibility']='4.5e-5 0.0 4.5e-5 0 0 0'
            mod_dict['deform']=f'0 {strain_vel:.3e} 0 0 0 0'
            params['output_deffnm'] = 'deform-y'
            params['traces']=['Box-Y','Pres-YY']
            params['scatter']=('Box-Y',['Pres-YY'],'tension_v_ylength.png')
        elif direction=='z':
            strain_vel=box[2][2]*edot
            mod_dict['ref_p']='1.0 1.0 0.0 0 0 0'
            mod_dict['compressibility']='4.5e-5 4.5e-5 0.0 0 0 0'
            mod_dict['deform']=f'0 0 {strain_vel:.3e} 0 0 0'
            params['output_deffnm'] = 'deform-z'
            params['traces']=['Box-Z','Pres-ZZ']
            params['scatter']=('Box-Z',['Pres-ZZ'],'tension_v_zlength.png')
        else:
            logger.error(f'Bad direction for uniaxial strain {direction}')
            return
        _strain_advice(edot,duration)

        mdp_modify(mdpname,mod_dict)

class PostSimShear(PostSimMD):
    """ a class to handle a constant-rate simple shear MD simulation

    The shear counterpart of :class:`PostSimDeform`.  Gromacs' ``deform`` moves one
    off-diagonal element of the box matrix at a constant rate; the shear modulus is then
    the slope of shear stress against engineering shear strain, exactly as Young's
    modulus is the slope of tensile stress against tensile strain.

    ``direction`` names the plane: ``xy`` shears the x face along y (box element YX),
    ``xz`` and ``yz`` likewise.  The normal pressures stay coupled at ``P`` while the
    off-diagonal is driven, because the anisotropic barostat is given zero
    compressibility off the diagonal and so does not fight the deformation.
    """
    # box element driven, the box length the strain is measured against, and the two
    # energy terms to trace, per shear plane.  The deform slot order is Gromacs'
    # own: XX YY ZZ YX ZX ZY.
    # 'zero' is the diagonal compressibility that must be switched off: driving an
    # off-diagonal element moves one Cartesian component of a box vector, and Gromacs
    # refuses to have the barostat acting on that same component of another vector
    # ("spurious periodicity effects").  Shearing YX or ZX moves x; ZY moves y.
    _planes={
        'xy': {'slot': 3, 'ref': 1, 'zero': 0, 'box': 'Box-YX', 'pres': 'Pres-XY'},
        'xz': {'slot': 4, 'ref': 2, 'zero': 0, 'box': 'Box-ZX', 'pres': 'Pres-XZ'},
        'yz': {'slot': 5, 'ref': 2, 'zero': 1, 'box': 'Box-ZY', 'pres': 'Pres-ZY'},
    }
    default_params={
        'subdir': f'{pfs.Dirs.postsim}/shear-xy',
        'input_top': f'{pfs.Dirs.systems_final}/final.top',
        'input_gro': f'{pfs.Dirs.postsim}/equilibrate/equilibrate.gro',
        'input_grx': f'{pfs.Dirs.systems_final}/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'output_deffnm':'shear-xy',
        'traces': ['Box-YX','Pres-XY'],
        'scatter': ('Box-YX',['Pres-XY'],'shearstress_v_yxbox.png'),
        'direction':'xy',
        'T':300.0,
        'P':1.0,
        'ps':1000,
        'edot': 0.001 # shear rate in ps^-1
    }

    def build_mdp(self,mdpname,**kwargs):
        """Builds the Gromacs mdp file required for a constant-rate simple shear.

        Args:
            mdpname (str): name of mdp file
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        nsteps=int(params['ps']/timestep)
        box=kwargs.get('box',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]))
        direction=params.get('direction','xy')
        plane=self._planes.get(direction)
        if plane is None:
            logger.error(f'Bad plane for simple shear {direction}; expected one of '
                         + ', '.join(sorted(self._planes)))
            return
        # the box element is driven at L*edot, so edot is a true strain rate in ps^-1
        rate=box[plane['ref']][plane['ref']]*params.get('edot',0.0)
        deform=['0']*6
        deform[plane['slot']]=f'{rate:.3e}'
        compress=['4.5e-5','4.5e-5','4.5e-5','0','0','0']
        compress[plane['zero']]='0'
        mod_dict={
            'ref_t':params['T'],
            'gen-temp':params['T'],
            'gen-vel':'yes',
            'tcoupl':'v-rescale',
            'nsteps': nsteps,
            'rlist': 1.2,
            'rcoulomb': 1.2,
            'rvdw': 1.2,
            'tau_t':0.5,
            'tau_p':1.0,
            'refcoord_scaling': 'com',
            'pcoupltype': 'anisotropic',
            # the Cartesian component the shear moves is uncoupled; the other two
            # normal directions stay at P.  Gromacs rejects the run otherwise.
            'ref_p':f'{params["P"]} {params["P"]} {params["P"]} 0 0 0',
            'compressibility':' '.join(compress),
            'deform':' '.join(deform),
            # as in PostSimDeform: the fitted modulus is only as good as the stress
            # sampling, and denser output is the cheapest way to improve it
            'nstcalcenergy': 50,
            'nstenergy': 50,
            # required from Gromacs 2025 whenever deform is combined with generated
            # velocities: without it the initial velocities carry no flow profile and
            # grompp refuses the run outright
            'deform-init-flow':'yes',
            }
        params['output_deffnm']=f'shear-{direction}'
        params['traces']=[plane['box'],plane['pres']]
        params['scatter']=(plane['box'],[plane['pres']],
                           f'shearstress_v_{direction}.png')
        _strain_advice(params.get('edot',0.0),params['ps'],what='shear strain')
        mdp_modify(mdpname,mod_dict)

class PostsimConfiguration:
    """ handles reading and parsing a postcure simulation input config file.
        Config file format
        
        - { key1: {<paramdict>}}
        - { key2: {<paramdict>}} 
        
        ...

        The config file is a list of single-element dictionaries, whose single keyword
        indicates the type of MD simulation to be run; simulations are run in the order
        they appear in the config file.

        Currently allowed simulation types:

        - 'equilibrate': simple NPT equilibration;
        - 'anneal': simple simulated annealing;
        - 'ladder': temperature ladder;
        - 'deform': constant strain-rate uniaxial deformation;
        - 'shear': constant-rate simple shear, for the shear modulus;
        
        """
    default_classes={'equilibrate':PostSimMD,'anneal':PostSimAnneal,'ladder':PostSimLadder,'deform':PostSimDeform,'shear':PostSimShear}
    def __init__(self):
        self.cfgFile=''
        self.baselist=[]
        self.stagelist=[]

    @classmethod
    def read(cls,filename,parse=True,**kwargs):
        """Generates a new PostsimConfiguration object by reading in the JSON or YAML file indicated by filename.

        Args:
            filename (str): name of file from which to read new PostsimConfiguration object
            parse (bool): if True, parse the input configuration file, defaults to True

        Raises:
            Exception: if extension of filename is not '.json' or '.yaml' or '.yml'

        Returns:
            PostsimConfiguration: a new PostsimConfiguration object
        """
        basename,extension=os.path.splitext(filename)
        if extension=='.json':
            return cls._read_json(filename,parse,**kwargs)
        elif extension=='.yaml' or extension=='.yml':
            return cls._read_yaml(filename,parse,**kwargs)
        else:
            raise Exception(f'Unknown config file extension {extension}')

    @classmethod
    def _read_json(cls,filename,parse=True,**kwargs):
        """Creates a new PostsimConfiguration object by reading from JSON input.

        Args:
            filename (str): name of JSON file
            parse (bool): if True, parse the JSON data, defaults to True

        Returns:
            PostsimConfiguration: a new PostsimConfiguration object
        """
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.baselist=json.load(f)
            assert type(inst.baselist)==list,f'Poorly formatted {filename}'
        if parse: inst.parse(**kwargs)
        return inst

    @classmethod
    def _read_yaml(cls,filename,parse=True,**kwargs):
        """Creates a new PostsimConfiguration object by reading from YAML input.

        Args:
            filename (str): name of YAML file
            parse (bool): if True, parse the YAML data, defaults to True

        Returns:
            PostsimConfiguration: a new PostsimConfiguration object
        """
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.baselist=yaml.safe_load(f)
            assert type(inst.baselist)==list,f'Poorly formatted {filename}'
        if parse: inst.parse(**kwargs)
        return inst

    def parse(self,**kwargs):
        """Parses a PostsimConfiguration file to build the list of stages to run."""
        for p in self.baselist:
            assert len(p)==1,f'Poorly formatted {self.cfgFile}; each stanza may have only one keyword'
            simtype=list(p.keys())[0]
            assert simtype in self.default_classes,f'Simulation type "{simtype}" in {self.cfgFile} not understood.'
            logger.info(f'passing in {p[simtype]}')
            self.stagelist.append(self.default_classes[simtype](p[simtype]))

def postsim(args):
    """Handles the postsim subcommand for managing post-cure production MD simulations.

    Args:
        args (argparse.Namespace): command-line arguments
    """
    setup_logging(args.loglevel, no_banner=args.no_banner)
    ess='y' if len(args.proj)==0 else 'ies'
    ogromacs={}
    if args.ocfg:
        ocfg=Configuration.read(args.ocfg)
        ogromacs=ocfg.gromacs
    cfg=PostsimConfiguration.read(args.cfg)
    logger.debug(f'{cfg.baselist}')
    logger.info(f'Project director{ess}: {args.proj}')
    software.sw_setup()
    logger.debug(f'ogromacs {ogromacs}')
    for d in args.proj:
        pfs.pfs_setup(root=os.getcwd(),topdirs=pfs.Dirs.postsim_topdirs,verbose=True,projdir=d,reProject=False,userlibrary=pfs.resolve_user_library(args.lib))
        pfs.go_to(pfs.Dirs.postsim)
        for stage in cfg.stagelist:
            stage.do(mdp_pfx='local',**ogromacs)
        pfs.go_root()

