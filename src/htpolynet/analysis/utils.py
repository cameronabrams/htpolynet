"""Various utility methods for plotting and postprocessing.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging
import os

import networkx as nx
import numpy as np
import pandas as pd

from scipy.constants import physical_constants
from scipy.optimize import curve_fit

from ..core.coordinates import Coordinates
from ..core.topocoord import TopoCoord
from ..external.gromacs import gmx_energy_trace

logger=logging.getLogger(__name__)

def density_from_gro(gro,mollib='./lib/molecules/parameterized',units='SI'):
    """Computes density from a Gromacs gro file.

    Args:
        gro (str): name of gro file
        mollib (str): location of parameterized molecular templates, defaults to './lib/molecules/parameterized'
        units (str): string indicating unit system, defaults to 'SI'

    Returns:
        float: density
    """
    C=Coordinates.read_gro(gro)
    resnames=list(set(C.A['resName']))
    templates={x:TopoCoord.from_top_gro(os.path.join(mollib,f'{x}.top'),os.path.join(mollib,f'{x}.gro')) for x in resnames}
    mass=0.0
    counts={}
    for nm,T in templates.items():
        if not nm in counts:
            counts[nm]={}
            for an in T.Topology.D['atoms']['atom']:
                counts[nm][an]=0
        these=C.A[C.A['resName']==nm]['atomName']
        for n in these:
            counts[nm][n]+=1
    for nm,adict in counts.items():
        T=templates[nm].Topology.D['atoms']
        for aname,c in adict.items():
            mass+=c*T[T['atom']==aname]['mass'].values[0]
    volume=np.prod(C.box.diagonal())
    fac=1.0
    if units=='SI':
        mfac=physical_constants['atomic mass constant'][0]
        lfac=1.e-9
        vfac=lfac**3
        fac=mfac/vfac
    return fac*mass/volume

""" These globals are used in traversing a project directory to extract data from edr files """
_system_dirs=['densification','precure',r'iter-{iter:d}','capping','postcure']
_postsim_dirs=['anneal','equilibrate']
_measurement_dirs=['Tg','E']
_md_ensembles={'nvt':['Temperature','Potential'],'npt':['Temperature','Potential','Density'],'default':['Temperature','Potential']}
_indir_pfx={}
_indir_pfx['densification']=[r'densified-{ens:s}',r'densified-repeat-{repeat:d}-{ens:s}']
_indir_pfx['precure']=[r'preequilibration-{ens:s}',r'annealed',r'postequilibration-{ens:s}']
_indir_pfx['iter-n']=[r'1-cure_drag-stage-{stage:d}-{ens:s}',r'3-cure_relax-stage-{stage:d}-{ens:s}',r'4-cure_equilibrate-{ens:s}']
_indir_pfx['capping']=[r'7-cap_relax-stage-{stage:d}-{ens:s}',r'8-cap_equilibrate-{ens:s}']
_indir_pfx['postcure']=[r'preequilibration-{ens:s}',r'annealed',r'postequilibration-{ens:s}']
def _concat_from_edr(df,edr,names,add=[],add_if_missing=[('Density',0.0)]):
    """Concatenates rows onto dataframe df by reading data from an edr file.

    Args:
        df (pandas.DataFrame): a dataframe with data read in from edr files
        edr (str): name of a new edr file to read from
        names (list of str): names of energy-like quantities to be read in; each must have a column in df already
        add (list): names and values to add to dataset, defaults to []
        add_if_missing (list): names and values to add if not found in edr file, defaults to [('Density',0.0)]

    Returns:
        tuple: the dataframe and the scalar value of last time before this data increment is read in
    """
    xshift=0.0
    if not df.empty: xshift=df.iloc[-1]['time(ps)']
    if len(names)==0: return df,xshift
    # print(f'{add_if_missing}')
    this_df=gmx_energy_trace(edr,names,xshift=xshift)
    for m in add:
        nm,val=m
        this_df[nm]=np.ones(this_df.shape[0],dtype=type(val))*val
    for m in add_if_missing:
        nm,val=m
        if not nm in this_df:
            this_df[nm]=np.ones(this_df.shape[0],dtype=type(val))*val
    df=pd.concat((df,this_df),ignore_index=True)
    return df,xshift

def postsim_density_evolution(proj_dir,append_dirname=False):
    """Returns a single dataframe that is a concatenation of the csv files in the 'postsim' subdirectories.

    Args:
        proj_dir (str): name of complete project directory

    Returns:
        pandas.DataFrame: the dataframe
    """
    if not os.path.exists(proj_dir): return
    sysd=os.path.join(proj_dir,'postsim')
    fn=[f'{sysd}/anneal/anneal.csv',f'{sysd}/equilibrate/equilibrate.csv']
    df=pd.DataFrame()
    lt=0.0
    for f in fn:
        t=pd.read_csv(f,header=0,index_col=None)
        logger.debug(f'shifting by {lt}')
        t['time(ps)']+=lt
        lt=t['time(ps)'].to_list()[-1]
        df=pd.concat((df,t))
    if append_dirname:
        df.rename(columns={x:f'{x}-{proj_dir}' for x in df.columns})
    return df

def density_evolution(proj_dir):
    """Returns a single dataframe containing density, temperature, number of bonds vs time by reading all edrs in the correct order from a complete project directory.

    Args:
        proj_dir (str): name of complete project directory

    Returns:
        pandas.DataFrame: the dataframe
    """
    if not os.path.exists(proj_dir): return
    sysd=os.path.join(proj_dir,'systems')
    df=pd.DataFrame()
    xshift=0.0
    transition_times=[xshift]
    interval_labels=[]
    markers=[]
    nbonds=0
    for subd in _system_dirs:
        # print(f'subd {subd}')
        if subd==r'iter-{iter:d}' or subd=='capping' or subd=='densification':
            # loooking for iterations
            # print(iter_mark0)
            iter=1
            iter_subd_rel=subd if subd in ['capping','densification'] else subd.format(iter=iter)
            iter_subd=os.path.join(sysd,iter_subd_rel)
            while os.path.exists(iter_subd):
                logger.info(f'Reading edrs in {iter_subd_rel}...')
                if r'iter' in subd:
                    markers.append(df.iloc[-1]['time(ps)'])
                dirkey='iter-n' if subd==r'iter-{iter:d}' else subd
                if subd==r'iter-{iter:d}':
                    if os.path.exists(os.path.join(iter_subd,'2-cure_update-bonds.csv')):
                        bdf=pd.read_csv(os.path.join(iter_subd,'2-cure_update-bonds.csv'),sep=r'\s+',header=0,index_col=None)
                        nbonds+=bdf.shape[0]
                    else:
                        if os.path.exists(os.path.join(iter_subd,'6-cap_update-bonds.csv')):
                            bdf=pd.read_csv(os.path.join(iter_subd,'6-cap_update-bonds.csv'),sep=r'\s+',header=0,index_col=None)
                            nbonds+=bdf.shape[0]
                # print(f'iter_subd {iter_subd} dirkey {dirkey}')
                for pfx in _indir_pfx[dirkey]:
                    # print(f'pfx {pfx}')
                    if r'stage' in pfx or r'repeat' in pfx:
                        stg=1
                        stg_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(stage=stg,repeat=stg,ens=x))+r'.edr') for x in _md_ensembles])
                        # print(f'{iter_subd} {stg} {pfx} stg_present {stg_present}')
                        while stg_present:
                            for ens,qtys in _md_ensembles.items():
                                edr_pfx=os.path.join(iter_subd,pfx.format(stage=stg,repeat=stg,ens=ens))
                                gro=edr_pfx+r'.gro'
                                # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                                if os.path.exists(edr_pfx+r'.edr'):
                                    density=0.0
                                    if ens!='npt':
                                        if os.path.exists(gro):
                                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                                            # print(f'density from gro {density}')
                                    df,xshift=_concat_from_edr(df,edr_pfx,qtys,add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                                    transition_times.append(xshift)
                                    interval_labels.append([edr_pfx])
                            stg+=1
                            stg_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(stage=stg,repeat=stg,ens=x))+r'.edr') for x in _md_ensembles])
                            # print(f'->{iter_subd} {stg} {pfx} stg_present {stg_present}')
                    elif r'equil' in pfx or r'densif' in pfx:
                        eq_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(ens=x))+r'.edr') for x in _md_ensembles])
                        # print(f'{iter_subd} eq_present? {eq_present}')
                        if eq_present:
                            for ens,qtys in _md_ensembles.items():
                                edr_pfx=os.path.join(iter_subd,pfx.format(ens=ens))
                                gro=edr_pfx+r'.gro'
                                # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'),ens,qtys)
                                if os.path.exists(edr_pfx+r'.edr'):
                                    density=0.0
                                    if ens!='npt':
                                        if os.path.exists(gro):
                                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                                    df,xshift=_concat_from_edr(df,edr_pfx,qtys,add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                                    transition_times.append(xshift)
                                    interval_labels.append([edr_pfx])
                iter+=1       
                iter_subd_rel=subd.format(iter=iter)
                iter_subd=os.path.join(sysd,iter_subd_rel) if r'iter' in subd else 'I_BET_THIS_FILE_DNE'
            # if r'iter' in subd: markers.append(df.iloc[-1]['time (ps)'])
        else:
            this_subd=os.path.join(sysd,subd)
            logger.info(f'Reading edrs in {subd}...')

            for pfx in _indir_pfx[subd]:
                # print(f'this_subd {this_subd} pfx {pfx}')
                if 'ens' in pfx:
                    for ens,qtys in _md_ensembles.items():
                        edr_pfx=os.path.join(this_subd,pfx.format(ens=ens))
                        gro=edr_pfx+r'.gro'
                        # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                        # print(gro,os.path.exists(gro))
                        if os.path.exists(edr_pfx+r'.edr'):
                            density=0.0
                            if ens=='nvt' and os.path.exists(gro):
                                density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                            df,xshift=_concat_from_edr(df,edr_pfx,qtys,add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                            transition_times.append(xshift)
                            interval_labels.append([edr_pfx])
                else:
                    edr_pfx=os.path.join(this_subd,pfx.format(ens=ens))
                    gro=edr_pfx+r'.gro'
                    # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                    if os.path.exists(edr_pfx+r'.edr'):
                        density=0.0
                        if os.path.exists(gro):
                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                        df,xshift=_concat_from_edr(df,edr_pfx,_md_ensembles['default'],add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                        transition_times.append(xshift)
                        interval_labels.append([edr_pfx])
    return df,transition_times,markers,interval_labels

def graph_from_bondsfile(bondsfile):
    """Generates a networkx Graph in which each node is a molecule (from 'mi' and 'mj' records in the bondsfile) and edges indicate two molecules are joined by a covalent bond.

    Args:
        bondsfile (str): name of bondsfile

    Returns:
        networkx.Graph: the graph
    """
    G=nx.Graph()
    df=pd.read_csv(bondsfile,header=0,index_col=None,sep=r'\s+')
    for i,r in df.iterrows():
        G.add_edge(r['mi'],r['mj'])
    return G

def mwbxl(G:nx.Graph,crosslinker='GMA',monomer='STY'):
    """Computes the histogram of monomer counts 'n' between crosslinking sites using a molecular connectivity graph; used mainly for vinyl-based polymerizations.

    Args:
        G (nx.Graph): molecular connectivity graph
        crosslinker (str): name of crosslinker molecule, defaults to 'GMA'
        monomer (str): name of monomer, defaults to 'STY'

    Returns:
        pd.DataFrame: a dataframe of 'n' and 'count'
    """
    xG=G.copy()
    # traverse edges and label hetero/homo participants
    for n in xG.nodes():
        xG.nodes[n]['neighbor_count']=0
    for u,v in xG.edges():
        n,m=xG.nodes[u],xG.nodes[v]
        n['neighbor_count']+=1
        m['neighbor_count']+=1
        if n['molecule_name']!=m['molecule_name']:
            n['monomer_type']='hetero'
            m['monomer_type']='hetero'
        else:
            n['monomer_type']='homo'
            m['monomer_type']='homo'
    wipe_us=[]
    for n in xG.nodes():
        if xG.nodes[n]['molecule_name']==crosslinker:
            wipe_us.append(n)
    xG.remove_nodes_from(wipe_us)
    chaintype=[]
    ch=nx.connected_components(xG)
    mwhist={}
    nb=0
    for c in ch:
        nhetero=0
        for n in c:
            node=xG.nodes[n]
            if node.get('monomer_type','none')=='hetero':
                nhetero+=1
        if nhetero==0:
            chaintype.append('isolated')
        elif nhetero==1:
            if len(c)==1:
                node=list(xG.nodes.values())[0]
                if node['neighbor_count']==1:
                    chaintype.append('dangling')
                elif node['neighbor_count']==2:
                    chaintype.append('bridging')
                    if len(c) not in mwhist:
                        mwhist[len(c)]=0
                    mwhist[len(c)]+=1
            else:
                chaintype.append('dangling')
        elif nhetero==2:
            chaintype.append('bridging')
            if len(c) not in mwhist:
                mwhist[len(c)]=0
            mwhist[len(c)]+=1
            nb+=1
    maxlen=max(list(mwhist.keys()))
    mw=[]
    counts=[]
    for n in range(maxlen):
        cnt=mwhist.get(n,0)
        mw.append(n)
        counts.append(cnt)
    df=pd.DataFrame({'n':mw,'counts':counts})
    return df

def clusters(G:nx.Graph):
    """Performs a clustering analysis and returns a histogram of cluster sizes (in numbers of molecules) as a pandas DataFrame.

    Args:
        G (nx.Graph): molecular connectivity graph

    Returns:
        pd.DataFrame: cluster size histogram
    """
    counts={}
    for c in sorted(nx.connected_components(G),key=len,reverse=True):
        size=len(c)
        if not size in counts:
            counts[size]=0
        counts[size]+=1
    sizes=[]
    numbers=[]
    for s,c in counts.items():
        sizes.append(s)
        numbers.append(c)
    df=pd.DataFrame({'sizes':sizes,'counts':numbers})
    df.sort_values(by='sizes',ascending=False,inplace=True)
    return df

def compute_tg(T,v,n_points=[10,20]):
    """Performs a Tg determination from (volume or density) vs temperature data by fitting lines to the low-T region and another to the high-T region, taking Tg as the temperature at which they intersect.

    Args:
        T (numpy.array): temperature values
        v (numpy.array): volume or density data
        n_points (list): number of points on the low and high side to fit lines to, defaults to [10,20]
    """
    def func(x,a,b):
        return x*a+b
    hot_par=[]
    cold_par=[]
    x=np.array(T)[:n_points[0]]
    y=np.array(v)[:n_points[0]]
    # logger.info(f'Tg: {len(x)} data points on cold side; begins at {x[0]}')
    cold_par,pcov=curve_fit(func,x,y)
    sse=(y-func(x,*cold_par))**2
    sst=(y-x.mean())**2
    r2=1-sse.sum()/sst.sum()
    # logger.info(f'cold fit: {cold_par} {r2}')

    x=np.array(T)[-n_points[1]:]
    y=np.array(v)[-n_points[1]:]
    # logger.info(f'Tg: {len(x)} data points on hot side; ends at {x[-1]}')
    hot_par,pcov=curve_fit(func,x,y)
    sse=(y-func(x,*hot_par))**2
    sst=(y-x.mean())**2
    r2=1-sse.sum()/sst.sum()
    # logger.info(f'hot fit: {hot_par} {r2}')
    Tg=-1
    if len(hot_par)>0 and len(cold_par)>0:
        Tg=-(hot_par[1]-cold_par[1])/(hot_par[0]-cold_par[0])
    return Tg,cold_par,hot_par

def compute_G(shear_strain,shear_stress,fit_domain=[10,100],fit_strain=None):
    """Computes the shear modulus by a linear fit to the elastic regime of a simple shear.

    The same fit as :func:`compute_E`, named separately because the inputs are not the
    same quantities: the engineering shear strain is the driven off-diagonal box element
    over the length of the sheared face (it starts at zero, so there is no -1), and the
    shear stress is the negated off-diagonal of the pressure tensor.
    ``htpolynet postsim`` writes both columns for a ``shear`` stage.

    Args:
        shear_strain (numpy.array): engineering shear strain, dimensionless
        shear_stress (numpy.array): shear stress, in the units of the pressure trace
        fit_domain (list): domain over which the fit is made, in rows, defaults to [10,100]
        fit_strain (list): [lo, hi] strain window to fit instead; overrides fit_domain

    Returns:
        tuple(float,float): G and R2 from the fit
    """
    return compute_E(shear_strain,shear_stress,fit_domain=fit_domain,fit_strain=fit_strain)

def compute_E(strain,stress,fit_domain=[10,100],fit_strain=None):
    """Computes the Young's modulus by performing a linear fit to an elastic regime in stress-vs-strain data.

    Args:
        strain (numpy.array): strain values
        stress (numpy.array): stress values
        fit_domain (list): domain over which fit is made, in rows, defaults to [10,100]
        fit_strain (list): [lo, hi] strain window to fit instead; overrides fit_domain.
            A window in rows covers a different strain whenever the rate, the run
            length or the energy output interval changes; a strain window does not.

    Returns:
        tuple(float,float): E and R2 from fit
    """
    if fit_strain is not None:
        lo,hi=fit_strain
        s=np.asarray(strain,dtype=float)
        mask=(s>=lo)&(s<=hi)
        if mask.sum()<2:
            raise ValueError(f'fewer than two points with strain in [{lo}, {hi}]; '
                             f'the data reach {s.max():.4g}')
        x=s[mask]
        y=np.asarray(stress,dtype=float)[mask]
    else:
        x=np.array(strain[fit_domain[0]:fit_domain[1]])
        y=np.array(stress[fit_domain[0]:fit_domain[1]])
    # logger.info(f'x: {x[0]} -> {x[-1]}')
    # logger.info(f'y: {y[0]} -> {y[-1]}')
    # The intercept is free.  A modulus is the slope; forcing the line through the
    # origin assumes the structure starts unstressed, and a glass equilibrated for
    # tens of ps does not -- example 6's reference structure started at +150-200 bar,
    # which through the origin adds offset/strain to E, about +1 GPa over a 2% window.
    a,b=np.polyfit(x,y,1)
    sse=(y-(a*x+b))**2
    sst=(y-y.mean())**2
    r2=1-sse.sum()/sst.sum()
    return a,r2