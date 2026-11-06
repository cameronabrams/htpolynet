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
from scipy.stats import f as f_dist

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

def tg_window_range(T,v,glassy=(10,15,20,25),rubbery=(10,15,20,25)):
    """Lowest and highest Tg over a grid of fit windows.

    A Tg is the crossing of two fitted lines, so when the lines are nearly parallel a
    small change in either window moves it a long way --- or off the temperature scale.
    Seen 2026-10-15 on a rebuild of example 6: the heating ladder's glassy slope was a
    third of the usual, the default windows gave 332.5 K against the reference build's
    529.2 K, and other windows gave anything from -1128 to 1640 K.  The spread over
    windows is what says whether a Tg can be quoted.

    Args:
        T (array-like): temperatures, sorted ascending
        v (array-like): density or volume
        glassy (tuple): glassy-side point counts to try
        rubbery (tuple): rubbery-side point counts to try

    Returns:
        tuple: (lowest, highest) Tg over the windows that fit in the data, or None if none do
    """
    n=len(T)
    tgs=[compute_tg(T,v,n_points=[c,h])[0] for c in glassy for h in rubbery if c+h<=n]
    return (min(tgs),max(tgs)) if tgs else None

def _tg_hyperbola(T,rho0,T0,a,b,c):
    """rho0 + a (T-T0) + b H(T-T0), H(t) = t/2 + sqrt(t^2/4 + exp(c)).

    Asymptotes: slope a well below T0 (glassy) and a+b well above (rubbery), crossing
    at T0, which is Tg; exp(c/2) sets the width of the transition between them."""
    t=np.asarray(T,dtype=float)-T0
    return rho0+a*t+b*(0.5*t+np.sqrt(0.25*t*t+np.exp(c)))

def fit_tg_hyperbola(T,v,width_bounds=(1.0,200.0),P_hat=0.9,n_synth=0,seed=0):
    """Tg from one hyperbola fitted to the whole density (or volume) curve.

    The two-line intersection of :func:`compute_tg` fits each side separately and
    intersects them, so how many points each side gets moves Tg -- by 100 K or more on
    some shipped examples.  Following Patrone et al., Polymer 87, 246 (2016), a
    hyperbola whose asymptotes are the glassy and rubbery lines is fitted to every point
    at once, weighted by sigma_j proportional to T_j, and Tg is read off it two ways:
    ``Tg`` is T0, the hyperbola's center, and ``Tg_tangents`` is where its tangents at
    the two ends of the ladder cross.  When the data reach both asymptotes the two
    agree; their gap is the paper's own check that they do.

    The screening criterion (their Eqs. 4-5): the slope reaches a fraction P_hat of the
    way between its asymptotes at T0 +/- DeltaT, DeltaT = e^(c/2) (2P-1)/sqrt(P(1-P)).
    A ladder that does not extend past that interval on both sides does not contain the
    asymptotic regimes, and its Tg is held out.  Also held out: a width on its upper
    bound, and a bend that does not beat a straight line (F-test, p > 0.01) -- a
    straight line otherwise fits as a sharp break at a noise feature.  A width on its lower
    bound is kept and flagged ``sharp``: the fit is then a two-line fit with its break
    where the misfit is least.

    Args:
        T (array-like): temperatures, K
        v (array-like): density or volume
        width_bounds (tuple): allowed e^(c/2), K
        P_hat (float): slope convergence that defines the transition interval
        n_synth (int): synthetic datasets for the within-build uncertainty (their
            section 5.1): residuals give a noise level, rising as T, from which data
            are redrawn about the fit and refitted; 0 skips it
        seed (int): random seed for the synthetic data

    Returns:
        dict: Tg, Tg_tangents, Tg_se (from the covariance), Tg_synth_sd, glassy and
            rubbery slopes, width, transition (T0-DeltaT, T0+DeltaT), sharp, held_out,
            reason, params; None if the fit does not converge
    """
    T=np.asarray(T,dtype=float); v=np.asarray(v,dtype=float)
    o=np.argsort(T); T,v=T[o],v[o]
    p,cov=_fit_hyperbola(T,v,width_bounds)
    if p is None:
        return None
    rho0,T0,a,b,c=p
    w_lo,w_hi=width_bounds
    width=float(np.exp(c/2))
    dT=width*(2*P_hat-1)/np.sqrt(P_hat*(1-P_hat))
    sharp=bool(width<=w_lo*1.001)
    # does the bend beat a straight line?  An F-test on the weighted residuals, not the
    # covariance of b: at the sharp-width bound that covariance is ill-conditioned and
    # would hold out a real break
    n=len(T)
    ss_h=float(np.sum(((v-_tg_hyperbola(T,*p))/T)**2))
    lin=np.polyfit(T,v,1,w=1.0/T)
    ss_l=float(np.sum(((v-np.polyval(lin,T))/T)**2))
    F=((ss_l-ss_h)/3.0)/(ss_h/(n-5)) if ss_h>0 and n>5 else float('inf')
    p_line=float(f_dist.sf(F,3,n-5)) if np.isfinite(F) else 0.0
    reason=None
    if T0<=T[0]+1.0 or T0>=T[-1]-1.0:
        reason='Tg at the edge of the data'
    elif width>=w_hi*0.999:
        reason=f'transition width at its {w_hi:g} K bound'
    elif p_line>0.01:
        reason=f'no significant change of slope (a straight line fits as well, p = {p_line:.2f})'
    elif T[0]>T0-dT or T[-1]<T0+dT:
        reason=(f'the ladder ({T[0]:.0f}-{T[-1]:.0f} K) does not reach past the transition, '
                f'{T0-dT:.0f}-{T0+dT:.0f} K at {P_hat:.0%} slope convergence')
    Tt=_tangent_crossing(p,T[0],T[-1])
    se=float(np.sqrt(cov[1,1])) if np.isfinite(cov[1,1]) and not sharp else float('nan')
    synth=None
    if n_synth>0:
        resid=v-_tg_hyperbola(T,*p)
        scale=np.sqrt(np.mean((resid/T)**2))
        rng=np.random.default_rng(seed)
        tgs=[]
        for _ in range(n_synth):
            q,_=_fit_hyperbola(T,_tg_hyperbola(T,*p)+rng.normal(0.0,scale*T),width_bounds,p0=p)
            if q is not None:
                tgs.append(q[1])
        synth=float(np.std(tgs,ddof=1)) if len(tgs)>1 else None
    return {'Tg':float(T0),'Tg_tangents':Tt,'Tg_se':se,'Tg_synth_sd':synth,
            'glassy':float(a),'rubbery':float(a+b),'width':width,
            'transition':(float(T0-dT),float(T0+dT)),'sharp':sharp,
            'held_out':reason is not None,'reason':reason,'params':p}

def _fit_hyperbola(T,v,width_bounds,p0=None):
    """Weighted least-squares fit of :func:`_tg_hyperbola`; (params, cov) or (None, None)."""
    w_lo,w_hi=width_bounds
    if p0 is None:
        k=max(3,len(T)//4)
        a0=np.polyfit(T[:k],v[:k],1)[0]
        a1=np.polyfit(T[-k:],v[-k:],1)[0]
        T0=0.5*(T[0]+T[-1])
        p0=[np.interp(T0,T,v),T0,a0,a1-a0,2*np.log(20.0)]
    lo=[-np.inf,T[0],-np.inf,-np.inf,2*np.log(w_lo)]
    hi=[np.inf,T[-1],np.inf,np.inf,2*np.log(w_hi)]
    p0=np.clip(np.asarray(p0,dtype=float),np.array(lo)+1e-9,np.array(hi)-1e-9)
    try:
        # sigma_j ~ T_j: Patrone et al. take sigma_j ~ T_j^p, 1/2 <= p <= 2, and find p
        # barely matters
        return curve_fit(_tg_hyperbola,T,v,p0=p0,sigma=T,bounds=(lo,hi),maxfev=40000)
    except (RuntimeError,ValueError):
        return None,None

def _tangent_crossing(p,T_lo,T_hi):
    """Where the fitted hyperbola's tangents at T_lo and T_hi cross."""
    def tangent(Tx,h=1e-3):
        y=_tg_hyperbola(Tx,*p)
        m=(_tg_hyperbola(Tx+h,*p)-_tg_hyperbola(Tx-h,*p))/(2*h)
        return m,y-m*Tx
    m1,b1=tangent(T_lo)
    m2,b2=tangent(T_hi)
    return float((b2-b1)/(m1-m2)) if abs(m1-m2)>1e-12 else float('nan')

def _align_builds(curves,width_bounds=(1.0,200.0)):
    """Each build's curve shifted onto a common density, by a joint fit.

    Replica builds of one system differ in density by a few kg/m^3 (3-11 across the
    v2.16.0 replicates), which is as large as the change of slope a Tg fit looks for,
    so pooling raw points turns the offsets into noise and hides the transition.  This
    fits one hyperbola -- shared T0, slopes and width -- with a separate rho0 for each
    build, and returns the curves with those offsets removed.
    """
    Ts=[np.asarray(c[0],dtype=float) for c in curves]
    Vs=[np.asarray(c[1],dtype=float) for c in curves]
    n=len(curves)
    T=np.concatenate(Ts); idx=np.concatenate([np.full(len(t),k) for k,t in enumerate(Ts)])
    V=np.concatenate(Vs)
    means=np.array([v.mean() for v in Vs])
    first=fit_tg_hyperbola(T,V-means[idx]+means.mean(),width_bounds)
    if first is None:
        return None
    _,T0,a0,b0,c0=first['params']
    def model(X,T0,a,b,c,*rho0):
        return np.asarray(rho0)[X[1].astype(int)]+_tg_hyperbola(X[0],0.0,T0,a,b,c)
    w_lo,w_hi=width_bounds
    lo=[T.min(),-np.inf,-np.inf,2*np.log(w_lo)]+[-np.inf]*n
    hi=[T.max(),np.inf,np.inf,2*np.log(w_hi)]+[np.inf]*n
    p0=[T0,a0,b0,c0]+[float(np.mean(v-_tg_hyperbola(t,0.0,T0,a0,b0,c0))) for t,v in zip(Ts,Vs)]
    p0=np.clip(np.asarray(p0,dtype=float),np.array(lo)+1e-9,np.array(hi)-1e-9)
    try:
        p,_=curve_fit(model,np.vstack([T,idx]),V,p0=p0,sigma=T,bounds=(lo,hi),maxfev=40000)
    except (RuntimeError,ValueError):
        return None
    rho0=np.asarray(p[4:])
    return [(t,v-(rho0[k]-rho0.mean())) for k,(t,v) in enumerate(zip(Ts,Vs))]

def pooled_tg(curves,n_boot=500,seed=0):
    """Tg from replica builds' curves fitted together, with a between-build error bar.

    Each build's curve is noisy on its own; fitting all builds at once uses them as data
    rather than averaging separately fitted answers.  Builds are first aligned on
    density (:func:`_align_builds`), since replicas differ by about as much as the
    signal.  The points within one build are correlated -- one structure, one
    trajectory -- so the pooled fit's own standard error understates the uncertainty,
    which is taken instead from a bootstrap over builds: resample whole builds with
    replacement, realign, refit.  That bootstrap is the uncertainty of the pooled Tg.
    Beside it, ``per_build_sd`` is the s.d. of the Tg of the builds that pass on their
    own, with ``n_kept`` saying how many.  It is conditioned on that filter, so it is
    not the scatter of an arbitrary build: on example 3 heating it is 7.7 K over the
    3 of 5 kept, and 73.6 K over all five, one of which the fit puts on the ladder's
    edge.  Quote it with k/n.  A pooled fit that is sharp when most builds are refused
    on their own is held out.
    Patrone et al. (2016) found exactly this: two
    simulations of one chemistry gave Tg distributions that did not overlap, because
    the structures crosslinking produces differ.

    Args:
        curves (list): one (T, v) pair of arrays per build
        n_boot (int): bootstrap resamples
        seed (int): random seed for the resampling

    Returns:
        dict: Tg and fit (the pooled :func:`fit_tg_hyperbola` result), Tg_boot_sd (None
            with one build or too few usable resamples), per_build (each build's own
            fit), n_builds, n_kept, per_build_sd (s.d. over the n_kept kept builds'
            Tg, None with fewer than two)
    """
    per=[fit_tg_hyperbola(t,y) for t,y in curves]
    def fit_pool(cs):
        aligned=_align_builds(cs) if len(cs)>1 else cs
        if aligned is None:
            return None
        return fit_tg_hyperbola(np.concatenate([c[0] for c in aligned]),
                                np.concatenate([c[1] for c in aligned]))
    pool=fit_pool(curves)
    kept=[f['Tg'] for f in per if f is not None and not f['held_out']]
    if (pool is not None and not pool['held_out'] and pool['sharp']
            and len(curves)>1 and 2*len(kept)<len(curves)):
        # pooling n builds buys significance with n: five straight lines can pool into
        # a "sharp" break that none of them shows (example 2 cooling, 2026-10-23: 4 of 5
        # builds refused, pooled kept at width 1 K, local slope showing no step)
        pool['held_out']=True
        pool['reason']=(f'a sharp break in the pool that {len(curves)-len(kept)} of '
                        f'{len(curves)} builds do not show on their own')
    sd=None
    if len(curves)>1 and pool is not None and not pool['held_out']:
        rng=np.random.default_rng(seed)
        tgs=[]
        for _ in range(n_boot):
            pick=rng.integers(0,len(curves),len(curves))
            if len(set(pick.tolist()))<2:
                continue
            r=fit_pool([curves[i] for i in pick])
            if r is not None and not r['held_out']:
                tgs.append(r['Tg'])
        sd=float(np.std(tgs,ddof=1)) if len(tgs)>1 else None
    return {'Tg':pool['Tg'] if pool and not pool['held_out'] else None,'fit':pool,
            'Tg_boot_sd':sd,'per_build':per,'n_builds':len(curves),'n_kept':len(kept),
            'per_build_sd':float(np.std(kept,ddof=1)) if len(kept)>1 else None}

def ceiling_drift(curves,cuts=(50.0,100.0)):
    """Pooled Tg refitted with the top of the ladder cut away.

    Passing every check in :func:`fit_tg_hyperbola` says a Tg is not obviously broken,
    not that the ladder was tall enough for it to have settled.  Example 1's heating
    Tg moved from 436.4 K on an 800 K ladder to 443.0 K at 700 K and 450.8 K at 650 K,
    kept by every check at each (htpolynet-sweep, plan C, 2026-10-25).  How far Tg moves
    as the top comes off is the part of its uncertainty the bootstrap cannot see.

    Args:
        curves (list): one (T, v) pair of arrays per build
        cuts (tuple): how many K to remove from the top of the ladder, one refit each

    Returns:
        list: one dict per cut, with cut, ceiling, Tg (None when the refit is held out)
            and reason
    """
    top=min(float(np.max(t)) for t,_ in curves)
    out=[]
    for cut in cuts:
        ceiling=top-cut
        trimmed=[(np.asarray(t)[np.asarray(t)<=ceiling],np.asarray(y)[np.asarray(t)<=ceiling])
                 for t,y in curves]
        r=pooled_tg(trimmed,n_boot=0)
        f=r['fit']
        out.append({'cut':cut,'ceiling':ceiling,'Tg':r['Tg'],
                    'reason':None if r['Tg'] is not None
                    else (f['reason'] if f is not None else 'did not converge')})
    return out

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