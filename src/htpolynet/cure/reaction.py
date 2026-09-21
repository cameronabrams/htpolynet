"""Handles Reaction objects.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

from copy import deepcopy
from enum import Enum

from ..analysis.plot import draw_reaction_dag
logger=logging.getLogger(__name__)

class reaction_stage(Enum):
    """Enumerated reaction stage
    """
    build=0   # only used for building molecules that use parameterized oligomers
    param=1   # used to generate parameterized oligomers that are not cure reactions
    cure=2    # used to generate parameterized oligormers that ARE cure reactions
    cap=3     # used to generate parameterized capped monomers
    repair=4  # parameterized template for postcure topology-repair operations
    unset=99
    def __str__(self):
        return self.name

class Reaction:
    default_directives={
        'name':'',
        'atoms':{},
        'bonds':[],
        'reactants':{},
        'product':'',
        'stage':reaction_stage.unset,
        'procession':{},
        'probability':1.0
    }
    def __init__(self,jsondict={}):
        """Generates a Reaction object using directives in jsondict.

        Args:
            jsondict (dict): directives for creating a new Reaction object, defaults to {}
        """
        self.jsondict=jsondict
        self.name=jsondict.get('name',self.default_directives['name'])
        self.atoms=jsondict.get('atoms',self.default_directives['atoms'])
        self.bonds=jsondict.get('bonds',self.default_directives['bonds'])
        self.reactants=jsondict.get('reactants',self.default_directives['reactants'])
        self.product=jsondict.get('product',self.default_directives['product'])
        self.stage=reaction_stage[jsondict.get('stage',str(self.default_directives['stage']))]
        self.procession=jsondict.get('procession',self.default_directives['procession'])
        self.probability=jsondict.get('probability',self.default_directives['probability'])
        for label in jsondict:
            if not label in self.default_directives:
                logging.debug(f'Ignoring unknown reaction directive "{label}"')
        self.symmetry_versions=[]
    
    def __str__(self):
        retstr=f'Reaction "{self.name}" ({str(self.stage)})\n'
        for i,r in self.reactants.items():
            retstr+=f'reactant {i}: {r}\n'
        retstr+=f'product {self.product}\n'
        for i,a in self.atoms.items():
            retstr+=f'atom {i}: {a}\n'
        for b in self.bonds:
            retstr+=f'bond {b}\n'
        return retstr

ReactionList = list[Reaction]

def parse_reaction_list(baselist:ReactionList):
    """Generates any new reactions inferred by procession (i.e., linear polymerization).

    Args:
        baselist (ReactionList): list of unprocessed reaction dicts read directly from configuration file

    Returns:
        ReactionList: entire set of all reactions, including new ones added here
    """
    processives=[]
    for i,R in enumerate(baselist):
        if R.procession:
            nReactions=[]
            logger.debug(f'Reaction {R.name} is processive: {R.procession}')
            final_product_name=R.product
            R.product=f'{final_product_name}_I0'
            original_name=R.name
            R.name=f'{R.name}-i0'
            base_reactant_key=R.procession['increment_resid']
            for c in range(R.procession['count']):
                nR=deepcopy(R)
                nR.name=f'{original_name}-i{c+1}'
                nR.product=f'{final_product_name}_I{c+1}' if c<(R.procession['count']-1) else final_product_name
                nR.reactants[base_reactant_key]=f'{final_product_name}_I{c}'
                for ak,arec in nR.atoms.items():
                    if arec['reactant']==base_reactant_key: arec['resid']+=(c+1)
                nR.procession={}
                nReactions.append(nR)
            processives.append((i,nReactions))
    lasti=0
    ret_reactions=baselist.copy()
    for rec in processives:
        atidx,sublist=rec
        atidx+=lasti
        ret_reactions=ret_reactions[0:atidx+1]+sublist+ret_reactions[atidx+1:]
        lasti+=len(sublist)
    return ret_reactions

def extract_molecule_reactions(rlist:ReactionList,plot=True):
    """Establishes the correct order of reactions so that molecular templates are created in a sensible order (reactants before products for all reactions).

    Args:
        rlist (ReactionList): list of all reactions
        plot (bool): flag to plot the reaction network, defaults to True

    Returns:
        list of tuples: ordered list of tuples, each of the form (name-of-product-molecule, Reaction)
    """
    working_rlist=rlist.copy()
    # Given an unsorted list of reactions (list of reactants +  one product), order a list of
    # of molecules so that no product comes before any of its reactants
    if not rlist: return []
    molecule_react_order=[]
    reactants=set()
    products=set()
    # first molecules to generate are the input reactants
    for R in working_rlist:
        for r in R.reactants.values():
            reactants.add(r)
        products.add(R.product)
    if plot:
        try:
            draw_reaction_dag(rlist, 'plots/reaction_network.png')
        except Exception as e:
            logger.warning(f'reaction_network.png render failed: {e}')
    input_reactants=reactants.intersection(reactants.symmetric_difference(products))
    logger.debug(f'Input reactants: {input_reactants}')
    for i in input_reactants:
        molecule_react_order.append((i,None))
    # once inputs are generated, any moleculed constructed only from inputs can be generated
    for R in rlist:
        if all([x in input_reactants for x in R.reactants.values()]):
            molecule_react_order.append((R.product,R))
            working_rlist.remove(R)
    
    while len(working_rlist)>0:
        molecules=[x[0] for x in molecule_react_order]
        for R in working_rlist:
            if all([x in molecules for x in R.reactants.values()]):
                molecule_react_order.append((R.product,R))
                molecules.append(R.product)
                working_rlist.remove(R)
            # else:
            #     logger.debug(f'Cannot place {R.product} since not all of {R.reactants.values()} are in the list yet')

    return molecule_react_order

def get_r(mname,RL:ReactionList):
    """Returns the Reaction that generates molecule mname as a product, if it exists, otherwise None.

    Args:
        mname (str): name of molecule to search
        RL (ReactionList): list of all Reactions

    Returns:
        Reaction or None: Reaction that generates mname, or None
    """
    if not RL: return None
    i=0
    while RL[i].product!=mname: 
        i+=1
        if i==len(RL):
            return None
    return RL[i]

def is_reactant(name:str,reaction_list:ReactionList,stage=reaction_stage.cure):
    """Tests to see if molecule 'name' appears as a reactant in any reaction in reaction_list.

    Args:
        name (str): name of molecule
        reaction_list (ReactionList): list of all Reactions
        stage (reaction_stage): stage of reactions to search, defaults to reaction_stage.cure

    Returns:
        bool: True if molecule is a reactant, False otherwise
    """
    reactants=[]
    for r in reaction_list:
        if r.stage==stage:
            for v in r.reactants.values():
                if not v in reactants:
                    reactants.append(v)
    return name in reactants
    
def product_sequence_resnames(R:Reaction,reactions:ReactionList):
    """Recursively generates sequence of product of R by traversing network of reactions; sequence is by definition an ordered list of monomer names.

    Args:
        R (Reaction): a Reaction
        reactions (ReactionList): list of all Reactions

    Returns:
        list: the sequence of R's product molecule as a list of monomer names
    """
    result=[]
    allProducts=[x.product for x in reactions]
    for rKey,rName in R.reactants.items():
        if rName in allProducts:
            result.extend(product_sequence_resnames(reactions[allProducts.index(rName)],reactions))
        else:
            result.extend([rName])
    return result

def molname_sequence_resnames(molname:str,reactions:ReactionList):
    """Determines monomer sequence of the molecule with name molname based on a traversal of the reaction network.

    Args:
        molname (str): name of molecule
        reactions (ReactionList): list of all Reactions

    Returns:
        list: monomer sequence of molecule as a list of monomer names
    """
    prodnames=[x.product for x in reactions]
    if molname in prodnames:
        R=reactions[prodnames.index(molname)]
        return product_sequence_resnames(R,reactions)
    else:
        return [molname]  # not a product in the reaction list?  must be a monomer

def reactant_resid_to_presid(R:Reaction,reactantName:str,resid:int,reactions:ReactionList):
    """Maps the resid of a reactant monomer to its inferred resid in the associated product in reaction R, by traversing the reaction network.

    Args:
        R (Reaction): a Reaction
        reactantName (str): the name of one of the reactants in R
        resid (int): the resid in that reactant we want mapped to the product
        reactions (ReactionList): list of all Reactions

    Returns:
        int: resid of this residue in the product of R
    """
    reactants=[x for x in R.reactants.values()]
    if reactantName in reactants:
        ridx=reactants.index(reactantName)
        logger.debug(f'reactantName {reactantName} {ridx} in {reactants} of {R.name}')
        S=0
        for i in range(ridx):
            seq=molname_sequence_resnames(reactants[i],reactions)
            logger.debug(f'{i} {seq} {S}')
            S+=len(seq)
        return S+resid
    else:
        return -1

def consumes_sacrificial_h(bondrec):
    """Returns True if forming this bond deletes one H from each of its atoms.

    True for a condensation, which is what htpolynet has always assumed, and false
    for an addition such as cyclotrimerization, where the ring closes with no atom
    lost.

    Args:
        bondrec (dict): one record from a Reaction's ``bonds``

    Returns:
        bool: whether to look for sacrificial hydrogens
    """
    return bool(bondrec.get('sacrificial_h', True))


def spanning_and_closing_bonds(R:Reaction):
    """Splits R's interreactant bonds into those that assemble the product and those that close a ring.

    A bond is *spanning* if it joins two reactants not yet connected by earlier bonds,
    and *closing* if they are already connected.  The distinction matters when a
    template is built: a spanning bond may position the incoming piece freely, while a
    closing bond's geometry is already fixed by the pieces it joins, so it cannot be
    placed and has to be brought together instead.

    Args:
        R (Reaction): a Reaction

    Returns:
        tuple: (spanning, closing) lists of indices into R.bonds
    """
    parent={}
    def find(x):
        parent.setdefault(x,x)
        while parent[x]!=x:
            parent[x]=parent[parent[x]]
            x=parent[x]
        return x
    spanning,closing=[],[]
    for i,bond in enumerate(R.bonds):
        a,b=[R.atoms[k]['reactant'] for k in bond['atoms']]
        if a==b:
            continue
        ra,rb=find(a),find(b)
        if ra==rb:
            closing.append(i)
        else:
            parent[ra]=rb
            spanning.append(i)
    return spanning,closing


def inter_reactant_bonds(R:Reaction):
    """Returns R's bonds that join two different reactants, as (reactant key, reactant key) pairs.

    Args:
        R (Reaction): a Reaction

    Returns:
        list: one (key, key) tuple per interreactant bond, in R.bonds order
    """
    pairs=[]
    for bond in R.bonds:
        i,j=[R.atoms[k]['reactant'] for k in bond['atoms']]
        if i!=j:
            pairs.append((i,j))
    return pairs


def is_ring_closing(R:Reaction):
    """Returns True if R's bonds close a cycle among its reactants.

    A cyclotrimerization does: three cyanate groups bond 1-2, 2-3 and 3-1, so the
    reactant graph is a triangle even though the new bonds alone form no cycle in the
    atom graph -- the ring closes through bonds each reactant already had.  A cure
    reaction joining two reactants with one bond is a tree and returns False.

    Args:
        R (Reaction): a Reaction

    Returns:
        bool: True if the reactant graph contains a cycle
    """
    edges=inter_reactant_bonds(R)
    nodes={k for e in edges for k in e}
    # a forest on n nodes has at most n-1 edges; more than that means a cycle
    return len(edges)>=len(nodes) and len(nodes)>0


def generate_product_name(R:Reaction):
    """Automatically generates the name of the product based on the names of the reactants and the bonds in R.

    Args:
        R (Reaction): a Reaction

    Returns:
        str: suggested name of product
    """
    pname=''
    for b in R.bonds:
        i,j=b['atoms']
        i_arec=R.atoms[i]
        j_arec=R.atoms[j]
        i_reKey=i_arec['reactant']
        i_reNm=R.reactants[i_reKey]
        i_aNm=i_arec['atom']
        j_reKey=j_arec['reactant']
        j_reNm=R.reactants[j_reKey]
        j_aNm=j_arec['atom']
        tbNm=f'{i_reNm}~{i_aNm}-{j_aNm}~{j_reNm}'
        if len(R.reactants)>1:
            if len(pname)==0:
                pname=tbNm
            else:
                pname+='---'+tbNm
    return pname

