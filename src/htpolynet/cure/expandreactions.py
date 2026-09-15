"""Handles expansion of reactions based on chains, stereoisomers, and symmetry.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

from copy import deepcopy
from itertools import combinations, product

from ..core.molecule import Molecule, MoleculeList, MoleculeDict
from ..cure.reaction import reaction_stage, Reaction, ReactionList, generate_product_name, reactant_resid_to_presid

logger=logging.getLogger(__name__)

def bondchain_expand_reactions(molecules:MoleculeDict):
    """Handles generation of new reactions and molecular templates implied by any C-C bondchains.

    Note:
        Must be called after all grx attributes are set for all molecules.

    Args:
        molecules (MoleculeDict): dictionary of molecular templates constructed from explicit declarations

    Returns:
        tuple(ReactionList, MoleculeDict): the list of new reactions and dictionary of new molecules
    """
    extra_reactions:ReactionList=[]
    extra_molecules:MoleculeDict={}
    monomers:MoleculeList=[]
    dimer_lefts:MoleculeList=[]
    dimer_rights:MoleculeList=[]
    for mname,M in molecules.items():
        if len(M.sequence)==1 and len(M.chain_manager.chains)>0 and M.generator==None and M.parentname==M.name:
            monomers.append(M)
        elif len(M.sequence)==2:
            A=molecules[M.sequence[0]]
            if len(A.chain_manager.chains)>0:
                dimer_lefts.append(M)
            A=molecules[M.sequence[1]]
            if len(A.chain_manager.chains)>0:
                dimer_rights.append(M)
    for mon in monomers:
        cnms=[]
        for c in mon.chain_manager.chains:
            cnms.append([mon.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in c.idx_list])
        logger.debug(f'Monomer {mon.name} has {len(mon.chain_manager.chains)} 2-chains: {[x for x in mon.chain_manager.chains]} {cnms}')

    for dim in dimer_lefts:
        logger.debug(f'Dimer_left {dim.name} has sequence {dim.sequence}')
        logger.debug(f'-> chains: {[x for x in dim.chain_manager.chains]}')
        for cl in dim.chain_manager.chains:
            nl=[dim.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in cl.idx_list]
            logger.debug(f'  -> {nl}')
    for dim in dimer_rights:
        logger.debug(f'Dimer_right {dim.name} has sequence {dim.sequence}')
        logger.debug(f'-> chains: {[x for x in dim.chain_manager.chains]}')
        for cl in dim.chain_manager.chains:
            nl=[dim.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in cl.idx_list]
            logger.debug(f'  -> {nl}')

    # monomer head attacks dimer tail
    MD=product(monomers,dimer_lefts)
    for m,d in MD:
        for mb in m.chain_manager.chains:
            h_idx=mb.idx_list[0]
            h_name=m.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':h_idx})
            # by definition, the dimer must have one C-C bondchain of length 4
            D4=[]
            for dc in d.chain_manager.chains:
                if len(dc.idx_list)==4:
                    D4.append(dc.idx_list)
            for DC in D4:
                t_idx=DC[-1]
                t_name=d.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':t_idx})
                new_mname=f'{m.name}~{h_name}={t_name}~{d.name}'
                '''construct reaction'''
                R=Reaction()
                R.reactants={1:m.name, 2:d.name}
                R.atoms={'A':{'reactant':1,'resid':1,'atom':h_name,'z':1},
                        'B':{'reactant':2,'resid':1,'atom':t_name,'z':1}}
                R.bonds=[{'atoms':['A','B'],'order':1}]
                R.stage=reaction_stage.param
                R.name=new_mname.lower()
                R.product=new_mname
                newP=Molecule.New(R.product,R).set_sequence_from_moldict(molecules)
                newP.origin='unparameterized'
                extra_molecules[R.product]=newP
                logger.debug(f'monomer atom {m.name}_{h_name} will attack dimer atom {d.name}[{d.sequence[0]}1_{t_name}] -> {new_mname}:')
                for ln in str(R).split('\n'):
                    logger.debug(ln)
                extra_reactions.append(R)
    # dimer head attacks monomer tail
    MD=product(monomers,dimer_rights)
    for m,d in MD:
        for mb in m.chain_manager.chains:
            assert len(mb.idx_list)==2,f'monomer {m.name} has a bondchain that is not length-2 -- this is IMPOSSIBLE'
            t_idx=mb.idx_list[-1]
            t_name=m.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':t_idx})
            D4=[]
            for dc in d.chain_manager.chains:
                if len(dc.idx_list)==4:
                    D4.append(dc.idx_list)
            for DC in D4:
                h_idx=DC[0]
                h_name=d.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':h_idx})
                new_mname=f'{d.name}~{h_name}={t_name}~{m.name}'
                '''construct reaction'''
                R=Reaction()
                R.reactants={1:d.name, 2:m.name}
                R.atoms={'A':{'reactant':1,'resid':2,'atom':h_name,'z':1},
                        'B':{'reactant':2,'resid':1,'atom':t_name,'z':1}}
                R.bonds=[{'atoms':['A','B'],'order':1}]
                R.stage=reaction_stage.param
                new_rxnname=new_mname.lower()
                R.name=new_rxnname
                R.product=new_mname
                newP=Molecule.New(R.product,R).set_sequence_from_moldict(molecules)
                newP.origin='unparameterized'
                extra_molecules[R.product]=newP
                logger.debug(f'dimer atom {d.name}[{d.sequence[1]}2_{h_name}] will attack monomer atom {m.name}_{t_name}-> {new_mname}:')
                for ln in str(R).split('\n'):
                    logger.debug(ln)
                extra_reactions.append(R)

    DD=product(dimer_rights,dimer_lefts)
    for dr,dl in DD:
        ''' head of dr attacks tail of dl '''
        for cr,cl in product(dr.chain_manager.chains,dl.chain_manager.chains):
            if len(cr.idx_list)==4 and len(cl.idx_list)==4:
                h_idx=cr.idx_list[0]
                h_name=dr.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':h_idx})
                t_idx=cl.idx_list[-1]
                t_name=dl.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':t_idx})
                new_mname=f'{dr.name}~{h_name}={t_name}~{dl.name}'
                '''construct reaction'''
                R=Reaction()
                R.reactants={1:dr.name, 2:dl.name}
                R.atoms={'A':{'reactant':1,'resid':2,'atom':h_name,'z':1},
                            'B':{'reactant':2,'resid':1,'atom':t_name,'z':1}}
                R.bonds=[{'atoms':['A','B'],'order':1}]
                R.stage=reaction_stage.param
                R.product=new_mname
                R.name=R.product.lower()
                # newP=Molecule(name=R.product,generator=R)
                newP=Molecule.New(R.product,R).set_sequence_from_moldict(molecules)
                newP.origin='unparameterized'
                extra_molecules[R.product]=newP
                logger.debug(f'dimer atom {dr.name}-{dr.sequence[1]}2_{h_name} will attack dimer atom {dl.name}-{dl.sequence[0]}1_{t_name} -> {new_mname}:')
                for ln in str(R).split('\n'):
                    logger.debug(ln)
                extra_reactions.append(R)

    return extra_reactions, extra_molecules


def generate_stereo_reactions(RL:ReactionList,MD:MoleculeDict):
    """Scans the list of reactions and generates any additional reactions in which all possible stereoisomers are reactants.

    Args:
        RL (ReactionList): list of Reactions
        MD (MoleculeDict): dictionary of available Molecules

    Returns:
        int: number of new reactions/molecular products created
    """
    adds=0
    terminal_reactions=[]
    for R in RL:
        if R.stage not in [reaction_stage.param,reaction_stage.build]: continue
        Prod=MD[R.product]
        logger.debug(f'Stereos for {R.name} ({str(R.stage)})')
        reactant_stereoisomers={k:[r]+list(MD[r].stereoisomers.keys()) for k,r in R.reactants.items()}
        for k,v in reactant_stereoisomers.items():
            assert all([m in MD for m in v])  # all stereoisomers must be in the dict of molecules
        logger.debug(reactant_stereoisomers)
        reactant_keys=list(R.reactants.keys())
        isomer_lists=list(reactant_stereoisomers.values())
        combos=product(*isomer_lists)
        next(combos)
        sidx=1
        for c in combos:
            nR=deepcopy(R)
            nR.name=R.name+f'S-{sidx}'
            nR.product=R.product+f'S-{sidx}'
            nR.stage=reaction_stage.build
            nR.reactants={k:v for k,v in zip(reactant_keys,c)}
            MD[nR.product]=Molecule.New(nR.product,nR)
            adds+=1
            MD[nR.product].sequence=MD[R.product].sequence
            MD[nR.product].parentname=R.product
            Prod.stereoisomers[nR.product]=MD[nR.product]
            logger.debug(c)
            sidx+=1
            terminal_reactions.append(nR)
    RL.extend(terminal_reactions)
    return adds


def generate_symmetry_reactions(RL:ReactionList,MD:MoleculeDict):
    """Scans reaction list to generate any new reactions implied by symmetry-equivalent atoms.

    Args:
        RL (ReactionList): list of Reactions
        MD (MoleculeDict): dict of available molecules

    Returns:
        int: number of new reactions/molecular products created
    """
    jdx=1
    terminal_reactions=[]
    tail_adds=0
    for R in RL:
        if R.stage not in [reaction_stage.param,reaction_stage.cure,reaction_stage.cap,reaction_stage.repair]: continue
        Prod=MD[R.product]
        logger.debug(f'Symmetry versions for {R.name} ({str(R.stage)})\n{str(R)}')
        sra_by_reactant={k:MD[rname].symmetry_relateds for k,rname in R.reactants.items()}
        logger.debug(f'  sra_by_reactant: {sra_by_reactant}')
        atom_options=[]
        for atom_key,atom_rec in R.atoms.items():
            this_atom_options=[]
            art=atom_rec['reactant']
            target_atom_name=atom_rec['atom']
            if art in sra_by_reactant:
                logger.debug(f'art {art} {sra_by_reactant[art]}')
                if len(sra_by_reactant[art])==0: sra_by_reactant[art]=[[target_atom_name]]
                for symm_set in sra_by_reactant[art]:
                    if target_atom_name in symm_set:
                        for atom_name in symm_set:
                            this_atom_options.append([atom_key,atom_name])
            atom_options.append(this_atom_options)
        logger.debug(f'  atom options: {atom_options}')
        if len(R.reactants)>1:
            olist=list(product(*atom_options))
        else:
            olist=list(zip(*atom_options))
        idx=1
        R.symmetry_versions=olist
        logger.debug(f'olist {olist}')
        if len(olist)==1: continue
        for P in olist[1:]:
            newR=deepcopy(R)
            newR.name=R.name+f'-S{idx}'
            logger.debug(f'Permutation {P}:')
            for pp in P:
                atomKey,atomName=pp
                newR.atoms[atomKey]['atom']=atomName
            pname=generate_product_name(newR)
            if len(pname)==0:
                pname=R.product+f'-{idx}'
            newR.product=pname
            newR.stage=R.stage
            logger.debug(f'Primary:')
            for ln in str(newR).split('\n'): logger.debug(ln)
            terminal_reactions.append(newR)
            MD[newR.product]=Molecule(name=newR.product,generator=newR)
            MD[newR.product].origin='symmetry_product'
            MD[newR.product].set_sequence_from_moldict(MD)
            for rR in [x for x in RL if R.product in x.reactants.values()]:
                reactantKey=list(rR.reactants.keys())[list(rR.reactants.values()).index(R.product)]
                logger.debug(f'  New product {newR.product} must replace reactant {reactantKey}:{R.product} in {rR.name}')
                nooR=deepcopy(rR)
                nooR.stage=rR.stage
                nooR.name=rR.name+f'-{reactantKey}:S{jdx}'
                nooR.reactants[reactantKey]=newR.product
                for naK,naRec in {k:v for k,v in nooR.atoms.items() if v['reactant']==reactantKey}.items():
                    na_resid=naRec['resid']
                    na_name=naRec['atom']
                    for p in P:
                        oaK,oa_name=p
                        oaRec=R.atoms[oaK]
                        oa_reactatnName=R.reactants[oaRec['reactant']]
                        oa_resid=oaRec['resid']
                        oa_resid_in_o_product=reactant_resid_to_presid(R,oa_reactatnName,oa_resid,RL)
                        if na_resid == oa_resid_in_o_product:
                            nooR.atoms[naK]['resid']=oa_resid_in_o_product
                            nooR.atoms[naK]['atom']=oa_name
                noor_pname=generate_product_name(nooR)
                if noor_pname in MD: continue
                if len(noor_pname)==0:
                    noor_pname=rR.product+f'-{jdx}'
                nooR.product=noor_pname
                logger.debug(f'Secondary:')
                for ln in str(nooR).split('\n'): logger.debug(ln)
                jdx+=1
                RL.append(nooR)
                tail_adds+=1
                MD[nooR.product]=Molecule.New(nooR.product,nooR)
                MD[nooR.product].origin='symmetry_product'
                MD[nooR.product].set_sequence_from_moldict(MD)
            idx+=1
        logger.debug(f'Symmetry expansion of reaction {R.name} ends')

    RL.extend(terminal_reactions)
    return len(terminal_reactions)+tail_adds


def _second_shell_symmetry_siblings(M:Molecule, atom_name:str):
    """Names of the atoms symmetry-equivalent to atom_name that sit exactly two bonds from it in M.

    Args:
        M (Molecule): a generated monomer
        atom_name (str): name of a cure-reactive atom of M

    Returns:
        list: sorted atom names; empty if atom_name has no such siblings
    """
    candidates = set()
    for symm_set in M.symmetry_relateds:
        if atom_name in symm_set:
            candidates.update(x for x in symm_set if x != atom_name)
    if not candidates:
        return []
    TC = M.TopoCoord
    A = TC.Coordinates.A
    idx_of = dict(zip(A['atomName'], A['globalIdx'].astype(int)))
    if atom_name not in idx_of:
        return []
    bl = TC.Topology.bondlist
    a = idx_of[atom_name]
    first = set(bl.partners_of(a))
    second = {j for i in first for j in bl.partners_of(i)} - first - {a}
    return sorted(n for n in candidates if n in idx_of and idx_of[n] in second)


def sibling_expand_reactions(molecules:MoleculeDict, reactions:ReactionList):
    """Generates templates for cure bonds that form next to atoms of the same residue that have already reacted.

    A cure template is the product of one bond, so it shows every other
    reactive atom of the residue unreacted.  Where such an atom is close
    enough to the bonding atom that the template splice overwrites its charge,
    as the ring carbons of a triazine are, each later bond on the residue
    resets the charges of the atoms that reacted before it.  For every cure
    reaction whose monomer reactant has symmetry-equivalent reactive atoms two
    bonds from its bonding atom, this builds one product per non-empty set of
    those atoms already bonded to the same partner, with the bond of interest
    formed last so that it is the product's template bond.  The templates'
    sibling context (see TopoCoord.get_siblings) lets find_template pick them.

    Note:
        Must be called after the monomers and the cure products are generated.

    Args:
        molecules (MoleculeDict): all molecules generated so far
        reactions (ReactionList): all reactions, symmetry-expanded

    Returns:
        tuple(ReactionList, MoleculeDict): the new param-stage reactions, and their products
            in the order they must be generated
    """
    extra_reactions:ReactionList = []
    extra_molecules:MoleculeDict = {}
    cure = [R for R in reactions if R.stage == reaction_stage.cure and len(R.reactants) == 2 and len(R.bonds) == 1]
    for R in cure:
        bond_keys = R.bonds[0]['atoms']
        for s_key, o_key in (bond_keys, bond_keys[::-1]):
            s_rec, o_rec = R.atoms[s_key], R.atoms[o_key]
            if s_rec['reactant'] == o_rec['reactant']:
                continue
            s_name = R.reactants[s_rec['reactant']]
            o_name = R.reactants[o_rec['reactant']]
            if s_name not in molecules or len(molecules[s_name].sequence) != 1:
                continue
            siblings = _second_shell_symmetry_siblings(molecules[s_name], s_rec['atom'])
            if not siblings:
                continue
            order = list(R.reactants.keys())
            s_resid = sum(len(molecules[R.reactants[k]].sequence) for k in order[:order.index(s_rec['reactant'])]) + s_rec['resid']
            for n in range(1, len(siblings) + 1):
                for S in combinations(siblings, n):
                    # start from the cure product that bonds the first sibling to the same partner atom
                    base = [B for B in cure if B.reactants == R.reactants
                            and B.atoms[s_key]['atom'] == S[0] and B.atoms[o_key]['atom'] == o_rec['atom']
                            and B.atoms[s_key]['resid'] == s_rec['resid'] and B.atoms[o_key]['resid'] == o_rec['resid']]
                    if not base:
                        logger.debug(f'{R.name}: no cure product bonds {s_name} {S[0]} to {o_name} {o_rec["atom"]}; '
                                     f'no sibling template for {S}')
                        continue
                    prev = base[0].product
                    prev_s_resid = s_resid
                    o_first = s_key == bond_keys[1]
                    for site in list(S[1:]) + [s_rec['atom']]:
                        # Molecule.prepare_new_bonds reads a bond's first atom from the first
                        # reactant, so reactant order follows the cure bond's atom order; when
                        # the partner comes first, it also shifts the bonding residue's resid
                        s_atom = {'resid': prev_s_resid, 'atom': site, 'z': 1}
                        o_atom = {'resid': o_rec['resid'], 'atom': o_rec['atom'], 'z': 1}
                        newR = Reaction()
                        if o_first:
                            newR.reactants = {1: o_name, 2: prev}
                            newR.atoms = {'A': {'reactant': 1, **o_atom}, 'B': {'reactant': 2, **s_atom}}
                            prev_s_resid += len(molecules[o_name].sequence)
                        else:
                            newR.reactants = {1: prev, 2: o_name}
                            newR.atoms = {'A': {'reactant': 1, **s_atom}, 'B': {'reactant': 2, **o_atom}}
                        newR.bonds = [{'atoms': ['A', 'B'], 'order': R.bonds[0].get('order', 1)}]
                        newR.stage = reaction_stage.param
                        newR.product = generate_product_name(newR)
                        newR.name = f'sibling:{newR.product}'
                        prev = newR.product
                        if newR.product in molecules or newR.product in extra_molecules:
                            continue
                        P = Molecule.New(newR.product, newR)
                        P.set_sequence_from_moldict({**molecules, **extra_molecules})
                        P.origin = 'unparameterized'
                        extra_molecules[newR.product] = P
                        extra_reactions.append(newR)
                        logger.debug(f'sibling template {newR.product}: {s_name} {site} bonds {o_name} {o_rec["atom"]}')
    return extra_reactions, extra_molecules
