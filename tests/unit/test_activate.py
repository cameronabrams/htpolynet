"""Tests deriving a monomer's active form from its real SMILES (`activate`).

The claim under test is strong: for every shipped example monomer that is drawn
in an active form, writing the real monomer plus `activate` gives the *same*
structure -- the same atoms, in the same order, embedded to the same
coordinates -- so a config can switch to the real form without changing a build.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import filecmp
import os
import shutil

import pytest
import yaml

from htpolynet.external.smiles_input import (activate_mol, activate_smiles, activated_bonds,
                                             activation_cap_reactions, materialize_smiles_inputs)

Chem = pytest.importorskip('rdkit.Chem')
from rdkit.Chem import AllChem  # noqa: E402

DEPOT = os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'htpolynet', 'resources', 'example_depot')

# Each example monomer the depot shipped in its active form, until the examples were
# switched to real form: (example file, constituent, the active SMILES it shipped, the
# real SMILES plus `activate` it ships now).  The real SMILES is written in the atom
# order of the active one, and labels on atoms the active SMILES left unmapped (the
# epoxide O of HIE, the isocyanate N of IPD) are ones it lacked.
REAL = [
    ('1-polystyrene.yaml', 'STY', 'c1ccccc1[CH2:1][CH3:2]',
     'c1ccccc1[CH:1]=[CH2:2]', [[1, 2]]),
    ('2-bisgma-styrene-thermoset.yaml', 'STY', 'c1ccccc1[CH2:1][CH3:2]',
     'c1ccccc1[CH:1]=[CH2:2]', [[1, 2]]),
    ('2-bisgma-styrene-thermoset.yaml', 'HIE', '[CH3:2][CH:1]([CH3])C(=O)OC[CH:3](O)[CH3:4]',
     '[CH2:2]=[C:1]([CH3])C(=O)OC[CH:3]1[O:5][CH2:4]1', [[2, 1], [5, 4]]),
    ('3-pacm-dgeba-epoxy-thermoset.yaml', 'DGE',
     'CC(C)(c1ccc(OC[CH:3]([OH:5])[CH3:1])cc1)c1ccc(OC[CH:4]([OH:6])[CH3:2])cc1',
     'CC(C)(c1ccc(OC[CH:3]2[O:5][CH2:1]2)cc1)c1ccc(OC[CH:4]2[O:6][CH2:2]2)cc1', [[5, 1], [6, 2]]),
    ('4-dfda-fde-epoxy-thermoset.yaml', 'FDE',
     'N(Cc1occc1)(C[CH:3]([OH:5])[CH3:1])C[CH:4]([OH:6])[CH3:2]',
     'N(Cc1occc1)(C[CH:3]1[O:5][CH2:1]1)C[CH:4]1[O:6][CH2:2]1', [[5, 1], [6, 2]]),
    ('5-htpb-ipdi.yaml', 'IPD', 'CC1(CN[CH:2]=O)C[CH:3](N[CH:1]=O)CC(C)(C)C1',
     'CC1(C[N:12]=[C:2]=O)C[CH:3]([N:11]=[C:1]=O)CC(C)(C)C1', [[12, 2], [11, 1]]),
]
IDS = [f'{e[:1]}-{n}' for e, n, _, _, _ in REAL]


def _constituent(example, name):
    with open(os.path.join(DEPOT, example)) as f:
        return yaml.safe_load(f)['constituents'][name]


def _structure(mol):
    """What smiles_input._rdkit_path embeds, as a molblock (it then hands this to obabel)."""
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    mol = Chem.AddHs(mol)
    assert AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) == 0
    AllChem.UFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)


@pytest.mark.parametrize('example,name,active,real,activate', REAL, ids=IDS)
def test_the_example_ships_the_real_form(example, name, active, real, activate):
    spec = _constituent(example, name)
    assert spec['smiles'] == real
    assert spec['activate'] == activate


@pytest.mark.parametrize('example,name,active,real,activate', REAL, ids=IDS)
def test_real_form_reproduces_the_former_active_form_exactly(example, name, active, real, activate):
    assert _structure(activate_mol(name, real, activate)) == _structure(active)


@pytest.mark.skipif(shutil.which('obabel') is None, reason='obabel not on PATH')
@pytest.mark.parametrize('example,name,active,real,activate', REAL, ids=IDS)
def test_real_form_writes_the_former_mol2_byte_for_byte(tmp_path, example, name, active, real, activate):
    spec = _constituent(example, name)
    former = {k: v for k, v in spec.items() if k != 'activate'}
    materialize_smiles_inputs({name: dict(former, smiles=active)}, inputs_dir=str(tmp_path / 'active'))
    materialize_smiles_inputs({name: spec}, inputs_dir=str(tmp_path / 'real'))
    assert filecmp.cmp(tmp_path / 'active' / f'{name}.mol2', tmp_path / 'real' / f'{name}.mol2', shallow=False)


def test_activated_smiles_keeps_the_input_atom_order():
    # the epoxide's ring bond is written last, so a canonical or DFS rewrite would move atoms
    real = '[CH2:1]1[O:3][CH:2]1CC'
    mol = activate_mol('X', real, [[3, 1]])
    assert [(a.GetSymbol(), a.GetAtomMapNum()) for a in mol.GetAtoms()][:3] == [('C', 1), ('O', 3), ('C', 2)]
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(Chem.MolFromSmiles('[CH3:1][CH:2]([OH:3])CC'))


@pytest.mark.parametrize('real,activate,expected', [
    ('[CH2:1]=[CH:2]C', [[1, 2]], 1),                 # vinyl: double -> single
    ('C[N:1]=[C:2]=O', [[1, 2]], 1),                  # isocyanate N=C -> single
    ('C[C:1]#[N:2]', [[1, 2]], 2),                    # nitrile: triple -> double
    ('[CH2:1]1[O:2][CH:3]1C', [[2, 1]], None),        # epoxide: ring opened
])
def test_each_pair_loses_one_bond_order_and_gains_two_hydrogens(real, activate, expected):
    before = Chem.AddHs(Chem.MolFromSmiles(real))
    after = Chem.AddHs(Chem.MolFromSmiles(activate_smiles('X', real, activate)))
    assert after.GetNumAtoms() == before.GetNumAtoms() + 2
    m2i = {a.GetAtomMapNum(): a.GetIdx() for a in after.GetAtoms() if a.GetAtomMapNum()}
    a, b = activate[0]
    bond = after.GetBondBetweenAtoms(m2i[a], m2i[b])
    if expected is None:
        assert bond is None
    else:
        assert bond.GetBondTypeAsDouble() == expected


@pytest.mark.parametrize('real,activate,match', [
    ('[CH:1]=[CH2:2]', [[1, 3]], 'not in its SMILES'),
    ('[CH2:1]=C[CH3:2]', [[1, 2]], 'not bonded'),
    ('[c:1]1[c:2]cccc1', [[1, 2]], 'aromatic'),
    ('[CH3:1][CH3:2]', [[1, 2]], 'acyclic single bond'),
    ('[CH:1]=[CH2:2]', [[1, 1]], 'with itself'),
    ('[CH:1]=[CH2:2]', [[1]], 'pair of atom-map labels'),
])
def test_impossible_activations_are_refused(real, activate, match):
    with pytest.raises(ValueError, match=match):
        activate_smiles('X', real, activate)


def test_opening_one_bond_twice_is_refused():
    with pytest.raises(ValueError, match='more than once'):
        activate_smiles('X', '[CH2:1]1[O:2][CH:3]1C', [[2, 1], [1, 2]])


def test_activated_bonds_reports_the_order_a_cap_restores():
    assert activated_bonds('X', 'C[N:1]=[C:2]=O.[CH2:3]1[O:4]C1', [[1, 2], [4, 3]]) == [(1, 2, 2), (4, 3, 1)]


# --- generated cap reactions ------------------------------------------------

def test_symmetric_epoxides_get_one_cap():
    spec = {'smiles': 'N(C[CH:3]1[O:5][CH2:1]1)C[CH:4]1[O:6][CH2:2]1',
            'reactive_atoms': {1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4', 5: 'O1', 6: 'O2'},
            'symmetry_equivalent_atoms': [['C1', 'C2'], ['C3', 'C4'], ['O1', 'O2']],
            'activate': [[5, 1], [6, 2]]}
    gen = activation_cap_reactions({'EPX': spec}, [])
    assert [{a['atom'] for a in g['atoms'].values()} for g in gen] == [{'O1', 'C1'}]


def test_a_hand_written_cap_on_the_symmetric_partner_suppresses_the_generated_one():
    spec = {'smiles': 'N(C[CH:3]1[O:5][CH2:1]1)C[CH:4]1[O:6][CH2:2]1',
            'reactive_atoms': {1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4', 5: 'O1', 6: 'O2'},
            'symmetry_equivalent_atoms': [['C1', 'C2'], ['C3', 'C4'], ['O1', 'O2']],
            'activate': [[5, 1], [6, 2]]}
    hand = {'stage': 'cap', 'reactants': {1: 'EPX'}, 'product': 'EPXC',
            'atoms': {'A': {'reactant': 1, 'resid': 1, 'atom': 'O2', 'z': 1},
                      'B': {'reactant': 1, 'resid': 1, 'atom': 'C2', 'z': 1}}}
    assert activation_cap_reactions({'EPX': spec}, [hand]) == []


# the hand-written caps the examples shipped before they were switched to real form,
# which the generated ones replace
SHIPPED_CAPS = {
    'STY': ('styCC', {'C1': 1, 'C2': 1}, 2),
    'HIE': ('hieCC', {'C1': 1, 'C2': 1}, 2),
    'DGE': ('Oxirane-formation', {'O1': 1, 'C1': 1}, 1),
    'FDE': ('Oxirane-formation', {'O1': 1, 'C1': 1}, 1),
}


@pytest.mark.parametrize('example,name', [
    ('1-polystyrene.yaml', 'STY'), ('2-bisgma-styrene-thermoset.yaml', 'STY'),
    ('2-bisgma-styrene-thermoset.yaml', 'HIE'), ('3-pacm-dgeba-epoxy-thermoset.yaml', 'DGE'),
    ('4-dfda-fde-epoxy-thermoset.yaml', 'FDE'),
])
def test_generated_cap_matches_the_one_the_example_used_to_write(example, name):
    """The configs no longer write these caps, so what they generate must be what they wrote."""
    spec = _constituent(example, name)
    with open(os.path.join(DEPOT, example)) as f:
        written = yaml.safe_load(f)['reactions']
    assert not [r for r in written if r['stage'] == 'cap' and list(r['reactants'].values()) == [name]]
    gen = activation_cap_reactions({name: spec}, written)
    _, atoms, order = SHIPPED_CAPS[name]
    assert len(gen) == 1
    g = gen[0]
    assert {a['atom']: a['z'] for a in g['atoms'].values()} == atoms
    assert g['bonds'][0]['order'] == order
    assert g['stage'] == 'cap' and g['reactants'] == {1: name}


def test_example_5_isocyanate_gets_no_cap_as_before():
    spec = _constituent('5-htpb-ipdi.yaml', 'IPD')
    assert activation_cap_reactions({'IPD': spec}, []) == []


def test_a_hand_written_cap_suppresses_the_generated_one():
    spec = _constituent('1-polystyrene.yaml', 'STY')
    hand = {'stage': 'cap', 'reactants': {1: 'STY'}, 'product': 'STYCC',
            'atoms': {'A': {'reactant': 1, 'resid': 1, 'atom': 'C1', 'z': 1},
                      'B': {'reactant': 1, 'resid': 1, 'atom': 'C2', 'z': 1}}}
    assert activation_cap_reactions({'STY': spec}, [hand]) == []


def test_cap_activated_false_writes_none():
    spec = {'smiles': '[CH2:1]=[CH:2]C', 'reactive_atoms': {1: 'C1', 2: 'C2'}, 'activate': [[1, 2]],
            'cap_activated': False}
    assert activation_cap_reactions({'PRO': spec}, []) == []


def test_an_unnamed_atom_gets_no_cap():
    # the isocyanate N is labeled only so `activate` can name it
    spec = {'smiles': 'C[N:9]=[C:1]=O', 'reactive_atoms': {1: 'C1'}, 'activate': [[9, 1]]}
    assert activation_cap_reactions({'ISO': spec}, []) == []


def test_symmetric_pairs_named_alike_get_one_cap():
    spec = {'smiles': '[CH2:1]=[CH:2]CC[CH:3]=[CH2:4]', 'reactive_atoms': {1: 'C1', 2: 'C2', 3: 'C2', 4: 'C1'},
            'activate': [[1, 2], [4, 3]]}
    assert len(activation_cap_reactions({'HDE': spec}, [])) == 1


def test_a_product_name_collision_is_refused():
    spec = {'smiles': '[CH2:1]=[CH:2]C', 'reactive_atoms': {1: 'C1', 2: 'C2'}, 'activate': [[1, 2]]}
    with pytest.raises(ValueError, match='already in use'):
        activation_cap_reactions({'PRO': spec, 'PRO-capped-C1-C2': {}}, [])


def test_configuration_appends_the_generated_caps():
    from htpolynet.core.configuration import Configuration
    C = Configuration.read(os.path.join(DEPOT, '1-polystyrene.yaml'))
    caps = [r for r in C.reaction_specs if str(r['stage']) == 'cap']
    assert [r['product'] for r in caps] == ['STY-capped-C1-C2']
    assert caps[0]['bonds'][0]['order'] == 2
