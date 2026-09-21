"""

.. module:: test_multibody_reactions
   :synopsis: reactions with more than two reactants, and the ring-closing case the
              pairwise bond search must refuse

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
from types import SimpleNamespace

import pandas as pd

from htpolynet.core.molecule import Molecule
from htpolynet.cure.curecontroller import CureController
from htpolynet.cure.reaction import Reaction, inter_reactant_bonds, is_ring_closing, reaction_stage
from htpolynet.geometry.bondlist import Bondlist

LOG = 'htpolynet.cure.curecontroller'


def reaction(reactants, atoms, bonds, stage=reaction_stage.cure, name='r', product='P'):
    """atoms: {key: (reactant, resid, atomname)}; bonds: [(keyA, keyB), ...]"""
    R = Reaction()
    R.name, R.stage, R.product = name, stage, product
    R.reactants = dict(reactants)
    R.atoms = {k: {'reactant': r, 'resid': i, 'atom': a, 'z': 1} for k, (r, i, a) in atoms.items()}
    R.bonds = [{'atoms': list(b), 'order': 1} for b in bonds]
    return R


def etherify():
    """The two-reactant, one-bond cure reaction of example 6."""
    return reaction({1: 'BPA', 2: 'TAZ'}, {'A': (1, 1, 'O1'), 'B': (2, 1, 'C1')}, [('A', 'B')])


def cyclotrimerize():
    """Three cyanate groups closing a triazine: 1-2, 2-3, 3-1."""
    return reaction({1: 'BCY', 2: 'BCY', 3: 'BCY'},
                    {'C1': (1, 1, 'C1'), 'N1': (1, 1, 'N1'),
                     'C2': (2, 1, 'C1'), 'N2': (2, 1, 'N1'),
                     'C3': (3, 1, 'C1'), 'N3': (3, 1, 'N1')},
                    [('N1', 'C2'), ('N2', 'C3'), ('N3', 'C1')],
                    name='cyclotrimerize', product='BCY3')


class TestRingClosingDetection(unittest.TestCase):
    def test_a_pairwise_cure_reaction_is_not_ring_closing(self):
        self.assertFalse(is_ring_closing(etherify()))
        self.assertEqual(inter_reactant_bonds(etherify()), [(1, 2)])

    def test_cyclotrimerization_is(self):
        R = cyclotrimerize()
        self.assertEqual(inter_reactant_bonds(R), [(1, 2), (2, 3), (3, 1)])
        self.assertTrue(is_ring_closing(R))

    def test_a_chain_of_three_reactants_is_not(self):
        # 1-2 and 2-3 is a path, not a ring: two bonds, three reactants
        R = reaction({1: 'A', 2: 'B', 3: 'C'},
                     {'a': (1, 1, 'X'), 'b1': (2, 1, 'Y'), 'b2': (2, 1, 'Z'), 'c': (3, 1, 'W')},
                     [('a', 'b1'), ('b2', 'c')])
        self.assertFalse(is_ring_closing(R))

    def test_an_intraresidue_bond_is_not_counted_as_connectivity(self):
        # the oxirane-formation capping reaction bonds two atoms of one reactant
        R = reaction({1: 'DGE'}, {'a': (1, 1, 'C1'), 'b': (1, 1, 'O1')}, [('a', 'b')],
                     stage=reaction_stage.cap)
        self.assertEqual(inter_reactant_bonds(R), [])
        self.assertFalse(is_ring_closing(R))

    def test_two_bonds_between_the_same_pair_close_a_ring(self):
        # both arms of one molecule into one partner: a ring through both molecules
        R = reaction({1: 'BCY', 2: 'TAZ'},
                     {'a1': (1, 1, 'C1'), 'a2': (1, 1, 'C2'), 'b1': (2, 1, 'N1'), 'b2': (2, 1, 'N2')},
                     [('a1', 'b1'), ('a2', 'b2')])
        self.assertTrue(is_ring_closing(R))


class TestPairwiseSearchRefusesRingClosing(unittest.TestCase):
    """The pairwise search would form the three bonds separately, in different
    iterations, which does not build a ring.  It must skip and say so."""

    def searcher(self, RL):
        cc = CureController({})
        adf = pd.DataFrame({'globalIdx': [1], 'resNum': [1], 'resName': ['BCY'], 'atomName': ['C1'],
                            'z': [1], 'nreactions': [0], 'reactantName': ['BCY'], 'molecule': [1]})
        TC = SimpleNamespace(
            gro_DataFrame=lambda name: adf if name == 'atoms' else None,
            files={'gro': 'nonexistent.gro'},
            linkcell_initialize=lambda *a, **k: None,
            linkcell_cleanup=lambda *a, **k: None,
            Topology=SimpleNamespace(bondlist=Bondlist()),
        )
        cc.state.current_radius = 0.5
        return cc, TC

    def test_it_warns_and_forms_nothing(self):
        cc, TC = self.searcher([cyclotrimerize()])
        with self.assertLogs(LOG, level='WARNING') as cm:
            bdf = cc._searchbonds(TC, [cyclotrimerize()], {})
        text = '\n'.join(cm.output)
        self.assertIn('closes a ring among 3 reactants', text)
        self.assertIn('skipping it', text)
        self.assertTrue(bdf.empty)

    def test_a_pairwise_reaction_in_the_same_list_is_not_skipped(self):
        # only the ring-closing one is named; the other is searched as usual, which
        # here finds nothing because the fake system holds a single atom
        cc, TC = self.searcher([etherify(), cyclotrimerize()])
        with self.assertLogs(LOG, level='WARNING') as cm:
            cc._searchbonds(TC, [etherify(), cyclotrimerize()],
                            {'BPA': FakeMol('BPA'), 'TAZ': FakeMol('TAZ')})
        text = '\n'.join(cm.output)
        self.assertIn('cyclotrimerize', text)
        self.assertNotIn('"r" closes', text)


class FakeMol:
    def __init__(self, *resnames):
        self.sequence = list(resnames)

    def get_resname(self, resid):
        return self.sequence[resid - 1]


class FakeCoordinates:
    def __init__(self, atoms):
        self.A = pd.DataFrame([{'globalIdx': i + 1, 'resNum': r, 'resName': rn, 'atomName': an}
                               for i, (r, rn, an) in enumerate(atoms)])


class FakeTC:
    def __init__(self, atoms, bonds):
        self.Coordinates = FakeCoordinates(atoms)
        self.Topology = type('T', (), {'bondlist': Bondlist.fromDataFrame(pd.DataFrame(bonds, columns=['ai', 'aj']))})()

    def get_gro_attribute_by_attributes(self, attr, d):
        A = self.Coordinates.A
        m = (A['resNum'] == d['resNum']) & (A['atomName'] == d['atomName'])
        return A.loc[m, attr].iloc[0]

    def get_bystanders(self, idx):
        return [[], []], [[], []], [[], []], [[], []]

    def get_oneaways(self, idx, chain_manager=None):
        return [None, None], [None, None], [None, None], [None, None]

    def get_siblings(self, idx):
        return [[], []]


class TestResidOffsetsForThreeReactants(unittest.TestCase):
    """A bond between the first and third reactant must count the second's residues."""

    def molecule(self, R, molecules):
        M = Molecule(name=R.product, generator=R)
        # one atom per residue is enough: the offsets are what is under test
        atoms, bonds = [], []
        resnum = 1
        names = ('C1', 'C2', 'N1', 'N2', 'O1', 'O2')
        for key in R.reactants:
            for rn in molecules[R.reactants[key]].sequence:
                first = len(atoms) + 1
                for an in names:
                    atoms.append((resnum, rn, an))
                for k in range(first + 1, len(atoms) + 1):
                    bonds.append((first, k))
                resnum += 1
        M.TopoCoord = FakeTC(atoms, bonds)
        M.prepare_new_bonds(available_molecules=molecules)
        return M

    def test_three_single_residue_reactants(self):
        molecules = {'BCY': FakeMol('BCY')}
        M = self.molecule(cyclotrimerize(), molecules)
        self.assertEqual([rb.resids for rb in M.reaction_bonds], [[1, 2], [2, 3], [3, 1]])
        self.assertEqual([bt.resnames for bt in M.bond_templates],
                         [['BCY', 'BCY'], ['BCY', 'BCY'], ['BCY', 'BCY']])

    def test_a_multi_residue_reactant_in_the_middle_shifts_the_third(self):
        # reactant 2 is a two-residue oligomer, so reactant 3's residue is 4, not 3
        R = reaction({1: 'A', 2: 'AB', 3: 'C'},
                     {'a': (1, 1, 'C1'), 'b1': (2, 1, 'N1'), 'b2': (2, 2, 'C1'), 'c': (3, 1, 'N1')},
                     [('a', 'b1'), ('b2', 'c')])
        molecules = {'A': FakeMol('A'), 'AB': FakeMol('A', 'B'), 'C': FakeMol('C')}
        M = self.molecule(R, molecules)
        self.assertEqual([rb.resids for rb in M.reaction_bonds], [[1, 2], [3, 4]])
        self.assertEqual([bt.resnames for bt in M.bond_templates], [['A', 'A'], ['B', 'C']])

    def test_two_reactants_are_unchanged(self):
        molecules = {'BPA': FakeMol('BPA'), 'TAZ': FakeMol('TAZ')}
        M = self.molecule(etherify(), molecules)
        self.assertEqual([rb.resids for rb in M.reaction_bonds], [[1, 2]])
        self.assertEqual([bt.resnames for bt in M.bond_templates], [['BPA', 'TAZ']])
