"""

.. module:: test_sibling_templates
   :synopsis: cure templates that know which neighbouring reactive atoms have already reacted

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import pandas as pd

from htpolynet.core.bondtemplate import BondTemplate
from htpolynet.core.molecule import Molecule
from htpolynet.core.topocoord import TopoCoord
from htpolynet.cure.expandreactions import _second_shell_symmetry_siblings, sibling_expand_reactions
from htpolynet.cure.reaction import Reaction, reaction_stage
from htpolynet.geometry.bondlist import Bondlist


def bt(siblings=None, names=('O1', 'C2')):
    return BondTemplate(list(names), ['BPA', 'TAZ'], False, 1, [[], []], [[], []],
                        [None, None], [None, None], siblings=siblings)


class TestBondTemplateSiblings(unittest.TestCase):
    def test_default_is_no_siblings(self):
        self.assertEqual(bt().siblings, [[], []])

    def test_bare_template_still_matches_a_ring_with_reacted_carbons(self):
        self.assertTrue(bt().matches(bt([[], [('C1', 'BPA')]])))

    def test_template_showing_a_reacted_carbon_needs_one_there(self):
        self.assertFalse(bt([[], [('C1', 'BPA')]]).matches(bt()))
        self.assertFalse(bt([[], [('C1', 'BPA')]]).matches(bt([[], [('C3', 'BPA')]])))
        self.assertTrue(bt([[], [('C1', 'BPA')]]).matches(bt([[], [('C1', 'BPA'), ('C3', 'BPA')]])))

    def test_more_reacted_context_is_more_specific(self):
        self.assertLess(bt().bystander_count(), bt([[], [('C1', 'BPA')]]).bystander_count())
        self.assertLess(bt([[], [('C1', 'BPA')]]).bystander_count(),
                        bt([[], [('C1', 'BPA'), ('C3', 'BPA')]]).bystander_count())

    def test_reverse_swaps_sides(self):
        t = bt([[], [('C1', 'BPA')]])
        t.reverse()
        self.assertEqual(t.siblings, [[('C1', 'BPA')], []])
        self.assertTrue(t.matches_reverse_of(bt([[], [('C1', 'BPA')]])))

    def test_order_of_siblings_does_not_matter(self):
        self.assertEqual(bt([[], [('C3', 'BPA'), ('C1', 'BPA')]]), bt([[], [('C1', 'BPA'), ('C3', 'BPA')]]))

    def test_siblings_are_part_of_equality(self):
        self.assertNotEqual(bt(), bt([[], [('C1', 'BPA')]]))


class FakeTopology:
    def __init__(self, bonds):
        self.bondlist = Bondlist.fromDataFrame(pd.DataFrame(bonds, columns=['ai', 'aj']))


class FakeCoordinates:
    def __init__(self, atoms):
        self.A = pd.DataFrame([{'globalIdx': i + 1, 'resNum': r, 'resName': rn, 'atomName': an}
                               for i, (r, rn, an) in enumerate(atoms)])


class FakeTC:
    def __init__(self, atoms, bonds):
        self.Coordinates = FakeCoordinates(atoms)
        self.Topology = FakeTopology(bonds)


RING = [(1, 'TAZ', 'C1'), (1, 'TAZ', 'N1'), (1, 'TAZ', 'C2'), (1, 'TAZ', 'N2'), (1, 'TAZ', 'C3'), (1, 'TAZ', 'N3')]
RING_BONDS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]


class TestGetSiblings(unittest.TestCase):
    def system(self):
        # BPA 2 already on C1 (atoms 7-8), BPA 3 bonding C2 now (atoms 9-10), and
        # a residue X bonded to N1, which is first-shell to C2
        atoms = RING + [(2, 'BPA', 'O1'), (2, 'BPA', 'C5'), (3, 'BPA', 'O1'), (3, 'BPA', 'C5'), (4, 'X', 'Q1')]
        bonds = RING_BONDS + [(7, 1), (7, 8), (9, 3), (9, 10), (11, 2)]
        return FakeTC(atoms, bonds)

    def test_reacted_ring_carbon_two_bonds_away_is_a_sibling(self):
        self.assertEqual(TopoCoord.get_siblings(self.system(), [9, 3]), [[], [('C1', 'BPA')]])

    def test_a_reacted_direct_neighbour_is_not(self):
        # N1 carries residue X but is bonded to C2 directly: the chain case
        self.assertNotIn(('N1', 'X'), TopoCoord.get_siblings(self.system(), [9, 3])[1])

    def test_nothing_on_an_unreacted_ring(self):
        TC = FakeTC(RING + [(3, 'BPA', 'O1')], RING_BONDS + [(7, 3)])
        self.assertEqual(TopoCoord.get_siblings(TC, [7, 3]), [[], []])


class FakeMonomer(Molecule):
    def __init__(self, name, atoms=None, bonds=None, symmetry=None):
        super().__init__(name=name)
        self.sequence = [name]
        self.symmetry_relateds = symmetry or []
        if atoms:
            self.TopoCoord = FakeTC(atoms, bonds)


def taz():
    return FakeMonomer('TAZ', RING + [(1, 'TAZ', 'H'), (1, 'TAZ', 'H1'), (1, 'TAZ', 'H2')],
                       RING_BONDS + [(1, 7), (3, 8), (5, 9)], [['C1', 'C2', 'C3'], ['N1', 'N2', 'N3']])


class TestSecondShellSymmetrySiblings(unittest.TestCase):
    def test_ring_carbons(self):
        self.assertEqual(_second_shell_symmetry_siblings(taz(), 'C1'), ['C2', 'C3'])

    def test_not_symmetry_related_is_not_a_sibling(self):
        self.assertEqual(_second_shell_symmetry_siblings(taz(), 'N1'), ['N2', 'N3'])
        M = taz(); M.symmetry_relateds = []
        self.assertEqual(_second_shell_symmetry_siblings(M, 'C1'), [])

    def test_adjacent_equivalents_are_not_siblings(self):
        # a vinyl-like pair: symmetry-equivalent but directly bonded
        M = FakeMonomer('VIN', [(1, 'VIN', 'C1'), (1, 'VIN', 'C2')], [(1, 2)], [['C1', 'C2']])
        self.assertEqual(_second_shell_symmetry_siblings(M, 'C1'), [])


class TestSiblingExpandReactions(unittest.TestCase):
    def setUp(self):
        self.molecules = {'BPA': FakeMonomer('BPA'), 'TAZ': taz()}
        self.reactions = []
        for c in ('C1', 'C2', 'C3'):
            R = Reaction()
            R.name = f'etherify-{c}'
            R.stage = reaction_stage.cure
            R.reactants = {1: 'BPA', 2: 'TAZ'}
            R.atoms = {'A': {'reactant': 1, 'resid': 1, 'atom': 'O1', 'z': 1},
                       'B': {'reactant': 2, 'resid': 1, 'atom': c, 'z': 1}}
            R.bonds = [{'atoms': ['A', 'B'], 'order': 1}]
            R.product = f'BPA~O1-{c}~TAZ'
            self.reactions.append(R)
            P = Molecule(name=R.product, generator=R)
            P.set_sequence_from_moldict(self.molecules)
            self.molecules[R.product] = P

    def test_one_template_per_bonding_carbon_and_set_of_reacted_siblings(self):
        new_reactions, new_molecules = sibling_expand_reactions(self.molecules, self.reactions)
        # 3 bonding carbons x {one of two others, both others}
        self.assertEqual(len(new_molecules), 9)
        self.assertEqual(len(new_reactions), 9)
        self.assertTrue(all(R.stage == reaction_stage.param for R in new_reactions))

    def test_the_bond_of_interest_is_formed_last(self):
        _, new_molecules = sibling_expand_reactions(self.molecules, self.reactions)
        full = new_molecules['BPA~O1-C3~BPA~O1-C2~BPA~O1-C1~TAZ']
        self.assertEqual(full.sequence, ['BPA', 'BPA', 'BPA', 'TAZ'])
        bond_atoms = full.generator.atoms
        self.assertEqual((bond_atoms['B']['resid'], bond_atoms['B']['atom']), (3, 'C3'))
        self.assertEqual(full.generator.reactants[2], 'BPA~O1-C2~BPA~O1-C1~TAZ')
        # the TAZ that C3 names really is residue 3 of that reactant
        self.assertEqual(new_molecules['BPA~O1-C2~BPA~O1-C1~TAZ'].sequence[2], 'TAZ')

    def test_the_template_bond_reads_the_same_way_round_as_the_cure_bond(self):
        # the cure bond is BPA O1 - TAZ C; a template listing C first would only
        # match through the reversed-bond path
        new_reactions, _ = sibling_expand_reactions(self.molecules, self.reactions)
        for R in new_reactions:
            first = R.atoms[R.bonds[0]['atoms'][0]]
            self.assertEqual(R.reactants[first['reactant']], 'BPA')
            self.assertEqual(first['atom'], 'O1')
            # and it comes from the first reactant, which prepare_new_bonds assumes
            self.assertEqual(first['reactant'], 1)

    def test_reactants_come_before_the_products_built_from_them(self):
        _, new_molecules = sibling_expand_reactions(self.molecules, self.reactions)
        order = list(new_molecules)
        for name, M in new_molecules.items():
            parent = M.generator.reactants[1]
            if parent in new_molecules:
                self.assertLess(order.index(parent), order.index(name))

    def test_nothing_for_a_residue_without_siblings(self):
        self.molecules['TAZ'].symmetry_relateds = []
        self.assertEqual(sibling_expand_reactions(self.molecules, self.reactions), ([], {}))
