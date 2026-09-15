"""

.. module:: test_repair_charge_neutrality
   :synopsis: postcure repair leaves every molecule it touched neutral, not just the system

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
import pandas as pd

from htpolynet.core.topology import Topology
from htpolynet.repair.topology_surgery import neutralize_touched_fragments


class FakeTC:
    def __init__(self, charges, bonds):
        self.Topology = Topology()
        self.Topology.D['atoms'] = pd.DataFrame({'nr': list(range(1, len(charges) + 1)),
                                                 'charge': [float(q) for q in charges]})
        self.Topology.D['bonds'] = pd.DataFrame(bonds, columns=['ai', 'aj'])

    def q(self, i):
        a = self.Topology.D['atoms']
        return float(a.loc[a['nr'] == i, 'charge'].iloc[0])

    def fragment_charge(self, idx):
        return sum(self.q(i) for i in idx)


class TestNeutralizeTouchedFragments(unittest.TestCase):
    def two_molecules(self):
        # molecule A (1-2-3) and molecule B (4-5-6), charged +0.24 and -0.24:
        # neutral in total, which is all a system-wide rebalance could see
        return FakeTC([0.10, 0.10, 0.04, -0.10, -0.10, -0.04], [(1, 2), (2, 3), (4, 5), (5, 6)])

    def test_each_touched_molecule_ends_neutral(self):
        TC = self.two_molecules()
        stats = neutralize_touched_fragments(TC, {2, 3, 5})
        self.assertAlmostEqual(TC.fragment_charge([1, 2, 3]), 0.0, places=9)
        self.assertAlmostEqual(TC.fragment_charge([4, 5, 6]), 0.0, places=9)
        self.assertEqual(stats['n_fragments'], 2)
        self.assertAlmostEqual(stats['max_excess'], 0.24)

    def test_only_touched_atoms_move(self):
        TC = self.two_molecules()
        neutralize_touched_fragments(TC, {2, 3, 5})
        self.assertAlmostEqual(TC.q(1), 0.10)
        self.assertAlmostEqual(TC.q(2), 0.10 - 0.12)
        self.assertAlmostEqual(TC.q(3), 0.04 - 0.12)
        self.assertAlmostEqual(TC.q(4), -0.10)
        self.assertAlmostEqual(TC.q(5), -0.10 + 0.24)

    def test_untouched_molecule_is_left_alone_and_system_stays_neutral(self):
        # B is charged but the repair never touched it; its charge is not
        # moved onto A's atoms by the per-molecule step, only by the final
        # system rebalance, which keeps the total at zero for Ewald
        TC = FakeTC([0.10, 0.10, 0.04, -0.10, -0.10, -0.04, 0.0], [(1, 2), (2, 3), (4, 5), (5, 6)])
        with self.assertLogs('htpolynet.repair.topology_surgery', level='INFO'):
            neutralize_touched_fragments(TC, {1, 2, 3})
        for i in (4, 5, 6):
            self.assertIn(TC.q(i), (-0.10, -0.04))
        self.assertAlmostEqual(TC.Topology.total_charge(), 0.0, places=9)

    def test_already_neutral_is_a_no_op(self):
        TC = FakeTC([0.2, -0.2, 0.1, -0.1], [(1, 2), (3, 4)])
        stats = neutralize_touched_fragments(TC, {1, 3})
        self.assertEqual(stats['n_fragments'], 0)
        self.assertEqual([TC.q(i) for i in (1, 2, 3, 4)], [0.2, -0.2, 0.1, -0.1])

    def test_nothing_touched(self):
        TC = self.two_molecules()
        self.assertEqual(neutralize_touched_fragments(TC, set())['n_fragments'], 0)
        self.assertAlmostEqual(TC.q(1), 0.10)

    def test_a_touched_atom_with_no_bonds_is_its_own_molecule(self):
        TC = FakeTC([0.3, -0.1, -0.2], [(2, 3)])
        neutralize_touched_fragments(TC, {1, 2})
        self.assertAlmostEqual(TC.q(1), 0.0, places=9)
        self.assertAlmostEqual(TC.fragment_charge([2, 3]), 0.0, places=9)

    def test_stale_indices_are_ignored(self):
        TC = FakeTC([0.10, 0.10, -0.20, 0.3], [(1, 2), (2, 3)])
        neutralize_touched_fragments(TC, {2, 4, 99})
        self.assertAlmostEqual(TC.fragment_charge([1, 2, 3]), 0.0, places=9)
        self.assertAlmostEqual(TC.q(4), 0.0, places=9)
