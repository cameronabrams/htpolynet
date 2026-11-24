"""

.. module:: test_mol2_partial_bonds
   :synopsis: a merged template keeps every reactant's mol2 bonds, and the bond check sees a short table

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

htpolynet-sweep, 2026-11-24: example 2's GMAS-4..15 were written with 24 of their
78 bonds.  Their reactant GM1S-k, a stereoisomer built on GM1's topology, had no mol2
bond table, and merging it in kept only the other reactant's bonds.  The check that
reported it compared the two lists with zip(), so a table short only at the end passed.
"""
import os
import shutil
import tempfile
import unittest

import pandas as pd

from htpolynet.core.topocoord import TopoCoord

_FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')


def _tc(with_mol2):
    TC = TopoCoord(topfilename=os.path.join(_FIX, 'config1.top'),
                   grofilename=os.path.join(_FIX, 'config1.gro'))
    if with_mol2:
        b = TC.Topology.D['bonds'][['ai', 'aj']].reset_index(drop=True)
        TC.Topology.D['mol2_bonds'] = pd.DataFrame({'bondIdx': range(1, len(b) + 1),
                                                    'ai': b['ai'].astype(int),
                                                    'aj': b['aj'].astype(int), 'order': '1'})
    return TC


def _pairs(df):
    return {(min(a, b), max(a, b)) for a, b in zip(df['ai'].astype(int), df['aj'].astype(int))}


class TestMergeKeepsEveryMol2Bond(unittest.TestCase):
    def check_merge(self, first_has, second_has):
        A, B = _tc(first_has), _tc(second_has)
        n = len(A.Topology.D['bonds'])
        with self.assertLogs('htpolynet.core.topology', level='WARNING') as logs:
            A.merge(B)
        T = A.Topology.D
        self.assertEqual(len(T['mol2_bonds']), 2 * n)
        self.assertEqual(_pairs(T['mol2_bonds']), _pairs(T['bonds']))
        self.assertEqual(list(T['mol2_bonds']['bondIdx']), list(range(1, 2 * n + 1)))
        self.assertIn('no mol2 bond table', '\n'.join(logs.output))

    def test_incoming_molecule_without_a_table(self):
        """The GMAS-4 case: the merged-in reactant had no table and was dropped."""
        self.check_merge(True, False)

    def test_molecule_being_built_without_a_table(self):
        """The reverse used to adopt the incoming table without shifting its indices."""
        self.check_merge(False, True)

    def test_both_with_tables_is_quiet(self):
        A, B = _tc(True), _tc(True)
        n = len(A.Topology.D['bonds'])
        with self.assertNoLogs('htpolynet.core.topology', level='WARNING'):
            A.merge(B)
        self.assertEqual(len(A.Topology.D['mol2_bonds']), 2 * n)

    def test_neither_with_a_table_invents_none(self):
        A, B = _tc(False), _tc(False)
        A.merge(B)
        self.assertNotIn('mol2_bonds', A.Topology.D)


class TestBondSourceCheck(unittest.TestCase):
    def test_complete_table_passes(self):
        self.assertTrue(_tc(True).Topology.bond_source_check())

    def test_table_short_only_at_the_end_fails(self):
        """zip() over a prefix used to call this consistent."""
        T = _tc(True).Topology
        mb = T.D['mol2_bonds'].sort_values(by=['ai', 'aj'])
        T.D['mol2_bonds'] = mb.iloc[:-3].reset_index(drop=True)
        with self.assertLogs('htpolynet.core.topology', level='ERROR') as logs:
            self.assertFalse(T.bond_source_check())
        self.assertIn('3 topology bonds are missing from the mol2', '\n'.join(logs.output))

    def test_reversed_pair_order_is_the_same_bond(self):
        T = _tc(True).Topology
        mb = T.D['mol2_bonds']
        T.D['mol2_bonds'] = mb.assign(ai=mb['aj'], aj=mb['ai'])
        self.assertTrue(T.bond_source_check())


class TestStereoisomerReadsItsOwnBonds(unittest.TestCase):
    """A topology loaded with the molecule's own mol2 carries that mol2's bond table."""

    def test_top_gro_plus_mol2_keeps_the_table(self):
        src = _tc(True)
        with tempfile.TemporaryDirectory() as td:
            m = os.path.join(td, 'X.mol2')
            src.write_mol2(m, molname='X')
            TC = TopoCoord(topfilename=os.path.join(_FIX, 'config1.top'),
                           grofilename=os.path.join(_FIX, 'config1.gro'), mol2filename=m)
        self.assertIn('mol2_bonds', TC.Topology.D)
        self.assertEqual(len(TC.Topology.D['mol2_bonds']), len(TC.Topology.D['bonds']))
        self.assertTrue(TC.Topology.bond_source_check())


if __name__ == '__main__':
    unittest.main()
