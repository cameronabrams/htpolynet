"""

.. module:: test_mol2_bond_fallback
   :synopsis: MOL2 writes for molecules loaded from a topology, with no mol2 bond table

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

from htpolynet.core.topocoord import TopoCoord

_FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')


def _bond_lines(path):
    with open(path) as f:
        txt = f.read()
    if '@<TRIPOS>BOND' not in txt:
        return []
    return txt.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0].strip().splitlines()


class TestMol2BondFallback(unittest.TestCase):
    """Symmetry siblings (e.g. GMAS-4) are built by loading the parent's top+gro and
    never get a mol2 bond table.  Writing one used to produce a bondless MOL2 and a
    warning; the 2026-09-13 example sweep logged twelve such writes in example 2."""

    def setUp(self):
        self.TC = TopoCoord(topfilename=os.path.join(_FIX, 'config1.top'),
                            grofilename=os.path.join(_FIX, 'config1.gro'))

    def test_fixture_really_has_no_mol2_bond_table(self):
        # guards the premise: if this fixture ever grows a mol2 table the test
        # below would pass without exercising the fallback at all
        self.assertNotIn('mol2_bonds', self.TC.Topology.D)
        self.assertGreater(len(self.TC.Topology.D['bonds']), 0)

    def test_topology_bonds_are_written_without_warning(self):
        n = len(self.TC.Topology.D['bonds'])
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, 'x.mol2')
            with self.assertNoLogs('htpolynet.io.mol2', level='WARNING'):
                self.TC.write_mol2(p, molname='X')
            lines = _bond_lines(p)
        self.assertEqual(len(lines), n)

    def test_fallback_bonds_are_single_and_reference_the_topology_pairs(self):
        b = self.TC.Topology.D['bonds'].iloc[0]
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, 'x.mol2')
            self.TC.write_mol2(p, molname='X')
            first = _bond_lines(p)[0].split()
        # bondIdx, ai, aj, order
        self.assertEqual((int(first[1]), int(first[2])), (int(b['ai']), int(b['aj'])))
        self.assertEqual(first[3], '1')
