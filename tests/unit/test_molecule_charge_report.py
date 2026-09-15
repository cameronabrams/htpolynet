"""

.. module:: test_molecule_charge_report
   :synopsis: build output reports molecules that carry net charge

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import pandas as pd

from htpolynet.core.runtime import Runtime
from htpolynet.core.topology import Topology

LOG = 'htpolynet.core.runtime'


def topology(charges, bonds):
    T = Topology()
    T.D['atoms'] = pd.DataFrame({'nr': list(range(1, len(charges) + 1)), 'charge': [float(q) for q in charges]})
    T.D['bonds'] = pd.DataFrame(bonds, columns=['ai', 'aj'])
    return T


class FakeRuntime:
    def __init__(self, T):
        self.TopoCoord = type('TC', (), {'Topology': T})()


class TestMoleculeCharges(unittest.TestCase):
    def test_one_entry_per_bonded_molecule_largest_first(self):
        T = topology([0.1, -0.1, 0.24, 0.0, -0.24], [(1, 2), (3, 4), (4, 5)])
        result = T.molecule_charges()
        self.assertEqual([n for n, _ in result], [3, 2])
        self.assertAlmostEqual(result[0][1], 0.0)
        self.assertAlmostEqual(result[1][1], 0.0)

    def test_unbonded_atom_is_its_own_molecule(self):
        T = topology([0.5, -0.5], [])
        self.assertEqual(sorted(q for _, q in T.molecule_charges()), [-0.5, 0.5])

    def test_empty(self):
        self.assertEqual(Topology().molecule_charges(), [])


class TestReportMoleculeCharges(unittest.TestCase):
    def test_neutral_system_is_an_info_line(self):
        R = FakeRuntime(topology([0.1, -0.1, 0.2, -0.2], [(1, 2), (3, 4)]))
        with self.assertLogs(LOG, level='INFO') as cm:
            charged = Runtime._report_molecule_charges(R, 'after cure')
        self.assertEqual(charged, [])
        self.assertIn('all 2 molecules neutral', '\n'.join(cm.output))
        self.assertNotIn('WARNING', '\n'.join(cm.output))

    def test_ions_that_sum_to_zero_are_a_warning(self):
        # the cyanate-cap repair bug: every molecule charged, system total zero
        R = FakeRuntime(topology([0.21, 0.0, -0.27, 0.06], [(1, 2), (3, 4)]))
        with self.assertLogs(LOG, level='WARNING') as cm:
            charged = Runtime._report_molecule_charges(R, 'after repair')
        text = '\n'.join(cm.output)
        self.assertEqual(len(charged), 2)
        self.assertIn('2 of 2 molecules', text)
        self.assertIn('-0.210 e', text)

    def test_below_tolerance_is_neutral(self):
        R = FakeRuntime(topology([0.005, 0.0, -0.005, 0.0], [(1, 2), (3, 4)]))
        with self.assertLogs(LOG, level='INFO'):
            self.assertEqual(Runtime._report_molecule_charges(R, 'in final'), [])
