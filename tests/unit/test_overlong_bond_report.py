"""

.. module:: test_overlong_bond_report
   :synopsis: a build says so when it hands on a bond that is not a bond any more

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import numpy as np
import pandas as pd

from htpolynet.core.runtime import Runtime

LOG = 'htpolynet.core.runtime'


def runtime(bonds, positions, names=None):
    """A Runtime stub carrying just the bonds and coordinates the report reads.

    `positions` is {globalIdx: (x, y, z)} in nm.  Lengths come from the same helper
    the real TopoCoord uses, so this exercises the report rather than a reimplementation
    of it.
    """
    idx = sorted(positions)
    A = pd.DataFrame({'globalIdx': idx,
                      'posX': [positions[i][0] for i in idx],
                      'posY': [positions[i][1] for i in idx],
                      'posZ': [positions[i][2] for i in idx],
                      'resNum': [names[i][0] if names else 1 for i in idx],
                      'resName': [names[i][1] if names else 'BDC' for i in idx],
                      'atomName': [names[i][2] if names else f'C{i}' for i in idx]})

    class TC:
        Coordinates = type('C', (), {'A': A})()
        Topology = type('T', (), {'D': {'bonds': pd.DataFrame(bonds, columns=['ai', 'aj'])}})()

        @staticmethod
        def add_length_attribute(bdf, attr_name='length'):
            p = A.set_index('globalIdx')
            bdf[attr_name] = [
                float(np.linalg.norm(p.loc[int(b.ai), ['posX', 'posY', 'posZ']].to_numpy(float)
                                     - p.loc[int(b.aj), ['posX', 'posY', 'posZ']].to_numpy(float)))
                for b in bdf.itertuples()]

    return type('R', (), {'TopoCoord': TC()})()


class TestOverlongBondReport(unittest.TestCase):
    def test_an_ordinary_structure_is_an_info_line(self):
        R = runtime([(1, 2), (2, 3)], {1: (0, 0, 0), 2: (0.15, 0, 0), 3: (0.30, 0, 0)})
        with self.assertLogs(LOG, level='INFO') as cm:
            out = Runtime._report_overlong_bonds(R, 'in final')
        self.assertTrue(out.empty)
        self.assertIn('no bond among 2', '\n'.join(cm.output))

    def test_a_stretched_bond_is_a_warning_that_names_it(self):
        # the measured case: a monomer's own backbone at 3.20 A, r0 1.516
        R = runtime([(1, 2)], {1: (0, 0, 0), 2: (0.320, 0, 0)},
                    names={1: (138, 'BDC', 'C5'), 2: (138, 'BDC', 'C8')})
        with self.assertLogs(LOG, level='WARNING') as cm:
            out = Runtime._report_overlong_bonds(R, 'in final')
        msg = '\n'.join(cm.output)
        self.assertEqual(out.shape[0], 1)
        self.assertIn('C5(1) of BDC138-C8(2) of BDC138 at 3.20 A', msg)

    def test_the_worst_come_first_and_the_list_is_capped(self):
        pos = {i: (0.0, 0.0, 0.0) for i in range(1, 9)}
        for n, i in enumerate(range(2, 9), start=1):
            pos[i] = (0.21 + 0.01 * n, 0.0, 0.0)
        R = runtime([(1, i) for i in range(2, 9)], pos)
        with self.assertLogs(LOG, level='WARNING') as cm:
            out = Runtime._report_overlong_bonds(R, 'in final', show=2)
        self.assertEqual(out.shape[0], 7)
        self.assertGreater(out.iloc[0]['length'], out.iloc[1]['length'])
        # capped at two named, so the third-worst is not in the message
        self.assertEqual('\n'.join(cm.output).count(' at '), 2)

    def test_a_system_with_no_bonds_is_silent(self):
        R = runtime([], {1: (0, 0, 0)})
        self.assertTrue(Runtime._report_overlong_bonds(R, 'in final').empty)
