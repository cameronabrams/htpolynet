"""

.. module:: test_plots_diag
   :synopsis: ``htpolynet plots diag`` reads both CURE and ring-cure diagnostic logs

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

import pandas as pd

from htpolynet.analysis.plot import diagnostics_graphs

# verbatim lines from the example 6 reference build (6ca6406), trimmed to three iterations
RING_LOG = """\
2026-10-03 11:19:19,971 htpolynet.core.runtime.my_logger INFO>  Ring cure begins: 720 reactive BDC group(s), (('N1', 'C1'), ('N2', 'C2')), template BDC3
2026-10-03 11:28:39,000 htpolynet.cure.ringcontroller.record INFO> Iteration 1: 97 ring(s), 291 of 720 groups consumed (conversion 0.404)
2026-10-03 11:33:00,000 htpolynet.cure.ringcontroller.record INFO> Iteration 2: 30 ring(s), 381 of 720 groups consumed (conversion 0.529)
2026-10-03 12:14:25,420 htpolynet.cure.ringcontroller.record INFO> Iteration 3: 2 ring(s), 702 of 720 groups consumed (conversion 0.975)
"""

CURE_LOG = """\
2026-05-29 04:26:48,821 htpolynet.core.runtime.my_logger INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
2026-05-29 04:29:10,323 htpolynet.cure.curecontroller.do_iter INFO> Iteration 1 current conversion 0.222 or 160 bonds
2026-05-29 04:31:10,323 htpolynet.cure.curecontroller.do_iter INFO> Iteration 2 current conversion 0.411 or 296 bonds
"""


class TestPlotsDiag(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as d:
            log = os.path.join(d, 'diagnostics.log')
            with open(log, 'w') as f:
                f.write(text)
            diagnostics_graphs([log], os.path.join(d, 'cure_info.png'))
            self.assertTrue(os.path.exists(os.path.join(d, 'cure_info.png')))
            return pd.read_csv(os.path.join(d, 'diagnostics.csv'), sep=' ')

    def test_ring_cure_log(self):
        df = self.parse(RING_LOG)
        self.assertEqual(df['iter'].tolist(), [0, 1, 2, 3])
        self.assertEqual(df['nbonds'].tolist(), [0, 291, 381, 702])
        self.assertEqual(df['conv'].tolist(), [0.0, 0.404, 0.529, 0.975])
        self.assertAlmostEqual(df['elapsed'].iloc[-1], (55 * 60 + 5.449) / 3600, places=4)

    def test_cure_log_still_parses(self):
        df = self.parse(CURE_LOG)
        self.assertEqual(df['iter'].tolist(), [0, 1, 2])
        self.assertEqual(df['nbonds'].tolist(), [0, 160, 296])

    def test_a_log_with_no_cure_says_so(self):
        with self.assertRaisesRegex(ValueError, 'no CURE or ring-cure progress lines'):
            self.parse('2026-05-29 04:26:48,821 htpolynet.core.runtime INFO> nothing here\n')
