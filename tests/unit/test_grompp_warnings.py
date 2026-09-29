"""

.. module:: test_grompp_warnings
   :synopsis: grompp's warnings reach the log instead of being discarded

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

from htpolynet.external.gromacs import report_grompp_warnings

LOG = 'htpolynet.external.gromacs'

BERENDSEN = """Generated 10 of the 10 non-bonded parameter combinations

WARNING 1 [file npt.mdp]:
  The Berendsen barostat does not generate any strictly correct ensemble,
  and should not be used for new production simulations (in our opinion).
  We recommend using the C-rescale barostat instead.


There was 1 warning
"""

TWO = """WARNING 1 [file shear.mdp, line 41]:
  An off-diagonal box element has deform set while compressibility > 0 for
  the same component of another box vector.

WARNING 2 [file shear.mdp]:
  The Berendsen barostat does not generate any strictly correct ensemble.
"""


class TestWarningsAreSurfaced(unittest.TestCase):
    """htpolynet passes -maxwarn 4 and `run` discards output on success, so a Gromacs
    warning on a run that proceeds reached neither the user nor the log.  That is how the
    Berendsen deprecation went unseen on every constant-pressure stage for as long as
    Gromacs has emitted it."""

    def test_a_warning_is_logged(self):
        with self.assertLogs(LOG, level='WARNING') as cm:
            found = report_grompp_warnings(BERENDSEN, '', 'npt')
        self.assertEqual(len(found), 1)
        self.assertIn('Berendsen barostat', '\n'.join(cm.output))
        self.assertIn('npt.mdp', '\n'.join(cm.output))

    def test_the_whole_body_is_kept_not_just_the_heading(self):
        found = report_grompp_warnings(BERENDSEN, '', 'npt')
        self.assertIn('C-rescale', found[0])

    def test_every_warning_is_reported(self):
        with self.assertLogs(LOG, level='WARNING') as cm:
            found = report_grompp_warnings(TWO, '', 'shear')
        self.assertEqual(len(found), 2)
        self.assertIn('off-diagonal box element', cm.output[0])

    def test_a_clean_run_says_nothing(self):
        import logging
        logger = logging.getLogger(LOG)
        with self.assertLogs(LOG, level='WARNING') as cm:
            self.assertEqual(report_grompp_warnings('all good\n', '', 'npt'), [])
            logger.warning('sentinel')
        self.assertEqual(len(cm.output), 1)

    def test_it_reads_stderr_too(self):
        self.assertEqual(len(report_grompp_warnings('', BERENDSEN, 'npt')), 1)
