"""

.. module:: test_checkout_shadowing
   :synopsis: a file of your own that replaces a packaged one says so

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
from unittest import mock

from htpolynet.core import projectfilesystem as pfs

LOG = 'htpolynet.core.projectfilesystem'


class TestShadowingIsAnnounced(unittest.TestCase):
    """A stale copy of a packaged file keeps its old version's settings for as long as
    it stays on the search path, and nothing used to say which file a run got.  An
    npt.mdp copied from 2.6.2 reverted LINCS to order 4 across four upgrades."""

    def note(self, packaged_exists, source='user library'):
        root = mock.MagicMock()
        root.joinpath.return_value.exists.return_value = packaged_exists
        lib = mock.MagicMock()
        lib._root = root
        with mock.patch.object(pfs, '_SYSTEM_LIBRARY_', lib):
            with self.assertLogs(LOG, level='WARNING') as cm:
                pfs._note_if_shadowing('mdp/npt.mdp', source)
                logger = __import__('logging').getLogger(LOG)
                logger.warning('sentinel')
            return '\n'.join(cm.output)

    def test_it_says_so_when_the_packaged_file_exists(self):
        msg = self.note(packaged_exists=True)
        self.assertIn('mdp/npt.mdp', msg)
        self.assertIn('replaces the copy', msg)
        self.assertIn('after any upgrade', msg)

    def test_it_is_silent_for_a_file_htpolynet_does_not_ship(self):
        msg = self.note(packaged_exists=False)
        self.assertNotIn('replaces the copy', msg)

    def test_it_names_which_source_the_file_came_from(self):
        self.assertIn('user cache', self.note(True, source='user cache'))

    def test_it_never_raises(self):
        # checkout must not fail because the warning could not be formed
        with mock.patch.object(pfs, '_SYSTEM_LIBRARY_', object()):
            pfs._note_if_shadowing('mdp/npt.mdp', 'user library')
