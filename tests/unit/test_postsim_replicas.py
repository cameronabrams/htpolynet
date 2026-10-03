"""

.. module:: test_postsim_replicas
   :synopsis: a ``deform`` stanza with ``replicas: N`` becomes N stages, each in its own directory

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

import yaml

from htpolynet.analysis.postsim import PostsimConfiguration, PostSimDeform, replica_subdirs


def read(stanzas):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'postsim.yaml')
        with open(f, 'w') as fh:
            yaml.safe_dump(stanzas, fh)
        return PostsimConfiguration.read(f)


class TestReplicas(unittest.TestCase):
    def test_one_replica_keeps_the_subdir(self):
        self.assertEqual(replica_subdirs({'subdir': 'postsim/deform-x'}), ['postsim/deform-x'])

    def test_replicas_expand_into_numbered_stages(self):
        cfg = read([{'equilibrate': {'ps': 10}},
                    {'deform': {'direction': 'x', 'subdir': 'postsim/deform-x', 'replicas': 3}}])
        self.assertEqual(len(cfg.stagelist), 4)
        deforms = [s for s in cfg.stagelist if isinstance(s, PostSimDeform)]
        self.assertEqual([s.params['subdir'] for s in deforms],
                         ['postsim/deform-x-r1', 'postsim/deform-x-r2', 'postsim/deform-x-r3'])

    def test_a_seed_gives_each_replica_its_own(self):
        cfg = read([{'deform': {'direction': 'y', 'subdir': 'postsim/deform-y', 'replicas': 2, 'seed': 100}}])
        self.assertEqual([s.params['seed'] for s in cfg.stagelist], [100, 101])

    def test_no_seed_is_left_to_gromacs(self):
        cfg = read([{'deform': {'direction': 'z', 'subdir': 'postsim/deform-z', 'replicas': 2}}])
        self.assertEqual([s.params['seed'] for s in cfg.stagelist], [None, None])

    def test_new_keys_are_accepted_not_ignored(self):
        cfg = read([{'deform': {'direction': 'x', 'fit_strain': [0.002, 0.03]}}])
        self.assertEqual(cfg.stagelist[0].params['fit_strain'], [0.002, 0.03])

    def test_zero_replicas_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'at least 1'):
            replica_subdirs({'subdir': 's', 'replicas': 0})


class TestSubdirsAreCheckedUpFront(unittest.TestCase):
    """A subdir go_to cannot reach is refused before any stage runs, not on its turn."""

    def setUp(self):
        from types import SimpleNamespace
        from unittest import mock
        import htpolynet.core.projectfilesystem as pfs
        patcher = mock.patch.object(pfs, '_PFS_', SimpleNamespace(projSubPaths={'postsim': '/p/postsim', 'systems': '/p/systems'}))
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_replica_subdirs_are_reachable(self):
        cfg = read([{'deform': {'direction': 'x', 'subdir': 'postsim/deform-x', 'replicas': 3}}])
        cfg.check_subdirs()

    def test_a_nested_subdir_is_refused_with_every_bad_one_named(self):
        # htpolynet-sweep, 2026-10-03: postsim/grid-A/deform-x-r1 died on its first pull
        cfg = read([{'equilibrate': {'ps': 10}},
                    {'deform': {'direction': 'x', 'subdir': 'postsim/grid-A/deform-x'}},
                    {'deform': {'direction': 'y', 'subdir': 'postsim/grid-A/deform-y'}}])
        with self.assertRaises(ValueError) as cm:
            cfg.check_subdirs()
        self.assertIn('postsim/grid-A/deform-x', str(cm.exception))
        self.assertIn('postsim/grid-A/deform-y', str(cm.exception))
