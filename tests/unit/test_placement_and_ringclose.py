"""

.. module:: test_placement_and_ringclose
   :synopsis: positioning a reactant for an addition, and the prescribed restraint
              ladder that pulls a ring shut

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import numpy as np
import pandas as pd

from htpolynet.core.topology import Topology
from htpolynet.geometry.placement import min_distance, place_group, sphere_directions


class TestSphereDirections(unittest.TestCase):
    def test_unit_vectors(self):
        d = sphere_directions(32)
        self.assertEqual(d.shape, (32, 3))
        self.assertTrue(np.allclose(np.linalg.norm(d, axis=1), 1.0))

    def test_deterministic(self):
        self.assertTrue(np.array_equal(sphere_directions(16), sphere_directions(16)))

    def test_spans_the_sphere(self):
        d = sphere_directions(64)
        for axis in range(3):
            self.assertLess(d[:, axis].min(), -0.5)
            self.assertGreater(d[:, axis].max(), 0.5)


class TestPlaceGroup(unittest.TestCase):
    """Positions an incoming reactant for an addition, where there are no sacrificial
    hydrogens to align on."""

    def group(self):
        # a three-atom piece: the attaching atom at the origin and two others along x
        return np.array([[0.0, 0.0, 0.0], [0.15, 0.0, 0.0], [0.30, 0.0, 0.0]])

    def test_the_attaching_atom_lands_a_bond_length_from_the_anchor(self):
        anchor = np.zeros(3)
        fixed = np.array([[0.0, 0.0, -0.5]])
        placed, _ = place_group(self.group(), 0, anchor, fixed, 0.147)
        self.assertAlmostEqual(float(np.linalg.norm(placed[0] - anchor)), 0.147, places=6)

    def test_the_group_keeps_its_shape(self):
        placed, _ = place_group(self.group(), 0, np.zeros(3), np.array([[0.0, 0.0, -0.5]]), 0.147)
        before = np.linalg.norm(self.group()[2] - self.group()[0])
        self.assertAlmostEqual(float(np.linalg.norm(placed[2] - placed[0])), float(before), places=6)

    def test_it_goes_away_from_the_crowd(self):
        # a wall of fixed atoms on -z: the piece should end up on the +z side
        anchor = np.zeros(3)
        wall = np.array([[x * 0.1, y * 0.1, -0.3] for x in range(-3, 4) for y in range(-3, 4)])
        placed, clearance = place_group(self.group(), 0, anchor, wall, 0.147)
        self.assertGreater(placed[:, 2].mean(), 0.0)
        self.assertGreater(clearance, 0.2)

    def test_clearance_ignores_the_attaching_atom(self):
        # the attaching atom is a bond length from the anchor by construction, so
        # counting it would cap every clearance at the bond length
        placed, clearance = place_group(self.group(), 0, np.zeros(3),
                                        np.array([[0.0, 0.0, 0.0]]), 0.147)
        self.assertGreater(clearance, 0.147)

    def test_a_single_atom_group_has_no_clearance_to_report(self):
        placed, clearance = place_group(np.zeros((1, 3)), 0, np.zeros(3),
                                        np.array([[1.0, 0.0, 0.0]]), 0.147)
        self.assertEqual(clearance, np.inf)


def topology_with_restraint(initial=0.4, kb=300000.0):
    T = Topology()
    T.D['atoms'] = pd.DataFrame({'nr': [1, 2, 3], 'type': ['ca', 'nb', 'c3'],
                                 'charge': [0.1, -0.1, 0.0], 'resnr': [1, 2, 3],
                                 'residue': ['A', 'B', 'C'], 'atom': ['C1', 'N1', 'C2']})
    T.D['bonds'] = pd.DataFrame({'ai': [1], 'aj': [3], 'funct': [1],
                                 'c0': [np.nan], 'c1': [np.nan]})
    # get_bond_parameters reads this table whether or not it needs it
    T.D['bondtypes'] = pd.DataFrame({'i': ['ca'], 'j': ['c3'], 'func': [1],
                                     'b0': [0.14], 'kb': [250000.0]})
    pairs = pd.DataFrame({'ai': [1], 'aj': [2], 'initial_distance': [initial]})
    T.add_restraints(pairs, typ=6, kb=kb)
    return T, pairs


class TestPrescribedLadder(unittest.TestCase):
    """attenuate_bond_parameters reads a bond's reference values from the row it
    overwrites.  On a restraint carrying explicit parameters that compounds, which is
    why the ring-closing ladder prescribes each stage instead."""

    def restraint(self, T):
        b = T.D['bonds']
        row = b[(b['funct'] == 6)]
        return float(row['c0'].iloc[0]), float(row['c1'].iloc[0])

    def test_the_restraint_goes_in_at_the_initial_distance(self):
        T, _ = topology_with_restraint(initial=0.4)
        self.assertEqual(self.restraint(T), (0.4, 300000.0))

    def test_a_prescribed_stage_sets_both_outright(self):
        T, pairs = topology_with_restraint()
        T.set_restraint_parameters(pairs, 0.25, 300000.0)
        self.assertEqual(self.restraint(T), (0.25, 300000.0))

    def test_stiffness_does_not_decay_over_a_ladder(self):
        T, pairs = topology_with_restraint(initial=0.4)
        n = 8
        for stage in range(n):
            frac = (stage + 1) / n
            T.set_restraint_parameters(pairs, 0.4 + frac * (0.15 - 0.4), 300000.0)
        b0, kb = self.restraint(T)
        self.assertAlmostEqual(b0, 0.15)
        self.assertEqual(kb, 300000.0)

    def test_attenuating_a_restraint_is_what_compounds(self):
        # the bug this replaced: eight stages left the spring 400x too weak
        T, pairs = topology_with_restraint(initial=0.4)
        n = 8
        for stage in range(n):
            T.attenuate_bond_parameters(pairs, stage, n, minimum_distance=0.15,
                                        init_colname='initial_distance')
        _, kb = self.restraint(T)
        self.assertLess(kb, 3000.0)

    def test_per_pair_lengths(self):
        T, _ = topology_with_restraint()
        pairs = pd.DataFrame({'ai': [1], 'aj': [2], 'initial_distance': [0.4]})
        T.set_restraint_parameters(pairs, [0.31], 1000.0)
        self.assertEqual(self.restraint(T), (0.31, 1000.0))
