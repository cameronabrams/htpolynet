"""

.. module:: test_ringcontroller
   :synopsis: the loop of a ring-closing cure: what it searches, when it widens, when
              it stops, and how it counts conversion

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import numpy as np
import pandas as pd

from htpolynet.cure.ringcontroller import RingController, RingCureState

from .test_triplesearch import BOX, triangle

LOG = 'htpolynet.cure.ringcontroller'


class TestConversionIsCountedInGroups(unittest.TestCase):
    """A cyclotrimerization consumes three groups per ring.  Conversion in groups is
    what experiments report; the pre-formed-triazine route has to infer it from a bond
    count instead."""

    def test_no_rings_is_zero(self):
        self.assertEqual(RingCureState(total_groups=1440).conversion, 0.0)

    def test_each_ring_consumes_three(self):
        s = RingCureState(total_groups=300)
        RingController(state=s).record([object(), object()])
        self.assertEqual(s.rings, 2)
        self.assertEqual(s.groups_consumed, 6)
        self.assertAlmostEqual(s.conversion, 0.02)

    def test_all_groups_consumed_is_one(self):
        self.assertEqual(RingCureState(groups_consumed=300, total_groups=300).conversion, 1.0)

    def test_an_empty_system_does_not_divide_by_zero(self):
        self.assertEqual(RingCureState().conversion, 0.0)


class TestStopping(unittest.TestCase):
    def controller(self, **d):
        c = RingController(d)
        c.setup(total_groups=300, max_radius=1.0)
        return c

    def test_runs_while_below_target(self):
        c = self.controller(desired_conversion=0.9)
        c.state.groups_consumed = 150
        self.assertFalse(c.is_cured())

    def test_stops_at_the_target(self):
        c = self.controller(desired_conversion=0.5)
        c.state.groups_consumed = 150
        self.assertTrue(c.is_cured())

    def test_stops_when_iterations_run_out(self):
        c = self.controller(desired_conversion=0.99, max_iterations=3)
        c.state.iter = 4
        self.assertTrue(c.is_cured())


class TestWidening(unittest.TestCase):
    def controller(self):
        c = RingController({'search_radius': 0.5, 'radial_increment': 0.1})
        c.setup(total_groups=300, max_radius=0.7)
        return c

    def test_starts_at_the_configured_radius(self):
        self.assertAlmostEqual(self.controller().radius, 0.5)

    def test_widens_by_one_increment(self):
        c = self.controller()
        self.assertFalse(c.widen())
        self.assertAlmostEqual(c.radius, 0.6)

    def test_stops_widening_at_the_maximum(self):
        c = self.controller()
        c.widen()
        c.widen()
        self.assertAlmostEqual(c.radius, 0.7)
        self.assertTrue(c.widen())

    def test_the_maximum_defaults_to_what_setup_was_given(self):
        c = RingController({})
        c.setup(total_groups=10, max_radius=1.25)
        self.assertEqual(c.dicts['max_radius'], 1.25)

    def test_a_configured_maximum_wins(self):
        c = RingController({'max_radius': 0.8})
        c.setup(total_groups=10, max_radius=1.25)
        self.assertEqual(c.dicts['max_radius'], 0.8)


class TestSearch(unittest.TestCase):
    def positions_from(self, adf, pos):
        return pos

    def test_it_finds_and_packs_a_ring(self):
        adf, pos = triangle(side=0.25)
        c = RingController({'search_radius': 0.6})
        c.setup(total_groups=3, max_radius=1.0)
        with self.assertLogs(LOG, level='INFO') as cm:
            sites, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(sites), 3)
        self.assertEqual(len(chosen), 1)
        self.assertIn('1 accepted', '\n'.join(cm.output))

    def test_nothing_within_the_radius(self):
        adf, pos = triangle(side=2.0)
        c = RingController({'search_radius': 0.3})
        c.setup(total_groups=3, max_radius=0.4)
        with self.assertLogs(LOG, level='INFO'):
            sites, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(chosen, [])

    def test_a_ring_per_iteration_limit_is_honored(self):
        # two separate triangles, but only one ring wanted this iteration
        a1, p1 = triangle(side=0.2, resnums=(1, 2, 3), molecules=(1, 2, 3))
        a2, p2 = triangle(side=0.2, resnums=(4, 5, 6), molecules=(4, 5, 6))
        a2 = a2.copy()
        a2['globalIdx'] += 100
        p2 = {k + 100: v + np.array([2.0, 0.0, 0.0]) for k, v in p2.items()}
        c = RingController({'search_radius': 0.6, 'max_rings_per_iteration': 1})
        c.setup(total_groups=6, max_radius=1.0)
        with self.assertLogs(LOG, level='INFO'):
            _, chosen = c.search(pd.concat([a1, a2], ignore_index=True), {**p1, **p2},
                                 BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(chosen), 1)

    def test_reacted_groups_are_not_offered_again(self):
        adf, pos = triangle(side=0.25)
        adf.loc[adf['resNum'] == 3, 'z'] = 0
        c = RingController({'search_radius': 0.6})
        c.setup(total_groups=3, max_radius=1.0)
        with self.assertLogs(LOG, level='INFO'):
            sites, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(sites), 2)
        self.assertEqual(chosen, [])


class TestRestartState(unittest.TestCase):
    def test_round_trips(self, ):
        import tempfile, os
        s = RingCureState(iteration=4, rings=12, groups_consumed=36, radius_index=2,
                          total_groups=1440)
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, 'ring_state.yaml')
            s.to_yaml(f)
            back = RingCureState.from_yaml(f)
        self.assertEqual((back.iter, back.rings, back.groups_consumed, back.radius_index,
                          back.total_groups), (4, 12, 36, 2, 1440))
        self.assertAlmostEqual(back.conversion, 36 / 1440)

    def test_a_resumed_controller_keeps_its_widened_radius(self):
        s = RingCureState(radius_index=3)
        c = RingController({'search_radius': 0.5, 'radial_increment': 0.05}, state=s)
        c.setup(total_groups=10, max_radius=1.0)
        self.assertAlmostEqual(c.radius, 0.65)


class TestClosureDefaults(unittest.TestCase):
    def test_the_ladder_has_defaults(self):
        c = RingController({})
        self.assertEqual(c.dicts['closure']['nstages'], 8)
        self.assertAlmostEqual(c.dicts['closure']['target'], 0.15)

    def test_partial_configuration_keeps_the_rest(self):
        c = RingController({'closure': {'nstages': 12}})
        self.assertEqual(c.dicts['closure']['nstages'], 12)
        self.assertAlmostEqual(c.dicts['closure']['target'], 0.15)
        self.assertTrue(c.dicts['closure']['equilibration'])
