"""

.. module:: test_ringcontroller
   :synopsis: the loop of a ring-closing cure: what it searches, when it widens, when
              it stops, and how it counts conversion

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest

import numpy as np
import pandas as pd

from htpolynet.cure.ringcontroller import RingController, RingCureState, _triple_key

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


class TestWideningWithinAnIteration(unittest.TestCase):
    """The pairwise cure grows its radius until an iteration has enough bonds.  Widening
    only when an iteration finds nothing leaves a rigid monomer trickling one or two
    rings at the starting radius forever -- 150 bisphenol A dicyanates reached 0.40
    rather than the 0.60 asked for, having never once widened."""

    def test_it_widens_until_it_finds_a_ring(self):
        adf, pos = triangle(side=0.45)
        c = RingController({'search_radius': 0.2, 'radial_increment': 0.1})
        c.setup(total_groups=3, max_radius=1.0)
        with self.assertLogs(LOG, level='INFO'):
            _, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(chosen), 1)
        self.assertGreater(c.radius, 0.2)

    def test_it_gives_up_at_the_maximum_radius(self):
        adf, pos = triangle(side=2.0)
        c = RingController({'search_radius': 0.2, 'radial_increment': 0.1})
        c.setup(total_groups=3, max_radius=0.5)
        with self.assertLogs(LOG, level='INFO'):
            _, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(chosen, [])
        self.assertAlmostEqual(c.radius, 0.5)

    def test_the_floor_is_the_configured_one_when_far_from_target(self):
        c = RingController({'min_rings_per_iteration': 4, 'desired_conversion': 0.6})
        c.setup(total_groups=300, max_radius=1.0)
        self.assertEqual(c.ring_floor(), 4)

    def test_the_floor_falls_to_what_is_left_near_the_target(self):
        # three groups short of the target is one ring, so do not widen for four
        c = RingController({'min_rings_per_iteration': 4, 'desired_conversion': 0.6})
        c.setup(total_groups=300, max_radius=1.0)
        c.state.groups_consumed = 177
        self.assertEqual(c.ring_floor(), 1)

    def test_a_per_iteration_cap_lowers_the_floor(self):
        # never widen in search of more rings than the iteration would accept
        c = RingController({'min_rings_per_iteration': 4, 'max_rings_per_iteration': 2})
        c.setup(total_groups=300, max_radius=1.0)
        self.assertEqual(c.ring_floor(), 2)

    def test_the_floor_is_never_zero(self):
        c = RingController({'desired_conversion': 0.5})
        c.setup(total_groups=300, max_radius=1.0)
        c.state.groups_consumed = 150
        self.assertEqual(c.ring_floor(), 1)


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


class TestDecliningRingsThatDidNotClose(unittest.TestCase):
    """Late in a cure the matrix is rigid enough that a ring sometimes stops well short.
    Bonding it anyway leaves the bond long for good: two survived relaxation and a 500 K
    anneal still 2.6 A apart in a conversion-0.95 build."""

    def bonds_and_work(self, worsts):
        """One ring per entry in `worsts`, three pairs each, that ring's worst last."""
        rows, final = [], []
        for n, w in enumerate(worsts):
            for k in range(3):
                rows.append({'ai': 10 * n + k, 'aj': 10 * n + k + 5, 'triple': n})
                final.append(0.15 if k < 2 else w)
        return (pd.DataFrame(rows),
                pd.DataFrame({'final_distance': final}),
                [object() for _ in worsts])

    def test_a_ring_that_closed_is_kept(self):
        c = RingController({})
        bdf, work, chosen = self.bonds_and_work([0.22])
        out_b, out_c = c.accept(bdf, work, chosen)
        self.assertEqual(len(out_c), 1)
        self.assertEqual(out_b.shape[0], 3)

    def test_a_ring_that_did_not_is_dropped(self):
        c = RingController({})
        bdf, work, chosen = self.bonds_and_work([0.327])
        with self.assertLogs(LOG, level='INFO'):
            out_b, out_c = c.accept(bdf, work, chosen)
        self.assertEqual(out_c, [])
        self.assertTrue(out_b.empty)

    def test_the_survivors_are_renumbered_contiguously(self):
        # form_rings zips its maps against the triple column, so a gap would misalign them
        c = RingController({})
        bdf, work, chosen = self.bonds_and_work([0.22, 0.40, 0.21])
        with self.assertLogs(LOG, level='INFO'):
            out_b, out_c = c.accept(bdf, work, chosen)
        self.assertEqual(len(out_c), 2)
        self.assertEqual(sorted(out_b['triple'].unique()), [0, 1])
        self.assertEqual(out_c, [chosen[0], chosen[2]])

    def test_a_zero_limit_accepts_anything(self):
        c = RingController({'closure': {'max_accept': 0}})
        bdf, work, chosen = self.bonds_and_work([0.9])
        out_b, out_c = c.accept(bdf, work, chosen)
        self.assertEqual(len(out_c), 1)


class TestATripleThatFailedIsNotTriedAgainWhileAnythingElseIsLeft(unittest.TestCase):
    """A declined ring leaves its groups unreacted, so the next search offered the same
    triples, ran the same ladder and declined them again: nine iterations running at
    conversion 0.97, 77 declines, before one finally closed and the cure finished at
    0.971.  The memory skips them; it does not forbid them, because that cure escaped
    only because the geometry between iterations had moved."""

    def bonds_and_work(self, worsts, donors=None):
        """One ring per entry in `worsts`; `donors` gives each ring's three donor atoms."""
        rows, final = [], []
        for n, w in enumerate(worsts):
            d = donors[n] if donors else [10 * n + k for k in range(3)]
            for k in range(3):
                rows.append({'ai': d[k], 'aj': 10 * n + k + 5, 'triple': n})
                final.append(0.15 if k < 2 else w)
        return (pd.DataFrame(rows),
                pd.DataFrame({'final_distance': final}),
                [object() for _ in worsts])

    def two_triangles(self):
        """Two rings far enough apart to be independent, and a controller that sees both."""
        a1, p1 = triangle(side=0.2, resnums=(1, 2, 3), molecules=(1, 2, 3))
        a2, p2 = triangle(side=0.2, resnums=(4, 5, 6), molecules=(4, 5, 6))
        a2 = a2.copy()
        a2['globalIdx'] += 100
        p2 = {k + 100: v + np.array([2.0, 0.0, 0.0]) for k, v in p2.items()}
        c = RingController({'search_radius': 0.6})
        c.setup(total_groups=6, max_radius=1.0)
        return c, pd.concat([a1, a2], ignore_index=True), {**p1, **p2}

    def test_a_declined_ring_is_remembered_by_its_donor_atoms(self):
        # row numbers into the site table name a different ring next iteration; the
        # donors do not move
        c = RingController({})
        bdf, work, chosen = self.bonds_and_work([0.40], donors=[[7, 21, 34]])
        with self.assertLogs(LOG, level='INFO'):
            c.accept(bdf, work, chosen)
        self.assertEqual(c.declined, {_triple_key([34, 7, 21])})

    def test_a_ring_that_closed_is_not_remembered(self):
        c = RingController({})
        bdf, work, chosen = self.bonds_and_work([0.22])
        c.accept(bdf, work, chosen)
        self.assertEqual(c.declined, set())

    def test_the_search_takes_the_other_ring_instead(self):
        c, adf, pos = self.two_triangles()
        with self.assertLogs(LOG, level='INFO'):
            sites, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(chosen), 2)
        first = _triple_key(sites.at[a, 'donor'] for a in chosen[0][1])
        c.declined.add(first)
        with self.assertLogs(LOG, level='INFO') as cm:
            sites, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(chosen), 1)
        self.assertNotEqual(_triple_key(sites.at[a, 'donor'] for a in chosen[0][1]), first)
        self.assertIn('1 passed over as already failed', '\n'.join(cm.output))

    def test_it_is_lifted_rather_than_leave_the_iteration_idle(self):
        # the relaxation between iterations moves the groups, so a triple that failed
        # can still close later -- one did, at iteration 41, and it finished that cure
        adf, pos = triangle(side=0.25)
        c = RingController({'search_radius': 0.6})
        c.setup(total_groups=3, max_radius=1.0)
        with self.assertLogs(LOG, level='INFO'):
            sites, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        c.declined.add(_triple_key(sites['donor']))
        with self.assertLogs(LOG, level='INFO') as cm:
            _, chosen = c.search(adf, pos, BOX, 'BCY', (('N1', 'C1'),))
        self.assertEqual(len(chosen), 1)
        self.assertEqual(c.declined, set())
        self.assertIn('offering them again', '\n'.join(cm.output))

    def test_both_orientations_share_one_key(self):
        # candidate_triples scores both and keeps the better, so the declined one is
        # the better one and the other has farther to go
        self.assertEqual(_triple_key([3, 1, 2]), _triple_key([1, 2, 3]))

    def test_a_ring_exactly_at_the_limit_is_kept(self):
        # the log rounds to three decimals, so a decline reported at the limit was over it
        c = RingController({'closure': {'max_accept': 0.30}})
        bdf, work, chosen = self.bonds_and_work([0.30])
        _, out_c = c.accept(bdf, work, chosen)
        self.assertEqual(len(out_c), 1)
        self.assertEqual(c.declined, set())
