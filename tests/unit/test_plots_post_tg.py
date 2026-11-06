"""

.. module:: test_plots_post_tg
   :synopsis: ``plots post`` says when a Tg depends on the fit window, and needs its config

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import types
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.plot import do_tg_plots, post_plots
from htpolynet.analysis.utils import (_tg_hyperbola, ceiling_drift, fit_tg_hyperbola,
                                      pooled_tg, tg_window_range)

TEMPS = np.arange(300.0, 605.0, 5.0)          # 61 rungs, as in example 6's ladders
LADDER = {'warmup_ps': 0, 'ps_per_rise': 1, 'ps_per_run': 10, 'deltaT': 5}


def kinked(T, Tg=450.0, glassy=-0.20, rubbery=-0.35, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    rho = np.where(T < Tg, 1200.0 + glassy * (T - 300.0),
                   1200.0 + glassy * (Tg - 300.0) + rubbery * (T - Tg))
    return rho + rng.normal(0.0, noise, T.shape)


def write_ladder(root, subdir, rho):
    """ladder.csv with one constant-T, constant-density block per rung."""
    rows = []
    step = LADDER['ps_per_rise'] + LADDER['ps_per_run']
    for k, (T, r) in enumerate(zip(TEMPS, rho)):
        for t in range(step):
            rows.append({'time(ps)': k * step + t, 'Temperature': T, 'Density': r,
                         'Volume': 1.0e5 / r})
    os.makedirs(os.path.join(root, subdir))
    pd.DataFrame(rows).to_csv(os.path.join(root, subdir, 'ladder.csv'), index=False)


class TestWindowRange(unittest.TestCase):
    def test_a_clean_kink_holds_its_tg_across_windows(self):
        lo, hi = tg_window_range(TEMPS, kinked(TEMPS))
        self.assertAlmostEqual(lo, 450.0, delta=1.0)
        self.assertAlmostEqual(hi, 450.0, delta=1.0)

    def test_nearly_parallel_lines_scatter_off_the_scale(self):
        # no visible transition: the two fitted lines are parallel up to noise
        lo, hi = tg_window_range(TEMPS, kinked(TEMPS, glassy=-0.30, rubbery=-0.30,
                                               noise=0.5, seed=0))
        self.assertTrue(lo < 300.0 or hi > 600.0, (lo, hi))

    def test_windows_that_do_not_fit_the_data_are_skipped(self):
        self.assertIsNone(tg_window_range(TEMPS[:15], kinked(TEMPS[:15])))


def smooth(T, Tg=450.0, width=25.0, noise=0.0, seed=0):
    """A transition of finite width, from the fitted form itself."""
    rng = np.random.default_rng(seed)
    rho = _tg_hyperbola(T, 1150.0, Tg, -0.20, -0.15, 2 * np.log(width))
    return rho + rng.normal(0.0, noise, T.shape)


class TestWholeCurveFit(unittest.TestCase):
    def test_a_smooth_transition_is_recovered_with_an_error_bar(self):
        f = fit_tg_hyperbola(TEMPS, smooth(TEMPS, noise=0.3, seed=1))
        self.assertAlmostEqual(f['Tg'], 450.0, delta=10.0)
        self.assertAlmostEqual(f['glassy'], -0.20, delta=0.03)
        self.assertAlmostEqual(f['rubbery'], -0.35, delta=0.03)
        self.assertFalse(f['held_out'])
        self.assertTrue(np.isfinite(f['Tg_se']))

    def test_a_sharp_kink_is_kept_and_flagged_sharp(self):
        f = fit_tg_hyperbola(TEMPS, kinked(TEMPS))
        self.assertAlmostEqual(f['Tg'], 450.0, delta=2.0)
        self.assertTrue(f['sharp'])
        self.assertFalse(f['held_out'])

    def test_no_transition_in_the_ladder_is_held_out(self):
        # the transition sits below the ladder: only one asymptote is in the data
        f = fit_tg_hyperbola(TEMPS, smooth(TEMPS, Tg=150.0, width=10.0, noise=0.3, seed=2))
        self.assertTrue(f['held_out'], f)
        self.assertIsNotNone(f['reason'])


    def test_both_tg_definitions_agree_when_the_ladder_reaches_both_asymptotes(self):
        f = fit_tg_hyperbola(TEMPS, smooth(TEMPS, width=15.0, noise=0.2, seed=4))
        self.assertAlmostEqual(f['Tg_tangents'], f['Tg'], delta=3.0)

    def test_a_ladder_too_narrow_for_its_transition_is_screened_out(self):
        # width 60 K centered at 450: the 90% interval is about 450 +/- 160 K, past both ends
        f = fit_tg_hyperbola(TEMPS, smooth(TEMPS, width=60.0, noise=0.05, seed=5))
        self.assertTrue(f['held_out'], f)
        self.assertIn('does not reach past the transition', f['reason'])
        lo, hi = f['transition']
        self.assertTrue(lo < TEMPS[0] or hi > TEMPS[-1])

    def test_synthetic_data_give_a_within_build_spread(self):
        f = fit_tg_hyperbola(TEMPS, smooth(TEMPS, noise=0.5, seed=6), n_synth=40)
        self.assertIsNotNone(f['Tg_synth_sd'])
        self.assertGreater(f['Tg_synth_sd'], 0.0)
        self.assertLess(f['Tg_synth_sd'], 30.0)


class TestPooledFit(unittest.TestCase):
    def curves(self, n, noise=1.5):
        return [(TEMPS, smooth(TEMPS, Tg=450.0 + 8.0 * (k - n / 2), noise=noise, seed=k))
                for k in range(n)]

    def test_pooling_builds_gives_a_between_build_error_bar(self):
        res = pooled_tg(self.curves(5), n_boot=100)
        self.assertEqual(res['n_builds'], 5)
        self.assertEqual(len(res['per_build']), 5)
        self.assertAlmostEqual(res['Tg'], 450.0, delta=15.0)
        self.assertIsNotNone(res['Tg_boot_sd'])
        self.assertGreater(res['Tg_boot_sd'], 0.0)

    def test_the_kept_builds_own_spread_is_reported(self):
        res = pooled_tg(self.curves(5), n_boot=40)
        own = [f['Tg'] for f in res['per_build'] if not f['held_out']]
        self.assertEqual(res['n_kept'], len(own))
        self.assertAlmostEqual(res['per_build_sd'], float(np.std(own, ddof=1)))

    def test_a_sharp_pool_most_builds_refuse_is_held_out(self):
        """Example 2 cooling, 2026-10-23: four of five builds refused on their own, the
        pool kept as a 1 K-wide break that the local slope does not show."""
        from unittest import mock
        import htpolynet.analysis.utils as U

        def fake(T, v, *a, **k):
            pooled = len(T) > len(TEMPS)
            refused = not pooled and v[0] != 0.0
            return {'Tg': 400.0, 'Tg_se': float('nan'), 'sharp': True,
                    'held_out': refused, 'reason': 'no significant change of slope' if refused else None}
        curves = [(TEMPS, np.full(TEMPS.shape, float(k))) for k in range(5)]   # build 0 kept
        with mock.patch.object(U, 'fit_tg_hyperbola', side_effect=fake), \
             mock.patch.object(U, '_align_builds', side_effect=lambda cs: cs):
            res = U.pooled_tg(curves, n_boot=5)
        self.assertTrue(res['fit']['held_out'])
        self.assertIn('4 of 5 builds do not show', res['fit']['reason'])
        self.assertIsNone(res['Tg'])

    def test_one_build_has_no_bootstrap(self):
        self.assertIsNone(pooled_tg(self.curves(1), n_boot=50)['Tg_boot_sd'])

    def test_the_bootstrap_is_reproducible(self):
        a = pooled_tg(self.curves(4), n_boot=60, seed=7)['Tg_boot_sd']
        b = pooled_tg(self.curves(4), n_boot=60, seed=7)['Tg_boot_sd']
        self.assertEqual(a, b)


def still_bending(T, noise=0.3, seed=0):
    """A rubbery branch that keeps steepening, so the top of the ladder pulls Tg along."""
    return smooth(T, noise=noise, seed=seed) - 0.0015 * np.clip(T - 450.0, 0, None) ** 2


class TestCeilingDrift(unittest.TestCase):
    def test_a_tall_enough_ladder_barely_moves(self):
        T = np.arange(300.0, 705.0, 5.0)
        curves = [(T, smooth(T, noise=1.0, seed=k)) for k in range(5)]
        full = pooled_tg(curves, n_boot=100)
        d = ceiling_drift(curves)
        self.assertEqual([x['cut'] for x in d], [50.0, 100.0])
        self.assertAlmostEqual(d[0]['ceiling'], 650.0)
        self.assertLess(abs(d[0]['Tg'] - full['Tg']), full['Tg_boot_sd'])

    def test_a_branch_still_bending_at_the_top_drifts(self):
        curves = [(TEMPS, still_bending(TEMPS, seed=k)) for k in range(5)]
        full = pooled_tg(curves, n_boot=20)
        self.assertIsNotNone(full['Tg'])
        d = ceiling_drift(curves)
        self.assertGreater(abs(d[0]['Tg'] - full['Tg']), 10.0)

    def test_a_cut_that_loses_the_transition_says_why(self):
        curves = [(TEMPS, smooth(TEMPS, noise=0.3, seed=k)) for k in range(3)]
        d = ceiling_drift(curves)
        self.assertIsNone(d[1]['Tg'])
        self.assertIn('does not reach past the transition', d[1]['reason'])


class TestDoTgPlots(unittest.TestCase):
    def run_tg(self, heat, cool, cooling_first=False, only=None):
        h = {'ladder': dict(LADDER, subdir='heat')}
        c = {'ladder': dict(LADDER, subdir='cool', deltaT=-LADDER['deltaT'])}
        phases = [c, h] if cooling_first else [h, c]
        if only is not None:
            phases = [h] if only == 'heating' else [c]
        with tempfile.TemporaryDirectory() as d:
            write_ladder(d, 'heat', heat)
            write_ladder(d, 'cool', cool)
            cwd = os.getcwd()
            os.chdir(d)
            try:
                with self.assertLogs('htpolynet.analysis.plot', level='INFO') as logs:
                    result = do_tg_plots(phases, [d])
            finally:
                os.chdir(cwd)
        return result, '\n'.join(logs.output)

    def test_two_sound_ladders_report_their_window_range_and_no_warning(self):
        result, log = self.run_tg(kinked(TEMPS), kinked(TEMPS, Tg=440.0))
        self.assertAlmostEqual(result['heating'], 450.0, delta=1.0)
        self.assertIn('over fit windows', log)
        self.assertIn('heating Tg, whole-curve fit: 450.0 K', log)
        self.assertNotIn('WARNING', log)
        self.assertEqual(result['unquotable'], [])

    def test_an_unfittable_heating_ladder_is_named(self):
        bad = kinked(TEMPS, glassy=-0.30, rubbery=-0.30, noise=0.5, seed=0)
        result, log = self.run_tg(bad, kinked(TEMPS))
        self.assertEqual(result['unquotable'], ['heating'])
        self.assertIn('no heating Tg can be read off this ladder', log)

    def test_direction_comes_from_the_rate_sign_not_the_order(self):
        """Plan B on example 5 listed cooling first and got its cooling ladder called heating."""
        result, log = self.run_tg(kinked(TEMPS, Tg=422.0), kinked(TEMPS, Tg=478.0), cooling_first=True)
        self.assertEqual(result['directions'], ['cooling', 'heating'])
        self.assertAlmostEqual(result['heating'], 422.0, delta=1.0)
        self.assertAlmostEqual(result['cooling'], 478.0, delta=1.0)
        self.assertIn('cooling Tg = 478', log)

    def test_a_single_ladder_runs_and_compares_nothing(self):
        result, log = self.run_tg(kinked(TEMPS), kinked(TEMPS), only='cooling')
        self.assertEqual(result['directions'], ['cooling'])
        self.assertIsNone(result['heating'])
        self.assertNotIn('differ by', log)

    def test_a_tg_that_follows_the_ladder_top_is_flagged(self):
        result, log = self.run_tg(still_bending(TEMPS), kinked(TEMPS))
        self.assertIn('heating Tg, whole-curve fit with the top of the ladder cut away', log)
        self.assertIn('the heating Tg moves', log)
        self.assertIn('raise its top', log)
        self.assertNotIn('the cooling Tg moves', log)

    def test_a_crashed_ladder_does_not_discard_the_finished_one(self):
        """Job 26725677: a heating ladder that died took its finished cooling ladder with it."""
        h = {'ladder': dict(LADDER, subdir='heat')}
        c = {'ladder': dict(LADDER, subdir='cool', deltaT=-LADDER['deltaT'])}
        with tempfile.TemporaryDirectory() as top:
            builds = [os.path.join(top, f'proj-{k}') for k in range(2)]
            for k, d in enumerate(builds):
                os.makedirs(d)
                write_ladder(d, 'cool', kinked(TEMPS, Tg=440.0, seed=k))
            write_ladder(builds[1], 'heat', kinked(TEMPS, seed=1))  # proj-0's heating died
            cwd = os.getcwd()
            os.chdir(top)
            try:
                with self.assertLogs('htpolynet.analysis.plot', level='INFO') as logs:
                    result = do_tg_plots([c, h], builds)
                written = sorted(f for f in os.listdir(top) if f.endswith('data.csv'))
            finally:
                os.chdir(cwd)
        log = '\n'.join(logs.output)
        self.assertEqual(written, ['ladder0-data.csv', 'ladder1-data.csv'])
        self.assertEqual(result['directions'], ['cooling', 'heating'])
        self.assertAlmostEqual(result['cooling'], 440.0, delta=2.0)
        self.assertAlmostEqual(result['heating'], 450.0, delta=2.0)
        self.assertIn('proj-0', log)
        self.assertIn('did this ladder finish', log)

    def test_a_ladder_no_build_finished_is_skipped(self):
        h = {'ladder': dict(LADDER, subdir='heat')}
        c = {'ladder': dict(LADDER, subdir='cool', deltaT=-LADDER['deltaT'])}
        with tempfile.TemporaryDirectory() as d:
            write_ladder(d, 'cool', kinked(TEMPS, Tg=440.0))
            cwd = os.getcwd()
            os.chdir(d)
            try:
                with self.assertLogs('htpolynet.analysis.plot', level='INFO') as logs:
                    result = do_tg_plots([c, h], [d])
                self.assertTrue(os.path.exists('ladder0-data.csv'))
                self.assertFalse(os.path.exists('ladder1-data.csv'))
            finally:
                os.chdir(cwd)
        self.assertEqual(result['directions'], ['cooling'])
        self.assertIsNone(result['heating'])
        self.assertIn('skipping it', '\n'.join(logs.output))

    def test_heating_and_cooling_far_apart_is_flagged(self):
        result, log = self.run_tg(kinked(TEMPS, Tg=422.0), kinked(TEMPS, Tg=478.0))
        self.assertIn('differ by 56 K', log)


class TestPostPlotsNeedsItsConfig(unittest.TestCase):
    def test_no_cfg_is_an_error_not_a_silent_success(self):
        args = types.SimpleNamespace(n_points=[10, 20], cfg=[], proj=['proj-0'])
        with self.assertRaisesRegex(ValueError, '--cfg'):
            post_plots(args)


if __name__ == '__main__':
    unittest.main()
