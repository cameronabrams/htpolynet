"""

.. module:: test_postsim_real_traces
   :synopsis: Tg and E on real trajectories, where the fit window matters

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.utils import compute_E, compute_tg

FIX = os.path.join(os.path.dirname(__file__), 'fixtures')


def ladder():
    return pd.read_csv(os.path.join(FIX, 'ladder-FDE-DFDA.csv'), sep=' ')


def deform():
    return pd.read_csv(os.path.join(FIX, 'deform-x-FDE-FCPDA.csv'))


class TestTgOnARealLadder(unittest.TestCase):
    """61 points from 300 to 599 K of a real FDE-DFDA build."""

    def test_it_lands_in_the_transition_region(self):
        d = ladder()
        Tg, cold, hot = compute_tg(d['Temperature'], d['Density'], n_points=[10, 20])
        self.assertAlmostEqual(Tg, 383.9, delta=1.0)
        self.assertGreater(Tg, d['Temperature'].min())
        self.assertLess(Tg, d['Temperature'].max())

    def test_the_glassy_side_is_the_stiffer_one(self):
        d = ladder()
        _, cold, hot = compute_tg(d['Temperature'], d['Density'], n_points=[10, 20])
        # density falls with temperature, faster above Tg
        self.assertLess(hot[0], cold[0])

    def test_it_reproduces_the_lines_stored_in_the_trace(self):
        # the file carries the glassy and rubbery lines fitted when it was made; their
        # intersection is an independent value to check against, and it is recovered
        # with the windows that produced it
        d = ladder()
        T = d['Temperature'].to_numpy()
        g = np.polyfit(T, d['glassy-line-Density'], 1)
        r = np.polyfit(T, d['rubbery-line-Density'], 1)
        stored = -(r[1] - g[1]) / (r[0] - g[0])
        Tg, _, _ = compute_tg(T, d['Density'].to_numpy(), n_points=[10, 30])
        self.assertAlmostEqual(Tg, stored, delta=0.5)

    def test_the_fit_window_moves_the_answer_by_tens_of_kelvin(self):
        """The reason this fixture exists.

        On synthetic data `compute_tg` recovers its input to a fraction of a degree.  On
        a real ladder the answer spans about 29 K across fit windows a reasonable person
        would choose, and nothing in htpolynet records which window produced a published
        number.  If this spread ever narrows, the fixture or the fit changed and someone
        should know which.
        """
        d = ladder()
        T, v = d['Temperature'].to_numpy(), d['Density'].to_numpy()
        got = [compute_tg(T, v, n_points=[c, h])[0]
               for c in (8, 10, 12, 15, 20) for h in (10, 15, 20, 25, 30)]
        self.assertAlmostEqual(min(got), 377.1, delta=1.0)
        self.assertAlmostEqual(max(got), 406.3, delta=1.0)
        self.assertGreater(max(got) - min(got), 25.0)


class TestYoungsModulusOnARealDeformation(unittest.TestCase):
    """1001 frames written by PostSimDeform itself, columns and all."""

    def test_the_stage_writes_the_columns_the_fit_needs(self):
        d = deform()
        for c in ('Box-X', 'Pres-XX', 'Box-X-strain', 'Pres-XX-stress'):
            self.assertIn(c, d.columns)

    def test_the_derived_columns_are_consistent_with_the_raw_ones(self):
        # stress is the negated pressure; strain starts at zero and rises
        d = deform()
        np.testing.assert_allclose(d['Pres-XX-stress'], -d['Pres-XX'], rtol=1e-9)
        self.assertAlmostEqual(d['Box-X-strain'].iloc[0], 0.0, places=9)
        self.assertGreater(d['Box-X-strain'].iloc[-1], d['Box-X-strain'].iloc[0])

    def test_it_gives_a_modulus_of_the_right_order(self):
        d = deform()
        E, _ = compute_E(d['Box-X-strain'], d['Pres-XX-stress'], fit_domain=[10, 200])
        self.assertAlmostEqual(E / 10.0, 4254.0, delta=50.0)   # bar -> MPa
        self.assertGreater(E / 10.0, 500.0)                    # a solid, not a liquid

    def test_the_fit_is_poor_even_where_the_modulus_is_sensible(self):
        """Instantaneous stress in MD is extremely noisy.

        R-squared here is well under 0.5, which on synthetic data would mean the fit had
        failed.  It does not mean that: it means single-frame pressure scatters hugely
        about the trend.  Anyone who adds an R-squared threshold as a quality gate should
        see this first.
        """
        d = deform()
        _, r2 = compute_E(d['Box-X-strain'], d['Pres-XX-stress'], fit_domain=[10, 200])
        self.assertLess(r2, 0.5)
        self.assertGreater(r2, 0.1)

    def test_the_fit_window_moves_the_modulus_by_nearly_a_factor_of_two(self):
        d = deform()
        got = [compute_E(d['Box-X-strain'], d['Pres-XX-stress'], fit_domain=w)[0] / 10.0
               for w in ([5, 50], [10, 100], [10, 200], [20, 300])]
        self.assertAlmostEqual(min(got), 2420.0, delta=50.0)
        self.assertAlmostEqual(max(got), 4398.0, delta=50.0)
        self.assertGreater(max(got) / min(got), 1.7)


class TestThereIsNoShearFixture(unittest.TestCase):
    def test_no_shear_trace_exists_yet(self):
        """`shear` has never been run, so G has no real-trajectory test.

        This asserts the absence deliberately: when someone adds a shear trace, this
        fails and the accompanying tests should be written at the same time.
        """
        self.assertEqual([f for f in os.listdir(FIX) if f.startswith('shear')], [])
