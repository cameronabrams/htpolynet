"""

.. module:: test_plots_post_E
   :synopsis: ``plots post`` fits E per pull over a strain window, reports the scatter between pulls, and reads replicas

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.plot import do_E_plots
from htpolynet.analysis.utils import compute_E


def write_deform(root, subdir, d, E_bar, noise_bar, seed, max_strain=0.03, npts=301):
    """One deform-<d>.csv: stress = E*strain + noise (bar), strain from 0 to max_strain."""
    rng = np.random.default_rng(seed)
    strain = np.linspace(0.0, max_strain, npts)
    stress = E_bar * strain + rng.normal(0.0, noise_bar, strain.shape)
    D = d.upper()
    sub = os.path.join(root, subdir)
    os.makedirs(sub)
    pd.DataFrame({f'Box-{D}-strain': strain, f'Pres-{D}{D}-stress': stress}).to_csv(
        os.path.join(sub, f'deform-{d}.csv'), index=False)


class TestComputeEStrainWindow(unittest.TestCase):
    def test_window_is_in_strain_not_rows(self):
        # linear to 2%, then a plateau: only a strain window independent of row
        # spacing recovers the linear slope for both samplings
        for npts in (101, 1001):
            strain = np.linspace(0.0, 0.05, npts)
            stress = np.where(strain < 0.02, 3000.0 * strain, 60.0)
            E, R2 = compute_E(strain, stress, fit_strain=[0.001, 0.02])
            self.assertAlmostEqual(E, 3000.0, delta=30.0)

    def test_too_few_points_in_window_raises(self):
        with self.assertRaisesRegex(ValueError, 'fewer than two points'):
            compute_E(np.array([0.0, 0.0001]), np.array([0.0, 1.0]), fit_strain=[0.001, 0.02])


class TestDoEPlots(unittest.TestCase):
    def run_E(self, phases, writer):
        with tempfile.TemporaryDirectory() as d:
            writer(d)
            cwd = os.getcwd()
            os.chdir(d)
            try:
                with self.assertLogs('htpolynet.analysis.plot', level='INFO') as logs:
                    result = do_E_plots(phases, [d])
                mean = pd.read_csv('E.csv', sep=' ')
                fits = pd.read_csv('E-fits.csv', sep=' ')
            finally:
                os.chdir(cwd)
        return result, mean, fits, logs.output

    @staticmethod
    def phases(**extra):
        return [{'deform': {'direction': d, 'subdir': f'postsim/deform-{d}', **extra}} for d in 'xyz']

    def test_three_directions_give_a_mean_and_an_error_bar(self):
        def w(root):
            for i, ax in enumerate('xyz'):
                write_deform(root, f'postsim/deform-{ax}', ax, 3.0e4, 1.0, seed=i)
        result, mean, fits, out = self.run_E(self.phases(), w)
        self.assertEqual(mean.shape[0], 301)
        self.assertEqual(len(fits), 3)
        self.assertAlmostEqual(result['E'] / 1000.0, 3.0, places=2)   # 3 GPa
        self.assertGreater(result['E_sem'], 0.0)
        self.assertFalse(any('WARNING' in line for line in out), out)

    def test_replicas_are_read_from_their_own_directories(self):
        def w(root):
            for i, ax in enumerate('xyz'):
                for k in (1, 2):
                    write_deform(root, f'postsim/deform-{ax}-r{k}', ax, 3.0e4, 1.0, seed=10 * i + k)
        result, _, fits, _ = self.run_E(self.phases(replicas=2), w)
        self.assertEqual(result['n'], 6)
        self.assertTrue(all(p.endswith(('-r1', '-r2')) for p in fits['pull']))

    def test_a_pull_short_of_the_window_is_reported(self):
        def w(root):
            for i, ax in enumerate('xyz'):
                write_deform(root, f'postsim/deform-{ax}', ax, 3.0e4, 1.0, seed=i, max_strain=0.01, npts=101)
        result, _, _, out = self.run_E(self.phases(), w)
        self.assertTrue(any('short of the fit window' in line for line in out), out)
        self.assertAlmostEqual(result['fit_strain'][1], 0.01)

    def test_a_noise_dominated_fit_is_flagged(self):
        # ex6 reference build: 1% strain against several hundred bar of noise
        def w(root):
            for i, ax in enumerate('xyz'):
                write_deform(root, f'postsim/deform-{ax}', ax, 3.0e4, 2000.0, seed=i, max_strain=0.01, npts=101)
        _, _, _, out = self.run_E(self.phases(fit_strain=[0.001, 0.01]), w)
        self.assertTrue(any('WARNING' in line and 'standard error is' in line for line in out), out)

    def test_a_residual_stress_does_not_bias_E(self):
        # ex6's reference structure started at +150-200 bar; through the origin that
        # offset/strain was added to E
        def w(root):
            for i, ax in enumerate('xyz'):
                write_deform(root, f'postsim/deform-{ax}', ax, 3.0e4, 1.0, seed=i)
                f = os.path.join(root, f'postsim/deform-{ax}', f'deform-{ax}.csv')
                df = pd.read_csv(f)
                df[f'Pres-{ax.upper()*2}-stress'] += 200.0
                df.to_csv(f, index=False)
        result, _, _, _ = self.run_E(self.phases(), w)
        self.assertAlmostEqual(result['E'] / 1000.0, 3.0, places=2)

    def test_pulls_that_disagree_are_flagged(self):
        # clean curves, but directions far apart: the fit of each is good and the
        # scatter between them is what says the number is not yet known
        def w(root):
            for i, (ax, E) in enumerate(zip('xyz', (1.0e4, 3.0e4, 6.0e4))):
                write_deform(root, f'postsim/deform-{ax}', ax, E, 1.0, seed=i)
        result, _, _, out = self.run_E(self.phases(), w)
        self.assertTrue(any('standard error is' in line for line in out), out)
