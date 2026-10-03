"""

.. module:: test_plots_post_E
   :synopsis: ``plots post`` averages deform curves without duplicate-column errors, and flags a noise fit

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.plot import do_E_plots


def write_deform(root, d, E_bar, noise_bar, seed):
    """One deform-<d>.csv: 101 points to 1% strain, stress = E*strain + noise (bar)."""
    rng = np.random.default_rng(seed)
    strain = np.arange(101) * 1.0e-4
    stress = E_bar * strain + rng.normal(0.0, noise_bar, strain.shape)
    D = d.upper()
    sub = os.path.join(root, 'postsim', f'deform-{d}')
    os.makedirs(sub)
    pd.DataFrame({f'Box-{D}-strain': strain, f'Pres-{D}{D}-stress': stress}).to_csv(
        os.path.join(sub, f'deform-{d}.csv'), index=False)


class TestDoEPlots(unittest.TestCase):
    phases = [{'deform': {'direction': d, 'subdir': f'postsim/deform-{d}'}} for d in 'xyz']

    def run_E(self, E_bar, noise_bar):
        with tempfile.TemporaryDirectory() as d:
            for i, ax in enumerate('xyz'):
                write_deform(d, ax, E_bar, noise_bar, seed=i)
            cwd = os.getcwd()
            os.chdir(d)
            try:
                with self.assertLogs('htpolynet.analysis.plot', level='INFO') as logs:
                    do_E_plots(self.phases, [d])
                mean = pd.read_csv('E.csv', sep=' ')
            finally:
                os.chdir(cwd)
        return mean, logs.output

    def test_three_directions_average_without_error(self):
        # pandas >= 3 refused the old side-by-side stack of same-named columns
        mean, out = self.run_E(E_bar=3.0e4, noise_bar=1.0)   # 3 GPa
        self.assertEqual(mean.shape[0], 101)
        self.assertIn('stress-std', mean.columns)
        self.assertTrue(any('E = 3.0' in line for line in out), out)
        self.assertFalse(any('WARNING' in line for line in out), out)

    def test_a_noise_dominated_fit_is_flagged(self):
        # ex6 reference build: 1% strain against several hundred bar of noise
        _, out = self.run_E(E_bar=3.0e4, noise_bar=500.0)
        self.assertTrue(any('WARNING' in line and 'not a modulus' in line for line in out), out)
