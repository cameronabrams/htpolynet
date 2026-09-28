"""

.. module:: test_postsim_moduli
   :synopsis: the fits behind the numbers postsim exists to produce -- Tg and E

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

import numpy as np

from htpolynet.analysis.postsim import PostSimDeform, PostsimConfiguration
from htpolynet.analysis.utils import compute_E, compute_tg

BOX = np.array([[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]])
MDP = 'dt = 0.002\nnsteps = 1\nref_t = 300\n'


class TestYoungsModulus(unittest.TestCase):
    """`compute_E` is the fit that turns a deformation run into the modulus a paper
    quotes, and it had no test."""

    def test_it_recovers_a_known_slope(self):
        strain = np.linspace(0.0, 0.1, 200)
        E, r2 = compute_E(strain, 2500.0 * strain)
        self.assertAlmostEqual(E, 2500.0, places=3)
        self.assertGreater(r2, 0.999)

    def test_it_fits_only_the_requested_window(self):
        # beyond the elastic regime the curve rolls over; the fit window is what keeps
        # that out of the modulus
        strain = np.linspace(0.0, 1.0, 400)
        stress = np.where(strain < 0.1, 2500.0 * strain, 250.0)
        E, _ = compute_E(strain, stress, fit_domain=[5, 35])
        self.assertAlmostEqual(E, 2500.0, delta=1.0)

    def test_a_rolled_over_curve_fitted_whole_understates_the_modulus(self):
        # the failure the window exists to prevent, stated as a test so the default
        # cannot quietly become "fit everything"
        strain = np.linspace(0.0, 1.0, 400)
        stress = np.where(strain < 0.1, 2500.0 * strain, 250.0)
        whole, _ = compute_E(strain, stress, fit_domain=[0, 400])
        self.assertLess(whole, 1000.0)

    def test_noise_does_not_bias_the_slope(self):
        rng = np.random.default_rng(0)
        strain = np.linspace(0.0, 0.1, 300)
        stress = 2500.0 * strain + rng.normal(0.0, 1.0, strain.size)
        E, r2 = compute_E(strain, stress, fit_domain=[0, 300])
        self.assertAlmostEqual(E, 2500.0, delta=25.0)
        self.assertGreater(r2, 0.9)


class TestGlassTransition(unittest.TestCase):
    """`compute_tg` intersects a cold-side and a hot-side line; Tg is where they cross."""

    def ladder(self, Tg=380.0, cold=1.0e-4, hot=5.0e-4, v0=1.0, n=60):
        T = np.linspace(250.0, 500.0, n)
        v = np.where(T < Tg, v0 + cold * (T - Tg), v0 + hot * (T - Tg))
        return T, v

    def test_it_finds_the_intersection(self):
        T, v = self.ladder(Tg=380.0)
        Tg, cold_par, hot_par = compute_tg(T, v, n_points=[15, 25])
        self.assertAlmostEqual(Tg, 380.0, delta=2.0)

    def test_the_two_slopes_differ_as_expected(self):
        T, v = self.ladder(cold=1.0e-4, hot=5.0e-4)
        _, cold_par, hot_par = compute_tg(T, v, n_points=[15, 25])
        self.assertLess(cold_par[0], hot_par[0])

    def test_it_tracks_a_shifted_transition(self):
        # windows chosen so neither straddles the transition; see the next test for what
        # happens when one does
        for want, n_points in ((320.0, [15, 25]), (420.0, [15, 15])):
            with self.subTest(want):
                T, v = self.ladder(Tg=want)
                got, _, _ = compute_tg(T, v, n_points=n_points)
                self.assertAlmostEqual(got, want, delta=3.0)

    def test_a_window_that_straddles_the_transition_biases_the_answer(self):
        # `compute_tg` fits the first n_points[0] and the last n_points[1] points and
        # takes it on faith that each window lies wholly on its own side of the
        # transition.  It does not check, and it does not warn.  With Tg at 420 K and a
        # hot window reaching down to about 400 K, the estimate comes back about 7 K
        # low -- small enough to look like a real Tg and to be quoted as one.
        T, v = self.ladder(Tg=420.0)
        straddled, _, _ = compute_tg(T, v, n_points=[15, 25])
        clean, _, _ = compute_tg(T, v, n_points=[15, 15])
        self.assertLess(straddled, clean - 4.0)
        self.assertAlmostEqual(clean, 420.0, delta=3.0)

    def test_a_single_straight_line_has_no_transition_to_find(self):
        # both fits land on the same line, so the intersection is degenerate; this
        # documents what comes back rather than asserting it is useful
        T = np.linspace(250.0, 500.0, 60)
        v = 1.0 + 1.0e-4 * T
        Tg, cold_par, hot_par = compute_tg(T, v, n_points=[15, 25])
        self.assertAlmostEqual(cold_par[0], hot_par[0], places=8)


class TestUniaxialDeformationMdp(unittest.TestCase):
    """The tensile counterpart of the shear stage, and the established one: its mdp had
    no test either."""

    def built(self, direction, edot=0.001):
        d = PostSimDeform({'direction': direction, 'edot': edot, 'ps': 1000})
        with tempfile.TemporaryDirectory() as t:
            f = os.path.join(t, 'deform.mdp')
            open(f, 'w').write(MDP)
            d.build_mdp(f, box=BOX)
            return open(f).read(), d.params

    def value(self, text, key):
        # match the key exactly: a prefix match picks up `deform-init-flow` for `deform`
        for l in text.split('\n'):
            if '=' in l and l.split('=')[0].strip() == key:
                return l.split('=')[1].split()
        raise AssertionError(f'{key} not in the mdp')

    def test_each_direction_drives_its_own_diagonal_element(self):
        for i, (d, L) in enumerate((('x', 5.0), ('y', 6.0), ('z', 7.0))):
            with self.subTest(d):
                text, _ = self.built(d)
                v = self.value(text, 'deform')
                self.assertAlmostEqual(float(v[i]), L * 0.001)
                self.assertEqual([v[j] for j in range(6) if j != i], ['0'] * 5)

    def test_the_driven_direction_is_not_pressure_coupled(self):
        # zero compressibility on the deformed axis; the others stay coupled
        for i, d in enumerate('xyz'):
            with self.subTest(d):
                text, _ = self.built(d)
                c = self.value(text, 'compressibility')
                self.assertEqual(float(c[i]), 0.0)
                self.assertTrue(all(float(c[j]) > 0 for j in range(3) if j != i))

    def test_the_traces_match_the_direction(self):
        for d, box_t, pres_t in (('x', 'Box-X', 'Pres-XX'), ('y', 'Box-Y', 'Pres-YY'),
                                 ('z', 'Box-Z', 'Pres-ZZ')):
            with self.subTest(d):
                _, p = self.built(d)
                self.assertEqual(p['traces'], [box_t, pres_t])
                self.assertEqual(p['output_deffnm'], f'deform-{d}')

    def test_the_flow_profile_is_initialized(self):
        # htpolynet shipped without this, so every deform run -- the Young's modulus
        # measurement -- was rejected by grompp on Gromacs 2025 before it started
        for d in 'xyz':
            with self.subTest(d):
                text, _ = self.built(d)
                self.assertEqual(self.value(text, 'deform-init-flow'), ['yes'])

    def test_an_unknown_direction_is_refused(self):
        d = PostSimDeform({'direction': 'q'})
        with tempfile.TemporaryDirectory() as t:
            f = os.path.join(t, 'deform.mdp')
            open(f, 'w').write(MDP)
            with self.assertLogs('htpolynet.analysis.postsim', level='ERROR'):
                d.build_mdp(f, box=BOX)


class TestEveryStageTypeIsRegistered(unittest.TestCase):
    def test_the_documented_stage_types_all_resolve(self):
        self.assertEqual(sorted(PostsimConfiguration.default_classes),
                         ['anneal', 'deform', 'equilibrate', 'ladder', 'shear'])
