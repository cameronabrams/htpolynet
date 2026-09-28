"""

.. module:: test_postsim_shear
   :synopsis: the simple-shear postsim stage that yields a shear modulus

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import unittest

import numpy as np

from htpolynet.analysis.postsim import PostSimShear, PostsimConfiguration
from htpolynet.analysis.utils import compute_E, compute_G

BOX = np.array([[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]])

MDP = 'dt = 0.002\nnsteps = 1\nref_t = 300\n'


def built(direction, edot=0.001, ps=1000):
    """Runs build_mdp on a throwaway mdp and returns (its contents, the stage params)."""
    s = PostSimShear({'direction': direction, 'edot': edot, 'ps': ps})
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'shear.mdp')
        open(f, 'w').write(MDP)
        s.build_mdp(f, box=BOX)
        return open(f).read(), s.params


class TestItIsRegistered(unittest.TestCase):
    def test_shear_is_a_postsim_stage(self):
        self.assertIs(PostsimConfiguration.default_classes['shear'], PostSimShear)


class TestTheDeformedBoxElement(unittest.TestCase):
    """Gromacs' deform takes six rates in the order XX YY ZZ YX ZX ZY, so which slot is
    driven decides which plane is sheared."""

    def slot(self, text):
        for l in text.split('\n'):
            if '=' in l and l.split('=')[0].strip() == 'deform':
                return l.split('=')[1].split()
        raise AssertionError('deform not in the mdp')

    def test_xy_drives_the_yx_element_at_the_y_length(self):
        text, _ = built('xy')
        v = self.slot(text)
        self.assertEqual([v[0], v[1], v[2], v[4], v[5]], ['0'] * 5)
        self.assertAlmostEqual(float(v[3]), 6.0 * 0.001)     # L_y * edot

    def test_xz_drives_the_zx_element_at_the_z_length(self):
        text, _ = built('xz')
        v = self.slot(text)
        self.assertAlmostEqual(float(v[4]), 7.0 * 0.001)
        self.assertEqual([v[0], v[1], v[2], v[3], v[5]], ['0'] * 5)

    def test_yz_drives_the_zy_element_at_the_z_length(self):
        text, _ = built('yz')
        v = self.slot(text)
        self.assertAlmostEqual(float(v[5]), 7.0 * 0.001)
        self.assertEqual(v[:5], ['0'] * 5)

    def test_the_rate_scales_with_edot(self):
        v = self.slot(built('xy', edot=0.004)[0])
        self.assertAlmostEqual(float(v[3]), 6.0 * 0.004)


class TestTheBarostatDoesNotFightTheDeformation(unittest.TestCase):
    """The driven element is off-diagonal, so the off-diagonal compressibilities must be
    zero while the normal directions stay coupled at P."""

    def value(self, text, key):
        # match the key exactly: a prefix match picks up `deform-init-flow` for `deform`
        for l in text.split('\n'):
            if '=' in l and l.split('=')[0].strip() == key:
                return l.split('=')[1].split()
        raise AssertionError(f'{key} not in the mdp')

    def test_off_diagonal_compressibility_is_zero(self):
        text, _ = built('xy')
        self.assertEqual(self.value(text, 'compressibility')[3:], ['0', '0', '0'])

    def test_the_sheared_cartesian_component_is_uncoupled(self):
        # grompp refuses a run where an off-diagonal element is deformed while the
        # barostat still acts on that same Cartesian component of another box vector:
        # "spurious periodicity effects".  YX and ZX move x; ZY moves y.
        for d, zeroed in (('xy', 0), ('xz', 0), ('yz', 1)):
            with self.subTest(d):
                c = self.value(built(d)[0], 'compressibility')
                self.assertEqual(float(c[zeroed]), 0.0)
                self.assertTrue(all(float(c[j]) > 0 for j in range(3) if j != zeroed))

    def test_the_flow_profile_is_initialized(self):
        # from Gromacs 2025 this is not optional: deform plus generated velocities is a
        # hard grompp ERROR without it, so the stage simply would not run
        text, _ = built('xy')
        self.assertEqual(self.value(text, 'deform-init-flow'), ['yes'])

    def test_stress_is_sampled_ten_times_per_ps(self):
        text, _ = built('xy')
        self.assertEqual(self.value(text, 'nstenergy'), ['50'])
        self.assertEqual(self.value(text, 'nstcalcenergy'), ['50'])

    def test_the_other_normal_directions_stay_coupled(self):
        text, _ = built('xy')
        self.assertEqual(self.value(text, 'ref_p')[3:], ['0', '0', '0'])

    def test_the_run_length_follows_ps_and_dt(self):
        text, _ = built('xy', ps=500)
        self.assertEqual(int(self.value(text, 'nsteps')[0]), 250000)


class TestTheTracesMatchThePlane(unittest.TestCase):
    def test_each_plane_traces_its_own_box_and_pressure_terms(self):
        for d, box_t, pres_t in (('xy', 'Box-YX', 'Pres-XY'),
                                 ('xz', 'Box-ZX', 'Pres-XZ'),
                                 ('yz', 'Box-ZY', 'Pres-ZY')):
            with self.subTest(d):
                _, p = built(d)
                self.assertEqual(p['traces'], [box_t, pres_t])
                self.assertEqual(p['output_deffnm'], f'shear-{d}')

    def test_an_unknown_plane_is_refused_rather_than_guessed(self):
        s = PostSimShear({'direction': 'zz'})
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, 'shear.mdp')
            open(f, 'w').write(MDP)
            with self.assertLogs('htpolynet.analysis.postsim', level='ERROR'):
                s.build_mdp(f, box=BOX)
            self.assertNotIn('deform', open(f).read())


class TestComputeG(unittest.TestCase):
    def test_it_recovers_a_known_slope(self):
        strain = np.linspace(0.0, 0.1, 200)
        G, r2 = compute_G(strain, 1200.0 * strain)
        self.assertAlmostEqual(G, 1200.0, places=3)
        self.assertGreater(r2, 0.999)

    def test_it_is_the_same_fit_as_youngs(self):
        strain = np.linspace(0.0, 0.1, 200)
        stress = 900.0 * strain
        self.assertAlmostEqual(compute_G(strain, stress)[0], compute_E(strain, stress)[0])
