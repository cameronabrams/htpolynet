"""

.. module:: test_mdp_ensembles
   :synopsis: the packaged mdps ask for couplings that sample the right ensemble

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import importlib.resources as ir
import unittest

DEPRECATED_BAROSTATS = {'berendsen'}
"""Gromacs 2025+: "The Berendsen barostat does not generate any strictly correct
ensemble, and should not be used for new production simulations." """


def mdps():
    d = ir.files('htpolynet').joinpath('resources/mdp')
    return {f.name: f.read_text() for f in d.iterdir() if f.name.endswith('.mdp')}


def setting(text, key):
    for line in text.split('\n'):
        if '=' in line and line.split('=')[0].strip() == key:
            return line.split('=')[1].split(';')[0].strip()
    return None


class TestNoPackagedMdpUsesADeprecatedBarostat(unittest.TestCase):
    """This went unnoticed for as long as Gromacs has warned about it, because
    htpolynet runs grompp with -maxwarn 4 and discarded its output on success.  A test
    is the only thing that would have caught it, so here is the test."""

    def test_every_pressure_coupled_mdp_uses_a_correct_ensemble(self):
        bad = {n: setting(t, 'pcoupl') for n, t in mdps().items()
               if (setting(t, 'pcoupl') or 'no').lower() in DEPRECATED_BAROSTATS}
        self.assertEqual(bad, {})

    def test_the_constant_pressure_mdps_are_actually_coupled(self):
        # guard against "fixing" the warning by switching pressure coupling off
        for name in ('npt.mdp', 'drag-npt.mdp', 'relax-npt.mdp'):
            with self.subTest(name):
                self.assertEqual((setting(mdps()[name], 'pcoupl') or 'no').lower(),
                                 'c-rescale')

    def test_c_rescale_is_only_used_where_it_is_supported(self):
        # C-rescale does not support anisotropic coupling; grompp errors outright.
        # The deformation stages override pcoupltype at runtime and set their own
        # barostat, so no packaged mdp may pair the two.
        for name, text in mdps().items():
            with self.subTest(name):
                if (setting(text, 'pcoupl') or '').lower() == 'c-rescale':
                    self.assertNotEqual((setting(text, 'pcoupltype') or '').lower(),
                                        'anisotropic')


class TestTheDeformationStagesCarryTheirOwnBarostat(unittest.TestCase):
    """They force anisotropic coupling, which C-rescale cannot do."""

    def built(self, cls, params):
        import os
        import shutil
        import tempfile

        import numpy as np
        d = tempfile.mkdtemp()
        f = os.path.join(d, 'm.mdp')
        shutil.copy(str(ir.files('htpolynet').joinpath('resources/mdp/npt.mdp')), f)
        cls(params).build_mdp(f, box=np.diag([5.6, 5.6, 5.6]))
        return open(f).read()

    def test_deform_uses_parrinello_rahman(self):
        from htpolynet.analysis.postsim import PostSimDeform
        text = self.built(PostSimDeform, {'direction': 'x', 'edot': 1e-4, 'ps': 300})
        self.assertEqual(setting(text, 'pcoupl'), 'Parrinello-Rahman')
        self.assertEqual(setting(text, 'pcoupltype'), 'anisotropic')

    def test_clamped_shear_uses_parrinello_rahman(self):
        from htpolynet.analysis.postsim import PostSimShear
        text = self.built(PostSimShear, {'direction': 'xy', 'edot': 1e-4, 'ps': 300})
        self.assertEqual(setting(text, 'pcoupl'), 'Parrinello-Rahman')

    def test_isochoric_shear_has_no_barostat_at_all(self):
        from htpolynet.analysis.postsim import PostSimShear
        text = self.built(PostSimShear, {'direction': 'xy', 'edot': 1e-4, 'ps': 300,
                                         'coupling': 'isochoric'})
        self.assertEqual(setting(text, 'pcoupl'), 'no')
