"""

.. module:: test_mdrun_gpu_options
   :synopsis: one set of mdrun GPU options must work for minimizations as well as dynamics

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from importlib.resources import files

import pytest

from htpolynet.external.gromacs import mdp_integrator, mdrun_options_for

GPU = {'nb': 'gpu', 'pme': 'gpu', 'bonded': 'gpu', 'update': 'gpu', 'ntomp': 8}


class TestMdrunOptionsFor:
    def test_dynamics_keeps_everything(self):
        assert mdrun_options_for(GPU, 'md') == GPU

    @pytest.mark.parametrize('integrator', ['steep', 'cg', 'l-bfgs'])
    def test_minimizer_downgrades_only_what_gromacs_rejects(self, integrator):
        got = mdrun_options_for(GPU, integrator)
        assert got['pme'] == 'auto' and got['update'] == 'auto' and got['bonded'] == 'auto'
        assert got['nb'] == 'gpu' and got['ntomp'] == 8

    def test_explicit_cpu_is_left_alone(self):
        assert mdrun_options_for({'pme': 'cpu', 'update': 'cpu'}, 'steep') == {'pme': 'cpu', 'update': 'cpu'}

    def test_the_configuration_is_not_mutated(self):
        opts = dict(GPU)
        mdrun_options_for(opts, 'steep')
        assert opts == GPU


class TestMdpIntegrator:
    mdp = files('htpolynet.resources') / 'mdp'

    @pytest.mark.parametrize('name,expected', [('min', 'steep'), ('drag-min', 'steep'), ('relax-min', 'steep'),
                                               ('npt', 'md'), ('relax-nvt', 'md')])
    def test_shipped_mdps(self, name, expected):
        # min.mdp writes "integrator = steep;" with the semicolon attached
        assert mdp_integrator(str(self.mdp / f'{name}.mdp')) == expected

    def test_absent_means_md(self, tmp_path):
        f = tmp_path / 'x.mdp'
        f.write_text('nsteps = 10\n')
        assert mdp_integrator(str(f)) == 'md'
