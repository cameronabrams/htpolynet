"""

.. module:: test_mdrun_gpu_options
   :synopsis: one set of mdrun GPU options must work for minimizations as well as dynamics

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from importlib.resources import files

import pytest

from htpolynet.external.gromacs import ensure_thread_mpi_ranks, mdp_integrator, mdrun_options_for

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


class TestThreadMpiRanks:
    """GROMACS refuses OpenMP threads plus a GPU without a rank count, which
    killed a GPU build at its first minimization (htpolynet-study, 2026-09-18)."""

    def test_gpu_id_and_ntomp_get_one_rank(self):
        got = ensure_thread_mpi_ranks({'gpu_id': 0, 'nb': 'gpu', 'ntomp': 12}, 'gmx mdrun')
        assert got['ntmpi'] == 1 and got['ntomp'] == 12

    @pytest.mark.parametrize('task', ['nb', 'pme', 'bonded', 'update'])
    def test_any_task_on_the_gpu_counts(self, task):
        assert ensure_thread_mpi_ranks({task: 'gpu', 'ntomp': 8}, 'gmx mdrun')['ntmpi'] == 1

    def test_no_gpu_asked_for_is_left_alone(self):
        opts = {'ntomp': 8, 'nb': 'cpu'}
        assert ensure_thread_mpi_ranks(opts, 'gmx mdrun') == opts

    def test_no_ntomp_is_left_alone(self):
        opts = {'gpu_id': 0, 'nb': 'gpu'}
        assert ensure_thread_mpi_ranks(opts, 'gmx mdrun') == opts

    @pytest.mark.parametrize('given', ['ntmpi', 'nt', 'npme'])
    def test_a_rank_count_already_given_wins(self, given):
        opts = {'gpu_id': 0, 'ntomp': 12, given: 4}
        assert ensure_thread_mpi_ranks(opts, 'gmx mdrun') == opts

    @pytest.mark.parametrize('cmd', ['mpirun -np 4 gmx_mpi mdrun', 'srun gmx_mpi mdrun',
                                     'gmx_mpi mdrun', 'mpiexec -n 2 mdrun_mpi'])
    def test_an_mpi_mdrun_is_left_alone(self, cmd):
        # gmx_mpi takes ranks from the launcher and rejects -ntmpi
        opts = {'gpu_id': 0, 'ntomp': 12}
        assert ensure_thread_mpi_ranks(opts, cmd) == opts

    def test_the_configuration_is_not_mutated(self):
        opts = {'gpu_id': 0, 'ntomp': 12}
        ensure_thread_mpi_ranks(opts, 'gmx mdrun')
        assert 'ntmpi' not in opts
