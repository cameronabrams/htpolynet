"""Manages the htpolynet application, provides the command-line interface entry point.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""

import logging
import os
import shutil
import textwrap

import argparse as ap
import yaml

from importlib.resources import files as pkg_files
from pathlib import Path

from .analysis.analyze import analyze
from .analysis.plot import plots
from .analysis.postsim import postsim
from .core import projectfilesystem as pfs
from .core.runtime import Runtime
from .external import software
from .external.slurm import generate_script, _normalize_keys
from .utils.banner import banner_message
from .utils.inputcheck import input_check
from .utils.logsetup import setup_logging
from .utils.stringthings import my_logger
from .analysis.piercings import piercings
from .utils.vmd_viz import make_viz

logger = logging.getLogger(__name__)


def info(args: ap.Namespace):
    """Handles the info subcommand.

    Args:
        args (argparse.Namespace): parsed arguments
    """
    print('This is some information on your installed version of htpolynet')
    l = pfs.lib_setup()
    software.sw_setup()
    print(l.info())
    print(software.to_string())
    system_mols, cached_mols = pfs.get_molecule_info()
    print('Available molecules in system library:')
    for m in system_mols:
        print(f'   {m}')
    if cached_mols:
        print('Cached parameterized molecules:')
        for m in cached_mols:
            print(f'   {m}')
    possibles = l.get_example_names()
    print('Available examples using htpolynet fetch-example')
    for p in possibles:
        if p[0].isdecimal():
            print(f'   {p}')

def run(args: ap.Namespace):
    """Handles the run subcommand.

    Args:
        args (argparse.Namespace): parsed arguments
    """
    setup_logging(args.loglevel, diag=args.diag, no_banner=args.no_banner)
    my_logger('htpolynet runtime begins', logger.info)
    if args.restart:
        logger.warning(
            'Restart (-restart) is experimental and currently broken at the '
            'cure stage: pre-cure stages skip correctly, but the first '
            'topology update of the resumed cure iteration fails with an '
            'internal "temp_idx ... already claimed" assertion (insufficient '
            'in-memory state reconstruction).  Re-running from scratch is the '
            'safe path.'
        )
    # Ensure the user-library directory tree exists so pfs.exists() / pfs.checkout()
    # can find any monomer inputs we generate from SMILES.  Without this, the
    # generated mol2 files would land in projPath/lib/... where pfs does not look.
    os.makedirs(os.path.join(args.lib, pfs.Dirs.molecules_inputs), exist_ok=True)
    os.makedirs(os.path.join(args.lib, pfs.Dirs.molecules_parameterized), exist_ok=True)
    userlib = args.lib
    software.sw_setup()
    my_logger(software.to_string(), logger.info)
    pfs.pfs_setup(root=os.getcwd(), topdirs=pfs.Dirs.run_topdirs, verbose=True, projdir=args.proj, reProject=args.restart, userlibrary=userlib)
    a = Runtime(cfgfile=args.config, restart=args.restart)
    if args.param_only:
        a.generate_molecules(force_checkin=args.force_checkin, force_parameterization=args.force_parameterization)
    else:
        a.do_workflow(force_checkin=args.force_checkin, force_parameterization=args.force_parameterization)
    my_logger('htpolynet runtime ends', logger.info)

def _fetch_one(depot, fullname):
    """Fetches one example by name (always a self-contained .yaml).

    Args:
        depot (str): path to the example depot directory
        fullname (str): example name without extension
    """
    import shutil
    yaml_path = os.path.join(depot, f'{fullname}.yaml')
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f'No example found for {fullname}')
    dest = os.path.basename(yaml_path)
    shutil.copy(yaml_path, dest)
    print(f'Fetched {dest}  (run with: htpolynet run {dest})')

def fetch_example(args):
    """Handles the fetch-example subcommand.

    Args:
        args (argparse.Namespace): parsed arguments
    """
    l = pfs.system()
    depot = l.get_example_depot_location()
    possibles = l.get_example_names()
    if args.n == 'all':
        for fullname in possibles:
            _fetch_one(depot, fullname)
        return
    if args.n.isdecimal():
        fullname = [x for x in possibles if x.startswith(args.n)][0]
    else:
        fullname = args.n
    _fetch_one(depot, fullname)


def setup_claude(args: ap.Namespace):
    """Handles the setup-claude subcommand.

    Copies the bundled skill into the user's skills directory.  This runs only
    when the user asks for it; installing the package never touches
    ``~/.claude/``.

    Args:
        args (argparse.Namespace): parsed arguments
    """
    # Resolve through htpolynet.resources, which is a real package; the claude/
    # subdirectory carries no __init__.py, matching every other resource
    # subdirectory.
    source = pkg_files('htpolynet.resources').joinpath('claude', 'SKILL.md')
    if not source.is_file():
        raise FileNotFoundError(f'Bundled skill not found at {source}')

    skill_dir = Path(args.skills_dir).expanduser() / 'htpolynet'
    target = skill_dir / 'SKILL.md'

    if target.exists() and not args.force:
        print(f'{target} already exists; re-run with --force to overwrite it')
        return

    skill_dir.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as fh, open(target, 'wb') as out:
        shutil.copyfileobj(fh, out)
    print(f'Written: {target}')
    print('Claude Code picks the skill up on its next session in any directory.')


def _add_run_options(p, loglevel='debug'):
    """Adds run options (no positional config) shared by run and gen-slurm-script."""
    p.add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    p.add_argument('-proj',type=str,default='next',help='project directory; "next" generates the next available; otherwise creates or resumes (with -restart) a named directory')
    p.add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    p.add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir (EXPERIMENTAL: broken at the cure stage; see usage docs)')
    p.add_argument('--no-banner',default=False,action='store_true',help='suppress the startup banner')
    p.add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    p.add_argument('--force-checkin',default=False,action='store_true',help='overwrite entries already in the user library (~/.htpolynet, or $HTPOLYNET_CACHE) with newly generated parameter files; a name the library does not yet hold is checked in either way')
    p.add_argument('--param-only',default=False,action='store_true',help='stop after parameterizing molecules without building the full system')
    p.add_argument('--loglevel',type=str,default=loglevel,help='log level for the diagnostic file (debug|info)')


def _add_runtime_args(p, loglevel='debug'):
    """Adds the common arguments shared by the run and parameterize subcommands."""
    p.add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    _add_run_options(p, loglevel=loglevel)

def gen_slurm_script(args: ap.Namespace):
    """Handles the gen-slurm-script subcommand.

    Args:
        args (argparse.Namespace): parsed arguments
    """
    # Load YAML config to pick up the optional slurm: section
    slurm_cfg = {}
    if os.path.exists(args.config):
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh) or {}
        slurm_cfg = _normalize_keys(cfg.get('slurm', {}))

    # CLI overrides take precedence over values from the YAML slurm: section
    cli_overrides = {
        'job_name':        args.job_name,
        'partition':       args.partition,
        'account':         args.account,
        'nodes':           args.nodes,
        'ntasks':          args.ntasks,
        'ntasks_per_node': args.ntasks_per_node,
        'cpus_per_task':   args.cpus_per_task,
        'time':            args.time,
        'mem':             args.mem,
        'gres':            args.gres,
    }
    for k, v in cli_overrides.items():
        if v is not None:
            slurm_cfg[k] = v

    # --sif / --apptainer-alias override YAML values
    sif = args.sif or slurm_cfg.get('sif')
    apptainer_alias = args.apptainer_alias or slurm_cfg.get('apptainer_alias', 'apptainer')

    if not args.native and not sif:
        raise SystemExit(
            'error: containerized mode requires --sif <path/to/image.sif> '
            '(or a sif: entry in the slurm: section of the config file)'
        )

    run_args = {
        'lib':                    args.lib,
        'proj':                   args.proj,
        'diag':                   args.diag,
        'restart':                args.restart,
        'no_banner':              args.no_banner,
        'force_parameterization': args.force_parameterization,
        'force_checkin':          args.force_checkin,
        'param_only':             args.param_only,
        'loglevel':               args.loglevel,
    }

    script = generate_script(
        config_file=args.config,
        slurm_cfg=slurm_cfg,
        run_args=run_args,
        native=args.native,
        sif=sif,
        apptainer_alias=apptainer_alias,
    )

    outfile = args.output or f'submit_{os.path.splitext(os.path.basename(args.config))[0]}.sh'
    with open(outfile, 'w') as fh:
        fh.write(script)
    print(f'Wrote SLURM script: {outfile}')


def _add_gen_slurm_args(p):
    """Adds arguments for the gen-slurm-script subcommand."""
    p.add_argument('config', type=str, help='htpolynet YAML config file')
    # run options (mirrors _add_run_options)
    _add_run_options(p, loglevel='info')
    # execution mode
    p.add_argument('--native', default=False, action='store_true',
                   help='run htpolynet natively instead of via Apptainer container')
    p.add_argument('--sif', type=str, default=None,
                   help='path to the Apptainer/Singularity .sif image (required unless --native)')
    p.add_argument('--apptainer-alias', type=str, default=None, dest='apptainer_alias',
                   help='container runtime command (default: apptainer)')
    # output
    p.add_argument('-o', '--output', type=str, default=None,
                   help='output script filename (default: submit_<config>.sh)')
    # common SLURM directives as CLI overrides
    p.add_argument('--job-name',        type=str, default=None, dest='job_name')
    p.add_argument('--partition',       type=str, default=None)
    p.add_argument('--account',         type=str, default=None)
    p.add_argument('--nodes',           type=int, default=None)
    p.add_argument('--ntasks',          type=int, default=None)
    p.add_argument('--ntasks-per-node', type=int, default=None, dest='ntasks_per_node')
    p.add_argument('--cpus-per-task',   type=int, default=None, dest='cpus_per_task')
    p.add_argument('--time',            type=str, default=None)
    p.add_argument('--mem',             type=str, default=None)
    p.add_argument('--gres',            type=str, default=None)


def _add_analysis_args(p):
    """Adds the common arguments shared by the postsim and analyze subcommands."""
    p.add_argument('-proj',type=str,default='',nargs='+',help='name of project directory')
    p.add_argument('-lib',type=str,default=None,help='local user library of molecular structures and parameterizations (default: ./lib if it exists, otherwise none)')
    p.add_argument('-ocfg',type=str,default='',help='original htpolynet config file used to generate project(s)')
    p.add_argument('-cfg',type=str,default='',help='config file specifying the operations to perform')
    p.add_argument('--no-banner',default=False,action='store_true',help='suppress the startup banner')
    p.add_argument('--loglevel',type=str,default='info',help='log level for the diagnostic file (debug|info)')

def cli():
    """Command-line interface."""
    l = pfs.lib_setup()
    example_names = l.get_example_names()
    example_ids = [x.split('-')[0] for x in example_names]

    subcommands = [
        ('run',              run,              'build a system using instructions in the config file and any required molecular structure inputs (use --param-only to stop after parameterization)'),
        ('info',             info,             'print version, software, and library information'),
        ('plots',            plots,            'generate plots summarizing aspects of the current completed build'),
        ('fetch-example',    fetch_example,    'fetch and unpack example(s): '+', '.join(f'"{x}"' for x in example_names)),
        ('input-check',      input_check,      'report the number of atoms that would be in the initial system based on config'),
        ('postsim',          postsim,          'perform specified post-cure MD simulations on final results in one or more project directories'),
        ('analyze',          analyze,          "perform 'gmx <command>' style analyses specified in the config file"),
        ('gen-slurm-script', gen_slurm_script, 'generate a SLURM submission script for running htpolynet on a cluster'),
        ('make-viz',         make_viz,         'regenerate VMD viz files (.viz.psf + .viz.tcl) from an existing gromacs top + gro pair'),
        ('piercings',        piercings,        'report any bond threaded through a ring in a gromacs top + gro pair; exits 1 if any are found'),
        ('setup-claude',     setup_claude,     "install htpolynet's Claude Code skill so an agent can drive htpolynet"),
    ]

    parser = ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers()
    subparsers.required = True
    cp = {}
    for name,handler,help_str in subcommands:
        cp[name] = subparsers.add_parser(name,help=help_str)
        cp[name].set_defaults(func=handler)

    _add_runtime_args(cp['run'])
    _add_gen_slurm_args(cp['gen-slurm-script'])

    cp['plots'].add_argument('source', type=str, choices=['diag', 'build', 'post'], default='build', help='source of data to plot')
    cp['plots'].add_argument('--diags', type=str, default=[], nargs='+', help='names of diagnostic log files (1 or more)')
    cp['plots'].add_argument('--proj', nargs='+', type=str, default=[], help='name of project director[y/ies]')
    cp['plots'].add_argument('--cfg', type=str, nargs='+', default=[], help='name of input config files')
    cp['plots'].add_argument('--buildplot', type=str, nargs='+', default=['t'], choices=['t', 'g', 'n', 'c'], help='build plot type: t=traces, g=2-D graph, n=homo-N between crosslinks, c=cluster-size distributions')
    cp['plots'].add_argument('--traces', type=str, nargs='+', default=['t', 'd', 'p'], choices=['t', 'd', 'p'], help='trace types: t=temperature, d=density, p=potential energy')
    cp['plots'].add_argument('--n_points', type=int, nargs=2, default=[10, 20], help='number of [cold-side,hot-side] data points for Tg line fits')
    cp['plots'].add_argument('--plotfile', type=str, default='', help='name of output plot file')
    cp['plots'].add_argument('--no-banner', default=False, action='store_true', help='suppress the startup banner')
    cp['plots'].add_argument('--loglevel', type=str, default='info', help='log level for the diagnostic file (debug|info)')

    cp['fetch-example'].add_argument('n', type=str, choices=example_ids+['all'], help='example to fetch: '+', '.join(example_names))

    cp['input-check'].add_argument('config', type=str, default=None, help='input configuration file in YAML format')
    cp['input-check'].add_argument('-lib', type=str, default='lib', help='local user library of molecular structures and parameterizations')

    cp['make-viz'].add_argument('-top', type=str, default='final.top', help='input gromacs topology file (default: final.top)')
    cp['make-viz'].add_argument('-gro', type=str, default='final.gro', help='input gromacs coordinate file (default: final.gro)')
    cp['make-viz'].add_argument('-grx', type=str, default=None, help='input htpolynet .grx (default: auto-detect <gro-stem>.grx; needed for the constituent-selection macros)')
    cp['make-viz'].add_argument('-prefix', type=str, default=None, help='output basename; the .viz.psf, .viz.tcl, and .viz.macros.tcl are written next to the input gro (default: stem of -gro)')

    cp['piercings'].add_argument('-top', type=str, default='final.top', help='input gromacs topology file (default: final.top)')
    cp['piercings'].add_argument('-gro', type=str, default='final.gro', help='input gromacs coordinate file (default: final.gro)')
    cp['piercings'].add_argument('-max-ring', dest='max_ring', type=int, default=8, help='largest ring to consider when the topology carries no ring list (default: 8)')
    cp['piercings'].add_argument('-json', type=str, default=None, help='also write the findings to this JSON file')

    cp['setup-claude'].add_argument('--skills-dir', type=str, default='~/.claude/skills',
                                    help='skills directory to install into (default: %(default)s); use ./.claude/skills to scope the skill to one project')
    cp['setup-claude'].add_argument('--force', default=False, action='store_true', help='overwrite an existing installed skill')

    _add_analysis_args(cp['postsim'])
    _add_analysis_args(cp['analyze'])

    args = parser.parse_args()
    args.func(args)
