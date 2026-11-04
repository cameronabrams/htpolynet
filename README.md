# htpolynet
> High-Throughput Polymer Network Atomistic Simulations

[![tests](https://github.com/cameronabrams/htpolynet/actions/workflows/test.yml/badge.svg)](https://github.com/cameronabrams/htpolynet/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/htpolynet.svg)](https://pypi.org/project/htpolynet/)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/htpolynet)](https://anaconda.org/conda-forge/htpolynet)
[![Python](https://img.shields.io/pypi/pyversions/htpolynet)](https://pypi.org/project/htpolynet/)
[![License: MIT](https://img.shields.io/pypi/l/htpolynet)](https://github.com/cameronabrams/htpolynet/blob/main/LICENSE)
[![Docs](https://readthedocs.org/projects/htpolynet/badge/?version=latest)](https://htpolynet.readthedocs.io/en/latest/)
[![PyPI Downloads](https://static.pepy.tech/badge/htpolynet)](https://pepy.tech/projects/htpolynet)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22070252-blue)](https://doi.org/10.5281/zenodo.22070252)

htpolynet is a Python utility for generating atomistic models of cross-linked polymer networks together with appropriate topology and parameter files required for molecular dynamics simulations using Gromacs.  It is intended as a fully automated system builder requiring as inputs only the molecular structures of any monomer species, a description of the polymerization chemistry, and a handful of options describing desired system size and composition.  htpolynet uses the Generalized Amber Force Field for atom-typing and parameter generation.

## Installation

From PyPI:
```bash
pip install htpolynet
```

From conda-forge:
```bash
conda install -c conda-forge htpolynet
```

From source:
```bash
git clone git@github.com:cameronabrams/htpolynet.git
cd htpolynet
pip install -e .
```

Once installed, the user has access to the main ``htpolynet`` command.

If you drive htpolynet with [Claude Code](https://claude.com/claude-code), install the bundled skill so the agent knows how to use it:
```bash
htpolynet setup-claude
```
This writes `~/.claude/skills/htpolynet/SKILL.md`; nothing is installed there unless you run it.

IMPORTANT NOTES: The programs ``antechamber``, ``parmchk2`` and ``tleap`` from AmberTools must be in your path.  These can be installed using the ``ambertools`` package from ``conda-forge`` or compiled from source.  You also need Gromacs installed so ``gmx`` is in your path.  The examples show how to build input monomer structures using OpenBabel, so to use them you need ``obabel`` in your path as well.

## Docker

As an alternative to a local installation, a prebuilt container image is published at ``ghcr.io/cameronabrams/htpolynet``.  It bundles htpolynet together with Gromacs, AmberTools, and OpenBabel, so no additional dependencies are required on the host beyond Docker.

Run htpolynet against a configuration file in the current directory:
```bash
docker run --rm -v $(pwd):/work ghcr.io/cameronabrams/htpolynet run config.yaml
```

**The image cannot use a GPU.**  Its Gromacs comes from conda-forge, built against OpenCL rather than CUDA, and Gromacs no longer drives NVIDIA devices through OpenCL.  Passing `--gpus all` starts the container and changes nothing about how it computes; on a cluster, target CPU partitions and do not request `--gres=gpu` or pass `--nv`.  If you need GPU-accelerated Gromacs, install htpolynet natively against a CUDA-enabled Gromacs.

A Docker Compose file is also provided in [docker/compose.yml](docker/compose.yml) for a shorter invocation (``docker compose run --rm htpolynet run config.yaml``).  See [docs/source/user-guide/container-usage.rst](docs/source/user-guide/container-usage.rst) for the full story, including Singularity/Apptainer use on HPC systems.

## Documentation

Please consult documentation at [htpolynet.readthedocs.io](https://htpolynet.readthedocs.io/).

## Repository relocation

This repository formerly lived at `AbramsGroup/HTPolyNet` and now lives at
[cameronabrams/htpolynet](https://github.com/cameronabrams/htpolynet).  GitHub redirects
the old URLs, so existing clones and forks continue to work; if you prefer, you can
update your remote explicitly:

```bash
git remote set-url origin git@github.com:cameronabrams/htpolynet.git
```

## Acknowledgments

htpolynet grew out of the original HTPolyNet prototype begun by Ming Huang in 2020.
Ketan S. Khare contributed early LAMMPS-related utilities, and S. Alexis Paz contributed
a bug fix.  The current package is a full rewrite, but the project owes its origins and
its published description to that earlier work.  When using htpolynet in published work,
please cite Huang and Abrams, *SoftwareX* **21**, 101303 (2023),
[doi:10.1016/j.softx.2022.101303](https://doi.org/10.1016/j.softx.2022.101303), along with
the GAFF and Gromacs papers listed in the [documentation](https://htpolynet.readthedocs.io/).
Papers by other groups that have used htpolynet are listed on the
[Publications Using HTPolyNet](https://htpolynet.readthedocs.io/en/latest/publications.html) page.

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/htpolynet/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

