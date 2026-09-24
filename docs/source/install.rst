##############################
Installation and Prerequisites
##############################

.. admonition:: Recommended: run htpolynet from the container
   :class: tip

   The easiest way to use ``htpolynet`` is the pre-built container image,
   ``ghcr.io/cameronabrams/htpolynet``.  It bundles everything ``htpolynet``
   needs --- Gromacs, AmberTools, OpenBabel, RDKit, and ``htpolynet`` itself
   --- so there is nothing else to install and no versions to line up.

   **On a desktop** (Docker), fetch the Compose file once, then run
   ``htpolynet`` through it:

   .. code-block:: console

      $ curl -O https://raw.githubusercontent.com/cameronabrams/htpolynet/main/docker/compose.yml
      $ docker compose run --rm htpolynet run config.yaml

   **On a cluster** (Singularity/Apptainer), pull the image once, then run it
   with your working directory bound in:

   .. code-block:: console

      $ singularity pull htpolynet.sif docker://ghcr.io/cameronabrams/htpolynet:latest
      $ singularity run --bind $(pwd):/work --pwd /work htpolynet.sif run config.yaml

   The default image runs Gromacs on CPUs only.  To use an NVIDIA GPU, use
   the ``:cuda`` tag instead.

   Nothing is installed on your system, and the image is always current ---
   which the ``conda-forge`` package is not; see
   :ref:`the note below <conda_forge_not_recommended>`.

   See :ref:`container_usage` for the details: persistent caches, GPU setup,
   and pinning a release.

The rest of this page covers installing ``htpolynet`` and its prerequisites
directly on your system, which you'll want for development or wherever
containers aren't available.

Software Prerequisites
----------------------

The following commands need to be on your ``PATH``:

1. ``antechamber``, ``parmchk2``, and ``tleap`` (`AmberTools
   <https://ambermd.org/GetAmber.php#ambertools>`_, version 22 or
   higher); the most convenient source is the ``conda-forge`` channel.
2. ``gmx`` or ``gmx_mpi`` (`Gromacs
   <https://manual.gromacs.org/documentation/current/index.html>`_,
   version 2022.1 or higher); available from ``conda-forge``, your
   distribution's package manager, or compiled from source.
3. ``obabel`` (`OpenBabel
   <https://github.com/openbabel/openbabel>`_); preferred installation
   via your Linux distribution's package manager.  OpenBabel is a
   required runtime dependency whenever you let ``htpolynet`` build
   monomer structures from SMILES strings (the recommended workflow —
   see :ref:`molecular_structure_inputs`).  RDKit on its own is not
   sufficient because it has no ``mol2`` writer, so ``htpolynet``
   always shells out to ``obabel`` for the final SDF→mol2 conversion.
   The only way to run ``htpolynet`` without ``obabel`` is to supply
   hand-prepared ``mol2``/``pdb`` files for every monomer.
4. ``dot`` (`Graphviz <https://graphviz.org/>`_); used to render the
   reaction-network plot ``plots/reaction_network.png`` at setup time.
   If ``dot`` is not on the ``PATH`` the build still proceeds; the
   plot is skipped with a warning.

`RDKit <https://www.rdkit.org/>`_ is a required Python dependency
(installed automatically with ``htpolynet``); it backs the
atom-mapping SMILES syntax that every depot example uses.  The
index-keyed ``rename_atoms`` form of the SMILES spec is still
supported if you don't need atom-map labels.

Installation
------------

Two supported workflows: **uv + Miniforge** (recommended for development
and for any system where you'd like to keep Python and the native MD
binaries cleanly separated), and **conda-only** (simpler one-stop
install if you already manage everything through ``conda-forge``).

uv + Miniforge (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Python side, use `uv <https://docs.astral.sh/uv/>`_ to manage a
per-project virtual environment.  For the native MD binaries
(AmberTools and Gromacs), use a `Miniforge
<https://github.com/conda-forge/miniforge>`_-managed conda-forge
environment.  Miniforge is community-maintained, uses only the
``conda-forge`` channel, and carries no Anaconda Inc. terms-of-service
encumbrance.

Install uv:

.. code-block:: console

    $ curl -LsSf https://astral.sh/uv/install.sh | sh

Install Miniforge (the installer at
https://github.com/conda-forge/miniforge#install includes one-liners
for every platform), then create the MD-tools environments:

.. code-block:: console

    $ mamba create -n gromacs    -c conda-forge gromacs
    $ mamba create -n ambertools -c conda-forge ambertools parmed

(Splitting AmberTools and Gromacs into separate environments avoids
solver conflicts between their MPI / CUDA-linkage variants.  Combine
them into one env if your platform doesn't trip on that.)

Add the env ``bin`` directories to your ``PATH`` so ``gmx``,
``antechamber``, ``tleap``, and ``parmchk2`` are always available
without needing to activate the envs:

.. code-block:: bash

    # ~/.bashrc (append-only so uv's python isn't shadowed)
    export PATH="$HOME/miniforge3/envs/gromacs/bin:$HOME/miniforge3/envs/ambertools/bin:$PATH"
    export AMBERHOME="$HOME/miniforge3/envs/ambertools"

Install ``htpolynet`` into a uv-managed virtualenv:

.. code-block:: console

    $ git clone git@github.com:cameronabrams/htpolynet.git
    $ cd htpolynet
    $ uv venv
    $ uv pip install -e .

For a global ``htpolynet`` command callable from any directory:

.. code-block:: console

    $ uv tool install --editable .

This installs ``htpolynet`` as a uv-managed tool with its own
dedicated environment and a shim in ``~/.local/bin/`` that's available
from any shell.

If your distribution doesn't already provide ``obabel`` and ``dot``,
install them from its package manager.  On Debian/Ubuntu:

.. code-block:: console

    $ sudo apt install openbabel graphviz

On openSUSE / Fedora / RHEL:

.. code-block:: console

    $ sudo zypper install openbabel graphviz    # openSUSE
    $ sudo dnf install openbabel graphviz       # Fedora/RHEL

Conda-only (simpler one-stop install)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _conda_forge_not_recommended:

.. admonition:: Do not install ``htpolynet`` itself from conda-forge
   :class: warning

   The ``htpolynet`` package on ``conda-forge`` is several releases behind
   and is **not currently recommended**.  What is published there installs
   and runs, so nothing warns you --- you simply get an old ``htpolynet``.

   It cannot be updated at the moment.  ``htpolynet`` now depends on
   `ycleptic <https://pypi.org/project/ycleptic/>`_, which is not yet
   available on ``conda-forge``; until it is, the feedstock's automatic
   version bumps produce a package that cannot import, so none has been
   published.

   You do not have to take this page's word for how far behind it is ---
   the page may itself be stale by the time you read it:

   .. code-block:: console

      $ conda search -c conda-forge htpolynet | tail -1   # what conda-forge has
      $ pip index versions htpolynet                      # what is current

   This applies only to the ``htpolynet`` package.  ``conda-forge`` remains
   the recommended source for AmberTools and Gromacs, as described above.

Using Miniforge (preferred) or any other ``conda``-style installer, you
can manage Python and the MD binaries in a single environment, then add
``htpolynet`` itself from PyPI:

.. code-block:: console

    $ mamba create -n htpolynet -c conda-forge python ambertools gromacs
    $ mamba activate htpolynet
    $ pip install htpolynet

Then install OpenBabel from your distribution's package manager as
described above (OpenBabel via ``conda-forge`` works too but the
distribution package is usually fresher and avoids a Python-extension
ABI link).

PyPI install (Python-only)
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you've installed the native MD binaries some other way, the
Python-only install of ``htpolynet`` from PyPI is:

.. code-block:: console

    $ pip install htpolynet

You're responsible for ensuring ``antechamber``, ``parmchk2``,
``tleap``, ``gmx``, and ``obabel`` are reachable on ``PATH``.

Compiling MD tools from source
------------------------------

If your distribution doesn't ship recent enough versions and you don't
want to use a conda-forge build, the native MD tools can be compiled
from source.  Below are reference build recipes.

AmberTools
^^^^^^^^^^

Requires ``csh``, ``flex``, and ``bison``:

.. code-block:: console

    $ tar jxf AmberTools24.tar.bz2
    $ cd amber_src
    $ ./configure --no-X11 --skip-python gnu
    $ source amber.sh
    $ make install

Gromacs
^^^^^^^

Reference CUDA-enabled (single-replica, no MPI) build:

.. code-block:: console

    $ tar xfz gromacs-2025.4.tar.gz
    $ cd gromacs-2025.4
    $ mkdir build
    $ cd build
    $ cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON \
               -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs
    $ make
    $ make check
    $ sudo make install

Add to your ``~/.bashrc``:

.. code-block:: bash

    source /usr/local/gromacs/bin/GMXRC

OpenBabel
^^^^^^^^^

Be sure to unpack `Eigen
<https://eigen.tuxfamily.org/index.php?title=Main_Page>`_ first so the
``conformer`` plug-in builds.  Example session where Eigen and
OpenBabel sources are in ``~/Downloads`` and the install prefix is
``~/opt/obabel``:

.. code-block:: console

    $ cd ~/build
    $ tar jxf ~/Downloads/eigen-3.4.0.tar.bz2
    $ tar jxf ~/Downloads/openbabel-3.1.1.tar.bz2
    $ cd openbabel-3.1.1 && mkdir build && cd build
    $ cmake .. -DEIGEN3_INCLUDE_DIR=${HOME}/build/eigen-3.4.0/ \
               -DCMAKE_INSTALL_PREFIX=${HOME}/opt/obabel
    $ make && make test && make install

Then set ``PATH``, ``LD_LIBRARY_PATH``, and ``BABEL_LIBDIR`` to point
at ``${HOME}/opt/obabel``.

Driving htpolynet with Claude Code
----------------------------------

``htpolynet`` ships a skill for `Claude Code
<https://claude.com/claude-code>`_ that teaches the agent how to use it:
start from the nearest example, describe monomers in their active form,
check a configuration before spending compute, and recognize the failure
modes that are known rather than mysterious.  Install it once, after
installing the package:

.. code-block:: console

    $ htpolynet setup-claude

Installing ``htpolynet`` never writes to ``~/.claude/`` on its own; the
skill is copied only when you run this.  See :doc:`/user-guide/usage` for the options,
including how to scope the skill to a single project.

Other Prerequisites
-------------------

To use ``htpolynet`` effectively, working knowledge of the following is
helpful:

1. MD simulation in general and Gromacs specifically.
2. The General Amber Force Field (GAFF), including

   a. how to use ``antechamber``, ``tleap``, and ``parmchk2`` to
      generate GAFF parameterizations; and
   b. how to use those parameterizations inside Gromacs.

3. Polymer chemistry, at least for the systems you intend to simulate.
