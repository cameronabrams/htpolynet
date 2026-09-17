# HTPolyNet container image
# Includes Gromacs, AmberTools (antechamber, tleap, parmchk2), and OpenBabel
# via conda-forge.
#
# The entrypoint auto-detects the host owner of the /work bind mount and drops
# privileges via gosu before invoking htpolynet, so the caller does not need
# --user, HOST_UID/HOST_GID env vars, or a matching /etc/passwd entry.
#
# Build context must be the repo root (this Dockerfile lives in docker/):
#
#   docker build -f docker/Dockerfile -t ghcr.io/cameronabrams/htpolynet:latest .
#
# Recommended usage is via Docker Compose; see docker/compose.yml and
# docs/source/user-guide/container-usage.rst for the full story (incl. GPU
# and Singularity/Apptainer on HPC).
#
#   docker compose run --rm htpolynet run config.yaml
#
# Raw `docker run` equivalent:
#
#   docker run --rm -v $(pwd):/work ghcr.io/cameronabrams/htpolynet run config.yaml
#
# GPUs: the DEFAULT image cannot use one.  Its Gromacs is the conda-forge
# OpenCL build, and Gromacs no longer drives NVIDIA devices through OpenCL, so
# --gpus all starts the container and changes nothing about how it computes.
#
# Build the CUDA variant instead, which installs gromacs=*=nompi_cuda*:
#
#   CONDA_OVERRIDE_CUDA=12.9 docker build -f docker/Dockerfile \
#       --build-arg GROMACS_BUILD=nompi_cuda \
#       --build-arg CONDA_OVERRIDE_CUDA=12.9 \
#       -t ghcr.io/cameronabrams/htpolynet:cuda .
#
# CONDA_OVERRIDE_CUDA is required to build on a machine with no NVIDIA driver:
# the CUDA package depends on the __cuda virtual package, which conda
# synthesizes only where a driver is present, so the solve fails on an
# ordinary CI runner without it.  Running the resulting image does need
# nvidia-container-toolkit and --gpus all.

FROM condaforge/miniforge3:latest

# Links the published GHCR package to its source repository.
LABEL org.opencontainers.image.source="https://github.com/cameronabrams/htpolynet"

# Install Gromacs (latest available on conda-forge); AmberTools follows below.
# Miniforge ships with mamba for fast solves; conda-forge is the only
# channel.  No Anaconda Inc. ToS dependency.
RUN apt-get update && apt-get install -y --no-install-recommends \
        openbabel \
        graphviz \
        gosu \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Which Gromacs build to install.  Empty (the default) takes whatever
# conda-forge considers best, which for linux-64 is the OpenCL build -- correct
# for a CPU image, and unable to drive an NVIDIA device.  Pass
# GROMACS_BUILD=nompi_cuda for the CUDA image; that package depends on the
# __cuda virtual package, so a build host without an NVIDIA driver also needs
# CONDA_OVERRIDE_CUDA set to a version to satisfy the solve.
ARG GROMACS_BUILD=
ARG CONDA_OVERRIDE_CUDA=

# Clean with conda, not mamba.  mamba 2.9.0 (miniforge3 26.7.2-0, the image
# :latest moved to on 2026-09-09) fails `mamba clean -afy` with "cannot get
# file size: No such file or directory" on the CPU AmberTools 26.0 package:
# its bin/amber.conda is a symlink that dangles inside the package cache, and
# mamba stats everything named *.conda.  The install itself succeeds; only the
# clean dies, taking the whole RUN with it.  conda clean handles the same cache
# correctly.  Reproduced locally in both directions before changing this.
RUN CONDA_OVERRIDE_CUDA="${CONDA_OVERRIDE_CUDA}" mamba install -y -n base -c conda-forge \
        "gromacs${GROMACS_BUILD:+=*=${GROMACS_BUILD}*}" \
        parmed \
        rdkit \
    && conda clean -afy

# AmberTools lives in its own environment, pinned, in both images.  Solved
# together with the CUDA Gromacs, conda-forge silently fell back to AmberTools
# 24.8 (the CUDA toolkit's nsight-compute needs an older krb5 than AmberTools
# 26 allows), so the CPU and CUDA images parameterized molecules with different
# AmberTools releases.  A separate environment keeps the two solves apart, and
# the pin keeps the images identical and makes a future conflict a build error
# instead of a quiet downgrade.  The environment is not activated: its
# wrappers (antechamber, tleap, parmchk2) set AMBERHOME from their own
# location, so putting its bin last on PATH is enough, and base's python and
# gmx still come first.
ARG AMBERTOOLS_VERSION=26.0
RUN mamba create -y -n amber -c conda-forge "ambertools=${AMBERTOOLS_VERSION}" \
    && conda clean -afy
ENV PATH=${PATH}:/opt/conda/envs/amber/bin

# Install HTPolyNet
COPY . /htpolynet
WORKDIR /htpolynet
RUN pip install --no-cache-dir .

# The commit this image was built from.  There is no git in the image and no
# .git beside the installed package, so without this an image cannot say what
# code it contains -- and its version string is actively misleading between
# releases, since the weekly rebuild builds from main HEAD.  Declared after the
# expensive layers so that a changed sha does not invalidate the mamba cache.
ARG HTPOLYNET_COMMIT=unknown
ENV HTPOLYNET_COMMIT=${HTPOLYNET_COMMIT}

# Pre-create HOME for the unprivileged user the entrypoint drops to.  When
# compose.yml mounts a named volume here, Docker initializes the fresh volume
# from this directory and preserves the 0777 mode.
RUN mkdir -p /home/htpolynet && chmod 0777 /home/htpolynet

COPY docker/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# User data is mounted here at runtime
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["--help"]
