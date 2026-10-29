#!/bin/bash
# Environment for the standalone PIC/electrostatics work in HET/pic (independent of bplsim).
# micromamba (user space, no sudo) + a conda-forge env with WarpX, gmsh and the Python stack.
#   bash pic/setup_env.sh          then:  source pic/activate.sh
set -e
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
mkdir -p "$HOME/.local/bin"
if [ ! -x "$HOME/.local/bin/micromamba" ]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME/.local" bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
micromamba --version
ENV="$MAMBA_ROOT_PREFIX/envs/hetpic"
if [ ! -d "$ENV" ]; then
  micromamba create -y -p "$ENV" -c conda-forge \
    python=3.12 warpx numpy scipy matplotlib h5py shapely gmsh python-gmsh pandas
fi
"$ENV/bin/python" -c "import pywarpx, gmsh, numpy, scipy, h5py; print('pywarpx', getattr(pywarpx, '__version__', '?'), 'gmsh', gmsh.__version__)"
ls "$ENV/bin" | grep -i warpx || true
