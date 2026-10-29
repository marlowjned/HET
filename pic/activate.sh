# source pic/activate.sh : puts the hetpic environment (WarpX, gmsh, numpy, ...) and GetDP on PATH
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$MAMBA_ROOT_PREFIX/envs/hetpic/bin:$HOME/tools/getdp/bin:$HOME/.local/bin:$PATH"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
