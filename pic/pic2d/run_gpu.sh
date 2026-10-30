#!/bin/bash
# het_breathing.py on the GPU build of WarpX (micromamba env warpx-gpu: CUDA, 2-D and RZ, double fields, single
# particles); arguments are passed through. Long runs are started detached from Windows so they outlive the shell:
#   Start-Process wsl.exe -WindowStyle Hidden -ArgumentList '-e','bash','-lc','bash ~/Projects/HET/pic/pic2d/run_gpu.sh --tag X ... > ~/Projects/HET/pic/pic2d/runs_X.log 2>&1'
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$MAMBA_ROOT_PREFIX/envs/warpx-gpu/bin:$PATH"
export CUDA_PATH="$MAMBA_ROOT_PREFIX/envs/warpx-gpu/targets/x86_64-linux"   # cupy JIT headers
cd "$(dirname "$0")" || exit 1
exec python het_breathing.py "$@"
