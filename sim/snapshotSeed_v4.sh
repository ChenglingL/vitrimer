#!/usr/bin/env bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES=0 ## 0 is fpr the TiTan / 1 is for Quadro
PY=python  # or full path to your env: /home/you/miniconda3/envs/hoomd/bin/python
OUT=/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/snapshotSeed
run_dir="${OUT}/"
mkdir -p "$run_dir"

$PY /home/cli428/vitrimer/sim/snapshotSeed_v4.py \
  --kT 0.2 \
  --rho 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 \
  --outdir "$run_dir" \
  --seed 0 


wait

echo "All jobs completed!"

