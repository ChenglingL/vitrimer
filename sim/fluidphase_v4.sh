#!/usr/bin/env bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES=0 ## 0 is fpr the TiTan / 1 is for Quadro
PY=python  # or full path to your env: /home/you/miniconda3/envs/hoomd/bin/python
OUT=/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/fluidPhase
run_dir="${OUT}/"
mkdir -p "$run_dir"

$PY /home/cli428/vitrimer/sim/fluidphase_v4.py \
  --kT 50 \
  --rho 0.85 0.95 1.05 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 \
  --outdir "$run_dir" \
  --dt 0.001 \
  --seed 0 


wait

echo "All jobs completed!"

