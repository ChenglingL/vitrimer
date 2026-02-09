#!/usr/bin/env bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES=0 ## 0 is fpr the TiTan / 1 is for Quadro
PY=python  # or full path to your env: /home/you/miniconda3/envs/hoomd/bin/python
OUT=/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/GPUfixed/fluidPhase
run_dir="${OUT}/"
fluid_dir=/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/snapshotSeed
mkdir -p "$run_dir"


declare -a KT_LIST=(10)
declare -a RHO_LIST=(1.0 1.1 1.2 1.3 1.4 1.45 1.5 1.6 1.7 2 2.2 2.5)
for kT in "${KT_LIST[@]}"; do
    for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/cli428/vitrimer/sim/fluidphase_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --outdir "$run_dir" \
      --fluid_dir "$fluid_dir" \
      --dt 0.001 \
      --seed 0 
    done
done
wait

echo "All jobs completed!"

