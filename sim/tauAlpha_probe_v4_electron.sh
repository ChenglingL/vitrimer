#!/usr/bin/env bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES=0 ## 0 is fpr the TiTan / 1 is for Quadro
PY=python  # or full path to your env: /home/you/miniconda3/envs/hoomd/bin/python
OUT=/home/chengling/Research/Project/vitrimer/data/test/vitrimerPaper/NVT/V4Test
run_dir="${OUT}/"
mkdir -p "$run_dir"

declare -a KT_LIST=(1.0 0.5 0.25 0.1)
declare -a RHO_LIST=(1.7 2.0)

for kT in "${KT_LIST[@]}"; do
  for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/chengling/Research/Project/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
      --per_decade 10 \
      --duration_after_wait 1000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --dt 0.001 \
      --seed 0 &
  done
done

declare -a KT_LIST=(1.0 0.1 0.01 0.001)
declare -a RHO_LIST=(0.5)

for kT in "${KT_LIST[@]}"; do
  for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/chengling/Research/Project/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
      --per_decade 10 \
      --duration_after_wait 1000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --dt 0.001 \
      --seed 0 &
  done
done

declare -a KT_LIST=(1.0 0.3 0.1 0.03)
declare -a RHO_LIST=(1.1 1.5)

for kT in "${KT_LIST[@]}"; do
  for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/chengling/Research/Project/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
      --per_decade 10 \
      --duration_after_wait 1000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --dt 0.001 \
      --seed 0 &
  done
done
wait

echo "All jobs completed!"

