#!/usr/bin/env bash
set -euo pipefail

PY=python  # or full path to your env: /home/you/miniconda3/envs/hoomd/bin/python
OUT=/home/chengling/Research/Project/vitrimer/data/test/vitrimerPaper/NVT/V5RevCrossDisable

mkdir -p "$OUT"

# Parameter sweeps
# declare -a KT_LIST=(0.0025 0.0022 0.001 0.00035)
# declare -a KT_LIST=(0.01 0.05 0.0005 0.0001 0.00005)
# declare -a KT_LIST=(0.1 0.05 0.025 0.01 0.005 0.0025 0.001 0.0007 0.0003 0.0002)
# declare -a KT_LIST=(1.0 0.1 0.01 0.001 0.0001 0.00001 0.000001)
# declare -a RHO_LIST=(0.923333)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/chengling/Research/Project/vitrimer/sim/test.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 20 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --dt 0.001 \
#       --seed 0 &
#   done
# done
declare -a KT_LIST=(1.0 0.1 0.03)
declare -a RHO_LIST=(0.92333)

for kT in "${KT_LIST[@]}"; do
  for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/chengling/Research/Project/vitrimer/sim/test.py \
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

# declare -a KT_LIST=(0.15 0.09 0.08 0.07 0.065 0.06)
# declare -a RHO_LIST=(2.5155)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     rho_formatted=$(printf "%.6f" "$rho")
#     run_dir="${OUT}/rho${rho_formatted}"
#     mkdir -p "$run_dir"

#     $PY /home/chengling/Research/Project/vitrimer/sim/test.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 20 \
#       --outdir "$run_dir" \
#       --dt 0.001 \
#       --seed 0 &
#   done
# done

wait

echo "All jobs completed!"

