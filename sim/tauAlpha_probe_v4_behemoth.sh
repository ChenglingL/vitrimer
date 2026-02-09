#!/usr/bin/env bash
set -euo pipefail
export CUDA_VISIBLE_DEVICES=0 ## 0 is fpr the TiTan / 1 is for Quadro
PY=python  # or full path to your env: /home/you/miniconda3/envs/hoomd/bin/python
OUT=/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/GPUfixed
run_dir="${OUT}/"
fluid_dir=/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/snapshotSeed
mkdir -p "$run_dir"


declare -a KT_LIST=(0.6 0.5 0.45)
declare -a RHO_LIST=(1.45)
for kT in "${KT_LIST[@]}"; do
    for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
      --per_decade 10 \
      --duration_after_wait 5000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --fluid_dir "$fluid_dir" \
      --dt 0.001 \
      --seed 0 
    done
done

declare -a KT_LIST=(10 3 1)
declare -a RHO_LIST=(1.45)
for kT in "${KT_LIST[@]}"; do
    for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
      --per_decade 10 \
      --duration_after_wait 5000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --fluid_dir "$fluid_dir" \
      --dt 0.001 \
      --seed 0 
    done
done

declare -a KT_LIST=(1.5 1.2 1.0 0.8)
declare -a RHO_LIST=(1.6)
for kT in "${KT_LIST[@]}"; do
    for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
      --per_decade 10 \
      --duration_after_wait 5000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --fluid_dir "$fluid_dir" \
      --dt 0.001 \
      --seed 0 
    done
done

declare -a KT_LIST=(10 5 2)
declare -a RHO_LIST=(1.6)
for kT in "${KT_LIST[@]}"; do
    for rho in "${RHO_LIST[@]}"; do
    run_dir="${OUT}/rho${rho}"
    mkdir -p "$run_dir"

    $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
      --kT "$kT" \
      --rho "$rho" \
      --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
      --per_decade 10 \
      --duration_after_wait 5000000 \
      --mttk_tau 1.0 \
      --outdir "$run_dir" \
      --fluid_dir "$fluid_dir" \
      --dt 0.001 \
      --seed 0 
    done
done

# declare -a KT_LIST=(0.0001)
# declare -a RHO_LIST=(1.0)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done
# declare -a KT_LIST=(0.0005)
# declare -a RHO_LIST=(1.0)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.045)
# declare -a RHO_LIST=(1.2)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.011)
# declare -a RHO_LIST=(0.5)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 5000000 10000000 15000000 20000000 25000000 30000000 35000000 40000000 45000000 50000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.45 0.4)
# declare -a RHO_LIST=(1.4)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.45)
# declare -a RHO_LIST=(1.2)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.001 0.0008)
# declare -a RHO_LIST=(1.0)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 3 1 0.3 0.1 0.03 0.01 0.005)
# declare -a RHO_LIST=(1.0)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.003)
# declare -a RHO_LIST=(1.0)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.008 0.005)
# declare -a RHO_LIST=(1.1)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.18 0.15)
# declare -a RHO_LIST=(1.3)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 3 1 0.3 0.1 0.05)
# declare -a RHO_LIST=(1.2)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 3 1 0.5 0.3)
# declare -a RHO_LIST=(1.4)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done



# declare -a KT_LIST=(0.0005 0.0001)
# declare -a RHO_LIST=(0.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(0.7)
# declare -a RHO_LIST=(1.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(1.6)
# declare -a RHO_LIST=(1.7)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(3.5)
# declare -a RHO_LIST=(2.0)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(4.5)
# declare -a RHO_LIST=(2.2)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(6.5)
# declare -a RHO_LIST=(2.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 500000 1000000 1500000 2000000 2500000 3000000 3500000 4000000 4500000 5000000\
#       --per_decade 10 \
#       --duration_after_wait 5000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done


# declare -a KT_LIST=(0.25 0.2)
# declare -a RHO_LIST=(1.3)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done
# declare -a KT_LIST=(0.001)
# declare -a RHO_LIST=(1.1)
# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 3 1 0.3 0.1)
# declare -a RHO_LIST=(1.3)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 3 1 0.1 0.01 0.001)
# declare -a RHO_LIST=(0.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 3 1 0.3 0.1 0.03 0.1)
# declare -a RHO_LIST=(1.1)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 5 2 1 0.9 0.8 0.7)
# declare -a RHO_LIST=(1.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 5 2 1.8 1.5)
# declare -a RHO_LIST=(1.7)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 6 5 4 3.8 3.5)
# declare -a RHO_LIST=(2.0)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 8 7 6.8 6.5)
# declare -a RHO_LIST=(2.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(10 7 6 5)
# declare -a RHO_LIST=(2.2)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(15 20 30)
# declare -a RHO_LIST=(1.7 2.0 2.2 2.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(1.0 0.5 0.25 0.1)
# declare -a RHO_LIST=(1.7 2.0)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --fluid_dir "$fluid_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(1.0 0.1 0.01 0.001)
# declare -a RHO_LIST=(0.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done

# declare -a KT_LIST=(1.0 0.3 0.1 0.03)
# declare -a RHO_LIST=(1.1 1.5)

# for kT in "${KT_LIST[@]}"; do
#   for rho in "${RHO_LIST[@]}"; do
#     run_dir="${OUT}/rho${rho}"
#     mkdir -p "$run_dir"

#     $PY /home/cli428/vitrimer/sim/tauAlpha_probe_v4.py \
#       --kT "$kT" \
#       --rho "$rho" \
#       --waits 0 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000\
#       --per_decade 10 \
#       --duration_after_wait 1000000 \
#       --mttk_tau 1.0 \
#       --outdir "$run_dir" \
#       --dt 0.001 \
#       --seed 0 
#   done
# done
wait

echo "All jobs completed!"

