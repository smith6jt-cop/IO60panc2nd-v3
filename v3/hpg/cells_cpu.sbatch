#!/bin/bash
# cells_cpu.sbatch -- pass C on HiPerGator, CPU half: step 2 (all measurements) + the two exports for
# ONE image. Submitted by submit_cells.sh with --dependency=afterok:<the image's GPU job>, so it
# starts only when step 0 passed and saved. No GPU is requested.
#
# Environment: IMAGE (exact image name, required), V3_OUT (output root; default $PROJ/results/v3/export).
#SBATCH --job-name=v3_cells_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=184gb
#SBATCH --time=04:30:00
#SBATCH --output=/blue/maigan/smith6jt/IO60panc2nd_v3/logs/%x_%j.log
set -euo pipefail
# shellcheck source=hpg_env.sh
source /home/smith6jt/IO60panc2nd/scripts/v3/hpg/hpg_env.sh

IMAGE=${IMAGE:?IMAGE (exact image name) is required}
OUT=${V3_OUT:-$PROJ/results/v3/export}
STEPS=2,C_cells,C_geojson
hpg_check_budget --steps "$STEPS" --margin-min 15
threads=$((${SLURM_CPUS_PER_TASK:-24} - 4))     # leave a few cores for GC and I/O
LOGDIR=$BLUE/logs/cells_cpu_${SLURM_JOB_ID:-0}
mkdir -p "$LOGDIR"
hpg_log "steps $STEPS for '$IMAGE' on $(hostname) (account ${SLURM_JOB_ACCOUNT:-?}, qos ${SLURM_JOB_QOS:-?}), step 2 N_THREADS=$threads"

rc=0
/usr/bin/time -v -o "$LOGDIR/time.txt" \
    "$V3DIR/run_v3.sh" --steps "$STEPS" --images "$IMAGE" --out "$OUT" --lanes-cpu 1 "${RUN_V3_HPG_OPTS[@]}" \
    --override "2:N_THREADS=$threads" \
    --run-id "hpg-cells-cpu-${SLURM_JOB_ID:-0}" || rc=$?
hpg_log "run_v3.sh exit $rc; wall $(grep -m1 'Elapsed' "$LOGDIR/time.txt" | awk '{print $NF}'), max RSS $(grep -m1 'Maximum resident' "$LOGDIR/time.txt" | awk '{printf "%.1f GB", $NF/1e6}')"
exit "$rc"
