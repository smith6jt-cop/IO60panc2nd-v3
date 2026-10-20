#!/bin/bash
# c2_cpu.sbatch -- pass C step 2 (all measurements) on HiPerGator CPU nodes, one image per array task.
#
# Step 2 is deterministic CPU arithmetic, so HiPerGator and workstation lanes may share the queue
# without changing a value. Needs step 0 done for the image (its provenance synced here), exactly as
# run_v3.sh demands locally. STEPS may add the pass C exports (C_cells needs duckdb in $PY).
#
#   IMAGES_FILE=... sbatch --account=clive --qos=clive --array=0-5%1 c2_cpu.sbatch     (hpg_submit.sh does this)
# Environment: IMAGES_FILE (required), STEPS (default 2; e.g. 2,C_cells,C_geojson).
# --time must exceed the sum of the STEPS' timeouts in steps.tsv plus a margin (step 2 alone: 1440 min):
# the job checks this first, and hpg_submit.sh computes --time for the STEPS it submits.
#SBATCH --job-name=v3_c2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=220gb
#SBATCH --time=25:00:00
#SBATCH --output=/blue/maigan/smith6jt/IO60panc2nd_v3/logs/%x_%A_%a.log
# TODO(HiPerGator): size --mem/--time from the pilot; steps.tsv caps step 2 at -Xmx160g, and the JVM's
# resident size can approach that cap, so --mem must stay well above it.
set -euo pipefail
# shellcheck source=hpg_env.sh
source /home/smith6jt/IO60panc2nd/scripts/v3/hpg/hpg_env.sh

IMAGE=$(hpg_image_for_task)
[[ -n $IMAGE ]] || { hpg_log "no image for task ${SLURM_ARRAY_TASK_ID:-?}"; exit 1; }
STEPS=${STEPS:-2}
hpg_check_budget --steps "$STEPS" --margin-min 30
threads=$((${SLURM_CPUS_PER_TASK:-32} - 4))     # leave a few cores for GC and I/O
hpg_log "steps $STEPS for '$IMAGE' on $(hostname), step 2 N_THREADS=$threads"

rc=0
"$V3DIR/run_v3.sh" --steps "$STEPS" --images "$IMAGE" --lanes-cpu 1 "${RUN_V3_HPG_OPTS[@]}" \
    --override "2:N_THREADS=$threads" \
    --run-id "hpg-c2-${SLURM_ARRAY_JOB_ID:-0}-${SLURM_ARRAY_TASK_ID:-0}" || rc=$?
hpg_log "run_v3.sh exit $rc"
exit "$rc"
