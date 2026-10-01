#!/bin/bash
# cells_gpu.sbatch -- pass C on HiPerGator, GPU half: preflight + step 0 (InstanSeg) for ONE image.
#
# One image per job, so a GPU is held only while it segments (the site kills a job whose GPU sits
# idle for an hour). Account, QOS, partition, GRES and --time come from submit_cells.sh on the
# command line (investment QOS only; GPU type in the order rtx_pro_6000, b200, l4).
#
# Environment: IMAGE (exact image name, required), GPU_TYPE (recorded), V3_OUT (output root; default
# $PROJ/results/v3/export), STAGE_TO_TMPDIR=1 to read the image from node-local disk.
#SBATCH --job-name=v3_cells_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=02:00:00
#SBATCH --output=/blue/maigan/smith6jt/IO60panc2nd_v3/logs/%x_%j.log
set -euo pipefail
# shellcheck source=hpg_env.sh
source /home/smith6jt/IO60panc2nd/scripts/v3/hpg/hpg_env.sh

IMAGE=${IMAGE:?IMAGE (exact image name) is required}
OUT=${V3_OUT:-$PROJ/results/v3/export}
STEPS=preflight,0
hpg_log "steps $STEPS for '$IMAGE' on $(hostname) (account ${SLURM_JOB_ACCOUNT:-?}, qos ${SLURM_JOB_QOS:-?}, partition ${SLURM_JOB_PARTITION:-?}, GPU_TYPE=${GPU_TYPE:-unset})"
hpg_check_budget --steps "$STEPS" --margin-min 15     # exits before touching the GPU if --time is too short
LOGDIR=$BLUE/logs/cells_gpu_${SLURM_JOB_ID:-0}
mkdir -p "$LOGDIR"

# ---- GPU model and utilisation log, stopped on exit ----
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv > "$LOGDIR/gpu.csv"
nvidia-smi dmon -s um -o DT > "$LOGDIR/nvidia_dmon.log" 2>&1 &
DMON=$!

LINK="" TARGET=""
restore() {
    if [[ -n $LINK && -n $TARGET ]]; then ln -sfn "$TARGET" "$LINK" && hpg_log "restored $LINK -> $TARGET"; fi
    kill "$DMON" 2> /dev/null || true
}
trap restore EXIT
trap 'exit 143' TERM

# ---- fast fail: DJL PyTorch must compute on CUDA before anything else happens ----
JAVA_TOOL_OPTIONS="-Xmx8g $V3_JAVA_OPTS" \
    "$V3_QUPATH" script -a "$LOGDIR/djl_cuda_check.json" "$V3DIR/hpg/djl_cuda_check.groovy" > "$LOGDIR/djl_cuda_check.log" 2>&1 || true
if ! python3 -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get("ok") else 1)' "$LOGDIR/djl_cuda_check.json" 2> /dev/null; then
    hpg_log "DJL/CUDA check FAILED -- stopping before the GPU can sit idle; see $LOGDIR/djl_cuda_check.log"
    exit 2
fi
hpg_log "DJL/CUDA ok on $(sed -n 2p "$LOGDIR/gpu.csv")"

if [[ ${STAGE_TO_TMPDIR:-0} == 1 ]]; then
    uri=$(python3 -c 'import json, sys
e = next(e for e in json.load(open(sys.argv[1]))["images"] if e["imageName"] == sys.argv[2])
print(e["serverBuilder"]["uri"])' "$PROJ/project.qpproj" "$IMAGE")
    LINK=${uri#file:}
    TARGET=$(readlink "$LINK")
    staged=$TMPDIR/${LINK##*/}
    t0=$SECONDS
    cp "$TARGET" "$staged"
    ln -sfn "$staged" "$LINK"
    hpg_log "staged ${LINK##*/} to $TMPDIR in $((SECONDS - t0)) s"
fi

rc=0
/usr/bin/time -v -o "$LOGDIR/time.txt" \
    "$V3DIR/run_v3.sh" --steps "$STEPS" --images "$IMAGE" --out "$OUT" --lanes-gpu 1 --lanes-light 1 "${RUN_V3_HPG_OPTS[@]}" \
    --run-id "hpg-cells-gpu-${SLURM_JOB_ID:-0}" || rc=$?
hpg_log "run_v3.sh exit $rc; wall $(grep -m1 'Elapsed' "$LOGDIR/time.txt" | awk '{print $NF}'), max RSS $(grep -m1 'Maximum resident' "$LOGDIR/time.txt" | awk '{printf "%.1f GB", $NF/1e6}')"
exit "$rc"
