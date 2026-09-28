#!/bin/bash
# c0_gpu.sbatch -- pass C step 0 (InstanSeg) on HiPerGator: ONE image per array task, step 0 ONLY.
#
# HiPerGator kills a job whose GPU is idle for an hour, so this task does nothing long without the
# GPU: check that DJL's PyTorch engine really computes on CUDA (seconds; stop otherwise), optionally
# stage the qptiff to node-local $TMPDIR (minutes), run step 0 through run_v3.sh (the same footer,
# provenance and checks as on the workstation), exit. nvidia-smi dmon is logged alongside.
#
# Submit with hpg_submit.sh (it picks account/QoS, the GPU type and the array throttle); by hand:
#   IMAGES_FILE=... GPU_TYPE=b200 sbatch --account=clive --qos=clive --array=0-2%1 \
#       --gres=gpu:b200:1 --partition=<gpu partition> c0_gpu.sbatch
# Environment: IMAGES_FILE (required; one exact image name per line), GPU_TYPE (recorded; pinning is done
# by --gres at submission), STAGE_TO_TMPDIR=1 to read the image from $TMPDIR.
# --time must exceed step 0's timeout in steps.tsv (480 min) plus staging and the save margin: the job
# checks this first (hpg_check_budget) and hpg_submit.sh computes --time from steps.tsv.
#SBATCH --job-name=v3_c0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=09:00:00
#SBATCH --output=/blue/maigan/smith6jt/IO60panc2nd_v3/logs/%x_%A_%a.log
# TODO(HiPerGator): confirm partition names per GPU type (`sinfo -o "%P %G"`), e.g. B200 may be its own
# partition; and size step 0's timeout (hence --time) from the pilot.
set -euo pipefail
# shellcheck source=hpg_env.sh
source /home/smith6jt/IO60panc2nd/scripts/v3/hpg/hpg_env.sh

IMAGE=$(hpg_image_for_task)
[[ -n $IMAGE ]] || { hpg_log "no image for task ${SLURM_ARRAY_TASK_ID:-?}"; exit 1; }
hpg_log "step 0 for '$IMAGE' on $(hostname) (GPU_TYPE=${GPU_TYPE:-unset}, CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset})"
hpg_check_budget --steps 0 --margin-min 45        # exits before touching the GPU if --time is too short
LOGDIR=$BLUE/logs/c0_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}
mkdir -p "$LOGDIR"

# ---- GPU utilization log, stopped on exit ----
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv > "$LOGDIR/gpu.csv"
nvidia-smi dmon -s um -o DT > "$LOGDIR/nvidia_dmon.log" 2>&1 &
DMON=$!

# ---- optional staging: repoint this image's symlink at a node-local copy; always restore ----
LINK="" TARGET=""
restore() {
    if [[ -n $LINK && -n $TARGET ]]; then ln -sfn "$TARGET" "$LINK" && hpg_log "restored $LINK -> $TARGET"; fi
    kill "$DMON" 2> /dev/null || true
}
trap restore EXIT
trap 'exit 143' TERM          # SLURM time limit / scancel: run the EXIT trap, then stop

# ---- fast fail: DJL PyTorch must compute on CUDA before anything else happens ----
"$V3_QUPATH" script -a "$LOGDIR/djl_cuda_check.json" "$V3DIR/hpg/djl_cuda_check.groovy" > "$LOGDIR/djl_cuda_check.log" 2>&1 || true
if ! python3 -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get("ok") else 1)' "$LOGDIR/djl_cuda_check.json" 2> /dev/null; then
    hpg_log "DJL/CUDA check FAILED -- stopping before the GPU can sit idle; see $LOGDIR/djl_cuda_check.log"
    exit 2
fi
hpg_log "DJL/CUDA ok: $(tr -d '\n ' < "$LOGDIR/djl_cuda_check.json" | cut -c1-200)"

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

# ---- step 0 through the driver: provenance/<tag>/hpg/0.json (+ 0.finish.json) is the success signal ----
rc=0
"$V3DIR/run_v3.sh" --steps 0 --images "$IMAGE" --lanes-gpu 1 "${RUN_V3_HPG_OPTS[@]}" \
    --run-id "hpg-c0-${SLURM_ARRAY_JOB_ID:-0}-${SLURM_ARRAY_TASK_ID:-0}" || rc=$?
hpg_log "run_v3.sh exit $rc (0 = step 0 passed and saved; see results/v3/provenance/<tag>/hpg/0.json)"
exit "$rc"
