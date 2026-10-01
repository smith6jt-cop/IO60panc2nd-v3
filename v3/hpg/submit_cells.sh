#!/usr/bin/env bash
# submit_cells.sh -- submit pass C for ONE image on HiPerGator: a GPU job (preflight + step 0) and a
# CPU job (step 2 + exports) chained afterok. Run it on a login node (ssh hpg 'bash -l -s').
#
# The maintainer's rules (2026-10-01):
#   - investment QOS only: --qos is the account name (clive or maigan); burst QOS is never used;
#   - GPU type in the order rtx_pro_6000, b200, l4: the first whose typed, single-partition request
#     the scheduler says would start now (sbatch --test-only) is taken; if none would start now the
#     job is queued on the first type and can be moved later with move_gpu_job.sh;
#   - no self-imposed cap on either investment (a job over the QOS limit is rejected: DenyOnLimit).
#
# Usage: submit_cells.sh --account clive|maigan --image "<exact image name>" [--out DIR] [--stage]
#                        [--cpu-heap 160] [--cpu-cpus 20] [--start-within-min 5] [--dry-run]
# --cpu-heap is step 2's JVM heap in GB (about 40 + 45 per million cells); the CPU job asks for
# heap + 24 GB, so several fit one investment (clive 312 GB, maigan 625 GB).
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=hpg_env.sh
source "$HERE/hpg_env.sh"

ACCOUNT="" IMAGE="" OUT="$PROJ/results/v3/export" STAGE=0 CPU_HEAP=160 CPU_CPUS=20 WITHIN_MIN=5 DRY=0
while (($#)); do
    case $1 in
        --account) ACCOUNT=$2; shift ;;
        --image) IMAGE=$2; shift ;;
        --out) OUT=$2; shift ;;
        --stage) STAGE=1 ;;
        --cpu-heap) CPU_HEAP=$2; shift ;;
        --cpu-cpus) CPU_CPUS=$2; shift ;;
        --start-within-min) WITHIN_MIN=$2; shift ;;
        --dry-run) DRY=1 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
[[ $ACCOUNT == clive || $ACCOUNT == maigan ]] || { echo "--account must be clive or maigan" >&2; exit 2; }
[[ -n $IMAGE ]] || { echo "--image is required" >&2; exit 2; }
QOS=$ACCOUNT                                   # the investment QOS has the account's name; never <account>-b
CPU_MEM=$((CPU_HEAP + 24))gb

# GPU types in the maintainer's order: partition:gres type
GPU_ORDER=(hpg-rtx6000:rtx_pro_6000 hpg-b200:b200 hpg-turin:l4)
GPU_COMMON=(--account="$ACCOUNT" --qos="$QOS" --cpus-per-task=16 --mem=128gb --time=02:15:00)

# seconds until the scheduler expects a typed single-partition GPU request to start; empty if refused
start_delay_s() {
    local part=$1 typ=$2 line when
    line=$(sbatch --test-only "${GPU_COMMON[@]}" --partition="$part" --gres="gpu:$typ:1" --wrap=true 2>&1 | tail -1)
    when=$(sed -n 's/.*to start at \([0-9T:-]*\).*/\1/p' <<< "$line")
    [[ -n $when ]] || { echo ""; return 0; }
    echo $(($(date -d "$when" +%s) - $(date +%s)))
}

# When this account's own GPU jobs already fill its QOS GPU limit, a new job waits for THEM, not for
# the partition, and --test-only then reports our own queue. Sending it to a lesser GPU type would
# only break the preference order, so it is queued on the first choice (found 2026-10-01: four clive
# jobs went to b200 / l4 behind one running rtx_pro_6000 job).
gpu_limit=$(sacctmgr -nP show qos "$QOS" format=grptres | tr ',' '\n' | sed -n 's#^gres/gpu=##p' | head -1)
own_gpu=$(squeue -h -A "$ACCOUNT" -u "$USER" -n v3_cells_gpu -t PENDING,RUNNING -o %i | wc -l)
pick="" pick_delay=""
if [[ -n $gpu_limit ]] && ((own_gpu >= gpu_limit)); then
    pick=${GPU_ORDER[0]}
    hpg_log "$ACCOUNT already has $own_gpu GPU job(s) for its limit of $gpu_limit: queueing on ${pick%%:*} (first choice)"
fi
for entry in "${GPU_ORDER[@]}"; do
    [[ -n $pick ]] && break
    part=${entry%%:*} typ=${entry##*:}
    d=$(start_delay_s "$part" "$typ")
    hpg_log "test-only $ACCOUNT $part gpu:$typ -> ${d:-refused}${d:+ s}"
    if [[ -n $d ]] && ((d <= WITHIN_MIN * 60)); then pick=$entry pick_delay=$d; break; fi
done
if [[ -z $pick ]]; then
    pick=${GPU_ORDER[0]}
    hpg_log "no GPU type would start within $WITHIN_MIN min; queueing on ${pick%%:*} (first choice)"
fi
PART=${pick%%:*} TYP=${pick##*:}

export_vars="ALL,IMAGE=$IMAGE,V3_OUT=$OUT,GPU_TYPE=$TYP,STAGE_TO_TMPDIR=$STAGE"
if ((DRY)); then
    echo "DRY RUN: $ACCOUNT '$IMAGE' -> $PART gpu:$TYP (expected start in ${pick_delay:-?} s); cpu mem $CPU_MEM; out $OUT"
    exit 0
fi
mkdir -p "$BLUE/logs" "$BLUE/work/submit"
gpu_job=$(sbatch --parsable "${GPU_COMMON[@]}" --partition="$PART" --gres="gpu:$TYP:1" --export="$export_vars" "$HERE/cells_gpu.sbatch")
cpu_job=$(sbatch --parsable --account="$ACCOUNT" --qos="$QOS" --cpus-per-task="$CPU_CPUS" --mem="$CPU_MEM" --time=04:45:00 \
    --dependency="afterok:$gpu_job" --export="ALL,IMAGE=$IMAGE,V3_OUT=$OUT,STEP2_HEAP_GB=$CPU_HEAP" "$HERE/cells_cpu.sbatch")
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -Is)" "$ACCOUNT" "$PART" "$TYP" "$gpu_job" "$cpu_job" "$IMAGE" >> "$BLUE/work/submit/cells_jobs.tsv"
echo "submitted: account=$ACCOUNT partition=$PART gpu=$TYP gpu_job=$gpu_job cpu_job=$cpu_job image='$IMAGE'"
