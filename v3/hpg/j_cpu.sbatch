#!/bin/bash
# j_cpu.sbatch -- pass J on HiPerGator for images whose script 2 ran here, one image per array task.
#
# The islets were reviewed on the workstation (gate 1) and exported there with
#   run_v3.sh --pass X --images <image>      -> results/v3/reviewed_islets/<tag>.geojson
# After that file (and results/v3/provenance/) is synced here, this job runs pass R then J:
#   R_import  replace every annotation except Tissue/pLN with the reviewed set (saved, verified by ID)
#   3, 8      signed distances, parenting
#   J_context, J_annotations, J_islets   the exports (J_context's Parquet post hook needs duckdb in $PY)
# run_v3.sh orders them (R_import before 3/8) and refuses to run pass J if R_import failed.
#
#   IMAGES_FILE=... sbatch --account=clive --qos=clive --array=0-5%1 j_cpu.sbatch     (hpg_submit.sh does this)
# --time must exceed the sum of the pass R,J timeouts in steps.tsv plus a margin (1020 + 30 min today):
# the job checks this first (hpg_check_budget); hpg_submit.sh computes --time.
#SBATCH --job-name=v3_j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180gb
#SBATCH --time=18:00:00
#SBATCH --output=/blue/maigan/smith6jt/IO60panc2nd_v3/logs/%x_%A_%a.log
set -euo pipefail
# shellcheck source=hpg_env.sh
source /home/smith6jt/IO60panc2nd/scripts/v3/hpg/hpg_env.sh

IMAGE=$(hpg_image_for_task)
[[ -n $IMAGE ]] || { hpg_log "no image for task ${SLURM_ARRAY_TASK_ID:-?}"; exit 1; }
hpg_check_budget --passes R,J --margin-min 30
tag=$(python3 -c 'import re, sys; print(re.sub(r"[^A-Za-z0-9._-]+", "_", sys.argv[1]))' "$IMAGE")
reviewed=$PROJ/results/v3/reviewed_islets/$tag.geojson
[[ -s $reviewed ]] || { hpg_log "missing $reviewed -- export it on the workstation (--pass X) and sync it"; exit 1; }
# run_v3 tracks step order, not input content: re-import when this file differs from the one imported
# (the current R_import record, as the merged provenance view sees it)
force=()
imported=$("$PY" "$V3DIR/v3ctl.py" --project "$PROJ/project.qpproj" show-record --image "$IMAGE" --step R_import \
    --field step_result.details.source_sha256 --require-pass 2> /dev/null || true)
if [[ $imported != "$(sha256sum "$reviewed" | cut -d' ' -f1)" ]]; then force=(--force R_import); fi
hpg_log "passes R,J for '$IMAGE' on $(hostname) with $reviewed ${force[*]:-(already imported)}"

rc=0
"$V3DIR/run_v3.sh" --pass R,J --images "$IMAGE" --lanes-cpu 1 --lanes-light 1 "${RUN_V3_HPG_OPTS[@]}" "${force[@]}" \
    --run-id "hpg-j-${SLURM_ARRAY_JOB_ID:-0}-${SLURM_ARRAY_TASK_ID:-0}" || rc=$?
hpg_log "run_v3.sh exit $rc"
exit "$rc"
