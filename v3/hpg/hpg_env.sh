# hpg_env.sh -- shared settings for the HiPerGator side of the v3 run. Source it; do not execute it.
#
# The point of the whole bundle: on HiPerGator the project lives at the SAME absolute path as on the
# workstation (/home/smith6jt/IO60panc2nd), with data/ and Images/ pointing into /blue. The .qpdata
# folders and results/v3/provenance then move both ways with no path rewriting, and run_v3.sh,
# steps.tsv and every script run unchanged.
# shellcheck shell=bash

export BLUE=${BLUE:-/blue/maigan/smith6jt/IO60panc2nd_v3}     # processing copy: data/ Images/ results/ work/
export PROJ=/home/smith6jt/IO60panc2nd                           # identical path to the workstation
export V3DIR=$PROJ/scripts/v3
export QUPATH_HOME=$HOME/QuPath                                  # portable QuPath 0.7 (Linux build) + v0.7/ user dir
export V3_QUPATH=$QUPATH_HOME/bin/QuPath
export DJL_CACHE_DIR=$HOME/.djl.ai                               # PyTorch 2.7.1-cu128 engine, copied from the workstation
export PY=${PY:-$BLUE/venv-v3/bin/python}                        # v3ctl.py needs only python3; post hooks need duckdb+pyarrow

# run_v3.sh working areas: locks node-local (flock is not trusted across nodes on shared filesystems),
# the rest on /blue so a failed job can be inspected.
export V3_LOCKDIR=${TMPDIR:-/tmp}/v3_locks
export V3_WORKROOT=$BLUE/work/driver
export V3_TMP=$BLUE/tmp/run_v3

# every run_v3.sh call on HiPerGator: no GUI can exist on a compute node, but other users' processes
# are visible to ps, so the busy check is skipped; the free-space floor applies to /blue, not '/'.
# (an array cannot be exported; every script that needs it sources this file)
RUN_V3_HPG_OPTS=(--allow-busy --root-path "$BLUE" --min-root-free-gb 200 --data-margin-gb 20)

hpg_log() { printf '%s [%s] %s\n' "$(date -Is)" "${SLURM_JOB_ID:-login}" "$*" >&2; }

# image name for this array task: line (SLURM_ARRAY_TASK_ID + 1) of IMAGES_FILE (one image name per line)
hpg_image_for_task() {
    local n=$((${SLURM_ARRAY_TASK_ID:-0} + 1))
    [[ -s ${IMAGES_FILE:-} ]] || { hpg_log "IMAGES_FILE not set or empty"; return 1; }
    sed -n "${n}p" "$IMAGES_FILE"
}
