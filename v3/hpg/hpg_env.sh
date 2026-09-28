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
# QuPath finds its extensions (InstanSeg, DJL) only through the 'userPath' Java preference; without it
# headless QuPath 0.7 loads NO extensions (tested 2026-09-28). The pipeline therefore reads its own
# preference tree, created by setup_hpg.sh, instead of the account's ~/.java (left untouched).
export V3_PREFS_ROOT=$QUPATH_HOME/pipeline_prefs
export V3_JAVA_OPTS="-Djava.util.prefs.userRoot=$V3_PREFS_ROOT"  # appended to every QuPath call by run_v3.sh
export PY=${PY:-$BLUE/venv-v3/bin/python}                        # v3ctl.py needs only python3; post hooks need duckdb+pyarrow

# Provenance site: every node of this cluster writes results/v3/provenance/<tag>/hpg/, so records
# synced from the workstation (<tag>/<workstation>/) are never overwritten; they count here only after
# `run_v3.sh --merge-provenance --images ...` adopted them (and the workstation merges ours likewise).
export V3_SITE=hpg

# run_v3.sh working areas. flock is node-local on a cluster filesystem, so image locks ALSO use mkdir
# on /blue (V3_SHARED_LOCKDIR): two jobs on different nodes can never process one image at once. A lock
# left by a killed job is broken once squeue no longer lists that job.
export V3_LOCKDIR=${TMPDIR:-/tmp}/v3_locks
export V3_SHARED_LOCKDIR=$BLUE/locks
export V3_WORKROOT=$BLUE/work/driver
export V3_TMP=$BLUE/tmp/run_v3

# every run_v3.sh call on HiPerGator: no GUI can exist on a compute node, but other users' processes
# are visible to ps, so the busy check is skipped; the free-space floor applies to /blue, not '/'.
# (an array cannot be exported; every script that needs it sources this file)
RUN_V3_HPG_OPTS=(--allow-busy --root-path "$BLUE" --min-root-free-gb 200 --data-margin-gb 20)

hpg_log() { printf '%s [%s] %s\n' "$(date -Is)" "${SLURM_JOB_ID:-login}" "$*" >&2; }

# the image lines of an images file: blank lines and '#' comments do not count (hpg_submit.sh writes
# normalized files, but the count and the lookup must agree on any file)
hpg_image_lines() { grep -v -E '^[[:space:]]*(#|$)' "$1" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//'; }

# image name for this array task: image line (SLURM_ARRAY_TASK_ID + 1) of IMAGES_FILE
hpg_image_for_task() {
    local n=$((${SLURM_ARRAY_TASK_ID:-0} + 1))
    [[ -s ${IMAGES_FILE:-} ]] || { hpg_log "IMAGES_FILE not set or empty"; return 1; }
    hpg_image_lines "$IMAGES_FILE" | sed -n "${n}p"
}

# SLURM durations: [D-]HH:MM:SS, MM:SS or SS -> seconds (fails on UNLIMITED/INVALID/NOT_SET)
hpg_duration_s() {
    local s=$1 d=0 h=0 m=0 x=0
    local -a p
    [[ $s =~ ^([0-9]+)-(.*)$ ]] && { d=${BASH_REMATCH[1]}; s=${BASH_REMATCH[2]}; }
    [[ $s =~ ^[0-9]+(:[0-9]+){0,2}$ ]] || return 1
    IFS=: read -r -a p <<< "$s"
    case ${#p[@]} in
        3) h=${p[0]} m=${p[1]} x=${p[2]} ;;
        2) m=${p[0]} x=${p[1]} ;;
        1) x=${p[0]} ;;
    esac
    echo $((10#$d * 86400 + 10#$h * 3600 + 10#$m * 60 + 10#$x))
}

# Refuse to start a job whose time limit cannot cover its steps: SLURM kills at the limit, and a step
# killed mid-save leaves data.qpdata half-written (.bkp behind). Needs: sum of the step timeouts in
# steps.tsv + a margin (v3ctl.py budget, the same numbers hpg_submit.sh uses for --time). Exports
# V3_DEADLINE_EPOCH, after which run_v3.sh starts no step and the footer starts no save.
#   hpg_check_budget --steps 0 | --passes R,J    [--margin-min N]
hpg_check_budget() {
    local need_min left_raw left_s
    need_min=$("$PY" "$V3DIR/v3ctl.py" budget "$@") || { hpg_log "budget: v3ctl.py budget $* failed"; exit 2; }
    left_raw=$(squeue -h -j "${SLURM_JOB_ID:?not inside a SLURM job}" -o %L 2> /dev/null | head -1 | tr -d '[:space:]')
    if ! left_s=$(hpg_duration_s "$left_raw"); then
        hpg_log "cannot read this job's time left ('$left_raw'); no deadline is enforced for the footer"
        return 0
    fi
    export V3_DEADLINE_EPOCH=$(($(date +%s) + left_s))
    if ((left_s < need_min * 60)); then
        hpg_log "time limit too short: ${left_s}s left, the steps need ${need_min} min (timeouts + margin); submit with hpg_submit.sh or a larger --time"
        exit 2
    fi
    hpg_log "time budget ok: ${left_s}s left >= ${need_min} min; job deadline $(date -d "@$V3_DEADLINE_EPOCH" -Is)"
}
