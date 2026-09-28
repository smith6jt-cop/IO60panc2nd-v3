#!/usr/bin/env bash
# run_v3.sh -- resumable scheduler for the v3 QuPath passes over the project images.
#
#   pass I : preflight -> 1 (islets) -> 4 (bands)
#   pass C : preflight -> 0 (InstanSeg over Tissue + pLN, GPU) -> 2 (measurements) -> C_cells (Parquet) + C_geojson
#   pass J : 3 (signed distance) -> 8 (parenting) -> J_context (Parquet), J_annotations, J_islets
#   pass X / R : export the reviewed islets locally / import them on HiPerGator before pass J there
# The steps, their scripts, heaps, lanes and checks live in steps.tsv; v3ctl.py makes every decision
# that needs JSON; this file runs processes, holds locks and fills lanes.
#
# Safety mechanics (see README.md):
#   - each step runs in its own headless QuPath call against a PRIVATE copy of project.qpproj
#     (data/ classifiers/ scripts/ resources/ symlinked), never the target project file itself;
#   - the step script is copied, its overrides applied, and a generated footer appended; the footer
#     checks the postcondition, re-checks that data.qpdata is still the file the step started from,
#     saves (footer-saves policy) and writes provenance/<tag>/<site>/<step>.json. No JSON = the step
#     failed, whatever QuPath's exit code says;
#   - a per-image flock means two steps never touch one image at once (plus a mkdir lock in
#     V3_SHARED_LOCKDIR across cluster nodes); lane slots (gpu/cpu/light) and a memory budget cap
#     concurrency; QuPath inherits the lock fds, so an orphaned QuPath still holds its image;
#   - free-space guard before every step; refuses to start while QuPath (GUI or headless) or a
#     Senior rebuild runs (lib.sh assert_project_idle), and stops starting steps if a GUI appears;
#   - a step whose provenance passed and was validated, whose outputs exist, whose script is still the
#     one the plan resolves, and whose upstream steps did not re-run (and no accepted external edit
#     came) after it is skipped; if data.qpdata changed outside run_v3 the image is blocked until
#     --accept-external records the change.
#   - Ctrl-C / TERM: running steps finish (QuPath runs in timeout's own process group), none start.
#
# Usage:
#   run_v3.sh --pass I|C|J|X|R[,..] | --steps S1,S2   [options]
#     --images LIST        comma list of EXACT image names, tags, or id:<entry id> (default: all);
#                          every token must match exactly one image
#     --images-file FILE   one image (name, tag or id:N) per line; blank lines and # comments ignored
#     --project FILE       target project (default: $PROJ/project.qpproj)
#     --out DIR            outputs + provenance (default: $V3 for the master, <project dir>/results/v3 otherwise)
#     --dry-run            print the plan (and busy/disk status); change nothing
#     --force S1,S2        re-run these steps even if done (their downstream steps re-run too); each must
#                          be among the selected steps
#     --rerun-from S       the same as --force S for a run; with --accept-external it is RECORDED: the
#                          edit invalidates S and everything downstream of it until they re-run
#     --lanes-gpu N  --lanes-cpu N  --lanes-light N     slot counts (default 1 each)
#     --mem-budget GB      sum of -Xmx (and DuckDB reservations) of concurrently running steps (default 300)
#     --duckdb-mem-gb N    DuckDB memory limit of the Parquet post hooks (default 48; reserved as N+8)
#     --save-margin-min N  no save starts later than N minutes before a step's deadline (default 10)
#     --script-set v3|orig|auto   v3 (default) requires the v3 copy of every step that has one
#     --script STEP=PATH   run PATH instead (testing);  --heap STEP=GB;  --override STEP:NAME=VALUE
#     --allow-busy         run although QuPath/a rebuild is running (forces --lanes-cpu 1)
#     --allow-drift        do not block an image whose step changed its server metadata
#     --accept-external --note TEXT [--rerun-from S]   record the current data.qpdata of --images as accepted
#     --merge-provenance [--from SITE] [--check]      adopt the records another site wrote for --images
#     --min-root-free-gb N (300)  --data-margin-gb N (50)  --root-path P (/)
#     --run-id ID
# Environment: V3_SITE (provenance site; default the short host name), V3_LOCKDIR, V3_SHARED_LOCKDIR
# (mkdir image locks on a shared filesystem), V3_DEADLINE_EPOCH (a job's end: no step saves after it),
# V3_WORKROOT, V3_TMP, V3_QUPATH, PY.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=lib.sh
source "$HERE/lib.sh"            # PROJ V3 SCRATCH QUPATH assert_project_idle

PY=${PY:-/scratch/io60/venv-v3/bin/python}
CTL=$HERE/v3ctl.py
QUPATH=${V3_QUPATH:-$QUPATH}
LOCKDIR=${V3_LOCKDIR:-$SCRATCH/locks}
SHARED_LOCKDIR=${V3_SHARED_LOCKDIR:-}
WORKROOT=${V3_WORKROOT:-$SCRATCH/work/driver}
TMPROOT=${V3_TMP:-$SCRATCH/tmp/run_v3}

PROJECT=$PROJ/project.qpproj OUT="" PASSES="" STEPS="" IMAGES="" IMAGES_FILE="" FORCE="" RERUN_FROM="" DRY=0
LANES_GPU=1 LANES_CPU=1 LANES_LIGHT=1 MEM_BUDGET_GB=300 SCRIPT_SET=v3 DUCKDB_MEM_GB=48 SAVE_MARGIN_MIN=10
ALLOW_BUSY=${SKIP_BUSY_CHECK:-0} ALLOW_DRIFT=0 ACCEPT=0 MERGE=0 MERGE_FROM="" CHECK_ONLY=0 NOTE="" RUN_ID=""
MIN_ROOT_FREE_GB=300 DATA_MARGIN_GB=50 ROOT_PATH=/
CTL_PASS=()                      # --script/--heap/--override, handed to every v3ctl call

usage() { sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit "${1:-0}"; }
while (($#)); do
    case $1 in
        --pass) PASSES=$2; shift ;;
        --steps) STEPS=$2; shift ;;
        --images) IMAGES=$2; shift ;;
        --images-file) IMAGES_FILE=$(realpath -e "$2"); shift ;;
        --project) PROJECT=$2; shift ;;
        --out) OUT=$2; shift ;;
        --dry-run) DRY=1 ;;
        --force) FORCE=$2; shift ;;
        --rerun-from) RERUN_FROM=$2; shift ;;
        --lanes-gpu) LANES_GPU=$2; shift ;;
        --lanes-cpu) LANES_CPU=$2; shift ;;
        --lanes-light) LANES_LIGHT=$2; shift ;;
        --mem-budget) MEM_BUDGET_GB=$2; shift ;;
        --duckdb-mem-gb) DUCKDB_MEM_GB=$2; shift ;;
        --save-margin-min) SAVE_MARGIN_MIN=$2; shift ;;
        --script-set) SCRIPT_SET=$2; shift ;;
        --script|--heap|--override) CTL_PASS+=("$1" "$2"); shift ;;
        --allow-busy) ALLOW_BUSY=1 ;;
        --allow-drift) ALLOW_DRIFT=1 ;;
        --accept-external) ACCEPT=1 ;;
        --merge-provenance) MERGE=1 ;;
        --from) MERGE_FROM=$2; shift ;;
        --check) CHECK_ONLY=1 ;;
        --note) NOTE=$2; shift ;;
        --min-root-free-gb) MIN_ROOT_FREE_GB=$2; shift ;;
        --data-margin-gb) DATA_MARGIN_GB=$2; shift ;;
        --root-path) ROOT_PATH=$2; shift ;;
        --run-id) RUN_ID=$2; shift ;;
        -h|--help) usage 0 ;;
        *) echo "unknown option: $1" >&2; usage 2 ;;
    esac
    shift
done
for v in LANES_GPU LANES_CPU LANES_LIGHT MEM_BUDGET_GB DUCKDB_MEM_GB SAVE_MARGIN_MIN; do
    [[ ${!v} =~ ^[0-9]+$ ]] || { echo "$v must be a whole number, not '${!v}'" >&2; exit 2; }
done

PROJECT=$(realpath -e "$PROJECT")
MASTER=$(realpath -e "$PROJ/project.qpproj")
[[ -n $OUT ]] || { if [[ $PROJECT == "$MASTER" ]]; then OUT=$V3; else OUT=$(dirname "$PROJECT")/results/v3; fi; }
RUN_ID=${RUN_ID:-$(date +%Y%m%dT%H%M%S)-$$}
RUNDIR=$WORKROOT/$RUN_ID
CTL_BASE=(--repo "$PROJ" --project "$PROJECT" --out "$OUT" --script-set "$SCRIPT_SET" --run-id "$RUN_ID" "${CTL_PASS[@]}")
ctl() { "$PY" "$CTL" "${CTL_BASE[@]}" "$@" < /dev/null; }

DRIVER_LOG=/dev/null
say() { printf '%s %s\n' "$(date +%T)" "$*" | tee -a "$DRIVER_LOG"; }
log() { say "$*"; }                                    # lib.sh's assert_project_idle reports through these
die() { say "FATAL: $*"; exit 1; }

# Real QuPath JVMs (GUI or headless) and Senior rebuild processes, matched by process name. lib.sh's
# assert_project_idle greps whole command lines, so it also flags shells (e.g. `bash -c "... QuPath ..."`
# wrappers of other sessions) that merely mention QuPath; this narrower list decides.
busy_procs() {
    ps -eo pid=,comm=,args= | awk '$2 == "QuPath" || ($2 == "java" && /qupath\.QuPath/) ||
        ($2 ~ /^(python[0-9.]*|make)$/ && /build_layers|make rebuild|refine\.py/)'
}
gui_open() {   # a QuPath GUI (not a headless 'script' run) could be editing any image
    ((ALLOW_BUSY)) && return 1
    ps -eo pid=,comm=,args= | awk '($2 == "QuPath" || ($2 == "java" && /qupath\.QuPath/)) && !/ script /' | grep -q .
}
project_idle() {   # 0 = idle. assert_project_idle first; if it objects, only real processes count.
    (assert_project_idle) > /dev/null 2>&1 && return 0
    local real
    real=$(busy_procs)
    if [[ -z $real ]]; then
        echo 'note: assert_project_idle matched only processes that mention QuPath/rebuild in their arguments (shell wrappers); no QuPath or rebuild process is running' >&2
        return 0
    fi
    printf '%s\n' "$real" | cut -c1-200 | sed 's/^/  busy: /' >&2
    return 1
}

IMG_OPTS=()
[[ -n $IMAGES ]] && IMG_OPTS+=(--images "$IMAGES")
[[ -n $IMAGES_FILE ]] && IMG_OPTS+=(--images-file "$IMAGES_FILE")

# ---- one-shot modes (v3ctl takes the image locks itself) -----------------------------------------
if ((ACCEPT)); then
    [[ ${#IMG_OPTS[@]} -gt 0 && -n $NOTE ]] || { echo '--accept-external needs --images (or --images-file) and --note' >&2; exit 2; }
    ACC_OPTS=()
    [[ -n $RERUN_FROM ]] && ACC_OPTS+=(--rerun-from "$RERUN_FROM")
    exec "$PY" "$CTL" "${CTL_BASE[@]}" accept-external "${IMG_OPTS[@]}" --note "$NOTE" "${ACC_OPTS[@]}" \
        --lockdir "$LOCKDIR" --shared-lockdir "$SHARED_LOCKDIR"
fi
if ((MERGE)); then
    [[ ${#IMG_OPTS[@]} -gt 0 ]] || { echo '--merge-provenance needs --images (or --images-file)' >&2; exit 2; }
    M_OPTS=()
    [[ -n $MERGE_FROM ]] && M_OPTS+=(--from "$MERGE_FROM")
    ((CHECK_ONLY)) && M_OPTS+=(--check)
    exec "$PY" "$CTL" "${CTL_BASE[@]}" merge-provenance "${IMG_OPTS[@]}" "${M_OPTS[@]}" \
        --lockdir "$LOCKDIR" --shared-lockdir "$SHARED_LOCKDIR"
fi
PLAN_OPTS=()
[[ -n $PASSES ]] && PLAN_OPTS+=(--passes "$PASSES")
[[ -n $STEPS ]] && PLAN_OPTS+=(--steps "$STEPS")
[[ -n $FORCE ]] && PLAN_OPTS+=(--force "$FORCE")
[[ -n $RERUN_FROM ]] && PLAN_OPTS+=(--rerun-from "$RERUN_FROM")
((ALLOW_DRIFT)) && PLAN_OPTS+=(--allow-drift)
[[ -n $PASSES || -n $STEPS ]] || { echo 'give --pass or --steps' >&2; usage 2; }
FIN_OPTS=()
((ALLOW_DRIFT)) && FIN_OPTS+=(--allow-drift)
DISK_OPTS=(--data-margin-gb "$DATA_MARGIN_GB" --min-root-free-gb "$MIN_ROOT_FREE_GB" --root-path "$ROOT_PATH")
# the global plan; a planning error (bad --images, unknown step, --force outside the selection,
# missing script) stops here. The resolved image list goes to the terminal (stderr).
plan_txt=$(ctl plan "${PLAN_OPTS[@]}" "${IMG_OPTS[@]}" --format tsv) || { echo 'planning failed (see above)' >&2; exit 2; }

if ((DRY)); then
    ctl plan "${PLAN_OPTS[@]}" "${IMG_OPTS[@]}" --format table 2> /dev/null
    echo
    if project_idle; then echo 'busy check : idle (a real run may start)'
    else echo 'busy check : BUSY -- a real run would refuse (QuPath or a Senior rebuild is running)'; fi
    while IFS=$'\t' read -r kind tag _; do
        [[ $kind == IMAGE ]] || continue
        printf 'disk guard : %s: %s\n' "$tag" "$(ctl diskcheck --image "$tag" "${DISK_OPTS[@]}" 2>&1 || true)"
    done <<< "$plan_txt"
    echo 'dry run: nothing was changed'
    exit 0
fi

# ---- refuse to start while something else could write the project ------------------------------
if ! project_idle; then
    if ((ALLOW_BUSY)); then
        echo 'WARNING: QuPath or a Senior rebuild is running; --allow-busy given, so continuing with one CPU lane.' >&2
        LANES_CPU=1
    else
        echo 'refusing to start (close QuPath / wait for the rebuild, or pass --allow-busy)' >&2
        exit 3
    fi
fi

# a reused --run-id starts clean: a STOP or status left by an earlier run must not steer this one
rm -rf -- "$RUNDIR/status" "$RUNDIR/STOP"
mkdir -p "$RUNDIR/status" "$LOCKDIR/mem" "$TMPROOT" "$OUT/logs/run_v3"
DRIVER_LOG=$OUT/logs/run_v3/$RUN_ID.driver.log
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || true)
say "run $RUN_ID: project $PROJECT -> $OUT; site ${V3_SITE:-$(hostname -s)}; passes '${PASSES}' steps '${STEPS}'; lanes gpu=$LANES_GPU cpu=$LANES_CPU light=$LANES_LIGHT; mem budget ${MEM_BUDGET_GB}g; scripts $SCRIPT_SET${SHARED_LOCKDIR:+; shared locks $SHARED_LOCKDIR}${V3_DEADLINE_EPOCH:+; job deadline $(date -d "@$V3_DEADLINE_EPOCH" -Is)}"

stop_requested() { [[ -e $RUNDIR/STOP ]]; }
runlog() { ctl runlog --image "$1" --step "$2" --action "$3" --status "$4" --message "${5:-}" > /dev/null 2>&1 || true; }
request_stop() { printf '%s %s\n' "$(date -Is)" "$*" >> "$RUNDIR/STOP"; say "STOP: $* -- running steps finish, no new step starts"; }
trap 'request_stop "interrupted by signal"' INT TERM

# ---- lane slots and the memory budget (flock files shared by every run_v3 on this host) --------
SLOT_FD="" SLOT_NAME="" MEM_FD="" MEM_FILE=""
acquire_slot() {   # $1 lane: take any free slot 1..N, waiting if all are busy
    local lane=$1 n i fd
    case $lane in gpu) n=$LANES_GPU ;; cpu) n=$LANES_CPU ;; *) n=$LANES_LIGHT ;; esac
    ((n >= 1)) || { say "no $lane lane on this run (--lanes-$lane $n)"; return 2; }
    while :; do
        for ((i = 1; i <= n; i++)); do
            exec {fd}>>"$LOCKDIR/lane-$lane-$i.lock"
            if flock -n "$fd"; then SLOT_FD=$fd SLOT_NAME=$lane-$i; return 0; fi
            exec {fd}>&-
        done
        stop_requested && return 1
        sleep 5
    done
}
release_slot() { [[ -n $SLOT_FD ]] && exec {SLOT_FD}>&-; SLOT_FD="" SLOT_NAME=""; }

mem_acquire() {   # $1 GB. A reservation is a locked file; a file nobody holds (its holder died) is stale.
    local gb=$1 f used glock
    if ((gb > MEM_BUDGET_GB)); then say "reservation ${gb}g exceeds --mem-budget ${MEM_BUDGET_GB}g; it could never start"; return 1; fi
    while :; do
        exec {glock}>>"$LOCKDIR/mem.lock"
        flock "$glock"
        used=0
        for f in "$LOCKDIR"/mem/*.gb; do
            [[ -e $f ]] || continue
            if flock -n "$f" true; then rm -f "$f"; continue; fi
            used=$((used + $(cat "$f" 2>/dev/null || echo 0)))
        done
        if ((used + gb <= MEM_BUDGET_GB)); then
            MEM_FILE=$LOCKDIR/mem/$RUN_ID.$BASHPID.$RANDOM.gb
            exec {MEM_FD}>"$MEM_FILE"
            flock "$MEM_FD"
            echo "$gb" > "$MEM_FILE"
            exec {glock}>&-
            return 0
        fi
        exec {glock}>&-
        stop_requested && return 1
        sleep 10
    done
}
mem_release() { [[ -n $MEM_FILE ]] && rm -f "$MEM_FILE"; [[ -n $MEM_FD ]] && exec {MEM_FD}>&-; MEM_FD="" MEM_FILE=""; }

# ---- one step of one image --------------------------------------------------------------------
run_step() {   # name tag step lane heap threads timeout_min post  -> 0 pass, 1 fail, 3 stopped, 4 drift
    local name=$1 tag=$2 step=$3 lane=$4 heap=$5 threads=$6 tmo=$7 post=$8
    local logf=$OUT/logs/run_v3/$tag/$step.log worker prepf preperr msg rc=0 frc=0 t0 wall resv deadline now line key
    mkdir -p "${logf%/*}" "$RUNDIR/$tag"
    if stop_requested; then runlog "$tag" "$step" stop STOPPED 'run stop requested earlier'; return 3; fi
    if gui_open; then request_stop 'a QuPath GUI is open'; runlog "$tag" "$step" stop STOPPED 'QuPath GUI open'; return 3; fi
    if ! msg=$(ctl diskcheck --image "$tag" "${DISK_OPTS[@]}" 2>&1); then
        request_stop "free-space guard before $tag/$step: $msg"
        runlog "$tag" "$step" stop STOPPED "$msg"
        return 3
    fi
    # the footer starts no save later than SAVE_MARGIN_MIN before the step's deadline: its timeout, or
    # the end of the SLURM job that runs this driver, whichever comes first
    now=$(date +%s)
    deadline=$((now + tmo * 60))
    if [[ -n ${V3_DEADLINE_EPOCH:-} ]] && ((V3_DEADLINE_EPOCH < deadline)); then
        deadline=$V3_DEADLINE_EPOCH
        if ((deadline - now < (SAVE_MARGIN_MIN + 5) * 60)); then
            request_stop "only $(((deadline - now) / 60)) min left before the job deadline; $tag/$step could not save"
            runlog "$tag" "$step" stop STOPPED 'job deadline too close'
            return 3
        fi
    fi
    if ! acquire_slot "$lane"; then
        if stop_requested; then runlog "$tag" "$step" stop STOPPED 'stopped while waiting for a slot'; return 3; fi
        runlog "$tag" "$step" slot FAIL "no $lane lane"; return 1
    fi
    # the JVM heap; a step with a Parquet post hook also needs DuckDB's memory after QuPath exits
    resv=$heap
    if [[ $post != - ]] && ((DUCKDB_MEM_GB + 8 > resv)); then resv=$((DUCKDB_MEM_GB + 8)); fi
    if ! mem_acquire "$resv"; then
        release_slot
        if stop_requested; then runlog "$tag" "$step" stop STOPPED 'stopped while waiting for memory'; return 3; fi
        runlog "$tag" "$step" memory FAIL "reservation ${resv}g > --mem-budget ${MEM_BUDGET_GB}g"; return 1
    fi
    worker=$RUNDIR/$tag/$step.$(date +%s%N)
    prepf=$worker.prepare.out preperr=$worker.prepare.err
    local ATTEMPT="" SEQ="" HEAP="" WSCRIPT="" WPROJECT="" PROVENANCE="" SCRIPT_SRC="" QARGS_FILE="" QARGS=()
    # prepare prints KEY=VALUE lines on stdout (parsed below, never eval'ed); its stderr stays separate
    if ! ctl prepare --image "$tag" --step "$step" --workdir "$worker" --tmp "$TMPROOT" \
            --slot "$SLOT_NAME" --gpu "$([[ $lane == gpu ]] && echo "$GPU_NAME")" \
            --deadline-epoch "$deadline" --save-margin-min "$SAVE_MARGIN_MIN" > "$prepf" 2> "$preperr"; then
        msg=$(cat "$preperr" "$prepf" 2> /dev/null | tail -c 1500)
        say "$tag: $step: PREPARE FAILED: $msg"
        runlog "$tag" "$step" prepare FAIL "$msg"
        mem_release; release_slot
        return 1
    fi
    while IFS= read -r line; do
        key=${line%%=*}
        case $key in
            ATTEMPT|SEQ|HEAP|WSCRIPT|WPROJECT|PROVENANCE|SCRIPT_SRC|QARGS_FILE) printf -v "$key" '%s' "${line#*=}" ;;
            *) msg="unexpected line from prepare: ${line:0:200}"; say "$tag: $step: PREPARE FAILED: $msg"
               runlog "$tag" "$step" prepare FAIL "$msg"; mem_release; release_slot; return 1 ;;
        esac
    done < "$prepf"
    if [[ -z $ATTEMPT || -z $WSCRIPT || -z $WPROJECT || ! -f $QARGS_FILE ]]; then
        say "$tag: $step: PREPARE FAILED: incomplete output ($(tr '\n' ' ' < "$prepf"))"
        runlog "$tag" "$step" prepare FAIL 'incomplete prepare output'; mem_release; release_slot; return 1
    fi
    mapfile -d '' -t QARGS < "$QARGS_FILE"
    [[ -s $preperr ]] && sed 's/^/prepare: /' "$preperr" >> "$logf"
    printf '\n===== %s run %s attempt %s seq %s slot %s heap %sg threads %s timeout %sm deadline %s\nscript %s\nsource %s\n' \
        "$(date -Is)" "$RUN_ID" "$ATTEMPT" "$SEQ" "$SLOT_NAME" "$HEAP" "$threads" "$tmo" "$(date -d "@$deadline" -Is)" \
        "$WSCRIPT" "$SCRIPT_SRC" >> "$logf"
    say "$tag: $step: start (slot $SLOT_NAME, heap ${HEAP}g, reserved ${resv}g, attempt $ATTEMPT)"
    t0=$SECONDS
    JAVA_TOOL_OPTIONS="-Xmx${HEAP}g -Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads${V3_JAVA_OPTS:+ $V3_JAVA_OPTS}" \
        /usr/bin/time -v -o "$worker/time.txt" \
        timeout --kill-after=120 "${tmo}m" \
        "$QUPATH" script -p "$WPROJECT" -i "$name" "${QARGS[@]}" "$WSCRIPT" < /dev/null >> "$logf" 2>&1 || rc=$?
    wall=$((SECONDS - t0))
    msg=$(ctl finish --image "$tag" --step "$step" --attempt "$ATTEMPT" --workdir "$worker" --exit-code "$rc" \
        --wall "$wall" --timefile "$worker/time.txt" --log "$logf" --slot "$SLOT_NAME" "${FIN_OPTS[@]}" 2>&1) || frc=$?
    say "$tag: $step: ${msg//$'\t'/ -- } [${wall}s, QuPath exit $rc]"
    if ((frc == 0)) && [[ $post != - ]]; then
        msg=$(ctl post --image "$tag" --step "$step" --attempt "$ATTEMPT" --tmp "$TMPROOT" \
            --duckdb-mem-gb "$DUCKDB_MEM_GB" 2>&1) || frc=$?
        say "$tag: $step: post $post: ${msg//$'\t'/ -- }"
    fi
    mem_release
    release_slot
    return "$frc"
}

# ---- one image: its chain of steps, under the image lock ---------------------------------------
run_image() {   # name tag  (runs in its own background subshell)
    local name=$1 tag=$2 status=$RUNDIR/status/$2 fd line kind step action lane heap threads tmo post reason rc
    local planout planerr=$RUNDIR/$2.plan.err lmsg n_img
    local -a lines F
    exec {fd}>>"$LOCKDIR/$tag.lock"
    if ! flock -n "$fd"; then
        echo LOCKED > "$status"; say "$tag: image lock held by another run -- skipped"
        runlog "$tag" - lock LOCKED 'image lock held by another run'; return 0
    fi
    if [[ -n $SHARED_LOCKDIR ]]; then   # cluster nodes do not share flocks: a mkdir lock on the shared filesystem
        if ! lmsg=$(ctl lock acquire --tag "$tag" --pid "$BASHPID" --shared-lockdir "$SHARED_LOCKDIR" 2>&1); then
            echo LOCKED > "$status"; say "$tag: shared image lock: $lmsg -- skipped"
            runlog "$tag" - lock LOCKED "shared lock: $lmsg"; return 0
        fi
        # shellcheck disable=SC2064   # expand now: this subshell's pid and paths
        trap "$(printf '%q ' "$PY" "$CTL" "${CTL_BASE[@]}" lock release --tag "$tag" --pid "$BASHPID" --shared-lockdir "$SHARED_LOCKDIR") > /dev/null 2>&1 < /dev/null || true" EXIT
    fi
    # re-plan under the lock: another run may have finished steps since the global plan. A failed plan
    # or one that does not describe exactly this image is an error, never "nothing to do".
    mkdir -p "$RUNDIR"
    if ! planout=$(ctl plan "${PLAN_OPTS[@]}" --images "$tag" --format tsv 2>> "$planerr"); then
        echo 'PLAN ERROR' > "$status"; say "$tag: re-plan failed: $(tail -3 "$planerr" | tr '\n' ' ')"
        runlog "$tag" - plan FAIL "re-plan failed: $(tail -3 "$planerr" | tr '\n' ' ')"; return 0
    fi
    n_img=$(awk -F'\t' -v t="$tag" '$1 == "IMAGE" && $2 == t' <<< "$planout" | wc -l)
    if [[ $n_img != 1 ]]; then
        echo 'PLAN ERROR' > "$status"; say "$tag: re-plan returned $n_img IMAGE lines for this image"
        runlog "$tag" - plan FAIL "re-plan returned $n_img IMAGE lines"; return 0
    fi
    mapfile -t lines <<< "$planout"
    for line in "${lines[@]}"; do
        IFS=$'\t' read -r -a F <<< "$line"
        kind=${F[0]:-}
        if [[ ${F[1]:-} != "$tag" ]]; then kind=BAD; fi
        case $kind in
            IMAGE)      # IMAGE tag name entry_id state message
                if [[ ${F[4]} != ok && ${F[4]} != fresh ]]; then
                    echo "BLOCKED (${F[4]})" > "$status"; say "$tag: BLOCKED: ${F[5]:-}"
                    runlog "$tag" - plan BLOCKED "${F[4]}: ${F[5]:-}"; return 0
                fi ;;
            STEP)       # STEP tag step action lane heap threads timeout post script reason
                IFS=$'\t' read -r _ _ step action lane heap threads tmo post _ reason <<< "$line"
                case $action in
                    skip) continue ;;
                    blocked) echo "BLOCKED at $step" > "$status"; say "$tag: $step: BLOCKED: $reason"
                             runlog "$tag" "$step" plan BLOCKED "$reason"; return 0 ;;
                    run) rc=0
                         run_step "$name" "$tag" "$step" "$lane" "$heap" "$threads" "$tmo" "$post" || rc=$?
                         case $rc in
                             0) ;;
                             3) echo "STOPPED before $step" > "$status"; return 0 ;;
                             4) echo "DRIFT at $step" > "$status"; return 0 ;;
                             *) echo "FAILED at $step" > "$status"; say "$tag: stopping this image (log: $OUT/logs/run_v3/$tag/$step.log)"; return 0 ;;
                         esac ;;
                    *) echo 'PLAN ERROR' > "$status"; say "$tag: plan error: $line"; return 0 ;;
                esac ;;
            *) echo 'PLAN ERROR' > "$status"; say "$tag: plan error: ${line:0:300}"; return 0 ;;
        esac
    done
    echo DONE > "$status"
}

# ---- main: global plan, then one chain runner per image, at most (sum of lanes) at a time -------
mapfile -t plan <<< "$plan_txt"
declare -A NAME_OF=() TODO=()
order=()
for line in "${plan[@]}"; do
    IFS=$'\t' read -r kind tag f3 f4 f5 _ <<< "$line"
    if [[ $kind == IMAGE ]]; then
        NAME_OF[$tag]=$f3
        order+=("$tag")
        say "image $tag = '$f3' (entry $f4): $f5"
        [[ $f5 != ok && $f5 != fresh ]] && TODO[$tag]=1
    fi
    [[ $kind == STEP && $f4 != skip ]] && TODO[$tag]=1
done
say "images: ${#order[@]} selected, ${#TODO[@]} with work to do"
MAX_INFLIGHT=$((LANES_GPU + LANES_CPU + LANES_LIGHT))
for tag in "${order[@]}"; do
    if [[ -z ${TODO[$tag]:-} ]]; then echo 'DONE (all skipped)' > "$RUNDIR/status/$tag"; continue; fi
    while (($(jobs -rp | wc -l) >= MAX_INFLIGHT)); do wait -n || true; done
    stop_requested && break
    run_image "${NAME_OF[$tag]}" "$tag" &
done
while [[ -n $(jobs -rp) ]]; do wait || true; done

# ---- summary ------------------------------------------------------------------------------------
say "== run $RUN_ID finished"
n_bad=0
for tag in "${order[@]}"; do
    st=$(cat "$RUNDIR/status/$tag" 2>/dev/null || echo 'NOT STARTED')
    say "  $(printf '%-24s' "$st") $tag"
    [[ $st == DONE* ]] || n_bad=$((n_bad + 1))
done
ctl report --workroot "$WORKROOT" 2>&1 | tee -a "$DRIVER_LOG"
[[ -e $RUNDIR/STOP ]] && { say "stop reason(s):"; sed 's/^/  /' "$RUNDIR/STOP" | tee -a "$DRIVER_LOG"; exit 3; }
((n_bad == 0))
