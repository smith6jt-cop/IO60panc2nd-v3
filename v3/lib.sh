# lib.sh -- shared settings for the v3 re-run helpers. Source it; do not execute it.
#
# Every helper writes what it did to results/v3/logs/ and what it produced to
# results/v3/manifests/, so a later step can check the earlier one instead of trusting it.

PROJ=/home/smith6jt/IO60panc2nd
V3=$PROJ/results/v3
SCRATCH=/scratch/io60
ARCHIVE=/data/IO60panc2nd_archive/v2_2026-09
QUPATH=/home/smith6jt/QuPath/bin/QuPath
SHARE=/mnt/ahc_share/PERSONAL/smith6jt

mkdir -p "$V3/logs" "$V3/manifests"

log() {
    printf '%s [%s] %s\n' "$(date -Is)" "${0##*/}" "$*" | tee -a "$V3/logs/phase0.log" >&2
}

die() {
    log "FATAL: $*"
    exit 1
}

# Refuse to run while something could be writing the project: QuPath (GUI or headless) or a
# Senior layer rebuild. Headless QuPath runs that this helper starts itself pass SKIP_BUSY_CHECK=1.
assert_project_idle() {
    [[ "${SKIP_BUSY_CHECK:-0}" == 1 ]] && return 0
    # Match on the process name (comm), not the whole command line: a shell that merely mentions
    # QuPath or refine.py in its arguments (another session's bash -c wrapper, a grep) is not busy.
    local busy
    busy=$(ps -eo pid=,comm=,args= | awk '$2 == "QuPath" || ($2 == "java" && /qupath\.QuPath/) ||
            ($2 ~ /^(python[0-9.]*|make)$/ && /build_layers|make rebuild|refine\.py/)' || true)
    [[ -z "$busy" ]] || die "project busy; close QuPath / wait for the rebuild first:
$(cut -c1-200 <<< "$busy")"
}
