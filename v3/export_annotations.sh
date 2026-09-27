#!/usr/bin/env bash
# export_annotations.sh -- run export_annotations.groovy on every project image, headless.
#
# Runs against a PRIVATE copy of project.qpproj (data/, classifiers/, scripts/ symlinked), so
# nothing QuPath does can touch the master project file. Sequential: one JVM at a time, heap
# capped, because the largest .qpdata is ~14 GB. Success is judged by the sentinel JSON the
# script writes, not by QuPath's exit code (headless QuPath exits 0 when a script throws).
set -euo pipefail
source "$(dirname "$0")/lib.sh"

assert_project_idle
work=$SCRATCH/work/export_ann
mkdir -p "$work" "$V3/logs/export_annotations"
cp -f "$PROJ/project.qpproj" "$work/project.qpproj"
for d in data classifiers scripts resources; do ln -sfn "$PROJ/$d" "$work/$d"; done

export JAVA_TOOL_OPTIONS="-Xmx160g"
status=$V3/manifests/export_annotations_status.tsv
printf 'image\ttag\tstatus\tseconds\n' > "$status"

mapfile -t names < <(python3 -c "
import json; p = json.load(open('$work/project.qpproj'))
print('\n'.join(e['imageName'] for e in p['images']))")

n_ok=0; n_bad=0
for name in "${names[@]}"; do
    tag=$(python3 -c "import re,sys; print(re.sub(r'[^A-Za-z0-9._-]+', '_', sys.argv[1]))" "$name")
    sentinel=$V3/manual_annotations/$tag.json
    t0=$SECONDS
    rm -f "$sentinel"
    SKIP_BUSY_CHECK=1 timeout 60m "$QUPATH" script -p "$work/project.qpproj" -i "$name" \
        -a "$V3" "$PROJ/scripts/v3/export_annotations.groovy" \
        > "$V3/logs/export_annotations/$tag.log" 2>&1 || true
    dt=$((SECONDS - t0))
    if [[ -s $sentinel ]]; then
        st=OK; n_ok=$((n_ok + 1))
        log "ok $name (${dt}s): $(python3 -c "import json; d=json.load(open('$sentinel')); print(d['manual_class_counts'], d['n_annotations'], 'annotations')")"
    else
        st=FAILED; n_bad=$((n_bad + 1))
        log "FAILED $name (${dt}s) -- see logs/export_annotations/$tag.log"
    fi
    printf '%s\t%s\t%s\t%s\n' "$name" "$tag" "$st" "$dt" >> "$status"
done

log "annotation export: $n_ok ok, $n_bad failed -> $status"
[[ $n_bad -eq 0 ]]
