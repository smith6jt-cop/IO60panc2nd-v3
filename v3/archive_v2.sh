#!/usr/bin/env bash
# archive_v2.sh -- compressed, verified copy of the v2 QuPath state before v3 overwrites it.
#
# One zstd tar per data/<entry>/ folder (all of them, orphans included), plus one tar of the
# project-level files (project.qpproj, classifiers, scripts, resources, notes, small results).
# Each archive is verified three ways before it is marked done:
#   1. zstd -t (the stream decodes),
#   2. the tar listing matches the source file list by path and byte size,
#   3. tar reported no "file changed as we read it" (exit 1) -- a live writer means a torn copy.
# Idempotent: an archive with a .ok marker is skipped, so an interrupted run just resumes.
# This COPIES; nothing under the project is modified or removed.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

assert_project_idle
mkdir -p "$ARCHIVE"
manifest=$ARCHIVE/MANIFEST.tsv
[[ -f $manifest ]] || printf 'archive\tsource\tn_files\tsource_bytes\tarchive_bytes\tsha256\tverified_at\n' > "$manifest"

NICE=(nice -n 10 ionice -c2 -n7)

# archive_one <name> <tar -C dir> <member>...
archive_one() {
    local name=$1 cdir=$2; shift 2
    local dest=$ARCHIVE/$name.tar.zst ok=$ARCHIVE/$name.tar.zst.ok
    if [[ -f $ok ]]; then log "skip $name (verified earlier)"; return 0; fi

    local srclist=$SCRATCH/tmp/$name.src.tsv tarlist=$SCRATCH/tmp/$name.tar.tsv
    (cd "$cdir" && find "$@" -type f -printf '%p\t%s\n' | sort) > "$srclist"
    local nfiles srcbytes
    nfiles=$(wc -l < "$srclist")
    srcbytes=$(awk -F'\t' '{s+=$2} END {print s+0}' "$srclist")
    log "archiving $name: $nfiles files, $srcbytes bytes"

    local rc=0
    "${NICE[@]}" tar -C "$cdir" -cf - "$@" | "${NICE[@]}" zstd -q -T16 -3 -f -o "$dest.part" || rc=$?
    [[ $rc -eq 0 ]] || die "$name: tar|zstd failed (rc=$rc; 1 = a file changed while being read)"

    "${NICE[@]}" zstd -q -t "$dest.part" || die "$name: zstd -t failed"
    # regular files only (dirs/symlinks are not in the find list); literal quoting so names with
    # spaces or non-ASCII compare byte-for-byte with find's output
    "${NICE[@]}" zstd -q -dc "$dest.part" | tar --quoting-style=literal -tvf - \
        | awk '$1 ~ /^-/ {size=$3; $1=$2=$3=$4=$5=""; sub(/^ +/,""); print $0"\t"size}' \
        | sort > "$tarlist"
    cmp -s "$srclist" "$tarlist" || die "$name: archive listing differs from source (see $srclist vs $tarlist)"

    mv "$dest.part" "$dest"
    local sha abytes
    sha=$("${NICE[@]}" sha256sum "$dest" | cut -d' ' -f1)
    abytes=$(stat -c %s "$dest")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name.tar.zst" "$cdir/$*" "$nfiles" "$srcbytes" "$abytes" "$sha" "$(date -Is)" >> "$manifest"
    printf 'files=%s src_bytes=%s archive_bytes=%s sha256=%s\n' "$nfiles" "$srcbytes" "$abytes" "$sha" > "$ok"
    rm -f "$srclist" "$tarlist"
    log "ok $name: $srcbytes -> $abytes bytes"
}

# 1. project-level files first (small; the irreplaceable pieces)
meta=()
# scripts/v3 is the new work in progress, not v2 state (and it changes while this runs)
for p in project.qpproj project.qpproj.backup classifiers resources transformations \
         donor_metadata_panc.xlsx ProcessingNotes.txt README.txt xenium_paths.csv \
         xenium_paths.csv.bak-2026-08-25 scripts.zip \
         results/norm_stats.csv results/AnnotationsMeasurements.csv results/AnnotationMeasurementsBatch2.csv; do
    [[ -e $PROJ/$p ]] && meta+=("$p")
done
for p in "$PROJ"/scripts/*; do
    [[ $(basename "$p") == v3 ]] || meta+=("scripts/$(basename "$p")")
done
archive_one project_meta "$PROJ" "${meta[@]}"

# 2. every data/<entry>/ folder, smallest first so the quick ones land early
while read -r _ d; do
    archive_one "data_${d//\//_}" "$PROJ/data" "$d"
done < <(cd "$PROJ/data" && du -s -- */ | sed 's#/$##' | sort -n)

cp "$manifest" "$V3/manifests/v2_archive_MANIFEST.tsv"
log "archive complete: $(($(wc -l < "$manifest") - 1)) archives in $ARCHIVE"
