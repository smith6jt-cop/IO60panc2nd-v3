#!/usr/bin/env bash
# globus_backup.sh -- off-machine copy of the irreplaceable pieces to HiPerGator /orange.
#
# What cannot be regenerated: the hand-drawn / curated annotations (Tissue, pLN), the project
# file, classifiers and scripts (all inside project_meta.tar.zst), plus the v2 archive manifest
# and the full v2 annotation export for comparison. Checksum sync, so re-running only sends changes.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

LOCAL_EP=f1b69b9e-f07a-11ef-8c40-0e26ca329435
HPG_EP=5dbaf795-8a7e-4dca-91aa-6e10d610c2b3
DEST=/orange/maigan/IO60panc2nd_backup/$(date +%F)

batch=$V3/manifests/globus_backup_batch.txt
cat > "$batch" <<B
"$ARCHIVE/project_meta.tar.zst" "$DEST/project_meta.tar.zst"
"$ARCHIVE/MANIFEST.tsv" "$DEST/v2_archive_MANIFEST.tsv"
--recursive "$V3/manual_annotations" "$DEST/manual_annotations"
--recursive "$V3/annotations_all" "$DEST/annotations_all"
"$V3/manifests/export_annotations_status.tsv" "$DEST/export_annotations_status.tsv"
B
globus mkdir "$HPG_EP:/orange/maigan/IO60panc2nd_backup/" 2>/dev/null || true
globus mkdir "$HPG_EP:$DEST/" 2>/dev/null || true
task=$(globus transfer "$LOCAL_EP" "$HPG_EP" --batch "$batch" \
        --label "IO60panc2nd irreplaceables $(date +%F)" --sync-level checksum --preserve-timestamp \
        --verify-checksum --notify failed,inactive,succeeded --format unix --jmespath task_id)
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$task" backup "$DEST" >> "$V3/manifests/globus_tasks.tsv"
log "Globus task $task submitted (irreplaceables -> $DEST)"
echo "$task"
