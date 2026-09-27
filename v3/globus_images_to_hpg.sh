#!/usr/bin/env bash
# globus_images_to_hpg.sh -- upload the project qptiffs to HiPerGator /blue for GPU/CPU processing.
#
# Sends from the local NVMe (Images/), never from the departmental share (~7 MB/s per task).
# Images already on /blue with an identical byte size are skipped (6533 and HDL112 were uploaded
# 2026-09-15 under Xenium_Analysis/data/phenocycler; setup_hpg.sh symlinks them into place).
# Must run BEFORE relocate_images.sh: Globus Connect Personal here only exports ~/ and /data.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

LOCAL_EP=f1b69b9e-f07a-11ef-8c40-0e26ca329435          # this workstation (Globus Connect Personal)
HPG_EP=5dbaf795-8a7e-4dca-91aa-6e10d610c2b3            # UFRC HiPerGator
DEST=/blue/maigan/smith6jt/IO60panc2nd_v3/Images
ALREADY=/blue/maigan/smith6jt/Xenium_Analysis/data/phenocycler

batch=$V3/manifests/globus_images_batch.txt
: > "$batch"
existing=$(globus ls -l "$HPG_EP:$ALREADY/" --format unix --jmespath 'DATA[].[name,size]' 2>/dev/null || true)
for f in "$PROJ"/Images/*.qptiff; do
    name=$(basename "$f"); size=$(stat -L -c %s "$f")
    if grep -qP "^\Q$name\E\t$size$" <<< "$existing"; then
        log "skip $name (already on /blue at $ALREADY, same size)"; continue
    fi
    printf '"%s" "%s/%s"\n' "$f" "$DEST" "$name" >> "$batch"
done
log "uploading $(wc -l < "$batch") images to $DEST"

globus mkdir "$HPG_EP:/blue/maigan/smith6jt/IO60panc2nd_v3/" 2>/dev/null || true
globus mkdir "$HPG_EP:$DEST/" 2>/dev/null || true
task=$(globus transfer "$LOCAL_EP" "$HPG_EP" --batch "$batch" \
        --label "IO60panc2nd v3 images $(date +%F)" --preserve-timestamp --verify-checksum \
        --notify failed,inactive,succeeded --format unix --jmespath task_id)
printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$task" images "$DEST" >> "$V3/manifests/globus_tasks.tsv"
log "Globus task $task submitted (images -> $DEST)"
echo "$task"
