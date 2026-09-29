#!/usr/bin/env bash
# relocate_images.sh -- move the largest qptiffs from /home to the second NVMe (/scratch/io60),
# leaving a per-file symlink at the original path.
#
# Why: /home is ~92% full and the v3 run needs room, but the images must stay on NVMe for the
# visual reviews (the /data HDD is too slow for QuPath on 30-170 GB images). QuPath URIs and the
# Senior manifests keep working because Images/<name> still resolves.
# Safety, per image:
#   - only if the share backup check matched it (results/v3/manifests/share_backup_check.tsv;
#     6551's renamed copy was verified by sampled chunks) -- /scratch is never the only copy;
#   - only if nothing has the file open and no QuPath/rebuild is running;
#   - rsync (verifies each transferred file), then an independent full-read checksum compare;
#   - the swap is an atomic rename of a symlink over the original, so a reader never sees a gap.
# Idempotent: an image whose Images/<name> is already a symlink to the destination is skipped.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

DEST=$SCRATCH/Images
MIN_ROOT_FREE_GB=300
# The 11 largest (~1.13 TB), chosen so / keeps >= 300 GB free for the OS and scratch.
IMAGES=(6505panc_Scan1.er.qptiff 6374_Scan1.er.qptiff 6551_Scan1.er.qptiff 6566panc_Scan1.er.qptiff
        6579_Scan2.er.qptiff 6476_Scan2.er.qptiff 6457_Scan1.er.qptiff 6623panc_Scan1.er.qptiff
        6593_Scan1.er.qptiff 6461_Scan1.er.qptiff 6442_Scan1.er.qptiff)

assert_project_idle
mkdir -p "$DEST"
check=$V3/manifests/share_backup_check.tsv
[[ -s $check ]] || die "run verify_share_backup.sh first"
manifest=$V3/manifests/relocated_images.tsv
[[ -f $manifest ]] || printf 'image\tbytes\tdest\tsha256\trelocated_at\n' > "$manifest"

for name in "${IMAGES[@]}"; do
    src=$PROJ/Images/$name dst=$DEST/$name
    if [[ -L $src && $(readlink -f "$src") == "$dst" ]]; then log "skip $name (already relocated)"; continue; fi
    [[ -f $src && ! -L $src ]] || die "$name: $src is not a regular file"
    status=$(awk -F'\t' -v n="$name" '$1==n {print $3}' "$check")
    [[ $status == MATCH || ( $name == 6551_Scan1.er.qptiff && $status == SIZE_ONLY_MATCH ) ]] \
        || die "$name: no verified share backup (status '$status')"
    if lsof -t -- "$src" >/dev/null 2>&1; then die "$name is open by another process"; fi

    bytes=$(stat -c %s "$src")
    free_gb=$(( $(df --output=avail -B1 "$SCRATCH" | tail -1) / 1024**3 ))
    need_gb=$(( bytes / 1024**3 + MIN_ROOT_FREE_GB ))
    [[ ! -e $dst ]] && (( free_gb < need_gb )) && die "$name: / has ${free_gb} GB free, need ${need_gb}"

    log "copying $name ($bytes bytes)"
    rsync -a --partial "$src" "$dst"
    [[ $(stat -c %s "$dst") == "$bytes" ]] || die "$name: size mismatch after copy"
    # independent verification: read both copies in full (in parallel: they are on different
    # NVMe drives) and compare digests
    sha256sum "$src" > "$SCRATCH/tmp/$name.src.sha" & p1=$!
    sha256sum "$dst" > "$SCRATCH/tmp/$name.dst.sha" & p2=$!
    wait "$p1" "$p2"
    s1=$(cut -d' ' -f1 < "$SCRATCH/tmp/$name.src.sha")
    s2=$(cut -d' ' -f1 < "$SCRATCH/tmp/$name.dst.sha")
    rm -f "$SCRATCH/tmp/$name.src.sha" "$SCRATCH/tmp/$name.dst.sha"
    [[ -n $s1 && $s1 == "$s2" ]] || die "$name: checksum mismatch ($s1 vs $s2)"

    if lsof -t -- "$src" >/dev/null 2>&1; then die "$name opened during the copy; not swapping"; fi
    ln -s "$dst" "$src.v3tmp"
    mv -T "$src.v3tmp" "$src"            # atomic replace of the file by the symlink
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$bytes" "$dst" "$s2" "$(date -Is)" >> "$manifest"
    log "relocated $name -> $dst (sha256 ${s2:0:16}...)"
done
log "relocation done: $(($(wc -l < "$manifest") - 1)) images on $DEST; /home free $(df -h --output=avail /home | tail -1)"
