#!/usr/bin/env bash
# verify_share_backup.sh -- is every project image also on the departmental share?
#
# The share (/mnt/ahc_share/PERSONAL/smith6jt) is the backup that lets a local image copy be
# relocated. This does not trust names alone: a copy counts only when its byte size matches.
# Read-only on both sides. Writes results/v3/manifests/share_backup_check.tsv and exits non-zero
# unless every image has a same-size copy.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

listing=$V3/manifests/share_qptiff_listing.tsv
out=$V3/manifests/share_backup_check.tsv

log "listing *.qptiff under $SHARE (skipping .temp scan intermediates)"
find "$SHARE" -type d -name .temp -prune -o -type f -iname '*.qptiff' -printf '%s\t%p\n' \
    2>/dev/null | sort -k2 > "$listing" || true
log "share listing: $(wc -l < "$listing") qptiff files"

printf 'image\tlocal_bytes\tstatus\tshare_bytes\tshare_path\n' > "$out"
n_ok=0; n_bad=0
for f in "$PROJ"/Images/*.qptiff; do
    name=$(basename "$f")
    size=$(stat -L -c %s "$f")
    # exact name + size first; then size alone, to catch a renamed copy
    hit=$(awk -F'\t' -v n="$name" -v s="$size" '{p=$2; sub(/.*\//,"",p)} p==n && $1==s {print $1"\t"$2; exit}' "$listing")
    if [[ -n "$hit" ]]; then
        status=MATCH
    else
        hit=$(awk -F'\t' -v s="$size" '$1==s {print $1"\t"$2; exit}' "$listing")
        if [[ -n "$hit" ]]; then
            status=SIZE_ONLY_MATCH
        else
            hit=$(awk -F'\t' -v n="$name" '{p=$2; sub(/.*\//,"",p)} p==n {print $1"\t"$2; exit}' "$listing")
            status=$([[ -n "$hit" ]] && echo SIZE_MISMATCH || echo MISSING)
            [[ -n "$hit" ]] || hit=$'-\t-'
        fi
    fi
    printf '%s\t%s\t%s\t%s\n' "$name" "$size" "$status" "$hit" >> "$out"
    if [[ $status == MATCH ]]; then n_ok=$((n_ok + 1)); else n_bad=$((n_bad + 1)); fi
done

log "share backup check: $n_ok MATCH, $n_bad not matched -> $out"
column -t -s $'\t' "$out" >&2
[[ $n_bad -eq 0 ]]
