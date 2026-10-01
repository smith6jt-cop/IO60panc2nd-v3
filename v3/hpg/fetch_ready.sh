#!/usr/bin/env bash
# fetch_ready.sh -- fetch every image whose cells pass finished on HiPerGator and is not here yet.
# Run on the workstation. An image is "here" when its local data file is no longer the stripped one.
# Usage: fetch_ready.sh [--out DIR] <entry id>...
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=../lib.sh
source "$HERE/../lib.sh"
OUT_ARGS=()
[[ ${1:-} == --out ]] && { OUT_ARGS=(--out "$2"); shift 2; }
for id in "$@"; do
    size=$(stat -c %s "$PROJ/data/$id/data.qpdata" 2> /dev/null || echo 0)
    if ((size > 50 * 1024 * 1024)); then echo "entry $id: already here ($size bytes)"; continue; fi
    if "$HERE/fetch_cells.sh" "$id" "${OUT_ARGS[@]}" --check > /dev/null 2>&1; then
        "$HERE/fetch_cells.sh" "$id" "${OUT_ARGS[@]}" 2>&1 | tail -2
    else
        echo "entry $id: not ready on HiPerGator"
    fi
done
