#!/usr/bin/env bash
# senior_v3.sh -- run one command inside the v3 Senior worktree, with the worktree's own code and data.
#
# The worktree shares the main checkout's .venv, whose editable install resolves `import phenocycler`
# to the MAIN checkout (measured 2026-10-01: phenocycler.roots.REPO = .../Islet-Explorer-Senior). A
# script that goes through the installed package would then read and write the v2 data. PYTHONPATH
# puts the worktree's package first; this wrapper refuses to run anything unless that took effect.
#
# Usage: senior_v3.sh [--cd SUBDIR] <command...>      e.g.  senior_v3.sh --cd core/phenocycler python -m phenocycler.pipeline run ...
set -euo pipefail
WT=/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3
SUB=""
[[ ${1:-} == --cd ]] && { SUB=$2; shift 2; }
(($#)) || { echo "usage: senior_v3.sh [--cd SUBDIR] <command...>" >&2; exit 2; }
export PYTHONPATH=$WT/core/phenocycler${PYTHONPATH:+:$PYTHONPATH}
export ISLET_LAYERS_DIR=/data/islet_layers_v3
export PATH=$WT/.venv/bin:$PATH
got=$("$WT/.venv/bin/python" -c 'import phenocycler.roots as r; print(r.REPO)')
[[ $got == "$WT" ]] || { echo "REFUSED: phenocycler resolves to $got, not the v3 worktree $WT" >&2; exit 3; }
cd "$WT/$SUB"
exec "$@"
