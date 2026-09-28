#!/usr/bin/env bash
# setup_hpg.sh -- one-time (and re-runnable) layout check/repair on HiPerGator. Login node only.
#
# Builds the identical absolute paths the workstation uses, so run_v3.sh runs unchanged:
#   /home/smith6jt/IO60panc2nd/project.qpproj   copied from the workstation (Globus)
#   /home/smith6jt/IO60panc2nd/classifiers/     copied
#   /home/smith6jt/IO60panc2nd/scripts/         copied (includes scripts/v3 and this bundle)
#   /home/smith6jt/IO60panc2nd/data     -> $BLUE/data        (symlink; the .qpdata folders)
#   /home/smith6jt/IO60panc2nd/results  -> $BLUE/results     (symlink; provenance + exports, too big for home)
#   /home/smith6jt/IO60panc2nd/Images/  a real directory of per-image symlinks -> $BLUE/Images/<file>
#       (per-file rather than one directory symlink, so c0_gpu.sbatch can repoint ONE image at a
#        staged copy in $TMPDIR and restore it afterwards)
#   ~/QuPath (portable 0.7 build, bin/QuPath), ~/QuPath/v0.7 (extensions, InstanSeg model),
#   ~/.djl.ai (PyTorch engine cache) -- copied from the workstation, never downloaded on a compute node
#   $BLUE/venv-v3 -- python3 with duckdb + pyarrow for the Parquet post hooks
#
# Usage: setup_hpg.sh [--check]     --check only reports; default also creates missing links/dirs.
# Nothing here copies data: the Globus transfers (globus_images_to_hpg.sh etc.) do that.
set -euo pipefail
# shellcheck source=hpg_env.sh
source "$(dirname "$0")/hpg_env.sh"
CHECK_ONLY=0
PHENO=${PHENO:-/blue/maigan/smith6jt/Xenium_Analysis/data/phenocycler}   # metabolic-panel qptiffs (Xenium lane)
[[ ${1:-} == --check ]] && CHECK_ONLY=1
problems=0
ok()  { printf '  ok       %s\n' "$*"; }
bad() { printf '  PROBLEM  %s\n' "$*"; problems=$((problems + 1)); }
fix() {   # fix <description> <command...>: run it, unless --check
    local what=$1; shift
    if ((CHECK_ONLY)); then bad "$what (missing; run without --check to create)"; return 0; fi
    "$@"
    printf '  created  %s\n' "$what"
}
link_dir() { mkdir -p "$BLUE/$1" && ln -s "$BLUE/$1" "$PROJ/$1"; }
make_venv() {
    # TODO(HiPerGator): confirm the module name with `module spider python` (e.g. python/3.12)
    if command -v module > /dev/null 2>&1; then module load python || true; fi
    python3 -m venv "$BLUE/venv-v3"
    "$BLUE/venv-v3/bin/pip" install -q duckdb pyarrow
}

echo "== project layout under $PROJ (identical to the workstation)"
if [[ -f $PROJ/project.qpproj ]]; then ok project.qpproj; else bad "$PROJ/project.qpproj (Globus it from the workstation)"; fi
if [[ -d $PROJ/classifiers ]]; then ok classifiers/; else bad "$PROJ/classifiers/"; fi
if [[ -f $V3DIR/run_v3.sh ]]; then ok scripts/v3/run_v3.sh; else bad "$V3DIR/run_v3.sh (copy scripts/)"; fi
for d in data results; do
    if [[ -L $PROJ/$d && $(readlink "$PROJ/$d") == "$BLUE/$d" ]]; then
        ok "$d -> $BLUE/$d"
    elif [[ -e $PROJ/$d ]]; then
        bad "$PROJ/$d exists but is not a symlink to $BLUE/$d -- fix by hand"
    else
        fix "$d -> $BLUE/$d" link_dir "$d"
    fi
done

echo "== Images/: one symlink per image named in project.qpproj"
if [[ -L $PROJ/Images ]]; then
    bad "$PROJ/Images is a directory symlink; replace it with a directory of per-image symlinks"
elif [[ ! -d $PROJ/Images ]]; then
    fix "$PROJ/Images/" mkdir -p "$PROJ/Images"
fi
if [[ -f $PROJ/project.qpproj && -d $PROJ/Images ]]; then
    while IFS= read -r uri; do
        f=${uri#file:}
        name=${f##*/}
        target=$BLUE/Images/$name
        if [[ $f != "$PROJ/Images/$name" ]]; then bad "image URI outside $PROJ/Images: $uri"; continue; fi
        if [[ ! -e $target && -s $PHENO/$name ]]; then
            # the metabolic-panel images (6533, HDL112, 6347, 6580) have their one /blue copy in the
            # Xenium project's phenocycler folder; link it rather than transferring a second copy
            fix "$target -> $PHENO/$name" ln -s "$PHENO/$name" "$target"
            ((CHECK_ONLY)) && continue
        fi
        if [[ ! -s $target ]]; then bad "$target (not transferred yet)"; continue; fi
        if [[ -L $f && $(readlink "$f") == "$target" ]]; then
            ok "Images/$name"
        elif [[ -L $f ]]; then
            # a staged repoint left behind by a killed c0 task; the image itself is untouched on /blue
            bad "Images/$name -> $(readlink "$f") (expected $target)"
            if ((!CHECK_ONLY)); then ln -sfn "$target" "$f"; printf '  repaired Images/%s\n' "$name"; fi
        else
            fix "Images/$name" ln -s "$target" "$f"
        fi
    done < <(python3 -c 'import json, sys
for e in json.load(open(sys.argv[1]))["images"]: print(e["serverBuilder"]["uri"])' "$PROJ/project.qpproj")
fi

echo "== QuPath, InstanSeg model, DJL cache"
if [[ -x $V3_QUPATH ]]; then ok "$V3_QUPATH"; else bad "$V3_QUPATH (portable QuPath 0.7 Linux build)"; fi
model=$QUPATH_HOME/v0.7/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1
if [[ -d $model ]]; then ok 'InstanSeg model'; else bad "$model (script 0 hard-codes this path)"; fi
# QuPath 0.7's DJL needs PyTorch 2.7.1; an older engine (e.g. 2.5.1-cu124) would trigger a download on a
# compute node without internet, and 2.7.1 needs the CUDA 12.8 build for Blackwell (B200)
if [[ -d $DJL_CACHE_DIR/pytorch/2.7.1-cu128-linux-x86_64 && -f $DJL_CACHE_DIR/pytorch/2.7.1.txt ]]; then
    ok "DJL PyTorch 2.7.1-cu128 engine in $DJL_CACHE_DIR"
else
    bad "$DJL_CACHE_DIR/pytorch/2.7.1-cu128-linux-x86_64 + 2.7.1.txt (copy them from the workstation's ~/.djl.ai)"
fi
# the pipeline's own Java preference tree: userPath -> $QUPATH_HOME/v0.7, so headless QuPath loads the
# InstanSeg + DJL extensions from there (see hpg_env.sh). Node names are Java's encodings of the
# preference path QuPath 0.7 uses (copied from a working workstation tree).
prefs_node=$V3_PREFS_ROOT/.java/.userPrefs/$'_!\'k!bw!u!\'c!a@"0!\'g!d@"i!#4!c@"1!(!!}@"0!\'g='/'_!$!!.g!3'
make_prefs() {
    mkdir -p "$prefs_node"
    printf '%s\n' '<?xml version="1.0" encoding="UTF-8" standalone="no"?>' \
        '<!DOCTYPE map SYSTEM "http://java.sun.com/dtd/preferences.dtd">' '<map MAP_XML_VERSION="1.0">' \
        "  <entry key=\"userPath\" value=\"$QUPATH_HOME/v0.7\"/>" '</map>' > "$prefs_node/prefs.xml"
}
# QuPath DELETES a userPath whose folder does not exist (and saves that on exit), so write it only after
# ~/QuPath/v0.7 has been copied over
if [[ ! -d $QUPATH_HOME/v0.7/extensions ]]; then
    bad "$QUPATH_HOME/v0.7/extensions (copy ~/QuPath from the workstation, then re-run setup_hpg.sh)"
elif grep -q "key=\"userPath\" value=\"$QUPATH_HOME/v0.7\"" "$prefs_node/prefs.xml" 2> /dev/null; then
    ok "pipeline QuPath preferences (userPath -> $QUPATH_HOME/v0.7)"
else
    fix "pipeline QuPath preferences in $V3_PREFS_ROOT" make_prefs
fi
# TODO(HiPerGator): run hpg/djl_cuda_check.groovy once in a short GPU session of each GPU type before the
# first array (c0_gpu.sbatch runs it anyway and stops on failure).

echo "== python for the v3ctl.py post hooks"
if [[ -x $PY ]] && "$PY" -c 'import duckdb, pyarrow' 2> /dev/null; then
    ok "$PY (duckdb, pyarrow)"
else
    fix "$BLUE/venv-v3 (duckdb, pyarrow)" make_venv
fi

echo "== working areas"
for d in "$BLUE/work/driver" "$BLUE/tmp/run_v3" "$BLUE/logs" "$BLUE/locks" "$BLUE/work/submit/jobs"; do
    if [[ -d $d ]]; then ok "$d"; else fix "$d" mkdir -p "$d"; fi
done

echo
if ((problems)); then
    echo "$problems problem(s) -- fix them before submitting"
    exit 1
fi
echo "layout ok (provenance site: $V3_SITE). Before every submission sync data/<id>/ and results/v3/provenance/"
echo 'from the workstation, then adopt its records: run_v3.sh --merge-provenance --images-file <file>.'
echo 'Until then run_v3.sh (and hpg_submit.sh) block those images as unmerged.'
