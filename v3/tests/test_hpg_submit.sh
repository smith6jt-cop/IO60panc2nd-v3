#!/usr/bin/env bash
# test_hpg_submit.sh -- hpg_submit.sh and hpg_env.sh against a mocked SLURM (no cluster, no QuPath).
#
# A synthetic 5-image project with provenance under site 'hpg' (as on HiPerGator):
#   A  preflight done                    -> c0 ready
#   B  preflight done, in our queued job -> refused (already queued/running)
#   C  preflight + step 0 done           -> c0: nothing to run; c2 ready
#   D  no record                         -> c0 blocked (preflight missing)
#   E  preflight + step 0 done           -> c2 ready
# Mocked: sacctmgr (GrpTRES per QoS), squeue (our jobs; the group's RUNNING and PENDING tres per QoS; a
# job's time left), sbatch (--parsable job ids, arguments logged), showQos, slurmInfo.
#   1  a token that is not an exact image          -> refused before anything else
#   2  c0: queued / nothing-to-run / blocked images -> refused (exit 2) and listed; --drop-unready submits
#      only A, with --time from steps.tsv, the GPU pinned, and the job registered for the next check
#   3  c2: PENDING jobs count against the shared QoS; maigan-b capped at --maigan-b-max; images split
#   4  hpg_image_for_task skips blank/comment lines; hpg_duration_s; hpg_check_budget refuses a short job
set -uo pipefail
V3DIR=$(cd "$(dirname "$0")/.." && pwd)
PY=/scratch/io60/venv-v3/bin/python
T=$(mktemp -d)
trap 'rm -rf -- "$T"' EXIT
FAILS=0
check() { local d=$1; shift; if "$@"; then printf '  PASS  %s\n' "$d"; else printf '  FAIL  %s\n' "$d"; FAILS=$((FAILS + 1)); fi; }

# ---- synthetic project + provenance (site hpg) ----------------------------------------------------
"$PY" - "$T" "$V3DIR" <<'EOF'
import json, sys
from pathlib import Path
T, V3DIR = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(V3DIR))
import v3ctl
names = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5}
proj = T / 'proj' / 'project.qpproj'
proj.parent.mkdir(parents=True)
proj.write_text(json.dumps({'images': [{'imageName': f'img {n}', 'entryID': i, 'serverBuilder': {}} for n, i in names.items()]}))
steps = v3ctl.load_steps(V3DIR / 'steps.tsv')
scripts = v3ctl.resolve_all_scripts(steps, Path('/home/smith6jt/IO60panc2nd'), 'v3', {})
for n, i in names.items():
    q = T / 'proj' / 'data' / str(i) / 'data.qpdata'
    q.parent.mkdir(parents=True)
    q.write_bytes(b'x' * (100 + i))
    done = {'A': ['preflight'], 'B': ['preflight'], 'C': ['preflight', '0'], 'D': [], 'E': ['preflight', '0']}[n]
    d = v3ctl.site_dir(T / 'out', v3ctl.tag_of(f'img {n}'), 'hpg')
    for seq, sid in enumerate(done, 1):
        s = v3ctl.stat_of(q)
        v3ctl.write_json(d / f'{sid}.json', {'schema': 'run_v3.provenance/1', 'step': sid, 'seq': seq, 'status': 'pass',
                                             'attempt_id': f'{n}{seq}', 'qpdata': {'before': s, 'after': s}, 'outputs': [],
                                             'script': {'source_path': scripts[sid]['path'], 'source_sha256': scripts[sid]['sha256']}})
        v3ctl.write_json(d / f'{sid}.finish.json', {'attempt_id': f'{n}{seq}', 'status': 'PASS'})
EOF

# ---- mocked SLURM -----------------------------------------------------------------------------------
mkdir -p "$T/bin"
cat > "$T/bin/sacctmgr" <<'EOF'
#!/bin/bash
case "$*" in
    *' clive '*)    echo 'cpu=40,mem=312G,gres/gpu=1' ;;
    *' maigan-b '*) echo 'cpu=720,mem=5625G' ;;
    *' maigan '*)   echo 'cpu=80,mem=625G,gres/gpu=2' ;;
esac
EOF
cat > "$T/bin/squeue" <<'EOF'
#!/bin/bash
args="$*"
if [[ $args == *'%F %K %T'* ]]; then echo '900 0 RUNNING'; exit 0; fi          # our job 900, task 0 = image B
if [[ $args == *'-o %L'* ]]; then echo "${MOCK_TIME_LEFT:-09:00:00}"; exit 0; fi
if [[ $args == *'tres-alloc'* ]]; then
    q=$(sed -E 's/.*-q ([^ ]+).*/\1/' <<< "$args")
    if [[ $q == maigan ]]; then
        echo 'cpu=32,mem=200G,node=1,billing=32'                                # RUNNING (another user)
        [[ $args == *PENDING* ]] && echo 'cpu=16,mem=100G,node=1,billing=16'    # PENDING: counts too
    fi
    exit 0
fi
echo '(squeue display)'
EOF
cat > "$T/bin/sbatch" <<'EOF'
#!/bin/bash
n=$(($(cat "$MOCK_DIR/jobs" 2>/dev/null || echo 1000) + 1)); echo "$n" > "$MOCK_DIR/jobs"
printf '%s\n' "$*" >> "$MOCK_DIR/sbatch.log"
echo "$n"
EOF
printf '#!/bin/bash\necho "(mock $0 $*)"\n' > "$T/bin/showQos"
cp "$T/bin/showQos" "$T/bin/slurmInfo"
chmod +x "$T/bin/"*
export PATH="$T/bin:$PATH" MOCK_DIR=$T BLUE=$T/blue PY V3_PROJECT_FILE=$T/proj/project.qpproj V3_OUT=$T/out USER=tester
unset SLURM_JOB_ID CUDA_VISIBLE_DEVICES SLURM_JOB_GPUS V3_SITE
mkdir -p "$BLUE/work/submit/jobs"
printf 'img B\n' > "$BLUE/work/submit/jobs/900.images"          # registered by an earlier submission
SUB=$V3DIR/hpg/hpg_submit.sh

echo '== 1. images must resolve exactly'
printf 'img A\nimg\n' > "$T/bad.txt"
bash "$SUB" c0 --images-file "$T/bad.txt" --gpu-type b200 > "$T/o1" 2>&1
rc=$?
check 'refused (exit 2), names the token' test "$rc" = 2 -a -n "$(grep -l "'img' matches no image" "$T/o1")"

echo '== 2. c0: queued, done and blocked images are refused; --drop-unready submits the rest'
printf '\n# pilot\nimg A\n\n  img B\nid:3\nimg D\n\n' > "$T/c0.txt"
bash "$SUB" c0 --images-file "$T/c0.txt" --gpu-type b200 > "$T/o2" 2>&1
rc=$?
check 'refused (exit 2)' test "$rc" = 2 -a ! -e "$T/sbatch.log"
check 'B: already in job 900' grep -q 'img B: already in job 900_0 RUNNING' "$T/o2"
check 'C: nothing to run' grep -q 'img C: nothing to run' "$T/o2"
check 'D: step 0 blocked (preflight missing)' grep -q 'img D: step 0 blocked: prerequisite not done: preflight' "$T/o2"
bash "$SUB" c0 --images-file "$T/c0.txt" --gpu-type b200 --drop-unready > "$T/o2b" 2>&1
rc=$?
check 'submitted (exit 0), one sbatch' test "$rc" = 0 -a "$(wc -l < "$T/sbatch.log")" = 1
line=$(cat "$T/sbatch.log" 2>/dev/null)
check '--time from steps.tsv (480 + 45 min), GPU pinned, array of 1' \
    test -n "$(grep -e '--time=08:45:00' <<< "$line" | grep -e '--gres=gpu:b200:1' | grep -e '--array=0-0%')"
part=$(sed -E 's/.*IMAGES_FILE=([^,]+),.*/\1/' <<< "$line")
check 'the job reads exactly "img A"' test "$(cat "$part")" = 'img A'
check 'the job is registered for the next check' test "$(cat "$BLUE/work/submit/jobs/1001.images" 2>/dev/null)" = 'img A'

echo '== 3. c2: pending jobs count, maigan-b is capped, the images are split'
rm -f "$T/sbatch.log"
printf 'img C\nimg E\n' > "$T/c2.txt"
bash "$SUB" c2 --images-file "$T/c2.txt" > "$T/o3" 2>&1
rc=$?
check 'submitted (exit 0)' test "$rc" = 0
check 'maigan: RUNNING + PENDING counted (32+16 cpu, 200+100 GB)' grep -q 'maigan limits cpu=80 mem=625G gpu=2; running+pending cpu=48 mem=300G' "$T/o3"
check 'maigan-b capped at --maigan-b-max 2' grep -q 'maigan-b limits .* -> 2 slot(s)' "$T/o3"
check 'two arrays: clive and maigan, one image each, --time 1440 + 30 min' \
    test "$(grep -c -e '--time=1-00:30:00' "$T/sbatch.log")" = 2 -a -n "$(grep -e '--qos=clive' "$T/sbatch.log")" -a -n "$(grep -e '--qos=maigan ' "$T/sbatch.log")"

echo '== 4. hpg_env.sh helpers'
# shellcheck source=../hpg/hpg_env.sh
source "$V3DIR/hpg/hpg_env.sh"
printf '\n# header\nX one\n\n   \nY two\n' > "$T/imgs"
check 'task 0 / 1 skip blank and comment lines' test "$(IMAGES_FILE=$T/imgs SLURM_ARRAY_TASK_ID=0 hpg_image_for_task)" = 'X one' -a \
    "$(IMAGES_FILE=$T/imgs SLURM_ARRAY_TASK_ID=1 hpg_image_for_task)" = 'Y two' -a -z "$(IMAGES_FILE=$T/imgs SLURM_ARRAY_TASK_ID=2 hpg_image_for_task)"
check 'durations: 1-02:03:04, 05:00, 42; UNLIMITED fails' \
    test "$(hpg_duration_s 1-02:03:04)" = 93784 -a "$(hpg_duration_s 05:00)" = 300 -a "$(hpg_duration_s 42)" = 42 -a -z "$(hpg_duration_s UNLIMITED)"
(export SLURM_JOB_ID=77 MOCK_TIME_LEFT=08:00:00; hpg_check_budget --steps 0 --margin-min 45) > "$T/o4" 2>&1
rc=$?
check 'a job with 8 h left refuses step 0 (needs 480 + 45 min)' test "$rc" = 2 -a -n "$(grep -l 'time limit too short' "$T/o4")"
dl=$(export SLURM_JOB_ID=77 MOCK_TIME_LEFT=09:00:00; hpg_check_budget --steps 0 --margin-min 45 2> /dev/null; echo "${V3_DEADLINE_EPOCH:-}")
check 'with 9 h left it passes and exports the job deadline' test -n "$dl" -a "$((dl - $(date +%s)))" -gt 32000
check 'V3_SITE=hpg and the shared lock dir on /blue' test "$V3_SITE" = hpg -a "$V3_SHARED_LOCKDIR" = "$BLUE/locks"

((FAILS == 0)) && echo 'ALL HPG TESTS PASSED' || echo "$FAILS HPG TEST(S) FAILED"
exit $((FAILS > 0))
