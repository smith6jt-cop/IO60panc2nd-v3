#!/usr/bin/env bash
# hpg_submit.sh -- check free capacity, choose accounts, submit one of the v3 arrays.
#
#   hpg_submit.sh c0|c2|j --images-file FILE [options]
#     c0  step 0 (GPU)  -> c0_gpu.sbatch        c2  step 2 (CPU) -> c2_cpu.sbatch        j  passes R,J -> j_cpu.sbatch
#   --account auto|clive|maigan|maigan-b   default auto: clive first (dedicated: 40 CPU / 312 GB / 1 GPU),
#                                          then maigan (80 CPU / 625 GB / 2 GPU, SHARED with the group),
#                                          then maigan-b burst (CPU arrays only)
#   --max-parallel N      cap on concurrently running tasks per account (default: what the capacity allows)
#   --maigan-share F      never hold more than this fraction of maigan's (or maigan-b's) CPU/memory (default 0.5)
#   --maigan-gpus N       GPUs we may take on maigan when they are free (default 1 of 2: never all of it)
#   --maigan-b-max N      at most N concurrent burst tasks (default 2)
#   --gpu-type T          pin every step-0 task to ONE GPU type (e.g. b200); --partition P to match
#   --stage               c0: copy the qptiff to $TMPDIR first (STAGE_TO_TMPDIR=1)
#   --steps S             c2: STEPS for c2_cpu.sbatch (default 2)
#   --drop-unready        submit the ready images and leave out the others (default: refuse the submission)
#   --dry-run             print the decision and the sbatch lines, submit nothing
#
# Before anything is submitted:
#   - every line of the images file must name exactly one project image (exact name, tag or id:N; blank
#     lines and # comments are ignored); the normalized list (one image name per line) is what the jobs read;
#   - an image already in one of YOUR queued or running v3 jobs is refused (jobs submitted here are
#     registered in $BLUE/work/submit/jobs/<jobid>.images); the mkdir image lock on /blue is the runtime
#     backstop for jobs submitted by hand;
#   - an image whose plan would not run anything (not merged, custody broken, prerequisite missing,
#     nothing to do) is refused: sync and merge provenance first (run_v3.sh --merge-provenance);
#   - --time is the sum of the job's step timeouts in steps.tsv plus a margin (v3ctl.py budget).
# Capacity per QoS = GrpTRES (sacctmgr) minus the RUNNING and PENDING jobs of the whole group (squeue).
# Every decision is appended to $BLUE/logs/submissions.log with the capacity it was based on.
# Run on a login node or a CPU-only session: never from an interactive GPU session, whose idle GPU
# would itself trip HiPerGator's one-hour idle rule.
set -euo pipefail
# shellcheck source=hpg_env.sh
source "$(dirname "$0")/hpg_env.sh"

KIND=${1:-}
shift || true
IMAGES_FILE="" ACCOUNT=auto MAX_PAR=0 MAIGAN_SHARE=0.5 MAIGAN_GPUS=1 MAIGAN_B_MAX=2 GPU_TYPE="" PARTITION=""
STAGE=0 STEPS=2 DRY=0 DROP=0
while (($#)); do
    case $1 in
        --images-file) IMAGES_FILE=$2; shift ;;
        --account) ACCOUNT=$2; shift ;;
        --max-parallel) MAX_PAR=$2; shift ;;
        --maigan-share) MAIGAN_SHARE=$2; shift ;;
        --maigan-gpus) MAIGAN_GPUS=$2; shift ;;
        --maigan-b-max) MAIGAN_B_MAX=$2; shift ;;
        --gpu-type) GPU_TYPE=$2; shift ;;
        --partition) PARTITION=$2; shift ;;
        --stage) STAGE=1 ;;
        --steps) STEPS=$2; shift ;;
        --drop-unready) DROP=1 ;;
        --dry-run) DRY=1 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
    shift
done
case $KIND in
    c0) SCRIPT=c0_gpu.sbatch; NEED_CPU=16; NEED_MEM_GB=128; NEED_GPU=1; SEL=(--steps 0); MARGIN=45 ;;
    c2) SCRIPT=c2_cpu.sbatch; NEED_CPU=32; NEED_MEM_GB=220; NEED_GPU=0; SEL=(--steps "$STEPS"); MARGIN=30 ;;
    j)  SCRIPT=j_cpu.sbatch;  NEED_CPU=24; NEED_MEM_GB=180; NEED_GPU=0; SEL=(--passes R,J); MARGIN=30 ;;
    *) echo 'usage: hpg_submit.sh c0|c2|j --images-file FILE [options]' >&2; exit 2 ;;
esac
[[ -s $IMAGES_FILE ]] || { echo "--images-file missing or empty" >&2; exit 2; }
IMAGES_FILE=$(realpath "$IMAGES_FILE")
if [[ -n ${SLURM_JOB_ID:-} && -n ${CUDA_VISIBLE_DEVICES:-}${SLURM_JOB_GPUS:-} ]]; then
    echo 'refusing: this shell is inside a GPU job; submit from a login node or a CPU-only session' >&2
    exit 2
fi
if [[ $KIND == c0 && -z $GPU_TYPE ]]; then
    echo 'c0 needs --gpu-type: every step-0 task must run on the ONE pinned GPU type' >&2
    exit 2
fi
# (V3_PROJECT_FILE / V3_OUT: tests point the checks at a fixture project; production uses the defaults)
CTL=("$PY" "$V3DIR/v3ctl.py" --project "${V3_PROJECT_FILE:-$PROJ/project.qpproj}" --out "${V3_OUT:-$PROJ/results/v3}")
SUBLOG=$BLUE/logs/submissions.log
REG=$BLUE/work/submit/jobs
mkdir -p "$BLUE/logs" "$BLUE/work/submit" "$REG"
STAMP=$(date +%Y%m%dT%H%M%S)

# ---- 0. the images: exact resolution, not queued already, ready to run ----------------------------
if ! names=$("${CTL[@]}" images --images-file "$IMAGES_FILE"); then
    echo "refusing: $IMAGES_FILE names an image that does not resolve exactly (see above)" >&2
    exit 2
fi
mapfile -t IMAGES <<< "$names"
[[ -n ${IMAGES[0]:-} ]] || { echo "no images in $IMAGES_FILE" >&2; exit 2; }

declare -A BUSY=()   # image name -> the queued/running job of ours that has it
while read -r base idx state; do
    [[ $idx =~ ^[0-9]+$ && -f $REG/$base.images ]] || continue
    img=$(sed -n "$((idx + 1))p" "$REG/$base.images")
    [[ -n $img ]] && BUSY[$img]="job ${base}_$idx $state"
done < <(squeue -u "$USER" -h -r -o '%F %K %T' 2> /dev/null || true)

plan=$("${CTL[@]}" plan "${SEL[@]}" --images-file "$IMAGES_FILE" --format tsv 2> /dev/null) \
    || { echo 'refusing: v3ctl.py plan failed for these images (run it by hand to see why)' >&2; exit 2; }
declare -A WHY=() NRUN=() NAMEOF=()
while IFS=$'\t' read -r kind tag f3 f4 f5 _ f7 _ _ _ f11; do
    case $kind in
        IMAGE) NAMEOF[$tag]=$f3; NRUN[$f3]=0
               [[ $f5 == ok || $f5 == fresh ]] || WHY[$f3]="image $f5" ;;
        STEP)  name=${NAMEOF[$tag]}
               [[ $f4 == run ]] && NRUN[$name]=$((NRUN[$name] + 1))
               [[ $f4 == blocked && -z ${WHY[$name]:-} ]] && WHY[$name]="step $f3 blocked: ${f11:0:120}" ;;
    esac
done <<< "$plan"
READY=() REFUSED=()
for img in "${IMAGES[@]}"; do
    why=${BUSY[$img]:+already in ${BUSY[$img]}}
    [[ -z $why ]] && why=${WHY[$img]:-}
    [[ -z $why && ${NRUN[$img]:-0} -eq 0 ]] && why='nothing to run (all steps done)'
    if [[ -n $why ]]; then REFUSED+=("$img: $why"); else READY+=("$img"); fi
done
if ((${#REFUSED[@]})); then
    printf 'not submittable:\n' >&2
    printf '  %s\n' "${REFUSED[@]}" >&2
    printf '%s %s %s refused=%s\n' "$(date -Is)" "$KIND" "$IMAGES_FILE" "${#REFUSED[@]}" >> "$SUBLOG"
    if ((!DROP)); then echo 'refusing the submission (fix these, or pass --drop-unready)' >&2; exit 2; fi
fi
((${#READY[@]})) || { echo 'no image left to submit' >&2; exit 3; }
TIME=$("${CTL[@]}" budget "${SEL[@]}" --margin-min "$MARGIN" --format slurm)

# ---- 1. show capacity (UF tools, for the human) --------------------------------------------------
echo "== $(date -Is) capacity before submitting $KIND for ${#READY[@]} image(s), --time=$TIME"
for q in maigan clive; do
    echo "-- showQos $q"; showQos "$q" 2>&1 || echo '(showQos unavailable)'
done
echo '-- slurmInfo'; slurmInfo 2>&1 | head -40 || echo '(slurmInfo unavailable)'
for a in maigan clive; do
    echo "-- squeue -A $a"; squeue -A "$a" -o '%.10i %.9P %.20j %.8u %.2t %.10M %.4C %.8m %.15b' 2>&1 | head -30 || true
done
# TODO(HiPerGator): paste one sample of `showQos maigan` and `slurmInfo` output into the repo and parse
# the group's free CPU/memory/GPU from it; until then the decision below uses standard SLURM only.

# ---- 2. free capacity per QoS: GrpTRES limit (sacctmgr) minus running AND pending jobs (squeue) ---
tres_field() {   # tres_field "<cpu=80,mem=625G,gres/gpu=2>" cpu|mem|gres/gpu -> integer (mem in GB); empty if absent
    local key=$2 v
    v=$(tr ',' '\n' <<< "$1" | awk -F= -v k="$key" '$1 == k {print $2; exit}')
    v=${v//[[:space:]]/}                                  # squeue -O pads its columns
    [[ -z $v ]] && return 0
    case $v in
        *T) echo $((${v%T} * 1024)) ;;
        *G) echo "${v%G}" ;;
        *M) echo $((${v%M} / 1024)) ;;
        *) if [[ $key == mem ]]; then echo $((v / 1024)); else echo "$v"; fi ;;   # sacctmgr mem without unit = MB
    esac
}
limit_of() {   # limit_of <qos> <cpu|mem|gres/gpu>
    local grp
    grp=$(sacctmgr -nP show qos "$1" format=GrpTRES 2> /dev/null | head -1)
    tres_field "$grp" "$2"
}
used_of() {    # used_of <qos> <cpu|mem|gres/gpu>: the whole group's RUNNING + PENDING jobs under that QoS
    # (-r lists every pending array task on its own line; a pending job holds nothing yet, but it will
    # start before ours, so counting it keeps us from over-committing the shared allocation)
    local total=0 line v
    while IFS= read -r line; do
        v=$(tres_field "$line" "$2")
        total=$((total + ${v:-0}))
    done < <(squeue -h -r -q "$1" -t RUNNING,PENDING -O tres-alloc:200 2> /dev/null)
    echo "$total"
}
slots_for() {   # slots_for <qos> -> how many of our tasks may run at once under it now
    local q=$1 lc lm lg uc um ug fc fm fg n
    lc=$(limit_of "$q" cpu); lm=$(limit_of "$q" mem); lg=$(limit_of "$q" gres/gpu)
    uc=$(used_of "$q" cpu); um=$(used_of "$q" mem); ug=$(used_of "$q" gres/gpu)
    [[ -n $lc && -n $lm ]] || { echo 0; return; }                  # unknown limit: do not guess
    fc=$((lc - uc)); fm=$((lm - um)); fg=$((${lg:-0} - ug))
    if [[ $q == maigan || $q == maigan-b ]]; then   # shared: never hold more than MAIGAN_SHARE of the group's allocation
        fc=$(awk -v l="$lc" -v f="$fc" -v s="$MAIGAN_SHARE" 'BEGIN {m = int(l * s); print (f < m ? f : m)}')
        fm=$(awk -v l="$lm" -v f="$fm" -v s="$MAIGAN_SHARE" 'BEGIN {m = int(l * s); print (f < m ? f : m)}')
        ((fg > MAIGAN_GPUS)) && fg=$MAIGAN_GPUS
    fi
    n=$((fc / NEED_CPU))
    ((fm / NEED_MEM_GB < n)) && n=$((fm / NEED_MEM_GB))
    if ((NEED_GPU)); then ((fg < n)) && n=$fg; fi
    ((n < 0)) && n=0
    [[ $q == maigan-b ]] && ((n > MAIGAN_B_MAX)) && n=$MAIGAN_B_MAX
    ((MAX_PAR > 0 && n > MAX_PAR)) && n=$MAX_PAR
    printf '%s limits cpu=%s mem=%sG gpu=%s; running+pending cpu=%s mem=%sG gpu=%s -> %s slot(s)\n' \
        "$q" "$lc" "$lm" "${lg:-0}" "$uc" "$um" "$ug" "$n" >&2
    echo "$n"
}

if [[ $ACCOUNT == auto ]]; then
    candidates=(clive maigan)
    ((NEED_GPU)) || candidates+=(maigan-b)
else
    candidates=("$ACCOUNT")
fi
if ((NEED_GPU)) && [[ " ${candidates[*]} " == *' maigan-b '* ]]; then
    echo 'maigan-b (burst) is for CPU arrays only' >&2; exit 2
fi
declare -A SLOTS=()
total=0
for q in "${candidates[@]}"; do
    SLOTS[$q]=$(slots_for "$q")
    total=$((total + SLOTS[$q]))
done
if ((total == 0)); then
    echo 'no free capacity on any candidate account now (or its limits could not be read); nothing submitted' >&2
    printf '%s %s %s images=%s decision=none\n' "$(date -Is)" "$KIND" "$IMAGES_FILE" "${#READY[@]}" >> "$SUBLOG"
    exit 3
fi

# ---- 3. split the images across the accounts in proportion to their slots, one array each --------
declare -A LIST=() COUNT=()
for img in "${READY[@]}"; do
    # weighted round robin: the account with the fewest images per slot so far takes the next one
    best="" best_load=""
    for q in "${candidates[@]}"; do
        ((SLOTS[$q] > 0)) || continue
        load=$(awk -v n="${COUNT[$q]:-0}" -v s="${SLOTS[$q]}" 'BEGIN {printf "%.6f", n / s}')
        if [[ -z $best ]] || awk -v a="$load" -v b="$best_load" 'BEGIN {exit !(a < b)}'; then best=$q best_load=$load; fi
    done
    LIST[$best]+="$img"$'\n'
    COUNT[$best]=$((${COUNT[$best]:-0} + 1))
done

for q in "${candidates[@]}"; do
    [[ -n ${LIST[$q]:-} ]] || continue
    part=$BLUE/work/submit/${KIND}_${q}_$STAMP.txt
    printf '%s' "${LIST[$q]}" > "$part"                  # exactly one image name per line, no blank lines
    n=${COUNT[$q]}
    [[ $(hpg_image_lines "$part" | wc -l) == "$n" ]] || { echo "internal error: $part does not hold $n images" >&2; exit 1; }
    # TODO(HiPerGator): confirm the account/QoS pairs (burst QoS maigan-b bills to account maigan)
    acct=${q%-b}
    args=(--account="$acct" --qos="$q" --array="0-$((n - 1))%${SLOTS[$q]}" --time="$TIME"
          --export="ALL,IMAGES_FILE=$part,GPU_TYPE=$GPU_TYPE,STAGE_TO_TMPDIR=$STAGE,STEPS=$STEPS")
    if ((NEED_GPU)); then args+=(--gres="gpu:${GPU_TYPE}:1"); fi
    if [[ -n $PARTITION ]]; then args+=(--partition="$PARTITION"); fi
    echo "sbatch ${args[*]} $V3DIR/hpg/$SCRIPT   # $n image(s)"
    if ((DRY)); then
        printf '%s %s qos=%s slots=%s images=%s file=%s time=%s dry-run\n' "$(date -Is)" "$KIND" "$q" "${SLOTS[$q]}" "$n" "$part" "$TIME" >> "$SUBLOG"
        continue
    fi
    jid=$(sbatch --parsable "${args[@]}" "$V3DIR/hpg/$SCRIPT")
    jid=${jid%%;*}
    cp "$part" "$REG/$jid.images"                        # the next submission refuses these images while the job lives
    echo "submitted job $jid"
    printf '%s %s qos=%s slots=%s images=%s file=%s time=%s job=%s\n' "$(date -Is)" "$KIND" "$q" "${SLOTS[$q]}" "$n" "$part" "$TIME" "$jid" >> "$SUBLOG"
done
((DRY)) && echo '(dry run: nothing submitted)'
exit 0
