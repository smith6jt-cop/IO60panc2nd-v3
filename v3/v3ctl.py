#!/usr/bin/env python3
"""v3ctl.py -- the decisions behind run_v3.sh: plan, instantiate, validate, convert, record.

run_v3.sh owns processes and lanes; everything that reads or writes JSON lives here.
Standard library only, except `post`, which imports duckdb lazily (venv-v3), so planning and
validation also work with a bare python3 (e.g. on HiPerGator).

Subcommands (run_v3.sh calls them; they are also usable by hand):
  steps             print the resolved step table
  images            resolve --images / --images-file exactly and print the image names
  plan              what would run / be skipped / is blocked, per image (tsv for the driver, or a table)
  budget            sum of the step timeouts of a selection (+ margin): the SLURM --time a job needs
  prepare           instantiate one step for one image: worker project + script copy + footer
  finish            validate the sentinel of an attempt, copy the script, check drift, log the row
  post              run a step's post hook (TSV -> Parquet + verification)
  diskcheck         the free-space guard
  accept-external   record that the current data.qpdata (changed outside run_v3) is accepted
  merge-provenance  adopt the records another site (e.g. HiPerGator) wrote for an image
  show-record       print the current record of a step (or one field of it)
  lock              acquire/release the shared (mkdir) image lock used across HiPerGator nodes
  runlog            append a run_log.csv row for an outcome without an attempt (blocked, locked...)
  report            summary of one run (statuses, drift, worker-vs-project entry differences)

Provenance layout, per image tag and per SITE (V3_SITE; default the short host name; 'hpg' on
HiPerGator), under <out>/provenance/<tag>/<site>/:
  <step>.json          written by the footer inside QuPath: the success sentinel + provenance
  <step>.finish.json   written by `finish`: the driver validated that sentinel (a sentinel without it
                       -- e.g. from an orphaned QuPath -- does not count as done)
  <step>.post.json     written by `post` (conversion + verification)
  <step>.groovy        the exact script that ran (step body + overrides + footer)
  <step>.drift.json    only if the step changed the image server metadata or the project entry
  external.<seq>.json  accepted modifications made outside run_v3 (e.g. gate-1 review in the GUI)
  history/             superseded records of earlier attempts
  merged/<site2>/      byte copies of the records site2 wrote, adopted by `merge-provenance`
  .running.<step>.json marker while an attempt is in flight
A site writes only into its own directory, so syncing provenance/ between the workstation and
HiPerGator never overwrites a record. Another site's records count only after merge-provenance
mirrored them into merged/<site2>/; until then the image is blocked.

Every record carries a per-image sequence number `seq` (unique across sites) and the size/mtime of
data.qpdata before and after it. Custody = the records, sorted by seq, form ONE chain (each record
starts from the file the previous one left) that ends at the current data.qpdata; two records with
the same seq (the image was processed on two sites from the same state) stop the image.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import errno
import fcntl
import filecmp
import getpass
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOOTER_TEMPLATE = HERE / 'footer_template.groovy'
CFG_PLACEHOLDER = '/*__V3_CFG__*/ null'
SAVE_POLICIES = {'self-saves', 'footer-saves', 'never'}
LANES = {'gpu', 'cpu', 'light'}
ON_EXTERNAL = {'stale', 'keep'}
POST_HOOKS = {'cells_parquet', 'context_parquet'}
POSTCONDITIONS = {'preflight', 'v3_result', 'detections_present', 'regions_detected', 'islets_present',
                  'islet_bands', 'cells_measured', 'signed_distance', 'parenting'}
# Custody: size must match exactly, mtime within 1 s -- Globus/rsync preserve whole seconds only,
# so a .qpdata that went to HiPerGator and back still matches its record.
MTIME_TOL_NS = 1_000_000_000
KEEP_TSV_DEFAULT = r'0080870_6533|0081366_HDL112'   # same-slide Xenium images: sameslide/ reads the TSV
LAST_DEFAULT = ('6457_',)                            # plan: 6457 is re-run last
# measurement-table columns that hold text; every other column is numeric (DOUBLE)
STRING_COLUMNS = {'Image', 'Object ID', 'Object type', 'Name', 'Class', 'Classification', 'Parent',
                  'Parent ID', 'ROI', 'TMA core'}
# cells_parquet: every row must carry these (non-null, non-NaN) -- step 2 measured every detection
CELLS_REQUIRED_COMPLETE = ('Cell: DAPI: Mean',)
# A cell outline smaller than this covers no pixel centre, so step 2 cannot measure it (NaN; 6533: 2 of 1,072,057).
# Such debris is allowed and counted here exactly as the footer's cells_measured check (debris_max_um2) allows it;
# any larger cell without a required value still fails. The geometry QC downstream drops < 20 um^2 anyway.
DEBRIS_MAX_UM2 = 5.0
DEBRIS_AREA_COLUMN = 'Cell: Area µm^2'
RUN_LOG_FIELDS = ['time_start', 'time_end', 'run_id', 'host', 'image', 'tag', 'step', 'action', 'status',
                  'seconds', 'exit_code', 'attempt_id', 'seq', 'lane', 'slot', 'heap_gb', 'peak_rss_kb',
                  'detections', 'qpdata_bytes', 'script', 'script_sha256', 'message', 'site']
PREPARE_KEYS = ('ATTEMPT', 'SEQ', 'HEAP', 'WSCRIPT', 'WPROJECT', 'PROVENANCE', 'SCRIPT_SRC', 'QARGS_FILE')
SITE_RE = r'[A-Za-z0-9][A-Za-z0-9._-]*'
DEFAULT_LOCKDIR = '/scratch/io60/locks'
SLURM_ENDED = {'COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT', 'NODE_FAIL', 'PREEMPTED', 'OUT_OF_MEMORY',
               'BOOT_FAIL', 'DEADLINE', 'REVOKED', 'SPECIAL_EXIT'}


class PlanError(Exception):
    pass


# ----------------------------------------------------------------------------------------------
# small helpers
def now_iso() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec='milliseconds')


def tag_of(name: str) -> str:
    # identical to export_annotations.sh, so every v3 output of one image shares its tag
    return re.sub(r'[^A-Za-z0-9._-]+', '_', name)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def stat_of(path: Path) -> dict:
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return {'exists': False, 'path': str(path)}
    return {'exists': True, 'path': str(path), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}


def same_stat(a: dict | None, b: dict | None) -> bool:
    if not a or not b:
        return False
    if not a.get('exists') and not b.get('exists'):
        return True
    if not (a.get('exists') and b.get('exists')):
        return False
    return a.get('size') == b.get('size') and abs(int(a['mtime_ns']) - int(b['mtime_ns'])) <= MTIME_TOL_NS


def fmt_stat(s: dict | None) -> str:
    if not s:
        return '<none>'
    if not s.get('exists'):
        return 'absent'
    t = dt.datetime.fromtimestamp(int(s['mtime_ns']) / 1e9).isoformat(timespec='seconds')
    return f"size={s.get('size')} mtime={t}"


def bkp_path(qpdata: Path) -> Path:
    return qpdata.with_name(qpdata.name + '.bkp')


def bkp_note(qpdata: Path) -> str:
    """QuPath's saveImageData() renames data.qpdata to data.qpdata.bkp, writes the new file and deletes
    the .bkp only after a complete write; readImageData() silently loads the .bkp when data.qpdata
    cannot be read. A leftover .bkp therefore means a save was interrupted or failed."""
    b = stat_of(bkp_path(qpdata))
    if not b['exists']:
        return ''
    return (f" {b['path']} exists ({fmt_stat(b)}): a QuPath save was interrupted or failed (QuPath moves "
            f"data.qpdata there before writing and deletes it only after a complete write; readImageData() "
            f"silently loads it if data.qpdata is unreadable). Check which file holds the image before re-running.")


def read_json(path: Path):
    try:
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write('\n')
    os.replace(tmp, path)


def groovy_literal(v) -> str:
    """A Groovy literal for a JSON-like value (single-quoted strings, so '$' is inert)."""
    if v is None:
        return 'null'
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, int):
        return f'{v}L' if abs(v) > 2 ** 31 - 1 else str(v)
    if isinstance(v, float):
        return repr(v) + 'd'
    if isinstance(v, (str, Path)):
        out = []
        for ch in str(v):
            if ch == '\\':
                out.append('\\\\')
            elif ch == "'":
                out.append("\\'")
            elif ch == '\n':
                out.append('\\n')
            elif ch == '\r':
                out.append('\\r')
            elif ch == '\t':
                out.append('\\t')
            elif ord(ch) < 32 or ord(ch) > 126:     # e.g. 'µ': \uXXXX, surrogate pairs for astral chars
                units = ch.encode('utf-16-be')
                out += ['\\u%04X' % int.from_bytes(units[i:i + 2], 'big') for i in range(0, len(units), 2)]
            else:
                out.append(ch)
        return "'" + ''.join(out) + "'"
    if isinstance(v, (list, tuple)):
        return '[' + ', '.join(groovy_literal(x) for x in v) + ']'
    if isinstance(v, dict):
        if not v:
            return '[:]'
        return '[' + ', '.join(f'{groovy_literal(str(k))}: {groovy_literal(x)}' for k, x in v.items()) + ']'
    raise TypeError(f'no Groovy literal for {type(v)}')


def resolve_site(explicit: str | None) -> str:
    """Which site this process writes provenance for. V3_SITE (the HiPerGator bundle sets 'hpg', so
    every compute node shares one site) or the short host name."""
    s = (explicit or os.environ.get('V3_SITE') or '').strip()
    if not s:
        if os.path.isdir('/blue'):   # HiPerGator: node host names vary, so they cannot name the site
            raise SystemExit('V3_SITE is not set on what looks like HiPerGator (/blue exists): source '
                             'scripts/v3/hpg/hpg_env.sh first (it sets V3_SITE=hpg)')
        s = socket.gethostname().split('.')[0]
    if not re.fullmatch(SITE_RE, s) or s in ('history', 'merged'):
        raise SystemExit(f'bad site name {s!r} (V3_SITE / --site)')
    return s


def git_info(scripts_repo: Path, source: Path | None) -> dict:
    """Commit of the scripts repository and whether the working tree (or the step script) differs."""
    def run(*args, timeout=60):
        return subprocess.run(['git', '-C', str(scripts_repo), *args], capture_output=True, text=True,
                              timeout=timeout)
    info = {'repo': str(scripts_repo), 'commit': None, 'dirty': None, 'source_clean': None}
    try:
        head = run('rev-parse', 'HEAD', timeout=20)
        if head.returncode != 0:
            info['error'] = (head.stderr.strip() or 'git rev-parse failed')[:200]
            return info
        info['commit'] = head.stdout.strip()
        st = run("--no-optional-locks", "status", "--porcelain")   # tracked changes and untracked files
        if st.returncode == 0:
            lines = [ln for ln in st.stdout.splitlines() if ln.strip()]
            info['dirty'] = bool(lines)
            info['dirty_files'] = lines[:20]
        if source is not None:
            try:
                rel = Path(source).resolve().relative_to(Path(scripts_repo).resolve())
            except ValueError:
                info['source_clean'] = None                  # a script outside the repo (--script)
            else:
                tracked = run('ls-files', '--error-unmatch', '--', str(rel), timeout=20).returncode == 0
                changed = run('--no-optional-locks', 'status', '--porcelain', '--', str(rel), timeout=20).stdout.strip()
                info['source_clean'] = tracked and not changed
    except (OSError, subprocess.TimeoutExpired) as e:
        info['error'] = str(e)[:200]
    return info


# ----------------------------------------------------------------------------------------------
# step table
@dataclass
class Step:
    id: str
    passes: list
    candidates: list
    args: list
    save: str
    lane: str
    heap_gb: int
    threads: int
    timeout_min: int
    after: list          # [(step id, optional)]
    needs: list
    on_external: str     # stale | keep
    overrides: dict
    postcondition: dict  # {'name': str, 'params': {key: [values]}}
    post: str | None
    description: str


def parse_postcondition(text: str) -> dict:
    m = re.fullmatch(r'\s*([A-Za-z_]\w*)\s*(?:\((.*)\))?\s*', text)
    if not m:
        raise ValueError(f'bad postcondition {text!r}')
    params = {}
    for part in (m.group(2) or '').split(';'):
        if part.strip():
            key, _, val = part.partition('=')
            params[key.strip()] = [x.strip() for x in val.split('|')]
    return {'name': m.group(1), 'params': params}


STEP_COLUMNS = ['step', 'passes', 'script', 'args', 'save', 'lane', 'heap_gb', 'threads', 'timeout_min', 'after',
                'needs', 'on_external', 'overrides', 'postcondition', 'post', 'description']


def load_steps(path: Path) -> dict:
    lines = [ln for ln in path.read_text(encoding='utf-8').splitlines() if ln.strip() and not ln.startswith('#')]
    header = lines[0].split('\t')
    if header != STEP_COLUMNS:
        raise SystemExit(f'{path}: header must be {STEP_COLUMNS}, got {header}')
    steps: dict[str, Step] = {}
    field = lambda v: '' if v.strip() == '-' else v.strip()
    for ln in lines[1:]:
        cols = ln.split('\t')
        if len(cols) != len(header):
            raise SystemExit(f'{path}: row {cols[0]!r} has {len(cols)} columns, header has {len(header)}')
        r = {k: field(v) for k, v in zip(header, cols)}
        overrides = {}
        for item in filter(None, (x.strip() for x in r['overrides'].split(';'))):
            k, _, v = item.partition('=')
            overrides[k.strip()] = v.strip()
        after = [(d.rstrip('?'), d.endswith('?')) for d in filter(None, (x.strip() for x in r['after'].split(',')))]
        st = Step(id=r['step'], passes=[p.strip() for p in r['passes'].split(',') if p.strip()],
                  candidates=[c.strip() for c in r['script'].split('|') if c.strip()],
                  args=[a for a in r['args'].split('|')] if r['args'] else [],
                  save=r['save'], lane=r['lane'], heap_gb=int(r['heap_gb']), threads=int(r['threads']),
                  timeout_min=int(r['timeout_min']), after=after,
                  needs=[d.strip() for d in r['needs'].split(',') if d.strip()], on_external=r['on_external'],
                  overrides=overrides, postcondition=parse_postcondition(r['postcondition']),
                  post=r['post'] or None, description=r.get('description', ''))
        if st.id in steps:
            raise SystemExit(f'{path}: duplicate step {st.id!r}')
        if not re.fullmatch(r'[A-Za-z0-9_]+', st.id) or st.id == 'external':
            raise SystemExit(f'{path}: step id {st.id!r} must be [A-Za-z0-9_]+ and not "external"')
        if st.save not in SAVE_POLICIES or st.lane not in LANES:
            raise SystemExit(f'{path}: step {st.id}: bad save policy {st.save!r} or lane {st.lane!r}')
        if st.on_external not in ON_EXTERNAL:
            raise SystemExit(f'{path}: step {st.id}: on_external must be one of {sorted(ON_EXTERNAL)}')
        if st.postcondition['name'] not in POSTCONDITIONS:
            raise SystemExit(f'{path}: step {st.id}: unknown postcondition {st.postcondition["name"]!r}')
        if st.post and st.post not in POST_HOOKS:
            raise SystemExit(f'{path}: step {st.id}: unknown post hook {st.post!r}')
        steps[st.id] = st
    for st in steps.values():
        for dep in [d for d, _ in st.after] + st.needs:
            if dep not in steps:
                raise SystemExit(f'{path}: step {st.id} depends on unknown step {dep!r}')
    return steps


def parse_kv_list(items, what) -> dict:
    out = {}
    for it in items or []:
        k, sep, v = it.partition('=')
        if not sep:
            raise SystemExit(f'{what}: expected KEY=VALUE, got {it!r}')
        out[k.strip()] = v.strip()
    return out


def resolve_script(st: Step, repo: Path, script_set: str, script_overrides: dict) -> tuple[Path, str]:
    if st.id in script_overrides:
        p = Path(script_overrides[st.id]).resolve()
        if not p.is_file():
            raise PlanError(f'--script {st.id}: {p} does not exist')
        return p, 'override'
    cands = [repo / c for c in st.candidates]
    if len(cands) == 1:
        chosen, how = cands[0], 'fixed'
    elif script_set == 'orig':
        chosen, how = cands[-1], 'orig'
    elif script_set == 'v3':
        chosen, how = cands[0], 'v3'
    else:
        chosen = next((c for c in cands if c.is_file()), cands[-1])
        how = 'auto:v3' if chosen == cands[0] else 'auto:orig'
    if not chosen.is_file():
        raise PlanError(f'step {st.id}: script {chosen} does not exist' +
                        (' (--script-set v3 never falls back to the original)' if how == 'v3' else ''))
    return chosen, how


def resolve_all_scripts(steps: dict, repo: Path, script_set: str, so: dict) -> dict:
    """step -> {'path', 'how', 'sha256'} or {'error'}; for every step (done() looks at dependencies too)."""
    out = {}
    for sid, st in steps.items():
        try:
            p, how = resolve_script(st, repo, script_set, so)
            out[sid] = {'path': str(p), 'how': how, 'sha256': sha256_file(p)}
        except PlanError as e:
            out[sid] = {'error': str(e)}
    return out


OVERRIDE_DECL = r'^(?P<pre>[ \t]*(?:final[ \t]+)?[A-Za-z_][\w.]*(?:<[^>\n]*>)?(?:\[\])*[ \t]+){name}(?P<mid>[ \t]*=[ \t]*)(?P<val>.*?)(?P<post>[ \t]*(?://.*)?)$'


def apply_overrides(text: str, overrides: dict) -> tuple[str, dict]:
    """Rewrite '<type> NAME = value' declarations; each must match exactly once."""
    applied = {}
    for name, value in overrides.items():
        pat = re.compile(OVERRIDE_DECL.format(name=re.escape(name)), re.M)
        hits = pat.findall(text)
        if len(hits) != 1:
            raise PlanError(f'override {name}={value}: the declaration matched {len(hits)} times (need exactly 1)')
        old = pat.search(text).group('val')
        text = pat.sub(lambda m: m.group('pre') + name + m.group('mid') + value + m.group('post'), text, count=1)
        applied[name] = {'from': old, 'to': value}
    return text, applied


def upstream_closure(steps: dict, sid: str) -> set:
    seen, todo = set(), [d for d, _ in steps[sid].after]
    while todo:
        d = todo.pop()
        if d not in seen:
            seen.add(d)
            todo += [x for x, _ in steps[d].after]
    return seen


def downstream_closure(steps: dict, roots) -> set:
    """roots plus every step that has one of them in its 'after' closure."""
    roots = set(roots)
    return roots | {s for s in steps if upstream_closure(steps, s) & roots}


# ----------------------------------------------------------------------------------------------
# project / images
@dataclass
class Image:
    name: str
    entry_id: int
    tag: str
    qpdata: Path
    entry: dict


def load_images(project: Path) -> list:
    proj = read_json(project)
    if proj is None:
        raise SystemExit(f'project not found: {project}')
    base = project.parent
    images = [Image(e['imageName'], e['entryID'], tag_of(e['imageName']),
                    base / 'data' / str(e['entryID']) / 'data.qpdata', e) for e in proj['images']]
    seen = {}
    for im in images:
        if im.tag in seen:
            raise SystemExit(f'{project}: images {seen[im.tag]!r} and {im.name!r} share the tag {im.tag!r}; '
                             f'outputs are keyed by tag, so rename one of them')
        seen[im.tag] = im.name
    return images


def image_tokens(spec: str | None = None, images_file: str | None = None) -> list:
    """--images: comma-separated tokens. --images-file: one token per line; blank lines and '#' comments
    are ignored (so a file's line count never disagrees with its images)."""
    toks = []
    if spec:
        toks += [t.strip() for t in spec.split(',') if t.strip()]
    if images_file:
        for ln in Path(images_file).read_text(encoding='utf-8').splitlines():
            ln = ln.strip()
            if ln and not ln.startswith('#'):
                toks.append(ln)
    return toks


def resolve_images(images: list, tokens: list, what: str = '--images') -> list:
    """Exact matching only: the image name or its tag, or 'id:<entry id>'. Every token must match
    exactly one image, and no image may be named twice (substring matching is gone: '6533' used to
    select every image whose name contains it)."""
    chosen = []
    for tok in tokens:
        if tok.startswith('id:'):
            try:
                eid = int(tok[3:])
            except ValueError:
                raise SystemExit(f"{what}: {tok!r}: 'id:' must be followed by an entry id number")
            hits = [im for im in images if int(im.entry_id) == eid]
        else:
            hits = [im for im in images if tok in (im.name, im.tag)]
        if len(hits) != 1:
            if hits:
                raise SystemExit(f'{what}: {tok!r} matches {len(hits)} images: {[h.name for h in hits]}')
            near = [f'{im.tag} (id:{im.entry_id})' for im in images if tok.lower() in im.name.lower()][:6]
            raise SystemExit(f'{what}: {tok!r} matches no image; give the exact image name, its tag, or '
                             f'id:<entry id>' + (f'. Names containing it: {", ".join(near)}' if near else ''))
        if hits[0] in chosen:
            raise SystemExit(f'{what}: {tok!r} names {hits[0].name!r} a second time')
        chosen.append(hits[0])
    return chosen


def select_images(images: list, spec: str | None = None, images_file: str | None = None,
                  what: str = '--images') -> list:
    toks = image_tokens(spec, images_file)
    if not toks:
        if images_file:
            raise SystemExit(f'--images-file {images_file}: no image names in it')
        first = [im for im in images if not im.name.startswith(LAST_DEFAULT)]
        return first + [im for im in images if im.name.startswith(LAST_DEFAULT)]
    return resolve_images(images, toks, what)


def print_resolved(chosen: list, stream=sys.stderr) -> None:
    for im in chosen:
        print(f'image  {im.tag}  = {im.name!r} (entry {im.entry_id})', file=stream)


# ----------------------------------------------------------------------------------------------
# records
@dataclass
class Rec:
    data: dict
    site: str
    slot: str       # canonical | history | external
    path: Path

    @property
    def seq(self):
        return self.data.get('seq')

    @property
    def step(self):
        return self.data.get('step')

    @property
    def attempt(self):
        return self.data.get('attempt_id')

    def label(self) -> str:
        return f"seq {self.seq} ({self.site}:{self.step})"


class Records:
    def __init__(self, site: str):
        self.site = site
        self.recs: list = []              # step sentinels + external records, local and merged
        self.post: dict = {}              # (site, step) -> post record
        self.finish: dict = {}            # (site, step) -> finish record
        self.drift: list = []             # drift record paths (local and merged)
        self.running: dict = {}           # step -> in-flight marker of THIS site
        self.foreign_running: list = []   # (site, step, marker) of another site, attempt not finished there
        self.unmerged: list = []          # (site, relpath, why)
        self.legacy: list = []            # records in the flat pre-site layout
        self.invalid: list = []           # unreadable record files

    def seq_records(self) -> list:
        return [r for r in self.recs if isinstance(r.data, dict) and isinstance(r.data.get('seq'), int)]

    def canonical(self, step: str):
        cands = [r for r in self.seq_records() if r.slot == 'canonical' and r.step == step]
        return max(cands, key=lambda r: (r.seq, r.site)) if cands else None

    def externals(self) -> list:
        return [r for r in self.seq_records() if r.slot == 'external']

    def max_seq(self) -> int:
        return max((r.seq for r in self.seq_records()), default=0)

    def ever_ran(self, step: str) -> bool:
        return any(r.step == step for r in self.seq_records())


def prov_dir(out: Path, tag: str) -> Path:
    return out / 'provenance' / tag


def site_dir(out: Path, tag: str, site: str) -> Path:
    return prov_dir(out, tag) / site


def native_files(sdir: Path) -> dict:
    """The record files a site wrote itself: top-level *.json (no markers) and history/*.json."""
    out = {}
    if not sdir.is_dir():
        return out
    for f in sdir.iterdir():
        if f.is_file() and f.name.endswith('.json') and not f.name.startswith('.'):
            out[f.name] = f
    h = sdir / 'history'
    if h.is_dir():
        for f in h.iterdir():
            if f.is_file() and f.name.endswith('.json') and not f.name.startswith('.'):
                out['history/' + f.name] = f
    return out


def _read_record(rec: Records, f: Path):
    try:
        return read_json(f)
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        rec.invalid.append(f'{f}: {e}')
        return None


def _load_site_dir(rec: Records, sdir: Path, site: str, local: bool) -> None:
    seen = set()   # (attempt_id or 'external', seq): an identical record reachable twice counts once
    for f in sorted(sdir.iterdir()) if sdir.is_dir() else []:
        n = f.name
        if not f.is_file() or not n.endswith('.json'):
            continue
        if m := re.fullmatch(r'\.running\.([A-Za-z0-9_]+)\.json', n):
            if local:
                rec.running[m.group(1)] = _read_record(rec, f)
            continue
        data = _read_record(rec, f)
        if m := re.fullmatch(r'([A-Za-z0-9_]+)\.json', n):
            if not isinstance(data, dict):
                rec.invalid.append(f'{f}: not a record')
                continue
            rec.recs.append(Rec(data, site, 'canonical', f))
            seen.add((data.get('attempt_id'), data.get('seq')))
        elif m := re.fullmatch(r'([A-Za-z0-9_]+)\.post\.json', n):
            rec.post[(site, m.group(1))] = data
        elif m := re.fullmatch(r'([A-Za-z0-9_]+)\.finish\.json', n):
            rec.finish[(site, m.group(1))] = data
        elif re.fullmatch(r'([A-Za-z0-9_]+)\.drift\.json', n):
            rec.drift.append(f)
        elif re.fullmatch(r'external\.\d+\.json', n):
            if isinstance(data, dict):
                rec.recs.append(Rec(data, site, 'external', f))
                seen.add(('external', data.get('seq')))
    hist = sdir / 'history'
    if hist.is_dir():
        for f in sorted(hist.glob('*.json')):
            data = _read_record(rec, f)
            if isinstance(data, dict) and str(data.get('schema', '')).startswith('run_v3.provenance'):
                key = (data.get('attempt_id') if data.get('kind') != 'external' else 'external', data.get('seq'))
                if key not in seen:
                    seen.add(key)
                    rec.recs.append(Rec(data, site, 'history', f))


def load_records(tagdir: Path, site: str, merged_override: dict | None = None) -> Records:
    """Records of one image as seen by `site`: its own directory, the other sites' records it merged
    (merged/<site2>/), and a list of what other sites wrote that is not merged yet.
    merged_override {site2: dir} substitutes a merged copy (merge-provenance validates a candidate)."""
    rec = Records(site)
    if not tagdir.is_dir():
        return rec
    rec.legacy = sorted(f for f in tagdir.iterdir() if f.is_file() and f.name.endswith('.json'))
    local = tagdir / site
    _load_site_dir(rec, local, site, local=True)
    merged = {}
    mroot = local / 'merged'
    if mroot.is_dir():
        merged = {d.name: d for d in mroot.iterdir() if d.is_dir() and not d.name.startswith('.')}
    merged.update(merged_override or {})
    for s2, d in sorted(merged.items()):
        _load_site_dir(rec, d, s2, local=False)
    for sd in sorted(tagdir.iterdir()):
        if not sd.is_dir() or sd.name == site or sd.name.startswith('.'):
            continue
        native = native_files(sd)
        mirror = native_files(merged[sd.name]) if sd.name in merged else {}
        for relpath, f in sorted(native.items()):
            m = mirror.get(relpath)
            if m is None:
                rec.unmerged.append((sd.name, relpath, 'not merged'))
            elif not filecmp.cmp(m, f, shallow=False):
                rec.unmerged.append((sd.name, relpath, 'changed since the merge'))
        for relpath in sorted(set(mirror) - set(native)):
            rec.unmerged.append((sd.name, relpath, 'removed at the source since the merge'))
        # an attempt of that site is in flight unless it wrote its finish record there
        for f in sorted(sd.glob('.running.*.json')):
            marker = read_json(f) if f.is_file() else None
            step = f.name[len('.running.'):-len('.json')]
            fin = read_json(sd / f'{step}.finish.json') if (sd / f'{step}.finish.json').is_file() else None
            if not (isinstance(marker, dict) and isinstance(fin, dict) and fin.get('attempt_id') == marker.get('attempt_id')):
                rec.foreign_running.append((sd.name, step, marker))
    return rec


def custody(recs: Records, current: dict) -> tuple[str, Rec | None, str]:
    """Is data.qpdata still the file the last recorded step left behind, and do the records form one
    chain (every record starts from the file the previous one left)?
    -> ('fresh'|'ok'|'broken'|'conflict', last record, message)"""
    chain = sorted((r for r in recs.seq_records() if (r.data.get('qpdata') or {}).get('after')),
                   key=lambda r: (r.seq, r.site))
    if not chain:
        return 'fresh', None, 'no run_v3 record yet'
    by_seq = {}
    for r in chain:
        by_seq.setdefault(r.seq, []).append(r)
    dups = {s: rs for s, rs in by_seq.items() if len(rs) > 1}
    if dups:
        top = max(by_seq)
        desc = '; '.join(f'seq {s}: ' + ', '.join(f'{r.site}:{r.step}' for r in rs) for s, rs in sorted(dups.items()))
        return 'conflict', None, (
            f'records share a sequence number ({desc}){" -- including the highest" if top in dups else ""}: the image '
            f'was processed on two sites from the same state, so neither record can be trusted to describe '
            f'data.qpdata. Keep one side: restore that side\'s data.qpdata, move the other side\'s records '
            f'out of provenance/<tag>/ (its own <site>/ directory and any merged/<site>/ copy), then plan again.')
    breaks = []
    for prev, nxt in zip(chain, chain[1:]):
        before = (nxt.data.get('qpdata') or {}).get('before')
        if not same_stat(before, prev.data['qpdata']['after']):
            breaks.append(f'{nxt.label()} started from {fmt_stat(before)} but {prev.label()} left '
                          f'{fmt_stat(prev.data["qpdata"]["after"])}')
    for r in chain:
        q = r.data['qpdata']
        if r.data.get('status') != 'pass' and not same_stat(q.get('before'), q.get('after')):
            breaks.append(f'{r.label()} failed but data.qpdata changed during it')
    if breaks:
        return 'conflict', None, ('custody chain broken: ' + '; '.join(breaks[:4]) +
                                  (f' (+{len(breaks) - 4} more)' if len(breaks) > 4 else '') +
                                  '. A record is missing, or records from different copies of the image were mixed.')
    ref = chain[-1]
    after = ref.data['qpdata']['after']
    if same_stat(after, current):
        return 'ok', ref, f'matches {ref.label()}'
    running = ', '.join(f"{s} (attempt {m.get('attempt_id')}, started {m.get('started')})"
                        for s, m in recs.running.items() if isinstance(m, dict))
    match = next((r for r in reversed(chain) if same_stat(r.data['qpdata']['after'], current)), None)
    msg = (f"data.qpdata changed outside run_v3: the last record {ref.label()} left {fmt_stat(after)}, "
           f"now {fmt_stat(current)}.")
    if match is not None:
        msg += f' The current file is the one {match.label()} left: an older copy (sync data/<id>/ from the site that ran {ref.label()}?).'
    if running:
        msg += f' An attempt was interrupted: {running}; every step is an idempotent re-run.'
    msg += (' If the change is intended (e.g. gate-1 review in the GUI) run: run_v3.sh --accept-external '
            '--images <image> --note "why" [--rerun-from STEP]; otherwise restore the file.')
    return 'broken', ref, msg


def image_state(recs: Records, current: dict, allow_drift: bool) -> tuple[str, str, Rec | None]:
    """(state, message, last record); a state other than fresh/ok blocks every step of the image."""
    if recs.legacy:
        return 'legacy_layout', (f'records in the old flat layout: {", ".join(str(p) for p in recs.legacy[:3])}; '
                                 f'move them into the {recs.site}/ subdirectory of that folder (they were '
                                 f'written by this site) and plan again'), None
    if recs.invalid:
        return 'records_invalid', 'unreadable record file(s): ' + '; '.join(recs.invalid[:3]), None
    if recs.unmerged:
        sites = sorted({s for s, _, _ in recs.unmerged})
        ex = '; '.join(f'{s}/{p} ({why})' for s, p, why in recs.unmerged[:4])
        return 'unmerged', (f'records from site(s) {", ".join(sites)} are not merged ({ex}{" ..." if len(recs.unmerged) > 4 else ""}). '
                            f'Sync the image\'s data/<id>/ from that site too, then run: run_v3.sh --merge-provenance '
                            f'--images <image>'), None
    if recs.foreign_running:
        ex = '; '.join(f"{s}: step {st} attempt {(m or {}).get('attempt_id')} started {(m or {}).get('started')}"
                       for s, st, m in recs.foreign_running)
        return 'foreign_in_flight', (f'another site has an attempt in flight or interrupted ({ex}). Wait for it; if '
                                     f'no job runs there any more, delete that .running.*.json marker on both sides'), None
    state, ref, msg = custody(recs, current)
    if state in ('broken', 'conflict'):
        return ('custody_broken' if state == 'broken' else 'custody_conflict'), msg, ref
    if recs.drift and not allow_drift:
        files = ', '.join(str(p) for p in recs.drift)
        return 'drift', (f'unresolved server-metadata / project-entry drift: {files}. Merge the change into the '
                         f'project in the GUI (or confirm it is unwanted), delete the drift file(s), then re-run '
                         f'(or pass --allow-drift).'), ref
    return state, msg, ref


# ----------------------------------------------------------------------------------------------
# planning
class Ctx:
    def __init__(self, steps: dict, recs: Records, scripts: dict | None = None):
        self.steps, self.recs, self.scripts, self.memo = steps, recs, scripts, {}
        self.externals = recs.externals()

    def _stale_by_external(self, sid: str, seq: int):
        st = self.steps[sid]
        for ext in sorted(self.externals, key=lambda r: r.seq):
            if ext.seq <= seq:
                continue
            note = ext.data.get('note') or ''
            if st.on_external == 'stale':
                return (f'stale: an external edit was accepted after it (seq {ext.seq} > {seq}: "{note}"); '
                        f'this step reads the edited image')
            rf = [s for s in (ext.data.get('rerun_from') or []) if s in self.steps]
            if rf and sid in downstream_closure(self.steps, rf):
                return f'stale: the external edit seq {ext.seq} ("{note}") declared --rerun-from {",".join(rf)}'
        return None

    def done(self, sid: str) -> tuple[bool, str]:
        if sid in self.memo:
            return self.memo[sid]
        self.memo[sid] = (False, 'dependency cycle')      # guard; replaced below
        st, r = self.steps[sid], self.recs.canonical(sid)
        res = self._done(sid, st, r)
        self.memo[sid] = res
        return res

    def _done(self, sid, st, r) -> tuple[bool, str]:
        if r is None:
            return False, 'never ran'
        rec = r.data
        if rec.get('status') != 'pass':
            pc = rec.get('postcondition') or {}
            return False, (f"last attempt failed: {pc.get('message') or ''} {pc.get('error') or ''} "
                           f"{(rec.get('save') or {}).get('policy_error') or ''}".strip())
        fin = self.recs.finish.get((r.site, sid))
        if not (isinstance(fin, dict) and fin.get('attempt_id') == rec.get('attempt_id')):
            return False, (f"the sentinel of attempt {rec.get('attempt_id')} was never validated by finish "
                           f"(the driver was interrupted, or an orphaned QuPath wrote it later); it re-runs")
        if fin.get('status') != 'PASS':
            return False, f"finish said {fin.get('status')}: {fin.get('reason') or ''}"
        if st.post:
            prec = self.recs.post.get((r.site, sid))
            if not (isinstance(prec, dict) and prec.get('status') == 'pass' and prec.get('attempt_id') == rec.get('attempt_id')):
                return False, f'post hook {st.post} not done'
        outs = list(rec.get('outputs') or [])
        if st.post:
            outs += list((self.recs.post.get((r.site, sid)) or {}).get('outputs') or [])
        for o in outs:
            s = stat_of(Path(o['path']))
            if not s['exists'] or s.get('size') != o.get('size'):
                return False, f"output missing or changed: {o['path']}"
        if self.scripts is not None:
            cur = self.scripts.get(sid) or {'error': 'not resolved'}
            ran = rec.get('script') or {}
            if 'error' in cur:
                return False, f"its script cannot be resolved now: {cur['error']}"
            if (ran.get('source_path'), ran.get('source_sha256')) != (cur['path'], cur['sha256']):
                return False, (f"script changed: it ran {ran.get('source_path')} (sha256 {str(ran.get('source_sha256'))[:12]}), "
                               f"the plan now uses {cur['path']} ({cur['sha256'][:12]}, {cur['how']})")
        why = self._stale_by_external(sid, rec['seq'])
        if why:
            return False, why
        for dep, optional in st.after:
            if optional and not self.recs.ever_ran(dep):
                continue
            ok, why = self.done(dep)
            if not ok:
                return False, f'upstream {dep} not done ({why})'
            dseq = self.recs.canonical(dep).seq
            if dseq > rec['seq']:
                return False, f"stale: {dep} re-ran after it (seq {dseq} > {rec['seq']})"
        # 'needs' only gates whether a step may START (plan_image checks it for steps that run); a step
        # done earlier stays done when a needed step is re-run or stale -- e.g. an accepted gate-1 edit
        # re-runs preflight, which must never make step 1 re-run over the reviewed islets
        return True, f"done (seq {rec.get('seq')}, {r.site})"


def selected_steps(steps: dict, passes: list, only: list) -> list:
    """The steps to consider, in dependency order (table order breaks ties)."""
    if only:
        bad = [s for s in only if s not in steps]
        if bad:
            raise SystemExit(f'--steps: unknown step(s) {bad}; known: {list(steps)}')
        chosen = [s for s in steps if s in only]
    elif passes:
        known = {p for st in steps.values() for p in st.passes}
        bad = [p for p in passes if p not in known]
        if bad:
            raise SystemExit(f'--pass: unknown pass(es) {bad}; known: {sorted(known)}')
        chosen = [s for s, st in steps.items() if set(st.passes) & set(passes)]
    else:
        raise SystemExit('give --pass (e.g. I or I,C,J) or --steps')
    order, left = [], list(chosen)
    while left:
        ready = next((s for s in left if not ({d for d, _ in steps[s].after} | set(steps[s].needs)) & (set(left) - {s})), None)
        if ready is None:
            raise SystemExit(f'dependency cycle among {left}')
        order.append(ready)
        left.remove(ready)
    return order


def plan_image(im: Image, steps: dict, sel: list, force: set, out: Path, allow_drift: bool,
               scripts: dict | None = None, site: str | None = None) -> dict:
    site = site or resolve_site(None)
    recs = load_records(prov_dir(out, im.tag), site)
    cur = stat_of(im.qpdata)
    state, msg, ref = image_state(recs, cur, allow_drift)
    if state not in ('fresh', 'ok'):
        return {'image': im, 'state': state, 'message': msg + bkp_note(im.qpdata), 'qpdata': cur,
                'steps': [(s, 'blocked', f'{state} (see image message)') for s in sel]}
    ctx = Ctx(steps, recs, scripts)
    run, why = {}, {}
    for s in sel:
        ok, reason = ctx.done(s)
        if s in force:
            run[s], why[s] = True, 'forced (--force / --rerun-from)'
        elif not ok:
            run[s], why[s] = True, reason
        else:
            run[s], why[s] = False, reason
    for s in sel:        # a step that re-runs makes everything downstream of it (in the selection) stale
        if not run[s]:
            ups = [u for u in upstream_closure(steps, s) if run.get(u)]
            if ups:
                run[s], why[s] = True, f'downstream of re-running {",".join(sorted(ups))}'
    rows, blocked = [], None
    for s in sel:
        if blocked:
            rows.append((s, 'blocked', f'after blocked step {blocked}'))
            continue
        if not run[s]:
            rows.append((s, 'skip', why[s]))
            continue
        st = steps[s]
        missing = []
        for dep, optional in st.after:
            if dep in sel and run.get(dep):
                continue
            if optional and not recs.ever_ran(dep):
                continue
            ok, r = ctx.done(dep)
            if not ok:
                missing.append(f'{dep} ({r})')
        for dep in st.needs:
            if dep in sel and run.get(dep):
                continue
            ok, r = ctx.done(dep)
            if not ok:
                missing.append(f'{dep} ({r})')
        if missing:
            blocked = s
            rows.append((s, 'blocked', 'prerequisite not done: ' + '; '.join(missing)))
        else:
            rows.append((s, 'run', why[s]))
    stale = [f'{s} ({m.get("attempt_id")})' for s, m in recs.running.items() if isinstance(m, dict)]
    message = msg + (f'; stale in-flight marker(s): {", ".join(stale)}' if stale else '') + bkp_note(im.qpdata)
    return {'image': im, 'state': state, 'message': message, 'qpdata': cur, 'steps': rows}


def static_checks(im: Image) -> list:
    """Cheap checks from project.qpproj alone (dry runs); the preflight step does the real ones."""
    md = im.entry.get('serverBuilder', {}).get('metadata', {})
    names = [c.get('name') for c in md.get('channels', [])]
    req = ['DAPI', 'CD68', 'CD44', 'Vimentin', 'CD99', 'E-cadherin', 'M2Gal3', 'HLA-A', 'Pan-Cytokeratin',
           'CD79a', 'CD3e', 'Beta-actin', 'Iba1', 'INS', 'GCG', 'SST']
    warn = []
    miss = [r for r in req if r not in names]
    if miss:
        warn.append(f'missing channels {miss}')
    pc = md.get('pixelCalibration', {})
    pw, ph = (pc.get('pixelWidth') or {}).get('value'), (pc.get('pixelHeight') or {}).get('value')
    if pw is None or abs(pw - 0.5078) > 1e-6 or abs(ph - 0.5078) > 1e-6:
        warn.append(f'pixel size {pw} x {ph} != 0.5078 (preflight will fail; fix the entry in the GUI)')
    return warn


# ----------------------------------------------------------------------------------------------
# locks: the per-host flock (run_v3.sh's LOCKDIR/<tag>.lock) and the shared mkdir lock that works
# across the nodes of a cluster filesystem (HiPerGator /blue: V3_SHARED_LOCKDIR)
def proc_start(pid: int):
    try:
        with open(f'/proc/{pid}/stat') as f:
            return int(f.read().rsplit(')', 1)[1].split()[19])
    except (OSError, ValueError, IndexError):
        return None


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def lock_owner(pid: int, run_id: str = '') -> dict:
    job = os.environ.get('SLURM_JOB_ID')
    if os.environ.get('SLURM_ARRAY_JOB_ID') and os.environ.get('SLURM_ARRAY_TASK_ID'):
        job = f"{os.environ['SLURM_ARRAY_JOB_ID']}_{os.environ['SLURM_ARRAY_TASK_ID']}"
    return {'host': socket.gethostname(), 'pid': pid, 'pid_start': proc_start(pid), 'slurm_job': job,
            'run_id': run_id, 'user': getpass.getuser(), 'time': now_iso()}


def owner_alive(info, dir_mtime: float) -> tuple[bool, str]:
    if not isinstance(info, dict):
        young = time.time() - dir_mtime < 120
        return young, ('its owner is still writing owner.json' if young else 'no owner.json (stale)')
    desc = (f"pid {info.get('pid')} on {info.get('host')} (run {info.get('run_id') or '-'}, "
            f"job {info.get('slurm_job') or '-'}, since {info.get('time')})")
    if info.get('host') == socket.gethostname():
        pid = info.get('pid')
        if isinstance(pid, int) and pid_alive(pid) and info.get('pid_start') in (None, proc_start(pid)):
            return True, desc
        return False, desc + ' -- that process is gone'
    job = info.get('slurm_job')
    if job and shutil.which('squeue'):
        try:
            r = subprocess.run(['squeue', '-h', '-j', str(job), '-o', '%T'], capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.TimeoutExpired) as e:
            return True, desc + f' (squeue failed: {e}; assuming it is alive)'
        state = (r.stdout.strip().splitlines() or [''])[0].strip()
        if r.returncode == 0 and state and state not in SLURM_ENDED:
            return True, desc + f' (job {state})'
        if (r.returncode == 0 and (not state or state in SLURM_ENDED)) or 'Invalid job id' in r.stderr:
            return False, desc + ' -- that job has ended'
        return True, desc + f' (squeue: {r.stderr.strip()[:100]}; assuming it is alive)'
    return True, desc + ' -- on another host, cannot check it from here (remove the lock directory by hand if it is gone)'


def shared_lock_acquire(root: Path, tag: str, owner: dict) -> tuple[bool, str]:
    root.mkdir(parents=True, exist_ok=True)
    d = root / f'{tag}.lock.d'
    for _ in range(5):
        try:
            d.mkdir()                                   # atomic on POSIX and on Lustre/GPFS
        except FileExistsError:
            try:
                mtime = d.stat().st_mtime
            except FileNotFoundError:
                continue
            info = read_json(d / 'owner.json') if (d / 'owner.json').is_file() else None
            alive, desc = owner_alive(info, mtime)
            if alive:
                return False, f'{d} held by {desc}'
            grave = root / f'.{tag}.lock.d.stale.{socket.gethostname()}.{os.getpid()}.{time.time_ns()}'
            try:
                os.rename(d, grave)                     # only one breaker wins
            except OSError:
                continue
            print(f'lock: broke stale {d} ({desc})', file=sys.stderr)
            shutil.rmtree(grave, ignore_errors=True)
            continue
        write_json(d / 'owner.json', owner)
        return True, f'acquired {d}'
    return False, f'{d}: could not acquire (contended)'


def shared_lock_release(root: Path, tag: str, pid: int) -> tuple[bool, str]:
    d = root / f'{tag}.lock.d'
    info = read_json(d / 'owner.json') if (d / 'owner.json').is_file() else None
    if not d.exists():
        return True, f'{d} not held'
    if not (isinstance(info, dict) and info.get('host') == socket.gethostname() and info.get('pid') == pid):
        return False, f'{d} is not ours (owner {info})'
    grave = root / f'.{tag}.lock.d.released.{os.getpid()}.{time.time_ns()}'
    os.rename(d, grave)
    shutil.rmtree(grave, ignore_errors=True)
    return True, f'released {d}'


@contextmanager
def image_locks(images: list, lockdir: Path, shared: Path | None, run_id: str):
    """Hold every image's host flock (non-blocking) and, with a shared lock dir, its mkdir lock."""
    held, taken = [], []
    try:
        lockdir.mkdir(parents=True, exist_ok=True)
        for im in images:
            f = open(lockdir / f'{im.tag}.lock', 'a')
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as e:
                f.close()
                if e.errno in (errno.EAGAIN, errno.EACCES):
                    raise SystemExit(f'{im.name}: image lock {lockdir / (im.tag + ".lock")} is held by a running '
                                     f'run_v3 (or an orphaned QuPath); try again when it has finished')
                raise
            held.append(f)
            if shared is not None:
                ok, msg = shared_lock_acquire(shared, im.tag, lock_owner(os.getpid(), run_id))
                if not ok:
                    raise SystemExit(f'{im.name}: shared image lock: {msg}')
                taken.append(im.tag)
        yield
    finally:
        for tag in taken:
            shared_lock_release(shared, tag, os.getpid())
        for f in held:
            f.close()


# ----------------------------------------------------------------------------------------------
# run log
def append_run_log(out: Path, row: dict) -> None:
    path = out / 'run_log.csv'
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {'site': row.get('site') or os.environ.get('V3_SITE', ''), **row}
    with open(path, 'a+', newline='', encoding='utf-8') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.seek(0)
            head = f.readline()
            fields = next(csv.reader([head])) if head.strip() else RUN_LOG_FIELDS   # keep an older file's columns
            f.seek(0, os.SEEK_END)
            if not head.strip():
                csv.writer(f).writerow(fields)
            csv.writer(f).writerow([row.get(k, '') for k in fields])
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


# ----------------------------------------------------------------------------------------------
# commands
def common(a):
    repo = Path(a.repo).resolve()
    project = Path(a.project).resolve()
    out = Path(a.out).resolve()
    steps = load_steps(Path(a.steps_file))
    return repo, project, out, steps


def find_image(project: Path, name: str) -> Image:
    return resolve_images(load_images(project), [name], '--image')[0]


def cmd_steps(a):
    repo, _, _, steps = common(a)
    so = parse_kv_list(a.script, '--script')
    print('\t'.join(['step', 'passes', 'lane', 'heap_gb', 'timeout_min', 'save', 'on_external', 'script', 'selection',
                     'postcondition', 'post']))
    for st in steps.values():
        try:
            p, how = resolve_script(st, repo, a.script_set, so)
        except PlanError as e:
            p, how = f'ERROR: {e}', '-'
        print('\t'.join([st.id, ','.join(st.passes), st.lane, str(st.heap_gb), str(st.timeout_min), st.save,
                         st.on_external, str(p), how, st.postcondition['name'], st.post or '-']))


def cmd_images(a):
    _, project, _, _ = common(a)
    chosen = select_images(load_images(project), a.images, a.images_file)
    print_resolved(chosen)
    for im in chosen:
        print(im.name)


def cmd_budget(a):
    _, _, _, steps = common(a)
    sel = selected_steps(steps, [p for p in (a.passes or '').split(',') if p], [s for s in (a.steps or '').split(',') if s])
    total = sum(steps[s].timeout_min for s in sel) + int(a.margin_min)
    print(f"steps {', '.join(f'{s}={steps[s].timeout_min}m' for s in sel)} + margin {a.margin_min}m = {total} min",
          file=sys.stderr)
    if a.format == 'slurm':
        d, rem = divmod(total, 24 * 60)
        print(f'{d}-{rem // 60:02d}:{rem % 60:02d}:00' if d else f'{rem // 60:02d}:{rem % 60:02d}:00')
    else:
        print(total)


def plan_force(a, steps: dict, sel: list) -> set:
    force = {s for s in (a.force or '').split(',') if s} | {s for s in (a.rerun_from or '').split(',') if s}
    unknown = force - set(steps)
    if unknown:
        raise SystemExit(f'--force/--rerun-from: unknown step(s) {sorted(unknown)}')
    outside = force - set(sel)
    if outside:
        raise SystemExit(f'--force/--rerun-from {",".join(sorted(outside))}: not among the selected steps '
                         f'({",".join(sel)}); add its pass with --pass or name it in --steps')
    return force


def cmd_plan(a):
    repo, project, out, steps = common(a)
    site = resolve_site(a.site)
    so = parse_kv_list(a.script, '--script')
    heaps = {k: int(v) for k, v in parse_kv_list(a.heap, '--heap').items()}
    sel = selected_steps(steps, [p for p in (a.passes or '').split(',') if p], [s for s in (a.steps or '').split(',') if s])
    force = plan_force(a, steps, sel)
    images = select_images(load_images(project), a.images, a.images_file)
    if a.images or a.images_file:
        print_resolved(images)
    scripts = resolve_all_scripts(steps, repo, a.script_set, so)
    plans = [plan_image(im, steps, sel, force, out, a.allow_drift, scripts, site) for im in images]
    if a.format == 'tsv':
        for p in plans:
            im = p['image']
            print('\t'.join(['IMAGE', im.tag, im.name, str(im.entry_id), p['state'],
                             p['message'].replace('\t', ' ').replace('\n', ' ')]))
            for s, action, reason in p['steps']:
                st = steps[s]
                spath = scripts[s].get('path') or f"ERROR: {scripts[s].get('error')}"
                if action == 'run' and 'error' in scripts[s]:
                    action, reason = 'blocked', spath
                print('\t'.join(['STEP', im.tag, s, action, st.lane, str(heaps.get(s, st.heap_gb)), str(st.threads),
                                 str(st.timeout_min), st.post or '-', spath, reason.replace('\t', ' ').replace('\n', ' ')]))
        return
    # human-readable table (dry runs)
    print(f'project : {project}\noutput  : {out}   (site {site})\nsteps   : {", ".join(sel)}  '
          f'(force: {", ".join(sorted(force)) or "-"})')
    print(f'scripts : (--script-set {a.script_set})')
    for s in sel:
        sc = scripts[s]
        ov = steps[s].overrides
        print(f"  {s:<14} {sc.get('how', 'ERROR'):<10} {sc.get('path') or 'ERROR: ' + sc.get('error', '')}" +
              (f'   overrides {ov}' if ov else ''))
    for p in plans:
        im = p['image']
        q = p['qpdata']
        size = f"{q['size'] / 1e9:.2f} GB" if q.get('exists') else 'missing'
        print(f"\n== {im.name}  [entry {im.entry_id}, tag {im.tag}]\n   data.qpdata {size}; custody: {p['state']} -- {p['message']}")
        for w in static_checks(im):
            print(f'   WARNING: {w}')
        for s, action, reason in p['steps']:
            st = steps[s]
            print(f"   {action.upper():<8} {s:<14} lane={st.lane:<5} heap={heaps.get(s, st.heap_gb):>3}g "
                  f"save={st.save:<12} {reason}")


def cmd_diskcheck(a):
    im = find_image(Path(a.project).resolve(), a.image)
    q = stat_of(im.qpdata)
    size = q.get('size', 0) if q.get('exists') else 0
    data_dir = im.qpdata.parent if im.qpdata.parent.exists() else im.qpdata.parent.parent
    free_data = shutil.disk_usage(data_dir).free
    need_data = 2 * size + a.data_margin_gb * 10 ** 9
    free_root = shutil.disk_usage(a.root_path).free
    need_root = a.min_root_free_gb * 10 ** 9
    problems = []
    if free_data < need_data:
        problems.append(f'{data_dir}: {free_data / 1e9:.0f} GB free < 2 x {size / 1e9:.1f} GB data.qpdata + '
                        f'{a.data_margin_gb} GB = {need_data / 1e9:.0f} GB (QuPath keeps the old file as .bkp while writing)')
    if free_root < need_root:
        problems.append(f'{a.root_path}: {free_root / 1e9:.0f} GB free < {a.min_root_free_gb} GB floor')
    if problems:
        print('DISK LOW: ' + '; '.join(problems))
        sys.exit(1)
    print(f'disk ok: data fs {free_data / 1e9:.0f} GB free (need {need_data / 1e9:.0f}), '
          f'{a.root_path} {free_root / 1e9:.0f} GB free (need {a.min_root_free_gb})')


def cmd_prepare(a):
    repo, project, out, steps = common(a)
    site = resolve_site(a.site)
    st = steps[a.step]
    im = find_image(project, a.image)
    so = parse_kv_list(a.script, '--script')
    extra = {}
    for item in a.override or []:
        sid, _, kv = item.partition(':')
        if sid == st.id:
            extra.update(parse_kv_list([kv], '--override'))
    heap = int(parse_kv_list(a.heap, '--heap').get(st.id, st.heap_gb))
    pdir = site_dir(out, im.tag, site)
    recs = load_records(prov_dir(out, im.tag), site)
    before = stat_of(im.qpdata)
    state, msg, _ = image_state(recs, before, allow_drift=True)
    if state not in ('fresh', 'ok'):      # the plan checked this; the file may have changed since
        raise PlanError(f'{im.name}: {state}: {msg}')
    seq = recs.max_seq() + 1
    attempt = uuid.uuid4().hex[:12]

    # resolve and instantiate first: a refused script/override must leave the old records alone
    src, how = resolve_script(st, repo, a.script_set, so)
    body = src.read_text(encoding='utf-8')
    body, applied = apply_overrides(body, {**st.overrides, **extra})
    pdir.mkdir(parents=True, exist_ok=True)

    # superseded records of this step go to history/ (custody still sees them there)
    hist = pdir / 'history'
    try:
        old = read_json(pdir / f'{st.id}.json')
    except (json.JSONDecodeError, UnicodeDecodeError):
        old = None
    old_seq = old.get('seq') if isinstance(old, dict) else 'x'
    for suffix in ('.json', '.finish.json', '.post.json', '.groovy', '.drift.json'):
        f = pdir / f'{st.id}{suffix}'
        if f.exists():
            hist.mkdir(exist_ok=True)
            os.replace(f, hist / f'{st.id}.seq{old_seq}.{time.time_ns()}{suffix}')

    worker = Path(a.workdir)
    worker.mkdir(parents=True, exist_ok=True)
    shutil.copy2(project, worker / 'project.qpproj')        # private copy: QuPath never sees the target file
    for d in ('data', 'classifiers', 'scripts', 'resources'):
        target = project.parent / d
        link = worker / d
        if target.exists() and not link.exists():
            link.symlink_to(target)
    step_cfg = worker / 'step_config.json'
    subst = {'{OUT}': str(out), '{TMP}': str(Path(a.tmp)), '{TAG}': im.tag, '{PROJ}': str(repo), '{STEP}': st.id,
             '{CFG}': str(step_cfg)}
    args = []
    for x in st.args:
        for k, v in subst.items():
            x = x.replace(k, v)
        args.append(x)
    script_path = worker / 'step.groovy'
    deadline_ms = int(float(a.deadline_epoch) * 1000) if a.deadline_epoch and float(a.deadline_epoch) > 0 else None
    cfg = {
        'run_id': a.run_id, 'attempt_id': attempt, 'seq': seq, 'step': st.id, 'passes': st.passes,
        'image': im.name, 'tag': im.tag, 'entry_id': im.entry_id, 'project': str(project), 'site': site,
        'worker_project': str(worker / 'project.qpproj'), 'out': str(out),
        'provenance': str(pdir / f'{st.id}.json'), 'save_policy': st.save, 'postcondition': st.postcondition,
        'lane': st.lane, 'slot': a.slot, 'heap_gb': heap, 'threads': st.threads,
        'host': socket.gethostname(), 'gpu': a.gpu or None,
        'source_script': str(src), 'source_sha256': sha256_file(src), 'script_selection': how,
        'repo_git': git_info(repo / 'scripts', src),
        'body_sha256': hashlib.sha256(body.encode('utf-8')).hexdigest(), 'instantiated_script': str(script_path),
        'overrides': applied, 'args': args, 'footer_template_sha256': sha256_file(FOOTER_TEMPLATE),
        'qpdata_before': before, 'driver_start': now_iso(), 'step_config': str(step_cfg),
        'deadline_epoch_ms': deadline_ms, 'save_margin_min': float(a.save_margin_min),
    }
    write_json(step_cfg, cfg)                  # scripts that save themselves (step 0) read it via {CFG}
    footer = FOOTER_TEMPLATE.read_text(encoding='utf-8')
    if footer.count(CFG_PLACEHOLDER) != 1:
        raise SystemExit(f'{FOOTER_TEMPLATE}: expected exactly one {CFG_PLACEHOLDER!r}')
    footer = footer.replace(CFG_PLACEHOLDER, groovy_literal(cfg))
    script_path.write_text(body.rstrip('\n') + '\n' + footer, encoding='utf-8')
    write_json(pdir / f'.running.{st.id}.json', {'attempt_id': attempt, 'seq': seq, 'run_id': a.run_id,
                                                  'started': cfg['driver_start'], 'host': cfg['host'], 'site': site,
                                                  'driver_pid': os.getppid(), 'worker': str(worker),
                                                  'qpdata_before': cfg['qpdata_before']})
    qargs = worker / 'qargs.bin'                # NUL-separated: the driver reads it with mapfile -d ''
    with open(qargs, 'wb') as f:
        for x in args:
            f.write(b'-a\0' + x.encode('utf-8') + b'\0')
    values = {'ATTEMPT': attempt, 'SEQ': seq, 'HEAP': heap, 'WSCRIPT': script_path, 'WPROJECT': worker / 'project.qpproj',
              'PROVENANCE': pdir / f'{st.id}.json', 'SCRIPT_SRC': src, 'QARGS_FILE': qargs}
    for k in PREPARE_KEYS:                      # KEY=VALUE lines, parsed (never eval'ed) by run_v3.sh
        v = str(values[k])
        if '\n' in v or '\0' in v:
            raise PlanError(f'{k} contains a newline: {v!r}')
        print(f'{k}={v}')


def entry_diff(a, b, path='') -> list:
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            out += entry_diff(a.get(k, '<absent>'), b.get(k, '<absent>'), f'{path}.{k}' if path else k)
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += entry_diff(x, y, f'{path}[{i}]')
        return out
    return [] if a == b else [{'path': path, 'project': a, 'now': b}]


#: Project-entry fields that are QuPath's own bookkeeping, not what a step loads. `randomizedName` is the
#: label for "mask image names"; the GUI re-rolls it for EVERY entry when it rewrites project.qpproj
#: (measured 2026-10-01 15:06: all 26 entries changed, nothing else did), so a step that was running
#: across a GUI save would otherwise be reported as having changed the image.
QUPATH_BOOKKEEPING = frozenset({'randomizedName'})


def drift_check(project: Path, im: Image, worker: Path, since_ns: int) -> dict | None:
    """Did the step change what QuPath will load next time? Headless runs without -s never write
    project.qpproj, but saveImageData() rewrites data/<id>/server.json whenever the image's server
    metadata (channel names, pixel size...) changed -- and on the next load QuPath takes the metadata
    from the project entry, so such a change silently reverts. Report both, never merge."""
    found = {}
    wp = read_json(worker / 'project.qpproj')
    if wp is not None:
        we = next((e for e in wp.get('images', []) if e.get('entryID') == im.entry_id), None)
        diffs = ([d for d in entry_diff(im.entry, we) if d['path'] not in QUPATH_BOOKKEEPING] if we is not None
                 else [{'path': '<entry>', 'project': 'present', 'now': 'absent'}])
        if diffs:
            found['worker_project_entry'] = diffs
    sj = im.qpdata.parent / 'server.json'
    s = stat_of(sj)
    if s['exists'] and s['mtime_ns'] >= since_ns:
        diffs = entry_diff(im.entry.get('serverBuilder'), read_json(sj))
        if diffs:
            found['server_json'] = {'path': str(sj), 'diffs': diffs}
    return found or None


def parse_time_file(path: Path) -> dict:
    out = {}
    try:
        for ln in path.read_text().splitlines():
            if 'Maximum resident set size' in ln:
                out['peak_rss_kb'] = int(ln.rsplit(':', 1)[1])
    except (FileNotFoundError, ValueError):
        pass
    return out


def refusal_lines(log: str, attempt: str) -> list:
    """The footer's / script 0's refusal messages of THIS attempt, from the per-step log."""
    try:
        text = Path(log).read_text(encoding='utf-8', errors='replace')
    except (OSError, TypeError):
        return []
    i = text.rfind(f'attempt {attempt} ')
    part = text[i:] if i >= 0 else text[-200_000:]
    keys = ('[run_v3] REFUSED', '[run_v3] SAVE FAILED', '[detect] REFUSED', 'REFUSED SAVE')
    out = []
    for ln in part.splitlines():
        if any(k in ln for k in keys):
            ln = ln.strip()
            if ln not in out:
                out.append(ln)
    return out[:3]


def cmd_finish(a):
    repo, project, out, steps = common(a)
    site = resolve_site(a.site)
    st = steps[a.step]
    im = find_image(project, a.image)
    pdir = site_dir(out, im.tag, site)
    marker = read_json(pdir / f'.running.{st.id}.json') or {}
    rec = read_json(pdir / f'{st.id}.json')
    worker = Path(a.workdir)
    script = worker / 'step.groovy'
    status, reason = 'FAIL', ''
    now = stat_of(im.qpdata)
    if rec is None:
        reason = ('no sentinel: the step threw, returned early, was killed or timed out, or refused to save '
                  f'(QuPath exit code {a.exit_code}; see {a.log})')
        refused = refusal_lines(a.log, a.attempt)
        if refused:
            reason += ' -- ' + ' | '.join(refused)
    elif rec.get('attempt_id') != a.attempt:
        reason = f"sentinel belongs to attempt {rec.get('attempt_id')}, not {a.attempt}"
    elif rec.get('status') != 'pass':
        pc = rec.get('postcondition') or {}
        reason = 'footer: ' + ' | '.join(x for x in [pc.get('message'), pc.get('error'),
                                                     (rec.get('save') or {}).get('policy_error')] if x)
    elif script.exists() and rec['script'].get('instantiated_sha256') != sha256_file(script):
        reason = 'the script QuPath ran is not the one prepared (sha256 mismatch)'
    elif not same_stat(rec['qpdata']['after'], now):
        reason = 'data.qpdata changed after the footer recorded it'
    else:
        status, reason = 'PASS', (rec.get('postcondition') or {}).get('message', '')
    note = bkp_note(im.qpdata)
    if note:
        reason = (reason + ' |' if reason else '') + note
    if script.exists():
        shutil.copy2(script, pdir / f'{st.id}.groovy')
    started_ns = int(dt.datetime.fromisoformat(marker['started']).timestamp() * 1e9) if marker.get('started') else 0
    drift = drift_check(project, im, worker, started_ns)
    if drift:
        write_json(pdir / f'{st.id}.drift.json', {'step': st.id, 'attempt_id': a.attempt, 'run_id': a.run_id,
                                                  'detected': now_iso(), **drift})
        if status == 'PASS' and not a.allow_drift:
            status = 'DRIFT'
            reason = f'step changed the image metadata or project entry: see {pdir / (st.id + ".drift.json")}'
    if rec is not None and rec.get('attempt_id') == a.attempt:
        # the driver's verdict on this sentinel; done() requires it (an orphan's sentinel has none)
        write_json(pdir / f'{st.id}.finish.json', {
            'schema': 'run_v3.finish/1', 'step': st.id, 'attempt_id': a.attempt, 'seq': rec.get('seq'),
            'status': status, 'reason': reason[:2000], 'run_id': a.run_id, 'site': site, 'time': now_iso(),
            'exit_code': a.exit_code, 'wall_seconds': a.wall, 'qpdata': now, 'qpdata_bkp': stat_of(bkp_path(im.qpdata)),
            'drift': bool(drift)})
    tinfo = parse_time_file(Path(a.timefile)) if a.timefile else {}
    append_run_log(out, {
        'time_start': marker.get('started', ''), 'time_end': now_iso(), 'run_id': a.run_id,
        'host': socket.gethostname(), 'image': im.name, 'tag': im.tag, 'step': st.id, 'action': 'run',
        'status': status, 'seconds': a.wall, 'exit_code': a.exit_code, 'attempt_id': a.attempt,
        'seq': marker.get('seq', ''), 'lane': st.lane, 'slot': a.slot, 'site': site,
        'heap_gb': (rec or {}).get('jvm', {}).get('heap_requested_gb', ''),
        'peak_rss_kb': ((rec or {}).get('jvm') or {}).get('vm_hwm_kb') or tinfo.get('peak_rss_kb', ''),
        'detections': ((rec or {}).get('counts') or {}).get('detections', ''),
        'qpdata_bytes': now.get('size', ''), 'script': (rec or {}).get('script', {}).get('source_path', ''),
        'script_sha256': (rec or {}).get('script', {}).get('source_sha256', ''), 'message': reason[:500]})
    (pdir / f'.running.{st.id}.json').unlink(missing_ok=True)
    print(f'{status}\t{reason}')
    sys.exit({'PASS': 0, 'DRIFT': 4}.get(status, 1))


def count_lines(path: Path) -> int:
    n = 0
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 24), b''):
            n += chunk.count(b'\n')
    return n


def tsv_to_parquet(tsv: Path, dest: Path, a, required_complete=()) -> dict:
    """TSV (QuPath measurement export) -> zstd Parquet with explicit types, then verified.
    required_complete: numeric columns that must be present, non-null and non-NaN on every row."""
    import duckdb   # venv-v3; imported here so the rest of v3ctl.py runs on a bare python3
    with open(tsv, encoding='utf-8') as f:
        header = f.readline().rstrip('\n').split('\t')
    if len(set(header)) != len(header):
        raise ValueError(f'duplicate column names in {tsv}')
    types = {c: ('VARCHAR' if c in STRING_COLUMNS else 'DOUBLE') for c in header}
    con = duckdb.connect()
    con.execute(f'SET threads={int(a.duckdb_threads)}')
    con.execute(f"SET memory_limit='{int(a.duckdb_mem_gb)}GB'")
    con.execute(f"SET temp_directory='{Path(a.tmp) / 'duckdb'}'")
    cols = '{' + ', '.join(f"'{c.replace(chr(39), chr(39) * 2)}': '{t}'" for c, t in types.items()) + '}'
    src = (f"read_csv('{tsv}', delim='\\t', header=true, auto_detect=false, columns={cols}, "
           f"quote='', escape='', nullstr='', strict_mode=true)")
    part = dest.with_name(dest.name + '.part')
    part.unlink(missing_ok=True)
    con.execute(f"COPY (SELECT * FROM {src}) TO '{part}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    n, n_id, n_distinct = con.execute(
        f"""SELECT count(*), count("Object ID"), count(DISTINCT "Object ID") FROM read_parquet('{part}')""").fetchone()
    complete, debris = {}, {}
    for c in required_complete:
        if c not in header:
            complete[c] = None
            continue
        q = c.replace('"', '""')
        complete[c] = con.execute(f"""SELECT count(*) FILTER (WHERE "{q}" IS NOT NULL AND NOT isnan("{q}"))
                                      FROM read_parquet('{part}')""").fetchone()[0]
        if DEBRIS_AREA_COLUMN in header:
            aq = DEBRIS_AREA_COLUMN.replace('"', '""')
            debris[c] = con.execute(f"""SELECT count(*) FILTER (WHERE ("{q}" IS NULL OR isnan("{q}"))
                                                                AND "{aq}" < {DEBRIS_MAX_UM2})
                                        FROM read_parquet('{part}')""").fetchone()[0]
    con.close()
    return {'columns': len(header), 'rows': n, 'object_id_non_null': n_id, 'object_id_distinct': n_distinct,
            'complete_counts': complete, 'debris_incomplete': debris, 'part': part, 'string_columns': [c for c in header if types[c] == 'VARCHAR']}


def cmd_post(a):
    repo, project, out, steps = common(a)
    site = resolve_site(a.site)
    st = steps[a.step]
    im = find_image(project, a.image)
    pdir = site_dir(out, im.tag, site)
    rec = read_json(pdir / f'{st.id}.json')
    t0 = time.time()
    result = {'schema': 'run_v3.post/1', 'step': st.id, 'hook': st.post, 'attempt_id': (rec or {}).get('attempt_id'),
              'seq': (rec or {}).get('seq'), 'run_id': a.run_id, 'image': im.name, 'tag': im.tag, 'site': site,
              'status': 'fail', 'started': now_iso(), 'outputs': [], 'checks': {}}
    try:
        if not (isinstance(rec, dict) and rec.get('status') == 'pass' and rec.get('attempt_id') == a.attempt):
            raise ValueError('the QuPath part of this attempt did not pass')
        tsv = Path(rec['transient_outputs'][0]['path'])
        expected = int(rec['counts']['detections'])
        sub = {'cells_parquet': 'cells', 'context_parquet': 'context'}[st.post]
        dest = out / sub / f'{im.tag}.parquet'
        dest.parent.mkdir(parents=True, exist_ok=True)
        tsv_rows = count_lines(tsv) - 1
        required = CELLS_REQUIRED_COMPLETE if st.post == 'cells_parquet' else ()
        conv = tsv_to_parquet(tsv, dest, a, required)
        part = conv.pop('part')
        checks = {'expected_rows (detections)': expected, 'tsv_rows': tsv_rows, **conv}
        problems = []
        if not (conv['rows'] == tsv_rows == expected):
            problems.append(f"row counts differ: parquet {conv['rows']}, tsv {tsv_rows}, detections {expected}")
        if conv['object_id_non_null'] != conv['rows']:
            problems.append(f"{conv['rows'] - conv['object_id_non_null']} rows without Object ID")
        if conv['object_id_distinct'] != conv['rows']:
            problems.append(f"Object ID not unique ({conv['object_id_distinct']} distinct of {conv['rows']})")
        for c, n in conv['complete_counts'].items():
            if n is None:
                problems.append(f"column '{c}' is missing (did step 2 run?)")
            elif n + conv['debris_incomplete'].get(c, 0) != conv['rows']:
                bad = conv['rows'] - n - conv['debris_incomplete'].get(c, 0)
                problems.append(f"'{c}' is null/NaN on {bad} of {conv['rows']} rows of cells >= {DEBRIS_MAX_UM2} um^2 "
                                f"(step 2 did not measure every detection)")
        if st.post == 'context_parquet':
            cells = out / 'cells' / f'{im.tag}.parquet'
            if cells.exists():
                import duckdb
                con = duckdb.connect()
                only_ctx, only_cells = con.execute(f"""
                    SELECT (SELECT count(*) FROM read_parquet('{part}') c ANTI JOIN read_parquet('{cells}') p USING ("Object ID")),
                           (SELECT count(*) FROM read_parquet('{cells}') p ANTI JOIN read_parquet('{part}') c USING ("Object ID"))
                """).fetchone()
                con.close()
                checks['object_ids_vs_cells_parquet'] = {'only_in_context': only_ctx, 'only_in_cells': only_cells}
                if only_ctx or only_cells:
                    problems.append(f'Object IDs differ from cells/{im.tag}.parquet (context-only {only_ctx}, cells-only {only_cells})')
        result['checks'] = checks
        if problems:
            part.unlink(missing_ok=True)
            raise ValueError('; '.join(problems))
        os.replace(part, dest)
        result['outputs'].append(stat_of(dest))
        if st.post == 'cells_parquet' and re.search(a.keep_tsv, im.tag):
            kept = out / 'cells_tsv' / f'{im.tag}.tsv'
            kept.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(tsv), kept)
            result['outputs'].append(stat_of(kept))
            result['tsv'] = f'kept at {kept} (same-slide lane reads TSV)'
        else:
            tsv.unlink()
            result['tsv'] = f'deleted {tsv}'
        result['status'] = 'pass'
        msg = f"{conv['rows']} rows x {conv['columns']} columns -> {dest}"
    except Exception as e:  # recorded, reported, and the step stays not-done
        result['error'] = f'{type(e).__name__}: {e}'
        msg = result['error']
    result['seconds'] = round(time.time() - t0, 1)
    result['finished'] = now_iso()
    write_json(pdir / f'{st.id}.post.json', result)
    append_run_log(out, {'time_start': result['started'], 'time_end': result['finished'], 'run_id': a.run_id,
                         'host': socket.gethostname(), 'image': im.name, 'tag': im.tag, 'step': f'{st.id}:post',
                         'action': 'post', 'status': result['status'].upper(), 'seconds': result['seconds'],
                         'attempt_id': a.attempt, 'seq': result['seq'], 'lane': st.lane, 'site': site,
                         'message': msg[:500]})
    print(f"{result['status'].upper()}\t{msg}")
    sys.exit(0 if result['status'] == 'pass' else 1)


def lock_dirs(a) -> tuple[Path, Path | None]:
    lockdir = Path(a.lockdir or os.environ.get('V3_LOCKDIR') or DEFAULT_LOCKDIR)
    shared = a.shared_lockdir if a.shared_lockdir is not None else os.environ.get('V3_SHARED_LOCKDIR', '')
    return lockdir, (Path(shared) if shared else None)


def cmd_accept_external(a):
    _, project, out, steps = common(a)
    site = resolve_site(a.site)
    rerun = [s for s in (a.rerun_from or '').split(',') if s]
    unknown = [s for s in rerun if s not in steps]
    if unknown:
        raise SystemExit(f'--rerun-from: unknown step(s) {unknown}; known: {list(steps)}')
    images = select_images(load_images(project), a.images, a.images_file, '--images')
    print_resolved(images, sys.stdout)
    stale_always = [s for s, st in steps.items() if st.on_external == 'stale']
    invalid = sorted(downstream_closure(steps, rerun)) if rerun else []
    lockdir, shared = lock_dirs(a)
    rc = 0
    with image_locks(images, lockdir, shared, a.run_id):
        for im in images:
            tagdir = prov_dir(out, im.tag)
            recs = load_records(tagdir, site)
            cur = stat_of(im.qpdata)
            state, msg, ref = image_state(recs, cur, allow_drift=True)
            if state != 'custody_broken':
                print(f'{im.name}: {state}; nothing to accept' + ('' if state in ('ok', 'fresh') else f' -- {msg}'))
                rc = rc or (0 if state in ('ok', 'fresh') else 1)
                continue
            seq = recs.max_seq() + 1
            write_json(site_dir(out, im.tag, site) / f'external.{seq}.json', {
                'schema': 'run_v3.provenance/1', 'kind': 'external', 'step': 'external', 'status': 'pass', 'seq': seq,
                'image': im.name, 'tag': im.tag, 'note': a.note, 'user': getpass.getuser(), 'host': socket.gethostname(),
                'site': site, 'time': now_iso(), 'rerun_from': rerun, 'invalidates': invalid,
                'stale_on_external': stale_always,
                'qpdata': {'before': ref.data['qpdata']['after'], 'after': cur}})
            print(f'{im.name}: accepted the current data.qpdata as seq {seq} ({a.note}). Now stale: '
                  f'{", ".join(sorted(set(stale_always) | set(invalid)))} (steps.tsv on_external=stale'
                  + (f', plus --rerun-from {",".join(rerun)} and its downstream' if rerun else '') + ')')
    sys.exit(rc)


def cmd_merge_provenance(a):
    _, project, out, _ = common(a)
    site = resolve_site(a.site)
    images = select_images(load_images(project), a.images, a.images_file, '--images')
    print_resolved(images, sys.stdout)
    lockdir, shared = lock_dirs(a)
    bad = 0
    with image_locks(images, lockdir, shared, a.run_id):
        for im in images:
            tagdir = prov_dir(out, im.tag)
            others = sorted(d.name for d in tagdir.iterdir() if d.is_dir() and d.name != site
                            and not d.name.startswith('.')) if tagdir.is_dir() else []
            if a.from_site:
                if a.from_site not in others:
                    print(f'{im.name}: no records from site {a.from_site!r} in {tagdir}')
                    continue
                others = [a.from_site]
            if not others:
                print(f'{im.name}: no other site has records; nothing to merge')
                continue
            staged, changes = {}, {}
            mroot = site_dir(out, im.tag, site) / 'merged'
            try:
                for s2 in others:
                    native = native_files(tagdir / s2)
                    old = native_files(mroot / s2)
                    changes[s2] = {'new': sorted(set(native) - set(old)),
                                   'changed': sorted(p for p in set(native) & set(old)
                                                     if not filecmp.cmp(native[p], old[p], shallow=False)),
                                   'removed': sorted(set(old) - set(native))}
                    tmp = mroot / f'.staging.{s2}.{os.getpid()}'
                    shutil.rmtree(tmp, ignore_errors=True)
                    for rel, f in native.items():
                        (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(f, tmp / rel)
                    tmp.mkdir(parents=True, exist_ok=True)
                    staged[s2] = tmp
                recs = load_records(tagdir, site, merged_override=staged)
                # judge the candidate on its own: other sites' unmerged records are not part of this merge
                recs.unmerged = []
                recs.foreign_running = [x for x in recs.foreign_running if x[0] in staged]
                cur = stat_of(im.qpdata)
                state, msg, ref = image_state(recs, cur, allow_drift=True)
                summary = '; '.join(f"{s2}: +{len(c['new'])} new, {len(c['changed'])} changed, {len(c['removed'])} removed"
                                    for s2, c in changes.items())
                if state not in ('ok', 'fresh'):
                    bad += 1
                    print(f'{im.name}: NOT merged ({summary}) -- {state}: {msg}')
                    continue
                if a.check:
                    print(f'{im.name}: would merge ({summary}); custody {state}: {msg}')
                    continue
                for s2, tmp in staged.items():
                    dest = mroot / s2
                    trash = mroot / f'.replaced.{s2}.{os.getpid()}'
                    if dest.exists():
                        os.rename(dest, trash)
                    os.rename(tmp, dest)
                    shutil.rmtree(trash, ignore_errors=True)
                staged = {}
                print(f'{im.name}: merged ({summary}); custody {state}: {msg}')
            finally:
                for tmp in staged.values():
                    shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(1 if bad else 0)


def cmd_show_record(a):
    _, project, out, steps = common(a)
    site = resolve_site(a.site)
    im = find_image(project, a.image)
    recs = load_records(prov_dir(out, im.tag), site)
    r = recs.canonical(a.step)
    if r is None or (a.require_pass and r.data.get('status') != 'pass'):
        sys.exit(1)
    val = r.data
    for k in [k for k in (a.field or '').split('.') if k]:
        val = val.get(k) if isinstance(val, dict) else None
    print(val if isinstance(val, str) else json.dumps(val, indent=2, ensure_ascii=False))


def cmd_lock(a):
    _, shared = lock_dirs(a)
    if shared is None:
        raise SystemExit('lock: no shared lock directory (--shared-lockdir / V3_SHARED_LOCKDIR)')
    if a.action == 'acquire':
        ok, msg = shared_lock_acquire(shared, a.tag, lock_owner(a.pid, a.run_id))
    else:
        ok, msg = shared_lock_release(shared, a.tag, a.pid)
    print(msg)
    sys.exit(0 if ok else 3)


def cmd_runlog(a):
    """A run_log row for an outcome without an attempt: blocked, locked, stopped, prepare failed."""
    _, project, out, _ = common(a)
    im = find_image(project, a.image)
    t = now_iso()
    append_run_log(out, {'time_start': t, 'time_end': t, 'run_id': a.run_id, 'host': socket.gethostname(),
                         'image': im.name, 'tag': im.tag, 'step': a.step, 'action': a.action,
                         'status': a.status, 'message': a.message[:500]})


def cmd_report(a):
    _, project, out, _ = common(a)
    rows = []
    path = out / 'run_log.csv'
    if path.exists():
        with open(path, encoding='utf-8') as f:
            rows = [r for r in csv.DictReader(f) if r['run_id'] == a.run_id]
    print(f'== run {a.run_id}: {len(rows)} logged step result(s)')
    for r in rows:
        print(f"  {r['status']:<7} {r['tag'][:44]:<44} {r['step']:<14} {r['seconds']:>7}s  "
              f"rss={r['peak_rss_kb'] or '-':>9}kB  {r['message'][:110]}")
    images = {im.tag: im for im in load_images(project)}
    work = Path(a.workroot) / a.run_id
    n_diff = 0
    if work.is_dir():
        for wp in sorted(work.glob('*/*/project.qpproj')):
            tag = wp.parent.parent.name
            im = images.get(tag)
            if im is None:
                continue
            we = next((e for e in (read_json(wp) or {}).get('images', []) if e.get('entryID') == im.entry_id), None)
            diffs = entry_diff(im.entry, we) if we else [{'path': '<entry missing>'}]
            if diffs:
                n_diff += 1
                print(f'  WORKER ENTRY DIFFERS: {wp}')
                for d in diffs[:10]:
                    print(f"      {d['path']}: project={str(d.get('project'))[:60]} worker={str(d.get('now'))[:60]}")
    drift = sorted(p for p in (out / 'provenance').rglob('*.drift.json') if 'history' not in p.parts) \
        if (out / 'provenance').is_dir() else []
    for d in drift:
        print(f'  DRIFT RECORD: {d}')
    print(f'  worker project entries differing from {project.name}: {n_diff}; drift records: {len(drift)} '
          f'(reported only -- nothing is merged into the project)')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--repo', default='/home/smith6jt/IO60panc2nd', help='root that steps.tsv script paths are relative to')
    ap.add_argument('--project', default='/home/smith6jt/IO60panc2nd/project.qpproj')
    ap.add_argument('--out', default='/home/smith6jt/IO60panc2nd/results/v3')
    ap.add_argument('--steps-file', default=str(HERE / 'steps.tsv'))
    ap.add_argument('--script-set', choices=['v3', 'orig', 'auto'], default='v3',
                    help='v3 (default) requires the v3 copy of every step that has one; orig forces the originals; '
                         'auto takes the v3 copy when it exists')
    ap.add_argument('--script', action='append', help='STEP=PATH: run PATH instead of the table script')
    ap.add_argument('--heap', action='append', help='STEP=GB')
    ap.add_argument('--override', action='append', help='STEP:NAME=VALUE, added to the steps.tsv overrides')
    ap.add_argument('--run-id', default='manual')
    ap.add_argument('--site', default=None, help='provenance site (default: $V3_SITE or the short host name)')
    sub = ap.add_subparsers(dest='cmd', required=True)

    sub.add_parser('steps')
    p = sub.add_parser('images')
    p.add_argument('--images')
    p.add_argument('--images-file')

    p = sub.add_parser('plan')
    p.add_argument('--passes')
    p.add_argument('--steps')
    p.add_argument('--images')
    p.add_argument('--images-file')
    p.add_argument('--force')
    p.add_argument('--rerun-from')
    p.add_argument('--allow-drift', action='store_true')
    p.add_argument('--format', choices=['tsv', 'table'], default='table')

    p = sub.add_parser('budget')
    p.add_argument('--passes')
    p.add_argument('--steps')
    p.add_argument('--margin-min', type=int, default=60)
    p.add_argument('--format', choices=['min', 'slurm'], default='min')

    p = sub.add_parser('diskcheck')
    p.add_argument('--image', required=True)
    p.add_argument('--data-margin-gb', type=float, default=50)
    p.add_argument('--min-root-free-gb', type=float, default=300)
    p.add_argument('--root-path', default='/')

    p = sub.add_parser('prepare')
    for x in ('--image', '--step', '--workdir', '--tmp'):
        p.add_argument(x, required=True)
    p.add_argument('--slot', default='')
    p.add_argument('--gpu', default='')
    p.add_argument('--deadline-epoch', default='', help='seconds since the epoch; the footer will not start a save '
                                                         'later than --save-margin-min before it')
    p.add_argument('--save-margin-min', type=float, default=10.0)

    p = sub.add_parser('finish')
    for x in ('--image', '--step', '--attempt', '--workdir'):
        p.add_argument(x, required=True)
    p.add_argument('--exit-code', default='')
    p.add_argument('--wall', default='')
    p.add_argument('--timefile')
    p.add_argument('--log', default='')
    p.add_argument('--slot', default='')
    p.add_argument('--allow-drift', action='store_true')

    p = sub.add_parser('post')
    for x in ('--image', '--step', '--attempt', '--tmp'):
        p.add_argument(x, required=True)
    p.add_argument('--keep-tsv', default=KEEP_TSV_DEFAULT)
    p.add_argument('--duckdb-threads', type=int, default=8)
    p.add_argument('--duckdb-mem-gb', type=int, default=48)

    for name in ('accept-external', 'merge-provenance'):
        p = sub.add_parser(name)
        p.add_argument('--images')
        p.add_argument('--images-file')
        p.add_argument('--lockdir', default=None)
        p.add_argument('--shared-lockdir', default=None)
        if name == 'accept-external':
            p.add_argument('--note', required=True)
            p.add_argument('--rerun-from', default='', help='STEP[,STEP]: the edit also invalidates these steps and '
                                                            'everything downstream of them')
        else:
            p.add_argument('--from', dest='from_site', default=None, help='merge only this site')
            p.add_argument('--check', action='store_true', help='validate only; change nothing')

    p = sub.add_parser('show-record')
    p.add_argument('--image', required=True)
    p.add_argument('--step', required=True)
    p.add_argument('--field', default='')
    p.add_argument('--require-pass', action='store_true')

    p = sub.add_parser('lock')
    p.add_argument('action', choices=['acquire', 'release'])
    p.add_argument('--tag', required=True)
    p.add_argument('--pid', type=int, required=True)
    p.add_argument('--lockdir', default=None)
    p.add_argument('--shared-lockdir', default=None)

    p = sub.add_parser('runlog')
    for x in ('--image', '--step', '--action', '--status'):
        p.add_argument(x, required=True)
    p.add_argument('--message', default='')

    p = sub.add_parser('report')
    p.add_argument('--workroot', required=True)

    a = ap.parse_args(argv)
    if a.cmd in ('accept-external', 'merge-provenance') and not (a.images or a.images_file):
        ap.error(f'{a.cmd} needs --images or --images-file')
    try:
        {'steps': cmd_steps, 'images': cmd_images, 'plan': cmd_plan, 'budget': cmd_budget, 'diskcheck': cmd_diskcheck,
         'prepare': cmd_prepare, 'finish': cmd_finish, 'post': cmd_post, 'accept-external': cmd_accept_external,
         'merge-provenance': cmd_merge_provenance, 'show-record': cmd_show_record, 'lock': cmd_lock,
         'runlog': cmd_runlog, 'report': cmd_report}[a.cmd](a)
    except PlanError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
