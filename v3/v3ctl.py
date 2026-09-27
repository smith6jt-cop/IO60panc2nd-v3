#!/usr/bin/env python3
"""v3ctl.py -- the decisions behind run_v3.sh: plan, instantiate, validate, convert, record.

run_v3.sh owns processes, locks and lanes; everything that reads or writes JSON lives here.
Standard library only, except `post`, which imports duckdb lazily (venv-v3), so planning and
validation also work with a bare python3 (e.g. on HiPerGator).

Subcommands (run_v3.sh calls them; they are also usable by hand):
  steps            print the resolved step table
  plan             what would run / be skipped / is blocked, per image (tsv for the driver, or a table)
  prepare          instantiate one step for one image: worker project + script copy + footer
  finish           validate the sentinel of an attempt, copy the script, check drift, log the row
  post             run a step's post hook (TSV -> Parquet + verification)
  diskcheck        the free-space guard
  accept-external  record that the current data.qpdata (changed outside run_v3) is accepted
  runlog           append a run_log.csv row for an outcome without an attempt (blocked, locked...)
  report           summary of one run (statuses, drift, worker-vs-project entry differences)

Provenance layout, per image tag, under <out>/provenance/<tag>/:
  <step>.json         written by the footer inside QuPath: the success sentinel + provenance
  <step>.post.json    written by `post` (conversion + verification)
  <step>.groovy       the exact script that ran (step body + overrides + footer)
  <step>.drift.json   only if the step changed the image server metadata or the project entry
  external.<seq>.json accepted modifications made outside run_v3 (e.g. gate-1 review in the GUI)
  history/            superseded records of earlier attempts
  .running.<step>.json marker while an attempt is in flight
Every record carries a per-image sequence number `seq`; a step is stale when a step it depends on
('after' in steps.tsv) has a larger seq.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
import getpass
import hashlib
import json
import os
import re
import shlex
import shutil
import socket
import sys
import time
import uuid
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOOTER_TEMPLATE = HERE / 'footer_template.groovy'
CFG_PLACEHOLDER = '/*__V3_CFG__*/ null'
SAVE_POLICIES = {'self-saves', 'footer-saves', 'never'}
LANES = {'gpu', 'cpu', 'light'}
POST_HOOKS = {'cells_parquet', 'context_parquet'}
POSTCONDITIONS = {'preflight', 'v3_result', 'detections_present', 'islets_present', 'islet_bands',
                  'cells_measured', 'signed_distance', 'parenting'}
# Custody: size must match exactly, mtime within 1 s -- Globus/rsync preserve whole seconds only,
# so a .qpdata that went to HiPerGator and back still matches its record.
MTIME_TOL_NS = 1_000_000_000
KEEP_TSV_DEFAULT = r'0080870_6533|0081366_HDL112'   # same-slide Xenium images: sameslide/ reads the TSV
LAST_DEFAULT = ('6457_',)                            # plan: 6457 is re-run last
# measurement-table columns that hold text; every other column is numeric (DOUBLE)
STRING_COLUMNS = {'Image', 'Object ID', 'Object type', 'Name', 'Class', 'Classification', 'Parent',
                  'Parent ID', 'ROI', 'TMA core'}
RUN_LOG_FIELDS = ['time_start', 'time_end', 'run_id', 'host', 'image', 'tag', 'step', 'action', 'status',
                  'seconds', 'exit_code', 'attempt_id', 'seq', 'lane', 'slot', 'heap_gb', 'peak_rss_kb',
                  'detections', 'qpdata_bytes', 'script', 'script_sha256', 'message']


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


def load_steps(path: Path) -> dict:
    lines = [ln for ln in path.read_text(encoding='utf-8').splitlines() if ln.strip() and not ln.startswith('#')]
    header = lines[0].split('\t')
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
                  needs=[d.strip() for d in r['needs'].split(',') if d.strip()], overrides=overrides,
                  postcondition=parse_postcondition(r['postcondition']), post=r['post'] or None,
                  description=r.get('description', ''))
        if st.id in steps:
            raise SystemExit(f'{path}: duplicate step {st.id!r}')
        if not re.fullmatch(r'[A-Za-z0-9_]+', st.id):
            raise SystemExit(f'{path}: step id {st.id!r} must be [A-Za-z0-9_]+')
        if st.save not in SAVE_POLICIES or st.lane not in LANES:
            raise SystemExit(f'{path}: step {st.id}: bad save policy {st.save!r} or lane {st.lane!r}')
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
        raise PlanError(f'step {st.id}: script {chosen} does not exist')
    return chosen, how


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
    return [Image(e['imageName'], e['entryID'], tag_of(e['imageName']),
                  base / 'data' / str(e['entryID']) / 'data.qpdata', e) for e in proj['images']]


def select_images(images: list, spec: str | None) -> list:
    if not spec:
        first = [im for im in images if not im.name.startswith(LAST_DEFAULT)]
        return first + [im for im in images if im.name.startswith(LAST_DEFAULT)]
    chosen = []
    for pat in (p.strip() for p in spec.split(',') if p.strip()):
        hits = [im for im in images if pat in (im.name, im.tag, str(im.entry_id))] or \
               [im for im in images if pat in im.name]
        if not hits:
            raise SystemExit(f'--images: nothing matches {pat!r}')
        chosen += [h for h in hits if h not in chosen]
    return chosen


# ----------------------------------------------------------------------------------------------
# records and planning
@dataclass
class Records:
    canonical: dict
    post: dict
    drift: dict
    running: dict
    history: list
    external: list

    def seq_records(self):
        return [r for r in list(self.canonical.values()) + self.history + self.external
                if isinstance(r, dict) and isinstance(r.get('seq'), int)]

    def max_seq(self) -> int:
        return max((r['seq'] for r in self.seq_records()), default=0)

    def ever_ran(self, step: str) -> bool:
        return step in self.canonical or any(r.get('step') == step for r in self.history if isinstance(r, dict))


def prov_dir(out: Path, tag: str) -> Path:
    return out / 'provenance' / tag


def load_records(pdir: Path) -> Records:
    rec = Records({}, {}, {}, {}, [], [])
    if not pdir.is_dir():
        return rec
    for f in pdir.iterdir():
        n = f.name
        if n.endswith('.tmp') or not n.endswith('.json'):
            continue
        if m := re.fullmatch(r'([A-Za-z0-9_]+)\.json', n):
            rec.canonical[m.group(1)] = read_json(f)
        elif m := re.fullmatch(r'([A-Za-z0-9_]+)\.post\.json', n):
            rec.post[m.group(1)] = read_json(f)
        elif m := re.fullmatch(r'([A-Za-z0-9_]+)\.drift\.json', n):
            rec.drift[m.group(1)] = f
        elif m := re.fullmatch(r'\.running\.([A-Za-z0-9_]+)\.json', n):
            rec.running[m.group(1)] = read_json(f)
        elif re.fullmatch(r'external\.\d+\.json', n):
            rec.external.append(read_json(f))
    hist = pdir / 'history'
    if hist.is_dir():
        for f in hist.glob('*.json'):
            r = read_json(f)
            if isinstance(r, dict) and r.get('schema', '').startswith('run_v3.provenance'):
                rec.history.append(r)
    return rec


def custody(recs: Records, current: dict) -> tuple[str, dict | None, str]:
    """Is data.qpdata still the file the last recorded step left behind?"""
    cands = [r for r in recs.seq_records() if (r.get('qpdata') or {}).get('after')]
    if not cands:
        return 'fresh', None, 'no run_v3 record yet'
    ref = max(cands, key=lambda r: r['seq'])
    after = ref['qpdata']['after']
    if same_stat(after, current):
        return 'ok', ref, f"matches seq {ref['seq']} ({ref.get('step', ref.get('kind'))})"
    running = ', '.join(f"{s} (attempt {m.get('attempt_id')}, started {m.get('started')})"
                        for s, m in recs.running.items() if isinstance(m, dict))
    msg = (f"data.qpdata changed outside run_v3: last record seq {ref['seq']} ({ref.get('step', ref.get('kind'))}) "
           f"left size={after.get('size')} mtime_ns={after.get('mtime_ns')}, now size={current.get('size')} "
           f"mtime_ns={current.get('mtime_ns')}.")
    if running:
        msg += f' An attempt was interrupted: {running}; every step is an idempotent re-run.'
    msg += (' If the change is intended (e.g. gate-1 review in the GUI) run: run_v3.sh --accept-external '
            '--images <image> --note "why"; otherwise restore the file (a data.qpdata.bkp may sit next to it).')
    return 'broken', ref, msg


class Ctx:
    def __init__(self, steps: dict, recs: Records):
        self.steps, self.recs, self.memo = steps, recs, {}

    def done(self, sid: str) -> tuple[bool, str]:
        if sid in self.memo:
            return self.memo[sid]
        self.memo[sid] = (False, 'dependency cycle')      # guard; replaced below
        st, rec = self.steps[sid], self.recs.canonical.get(sid)
        res = (True, f"done (seq {rec.get('seq')})") if isinstance(rec, dict) else (False, 'never ran')
        if isinstance(rec, dict) and rec.get('status') != 'pass':
            pc = rec.get('postcondition') or {}
            res = (False, f"last attempt failed: {pc.get('message') or ''} {pc.get('error') or ''} "
                          f"{(rec.get('save') or {}).get('policy_error') or ''}".strip())
        if res[0] and st.post:
            prec = self.recs.post.get(sid)
            if not (isinstance(prec, dict) and prec.get('status') == 'pass' and prec.get('attempt_id') == rec.get('attempt_id')):
                res = (False, f'post hook {st.post} not done')
        if res[0]:
            outs = list(rec.get('outputs') or [])
            if st.post:
                outs += list((self.recs.post.get(sid) or {}).get('outputs') or [])
            for o in outs:
                s = stat_of(Path(o['path']))
                if not s['exists'] or s.get('size') != o.get('size'):
                    res = (False, f"output missing or changed: {o['path']}")
                    break
        if res[0]:
            for dep, optional in st.after:
                if optional and not self.recs.ever_ran(dep):
                    continue
                ok, why = self.done(dep)
                if not ok:
                    res = (False, f'upstream {dep} not done ({why})')
                    break
                if self.recs.canonical[dep]['seq'] > rec['seq']:
                    res = (False, f"stale: {dep} re-ran after it (seq {self.recs.canonical[dep]['seq']} > {rec['seq']})")
                    break
        if res[0]:
            for dep in st.needs:
                ok, why = self.done(dep)
                if not ok:
                    res = (False, f'needs {dep} ({why})')
                    break
        self.memo[sid] = res
        return res


def upstream_closure(steps: dict, sid: str) -> set:
    seen, todo = set(), [d for d, _ in steps[sid].after]
    while todo:
        d = todo.pop()
        if d not in seen:
            seen.add(d)
            todo += [x for x, _ in steps[d].after]
    return seen


def selected_steps(steps: dict, passes: list, only: list) -> list:
    """The steps to consider, in dependency order (table order breaks ties)."""
    if only:
        bad = [s for s in only if s not in steps]
        if bad:
            raise SystemExit(f'--steps: unknown step(s) {bad}; known: {list(steps)}')
        chosen = [s for s in steps if s in only]
    elif passes:
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


def plan_image(im: Image, steps: dict, sel: list, force: set, out: Path, allow_drift: bool) -> dict:
    recs = load_records(prov_dir(out, im.tag))
    cur = stat_of(im.qpdata)
    state, ref, msg = custody(recs, cur)
    rows = []
    if state == 'broken':
        return {'image': im, 'state': 'custody_broken', 'message': msg, 'qpdata': cur,
                'steps': [(s, 'blocked', 'custody broken (see image message)') for s in sel]}
    if recs.drift and not allow_drift:
        files = ', '.join(str(p) for p in recs.drift.values())
        return {'image': im, 'state': 'drift', 'qpdata': cur,
                'message': f'unresolved server-metadata / project-entry drift: {files}. Merge the change into the '
                           f'project in the GUI (or confirm it is unwanted), delete the drift file(s), then re-run '
                           f'(or pass --allow-drift).',
                'steps': [(s, 'blocked', 'unresolved drift') for s in sel]}
    ctx = Ctx(steps, recs)
    run, why = {}, {}
    for s in sel:
        ok, reason = ctx.done(s)
        if s in force:
            run[s], why[s] = True, 'forced (--force)'
        elif not ok:
            run[s], why[s] = True, reason
        else:
            run[s], why[s] = False, reason
    for s in sel:        # a step that re-runs makes everything downstream of it (in the selection) stale
        if not run[s]:
            ups = [u for u in upstream_closure(steps, s) if run.get(u)]
            if ups:
                run[s], why[s] = True, f'downstream of re-running {",".join(sorted(ups))}'
    blocked = None
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
    message = msg + (f'; stale in-flight marker(s): {", ".join(stale)}' if stale else '')
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
# run log
def append_run_log(out: Path, row: dict) -> None:
    path = out / 'run_log.csv'
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'a', newline='', encoding='utf-8') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            if f.tell() == 0:
                csv.writer(f).writerow(RUN_LOG_FIELDS)
            csv.writer(f).writerow([row.get(k, '') for k in RUN_LOG_FIELDS])
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
    for im in load_images(project):
        if name in (im.name, im.tag):
            return im
    raise SystemExit(f'image {name!r} not in {project}')


def cmd_steps(a):
    repo, _, _, steps = common(a)
    so = parse_kv_list(a.script, '--script')
    print('\t'.join(['step', 'passes', 'lane', 'heap_gb', 'save', 'script', 'selection', 'postcondition', 'post']))
    for st in steps.values():
        try:
            p, how = resolve_script(st, repo, a.script_set, so)
        except PlanError as e:
            p, how = f'ERROR: {e}', '-'
        print('\t'.join([st.id, ','.join(st.passes), st.lane, str(st.heap_gb), st.save, str(p), how,
                         st.postcondition['name'], st.post or '-']))


def cmd_plan(a):
    repo, project, out, steps = common(a)
    so = parse_kv_list(a.script, '--script')
    heaps = {k: int(v) for k, v in parse_kv_list(a.heap, '--heap').items()}
    sel = selected_steps(steps, [p for p in (a.passes or '').split(',') if p], [s for s in (a.steps or '').split(',') if s])
    force = {s for s in (a.force or '').split(',') if s}
    if force - set(steps):
        raise SystemExit(f'--force: unknown step(s) {sorted(force - set(steps))}')
    images = select_images(load_images(project), a.images)
    scripts = {}
    for s in sel:
        try:
            scripts[s] = resolve_script(steps[s], repo, a.script_set, so)
        except PlanError as e:
            scripts[s] = (f'ERROR: {e}', 'error')
    plans = [plan_image(im, steps, sel, force, out, a.allow_drift) for im in images]
    if a.format == 'tsv':
        for p in plans:
            im = p['image']
            print('\t'.join(['IMAGE', im.tag, im.name, str(im.entry_id), p['state'], p['message'].replace('\t', ' ')]))
            for s, action, reason in p['steps']:
                st = steps[s]
                spath, how = scripts[s]
                if action == 'run' and how == 'error':
                    action, reason = 'blocked', spath
                print('\t'.join(['STEP', im.tag, s, action, st.lane, str(heaps.get(s, st.heap_gb)), str(st.threads),
                                 str(st.timeout_min), st.post or '-', str(spath), reason.replace('\t', ' ')]))
        return
    # human-readable table (dry runs)
    print(f'project : {project}\noutput  : {out}\nsteps   : {", ".join(sel)}  (force: {", ".join(sorted(force)) or "-"})')
    print('scripts :')
    for s in sel:
        spath, how = scripts[s]
        ov = steps[s].overrides
        print(f'  {s:<14} {how:<10} {spath}' + (f'   overrides {ov}' if ov else ''))
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
    st = steps[a.step]
    im = find_image(project, a.image)
    so = parse_kv_list(a.script, '--script')
    extra = {}
    for item in a.override or []:
        sid, _, kv = item.partition(':')
        if sid == st.id:
            extra.update(parse_kv_list([kv], '--override'))
    heap = int(parse_kv_list(a.heap, '--heap').get(st.id, st.heap_gb))
    pdir = prov_dir(out, im.tag)
    pdir.mkdir(parents=True, exist_ok=True)
    recs = load_records(pdir)
    seq = recs.max_seq() + 1
    attempt = uuid.uuid4().hex[:12]

    # resolve and instantiate first: a refused script/override must leave the old records alone
    src, how = resolve_script(st, repo, a.script_set, so)
    body = src.read_text(encoding='utf-8')
    body, applied = apply_overrides(body, {**st.overrides, **extra})

    # superseded records of this step go to history/ (custody still sees them there)
    hist = pdir / 'history'
    old = recs.canonical.get(st.id)
    old_seq = old.get('seq') if isinstance(old, dict) else 'x'
    for suffix in ('.json', '.post.json', '.groovy', '.drift.json'):
        f = pdir / f'{st.id}{suffix}'
        if f.exists():
            hist.mkdir(exist_ok=True)
            os.replace(f, hist / f'{st.id}.seq{old_seq}.{int(time.time())}{suffix}')

    worker = Path(a.workdir)
    worker.mkdir(parents=True, exist_ok=True)
    shutil.copy2(project, worker / 'project.qpproj')        # private copy: QuPath never sees the target file
    for d in ('data', 'classifiers', 'scripts', 'resources'):
        target = project.parent / d
        link = worker / d
        if target.exists() and not link.exists():
            link.symlink_to(target)
    subst = {'{OUT}': str(out), '{TMP}': str(Path(a.tmp)), '{TAG}': im.tag, '{PROJ}': str(repo), '{STEP}': st.id}
    args = []
    for x in st.args:
        for k, v in subst.items():
            x = x.replace(k, v)
        args.append(x)
    script_path = worker / 'step.groovy'
    cfg = {
        'run_id': a.run_id, 'attempt_id': attempt, 'seq': seq, 'step': st.id, 'passes': st.passes,
        'image': im.name, 'tag': im.tag, 'entry_id': im.entry_id, 'project': str(project),
        'worker_project': str(worker / 'project.qpproj'), 'out': str(out),
        'provenance': str(pdir / f'{st.id}.json'), 'save_policy': st.save, 'postcondition': st.postcondition,
        'lane': st.lane, 'slot': a.slot, 'heap_gb': heap, 'threads': st.threads,
        'host': socket.gethostname(), 'gpu': a.gpu or None,
        'source_script': str(src), 'source_sha256': sha256_file(src), 'script_selection': how,
        'body_sha256': hashlib.sha256(body.encode('utf-8')).hexdigest(), 'instantiated_script': str(script_path),
        'overrides': applied, 'args': args, 'footer_template_sha256': sha256_file(FOOTER_TEMPLATE),
        'qpdata_before': stat_of(im.qpdata), 'driver_start': now_iso(),
    }
    footer = FOOTER_TEMPLATE.read_text(encoding='utf-8')
    if footer.count(CFG_PLACEHOLDER) != 1:
        raise SystemExit(f'{FOOTER_TEMPLATE}: expected exactly one {CFG_PLACEHOLDER!r}')
    footer = footer.replace(CFG_PLACEHOLDER, groovy_literal(cfg))
    script_path.write_text(body.rstrip('\n') + '\n' + footer, encoding='utf-8')
    write_json(pdir / f'.running.{st.id}.json', {'attempt_id': attempt, 'seq': seq, 'run_id': a.run_id,
                                                  'started': cfg['driver_start'], 'host': cfg['host'],
                                                  'driver_pid': os.getppid(), 'worker': str(worker),
                                                  'qpdata_before': cfg['qpdata_before']})
    qargs = []
    for x in args:
        qargs += ['-a', x]
    print(f'ATTEMPT={attempt}')
    print(f'SEQ={seq}')
    print(f'HEAP={heap}')
    print(f'WSCRIPT={shlex.quote(str(script_path))}')
    print(f'WPROJECT={shlex.quote(str(worker / "project.qpproj"))}')
    print(f'PROVENANCE={shlex.quote(str(pdir / (st.id + ".json")))}')
    print(f'SCRIPT_SRC={shlex.quote(str(src))}')
    print('QARGS=(' + ' '.join(shlex.quote(x) for x in qargs) + ')')


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


def drift_check(project: Path, im: Image, worker: Path, since_ns: int) -> dict | None:
    """Did the step change what QuPath will load next time? Headless runs without -s never write
    project.qpproj, but saveImageData() rewrites data/<id>/server.json whenever the image's server
    metadata (channel names, pixel size...) changed -- and on the next load QuPath takes the metadata
    from the project entry, so such a change silently reverts. Report both, never merge."""
    found = {}
    wp = read_json(worker / 'project.qpproj')
    if wp is not None:
        we = next((e for e in wp.get('images', []) if e.get('entryID') == im.entry_id), None)
        diffs = entry_diff(im.entry, we) if we is not None else [{'path': '<entry>', 'project': 'present', 'now': 'absent'}]
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


def cmd_finish(a):
    repo, project, out, steps = common(a)
    st = steps[a.step]
    im = find_image(project, a.image)
    pdir = prov_dir(out, im.tag)
    marker = read_json(pdir / f'.running.{st.id}.json') or {}
    rec = read_json(pdir / f'{st.id}.json')
    worker = Path(a.workdir)
    script = worker / 'step.groovy'
    status, reason = 'FAIL', ''
    if rec is None:
        reason = ('no sentinel: the step threw, returned early, was killed or timed out '
                  f'(QuPath exit code {a.exit_code}; see {a.log})')
    elif rec.get('attempt_id') != a.attempt:
        reason = f"sentinel belongs to attempt {rec.get('attempt_id')}, not {a.attempt}"
    elif rec.get('status') != 'pass':
        pc = rec.get('postcondition') or {}
        reason = 'footer: ' + ' | '.join(x for x in [pc.get('message'), pc.get('error'),
                                                     (rec.get('save') or {}).get('policy_error')] if x)
    elif script.exists() and rec['script'].get('instantiated_sha256') != sha256_file(script):
        reason = 'the script QuPath ran is not the one prepared (sha256 mismatch)'
    elif not same_stat(rec['qpdata']['after'], stat_of(im.qpdata)):
        reason = 'data.qpdata changed after the footer recorded it'
    else:
        status, reason = 'PASS', (rec.get('postcondition') or {}).get('message', '')
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
    tinfo = parse_time_file(Path(a.timefile)) if a.timefile else {}
    qp = stat_of(im.qpdata)
    append_run_log(out, {
        'time_start': marker.get('started', ''), 'time_end': now_iso(), 'run_id': a.run_id,
        'host': socket.gethostname(), 'image': im.name, 'tag': im.tag, 'step': st.id, 'action': 'run',
        'status': status, 'seconds': a.wall, 'exit_code': a.exit_code, 'attempt_id': a.attempt,
        'seq': marker.get('seq', ''), 'lane': st.lane, 'slot': a.slot,
        'heap_gb': (rec or {}).get('jvm', {}).get('heap_requested_gb', ''),
        'peak_rss_kb': ((rec or {}).get('jvm') or {}).get('vm_hwm_kb') or tinfo.get('peak_rss_kb', ''),
        'detections': ((rec or {}).get('counts') or {}).get('detections', ''),
        'qpdata_bytes': qp.get('size', ''), 'script': (rec or {}).get('script', {}).get('source_path', ''),
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


def tsv_to_parquet(tsv: Path, dest: Path, a) -> dict:
    """TSV (QuPath measurement export) -> zstd Parquet with explicit types, then verified."""
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
    con.close()
    return {'columns': len(header), 'rows': n, 'object_id_non_null': n_id, 'object_id_distinct': n_distinct,
            'part': part, 'string_columns': [c for c in header if types[c] == 'VARCHAR']}


def cmd_post(a):
    repo, project, out, steps = common(a)
    st = steps[a.step]
    im = find_image(project, a.image)
    pdir = prov_dir(out, im.tag)
    rec = read_json(pdir / f'{st.id}.json')
    t0 = time.time()
    result = {'schema': 'run_v3.post/1', 'step': st.id, 'hook': st.post, 'attempt_id': (rec or {}).get('attempt_id'),
              'seq': (rec or {}).get('seq'), 'run_id': a.run_id, 'image': im.name, 'tag': im.tag,
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
        conv = tsv_to_parquet(tsv, dest, a)
        part = conv.pop('part')
        checks = {'expected_rows (detections)': expected, 'tsv_rows': tsv_rows, **conv}
        problems = []
        if not (conv['rows'] == tsv_rows == expected):
            problems.append(f"row counts differ: parquet {conv['rows']}, tsv {tsv_rows}, detections {expected}")
        if conv['object_id_non_null'] != conv['rows']:
            problems.append(f"{conv['rows'] - conv['object_id_non_null']} rows without Object ID")
        if conv['object_id_distinct'] != conv['rows']:
            problems.append(f"Object ID not unique ({conv['object_id_distinct']} distinct of {conv['rows']})")
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
                         'attempt_id': a.attempt, 'seq': result['seq'], 'lane': st.lane, 'message': msg[:500]})
    print(f"{result['status'].upper()}\t{msg}")
    sys.exit(0 if result['status'] == 'pass' else 1)


def cmd_accept_external(a):
    _, project, out, _ = common(a)
    for im in select_images(load_images(project), a.images):
        pdir = prov_dir(out, im.tag)
        recs = load_records(pdir)
        cur = stat_of(im.qpdata)
        state, ref, msg = custody(recs, cur)
        if state != 'broken':
            print(f'{im.name}: custody {state}; nothing to accept')
            continue
        seq = recs.max_seq() + 1
        write_json(pdir / f'external.{seq}.json', {
            'schema': 'run_v3.provenance/1', 'kind': 'external', 'step': 'external', 'status': 'pass', 'seq': seq,
            'image': im.name, 'tag': im.tag, 'note': a.note, 'user': getpass.getuser(), 'host': socket.gethostname(),
            'time': now_iso(), 'qpdata': {'before': ref['qpdata']['after'], 'after': cur}})
        print(f'{im.name}: accepted the current data.qpdata as seq {seq} ({a.note})')


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
    drift = sorted(out.glob('provenance/*/*.drift.json'))
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
    ap.add_argument('--script-set', choices=['auto', 'orig', 'v3'], default='auto')
    ap.add_argument('--script', action='append', help='STEP=PATH: run PATH instead of the table script')
    ap.add_argument('--heap', action='append', help='STEP=GB')
    ap.add_argument('--override', action='append', help='STEP:NAME=VALUE, added to the steps.tsv overrides')
    ap.add_argument('--run-id', default='manual')
    sub = ap.add_subparsers(dest='cmd', required=True)

    sub.add_parser('steps')
    p = sub.add_parser('plan')
    p.add_argument('--passes')
    p.add_argument('--steps')
    p.add_argument('--images')
    p.add_argument('--force')
    p.add_argument('--allow-drift', action='store_true')
    p.add_argument('--format', choices=['tsv', 'table'], default='table')

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

    p = sub.add_parser('accept-external')
    p.add_argument('--images', required=True)
    p.add_argument('--note', required=True)

    p = sub.add_parser('runlog')
    for x in ('--image', '--step', '--action', '--status'):
        p.add_argument(x, required=True)
    p.add_argument('--message', default='')

    p = sub.add_parser('report')
    p.add_argument('--workroot', required=True)

    a = ap.parse_args(argv)
    try:
        {'steps': cmd_steps, 'plan': cmd_plan, 'diskcheck': cmd_diskcheck, 'prepare': cmd_prepare,
         'finish': cmd_finish, 'post': cmd_post, 'accept-external': cmd_accept_external,
         'runlog': cmd_runlog, 'report': cmd_report}[a.cmd](a)
    except PlanError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
