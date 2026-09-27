"""Unit tests for v3ctl.py (planning, custody, overrides, literals, TSV->Parquet).

    /scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_v3ctl.py

No QuPath and no real data: records and projects are synthesized in tmp_path.
"""
import json
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import v3ctl  # noqa: E402

REPO = Path('/home/smith6jt/IO60panc2nd')
STEPS = v3ctl.load_steps(HERE.parent / 'steps.tsv')


# ---- literals and overrides -------------------------------------------------------------------
def test_groovy_literal_escapes():
    lit = v3ctl.groovy_literal({'a': "it's a \\ path $HOME\nµm", 'n': 2 ** 40, 'x': [True, None, 1.5]})
    assert lit == ("['a': 'it\\'s a \\\\ path $HOME\\n\\u00B5m', 'n': 1099511627776L, "
                   "'x': [true, null, 1.5d]]")
    assert v3ctl.groovy_literal({}) == '[:]'
    assert v3ctl.groovy_literal('\U0001F600') == "'\\uD83D\\uDE00'"     # astral char -> surrogate pair


def test_override_real_step2_script():
    text = (REPO / 'scripts/2_fast_cell_measurements.groovy').read_text()
    new, applied = v3ctl.apply_overrides(text, {'N_THREADS': '24'})
    assert applied == {'N_THREADS': {'from': '18', 'to': '24'}}
    assert 'int N_THREADS        = 24' in new and new.count('= 24') == 1
    assert len(new.splitlines()) == len(text.splitlines())


def test_override_aligned_declaration_and_comment_kept():
    text = 'int    N_THREADS    = 12\nString DEVICE       = "gpu0"      // "cpu" if no CUDA\n'
    new, _ = v3ctl.apply_overrides(text, {'DEVICE': '"cpu"'})
    assert 'String DEVICE       = "cpu"      // "cpu" if no CUDA' in new


@pytest.mark.parametrize('text', ['int X = 1\nint X = 2\n', 'int Y = 1\n', '// X = 3\n'])
def test_override_must_match_exactly_once(text):
    with pytest.raises(v3ctl.PlanError):
        v3ctl.apply_overrides(text, {'X': '5'})


def test_parse_postcondition():
    pc = v3ctl.parse_postcondition('parenting(forbid=Other: Gal3Norm|Other: X; k=v)')
    assert pc == {'name': 'parenting', 'params': {'forbid': ['Other: Gal3Norm', 'Other: X'], 'k': ['v']}}
    assert v3ctl.parse_postcondition('v3_result') == {'name': 'v3_result', 'params': {}}


# ---- step table ---------------------------------------------------------------------------------
def test_steps_table_is_consistent():
    assert list(STEPS)[:3] == ['preflight', '1', '4']
    assert STEPS['2'].overrides == {'N_THREADS': '24'}
    assert STEPS['0'].save == 'self-saves' and STEPS['0'].lane == 'gpu'
    assert ('R_import', True) in STEPS['3'].after


def test_selection_is_topological():
    order = v3ctl.selected_steps(STEPS, ['R', 'J'], [])
    assert order.index('R_import') < order.index('3') < order.index('J_context')
    assert v3ctl.selected_steps(STEPS, ['I', 'C'], []) == ['preflight', '1', '4', '0', '2', 'C_cells', 'C_geojson']


def test_script_selection(tmp_path):
    st = STEPS['1']
    p, how = v3ctl.resolve_script(st, REPO, 'orig', {})
    assert p.name == '1_panc_objects.groovy' and how == 'orig'
    fake = tmp_path / 'x.groovy'
    fake.write_text('println 1')
    assert v3ctl.resolve_script(st, REPO, 'auto', {'1': str(fake)}) == (fake, 'override')


# ---- custody and planning ------------------------------------------------------------------------
def make_image(tmp_path, size=100):
    q = tmp_path / 'proj' / 'data' / '1' / 'data.qpdata'
    q.parent.mkdir(parents=True)
    q.write_bytes(b'x' * size)
    return v3ctl.Image('img', 1, 'img', q, {})


def rec(step, seq, qp, status='pass', **kw):
    s = v3ctl.stat_of(qp)
    r = {'schema': 'run_v3.provenance/1', 'step': step, 'seq': seq, 'status': status, 'attempt_id': f'a{seq}',
         'qpdata': {'before': s, 'after': s}, 'outputs': [], 'postcondition': {'message': 'm'}}
    r.update(kw)
    return r


def write_recs(out, *recs):
    d = v3ctl.prov_dir(out, 'img')
    d.mkdir(parents=True, exist_ok=True)
    for r in recs:
        v3ctl.write_json(d / f"{r['step']}.json", r)


def actions(plan):
    return {s: a for s, a, _ in plan['steps']}


def test_custody_fresh_ok_broken(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    sel = v3ctl.selected_steps(STEPS, ['I'], [])
    assert v3ctl.plan_image(im, STEPS, sel, set(), out, False)['state'] == 'fresh'
    write_recs(out, rec('preflight', 1, im.qpdata))
    assert v3ctl.plan_image(im, STEPS, sel, set(), out, False)['state'] == 'ok'
    ns = os.stat(im.qpdata).st_mtime_ns
    os.utime(im.qpdata, ns=(ns, ns + 500_000_000))            # within the 1 s tolerance (Globus/rsync)
    assert v3ctl.plan_image(im, STEPS, sel, set(), out, False)['state'] == 'ok'
    im.qpdata.write_bytes(b'y' * 101)                          # someone saved the image
    p = v3ctl.plan_image(im, STEPS, sel, set(), out, False)
    assert p['state'] == 'custody_broken' and set(actions(p).values()) == {'blocked'}


def test_done_skip_and_downstream_of_forced(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata))
    sel = v3ctl.selected_steps(STEPS, ['I'], [])
    assert actions(v3ctl.plan_image(im, STEPS, sel, set(), out, False)) == {'preflight': 'skip', '1': 'skip', '4': 'skip'}
    forced = actions(v3ctl.plan_image(im, STEPS, sel, {'1'}, out, False))
    assert forced == {'preflight': 'skip', '1': 'run', '4': 'run'}    # 4 is downstream of 1


def test_stale_when_upstream_reran(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('4', 2, im.qpdata), rec('1', 3, im.qpdata))
    p = v3ctl.plan_image(im, STEPS, v3ctl.selected_steps(STEPS, ['I'], []), set(), out, False)
    assert actions(p)['4'] == 'run' and 'stale' in dict((s, r) for s, _, r in p['steps'])['4']


def test_failed_attempt_reruns_and_blocks_missing_prerequisite(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('0', 2, im.qpdata, status='fail'))
    pc = v3ctl.plan_image(im, STEPS, v3ctl.selected_steps(STEPS, ['C'], []), set(), out, False)
    assert actions(pc)['0'] == 'run' and actions(pc)['2'] == 'run'
    pj = v3ctl.plan_image(im, STEPS, v3ctl.selected_steps(STEPS, ['J'], []), set(), out, False)
    assert actions(pj)['3'] == 'blocked'                        # needs 0 (pass C) and 1 (pass I)
    assert all(a == 'blocked' for s, a in actions(pj).items() if s != '3')


def test_post_hook_and_outputs_are_part_of_done(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    parquet = tmp_path / 'cells.parquet'
    parquet.write_bytes(b'p' * 10)
    base = [rec('preflight', 1, im.qpdata), rec('0', 2, im.qpdata), rec('2', 3, im.qpdata),
            rec('C_cells', 4, im.qpdata), rec('C_geojson', 5, im.qpdata)]
    write_recs(out, *base)
    sel = ['C_cells']
    assert actions(v3ctl.plan_image(im, STEPS, sel, set(), out, False)) == {'C_cells': 'run'}   # no post record
    v3ctl.write_json(v3ctl.prov_dir(out, 'img') / 'C_cells.post.json',
                     {'status': 'pass', 'attempt_id': 'a4', 'outputs': [v3ctl.stat_of(parquet)]})
    assert actions(v3ctl.plan_image(im, STEPS, sel, set(), out, False)) == {'C_cells': 'skip'}
    parquet.unlink()
    assert actions(v3ctl.plan_image(im, STEPS, sel, set(), out, False)) == {'C_cells': 'run'}   # output gone


def test_optional_dependency_only_counts_once_it_ran(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    recs = [rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata),
            rec('0', 4, im.qpdata), rec('3', 5, im.qpdata)]
    write_recs(out, *recs)
    assert actions(v3ctl.plan_image(im, STEPS, ['3'], set(), out, False)) == {'3': 'skip'}
    write_recs(out, rec('R_import', 6, im.qpdata))              # reviewed islets imported after 3 ran
    assert actions(v3ctl.plan_image(im, STEPS, ['3'], set(), out, False)) == {'3': 'run'}


def test_drift_record_blocks_unless_allowed(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata))
    v3ctl.write_json(v3ctl.prov_dir(out, 'img') / '1.drift.json', {'server_json': {}})
    assert v3ctl.plan_image(im, STEPS, ['1'], set(), out, False)['state'] == 'drift'
    assert v3ctl.plan_image(im, STEPS, ['1'], set(), out, True)['state'] == 'ok'


# ---- TSV -> Parquet ---------------------------------------------------------------------------------
class A:
    duckdb_threads, duckdb_mem_gb = 2, 2

    def __init__(self, tmp):
        self.tmp = str(tmp)


def test_tsv_to_parquet_types_nan_and_ids(tmp_path):
    pytest.importorskip('duckdb')
    tsv = tmp_path / 'x.tsv'
    tsv.write_text('Image\tObject ID\tParent\tCentroid X µm\tCell: DAPI: Mean\n'
                   'im\tid-1\tIslet_1 (Islet)\t1.5\tNaN\n'
                   'im\tid-2\tAnnotation (Tissue)\t2.5\t\n', encoding='utf-8')
    res = v3ctl.tsv_to_parquet(tsv, tmp_path / 'x.parquet', A(tmp_path))
    assert (res['rows'], res['object_id_distinct'], res['columns']) == (2, 2, 5)
    import pyarrow.parquet as pq
    t = pq.read_table(res['part']).to_pydict()
    assert t['Centroid X µm'] == [1.5, 2.5]
    assert str(t['Cell: DAPI: Mean'][0]) == 'nan' and t['Cell: DAPI: Mean'][1] is None
    assert res['string_columns'] == ['Image', 'Object ID', 'Parent']


def test_tsv_with_text_in_numeric_column_fails(tmp_path):
    pytest.importorskip('duckdb')
    tsv = tmp_path / 'bad.tsv'
    tsv.write_text('Object ID\tSomething new\nid-1\thello\n', encoding='utf-8')
    with pytest.raises(Exception):
        v3ctl.tsv_to_parquet(tsv, tmp_path / 'bad.parquet', A(tmp_path))
