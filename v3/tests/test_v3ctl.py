"""Unit tests for v3ctl.py (planning, custody, provenance sites, staleness, locks, literals, Parquet).

    /scratch/io60/venv-v3/bin/python -m pytest -q scripts/v3/tests/test_v3ctl.py

No QuPath and no real data: records and projects are synthesized in tmp_path. The only real files
read are steps.tsv, the step scripts (hashing), footer_template.groovy and the scripts git repo
(read-only `git rev-parse` / `git --no-optional-locks status`).
"""
import csv
import fcntl
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import v3ctl  # noqa: E402

REPO = Path('/home/smith6jt/IO60panc2nd')
STEPS = v3ctl.load_steps(HERE.parent / 'steps.tsv')
SITE = 'ws'


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
    pc = v3ctl.parse_postcondition('parenting(forbid=Other: Gal3Norm|ROOT; k=v)')
    assert pc == {'name': 'parenting', 'params': {'forbid': ['Other: Gal3Norm', 'ROOT'], 'k': ['v']}}
    assert v3ctl.parse_postcondition('v3_result') == {'name': 'v3_result', 'params': {}}


# ---- step table ---------------------------------------------------------------------------------
def test_steps_table_is_consistent():
    assert list(STEPS)[:3] == ['preflight', '1', '4']
    assert STEPS['2'].overrides == {'N_THREADS': '24'}
    s0 = STEPS['0']
    assert s0.save == 'self-saves' and s0.lane == 'gpu' and s0.args == ['{CFG}']
    assert s0.candidates == ['scripts/v3/qupath/0_cell_detection_panc_v3.groovy']       # v3 copy only
    assert s0.postcondition == {'name': 'regions_detected', 'params': {'classes': ['Tissue', 'pLN']}}
    assert ('R_import', True) in STEPS['3'].after
    assert '2' in STEPS['3'].needs and '2' in STEPS['J_context'].needs                      # finding 5
    assert set(STEPS['8'].postcondition['params']['forbid']) >= {'Other: Gal3Norm', 'ROOT'}
    stale = {s for s, st in STEPS.items() if st.on_external == 'stale'}
    assert stale == {'preflight', 'C_cells', 'C_geojson', '3', '8', 'J_context', 'J_annotations', 'J_islets', 'R_export'}


def test_steps_header_must_include_on_external(tmp_path):
    f = tmp_path / 'steps.tsv'
    lines = (HERE.parent / 'steps.tsv').read_text().splitlines()
    f.write_text('\n'.join(ln.replace('\ton_external', '') if ln.startswith('step\t') else ln for ln in lines))
    with pytest.raises(SystemExit):
        v3ctl.load_steps(f)


def test_selection_is_topological():
    order = v3ctl.selected_steps(STEPS, ['R', 'J'], [])
    assert order.index('R_import') < order.index('3') < order.index('J_context')
    assert v3ctl.selected_steps(STEPS, ['I', 'C'], []) == ['preflight', '1', '4', '0', '2', 'C_cells', 'C_geojson']
    with pytest.raises(SystemExit):
        v3ctl.selected_steps(STEPS, ['Q'], [])


def test_script_selection(tmp_path):
    st = STEPS['1']
    p, how = v3ctl.resolve_script(st, REPO, 'orig', {})
    assert p.name == '1_panc_objects.groovy' and how == 'orig'
    fake = tmp_path / 'x.groovy'
    fake.write_text('println 1')
    assert v3ctl.resolve_script(st, REPO, 'auto', {'1': str(fake)}) == (fake, 'override')


def test_script_set_v3_never_falls_back(tmp_path):
    (tmp_path / 'orig.groovy').write_text('println 1')
    st = v3ctl.Step('s', ['I'], ['missing_v3.groovy', 'orig.groovy'], [], 'never', 'light', 1, 1, 1, [], [],
                    'keep', {}, {'name': 'v3_result', 'params': {}}, None, '')
    with pytest.raises(v3ctl.PlanError, match='never falls back'):
        v3ctl.resolve_script(st, tmp_path, 'v3', {})
    assert v3ctl.resolve_script(st, tmp_path, 'auto', {})[1] == 'auto:orig'


def test_default_script_set_is_v3():
    """finding 6: production default --script-set v3 (run_v3.sh passes its own default, v3, too)"""
    out = subprocess.run([sys.executable, str(HERE.parent / 'v3ctl.py'), '--help'], capture_output=True, text=True).stdout
    assert 'v3 (default)' in out
    assert 'SCRIPT_SET=v3' in (HERE.parent / 'run_v3.sh').read_text()


def test_budget_sums_timeouts(capsys):
    v3ctl.main(['budget', '--steps', '0', '--margin-min', '45', '--format', 'slurm'])
    assert capsys.readouterr().out.strip() == '08:45:00'
    v3ctl.main(['budget', '--passes', 'R,J', '--margin-min', '30'])
    assert capsys.readouterr().out.strip() == str(60 + 240 + 240 + 240 + 120 + 120 + 30)


# ---- image selection (finding 3) --------------------------------------------------------------
def fake_images():
    names = ['0080870_6533_pLNpanc_Scan3.er.qptiff - resolution #1', '6533_Scan1.er.qptiff - resolution #1',
             '6450panc_Scan1.er.qptiff - resolution #1']
    return [v3ctl.Image(n, i + 1, v3ctl.tag_of(n), Path(f'/x/{i + 1}/data.qpdata'), {}) for i, n in enumerate(names)]


def test_images_exact_name_tag_or_id_only():
    ims = fake_images()
    assert [im.entry_id for im in v3ctl.resolve_images(ims, ['6533_Scan1.er.qptiff - resolution #1'])] == [2]
    assert [im.entry_id for im in v3ctl.resolve_images(ims, ['6450panc_Scan1.er.qptiff_-_resolution_1', 'id:1'])] == [3, 1]
    with pytest.raises(SystemExit, match='matches no image'):
        v3ctl.resolve_images(ims, ['6533'])                     # the old substring match selected two images
    with pytest.raises(SystemExit, match='matches no image'):
        v3ctl.resolve_images(ims, ['1'])                        # an entry id needs the explicit id: prefix
    with pytest.raises(SystemExit, match='second time'):
        v3ctl.resolve_images(ims, ['id:3', '6450panc_Scan1.er.qptiff - resolution #1'])
    with pytest.raises(SystemExit, match="'id:' must be followed"):
        v3ctl.resolve_images(ims, ['id:x'])


def test_images_file_ignores_blank_lines_and_comments(tmp_path):
    f = tmp_path / 'imgs.txt'
    f.write_text('\n# pilot\n  id:2  \n\n6450panc_Scan1.er.qptiff - resolution #1\n\n')
    ims = v3ctl.select_images(fake_images(), None, str(f))
    assert [im.entry_id for im in ims] == [2, 3]
    (tmp_path / 'empty.txt').write_text('\n# nothing\n')
    with pytest.raises(SystemExit, match='no image names'):
        v3ctl.select_images(fake_images(), None, str(tmp_path / 'empty.txt'))


def test_duplicate_tags_in_project_are_refused(tmp_path):
    proj = tmp_path / 'project.qpproj'
    proj.write_text(json.dumps({'images': [{'imageName': 'a b', 'entryID': 1}, {'imageName': 'a_b', 'entryID': 2}]}))
    with pytest.raises(SystemExit, match='share the tag'):
        v3ctl.load_images(proj)


# ---- records, custody and planning ----------------------------------------------------------------
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


def write_recs(out, *recs, site=SITE, finish=True):
    d = v3ctl.site_dir(out, 'img', site)
    d.mkdir(parents=True, exist_ok=True)
    for r in recs:
        v3ctl.write_json(d / f"{r['step']}.json", r)
        if finish:
            v3ctl.write_json(d / f"{r['step']}.finish.json", {'schema': 'run_v3.finish/1', 'step': r['step'],
                                                              'attempt_id': r['attempt_id'], 'seq': r['seq'],
                                                              'status': 'PASS' if r['status'] == 'pass' else 'FAIL'})
    return d


def plan(im, sel, out, force=(), allow_drift=False, scripts=None, site=SITE):
    return v3ctl.plan_image(im, STEPS, sel, set(force), out, allow_drift, scripts, site)


def actions(p):
    return {s: a for s, a, _ in p['steps']}


def reasons(p):
    return {s: r for s, _, r in p['steps']}


def test_custody_fresh_ok_broken(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    sel = v3ctl.selected_steps(STEPS, ['I'], [])
    assert plan(im, sel, out)['state'] == 'fresh'
    write_recs(out, rec('preflight', 1, im.qpdata))
    assert plan(im, sel, out)['state'] == 'ok'
    ns = os.stat(im.qpdata).st_mtime_ns
    os.utime(im.qpdata, ns=(ns, ns + 500_000_000))            # within the 1 s tolerance (Globus/rsync)
    assert plan(im, sel, out)['state'] == 'ok'
    im.qpdata.write_bytes(b'y' * 101)                          # someone saved the image
    p = plan(im, sel, out)
    assert p['state'] == 'custody_broken' and set(actions(p).values()) == {'blocked'}


def test_leftover_bkp_is_named(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata))
    v3ctl.bkp_path(im.qpdata).write_bytes(b'old')
    p = plan(im, ['1'], out)
    assert 'data.qpdata.bkp exists' in p['message'] and 'readImageData' in p['message']


def test_done_skip_and_downstream_of_forced(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata))
    sel = v3ctl.selected_steps(STEPS, ['I'], [])
    assert actions(plan(im, sel, out)) == {'preflight': 'skip', '1': 'skip', '4': 'skip'}
    assert actions(plan(im, sel, out, force={'1'})) == {'preflight': 'skip', '1': 'run', '4': 'run'}


def test_stale_when_upstream_reran(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('4', 2, im.qpdata), rec('1', 3, im.qpdata))
    p = plan(im, v3ctl.selected_steps(STEPS, ['I'], []), out)
    assert actions(p)['4'] == 'run' and 'stale' in reasons(p)['4']


def test_failed_attempt_reruns_and_blocks_missing_prerequisite(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('0', 2, im.qpdata, status='fail'))
    pc = plan(im, v3ctl.selected_steps(STEPS, ['C'], []), out)
    assert actions(pc)['0'] == 'run' and actions(pc)['2'] == 'run'
    pj = plan(im, v3ctl.selected_steps(STEPS, ['J'], []), out)
    assert actions(pj)['3'] == 'blocked'                        # needs 0 and 2 (pass C) and 1 (pass I)
    assert all(a == 'blocked' for s, a in actions(pj).items() if s != '3')


def test_step3_needs_step2(tmp_path):
    """finding 5: pass J may not start on an image whose step 2 is not done (e.g. running in another job)."""
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata), rec('0', 4, im.qpdata))
    p = plan(im, ['3', 'J_context'], out)
    assert actions(p)['3'] == 'blocked' and 'needs 2' not in reasons(p)['3'] and '2 (never ran)' in reasons(p)['3']
    write_recs(out, rec('2', 5, im.qpdata))
    assert actions(plan(im, ['3'], out))['3'] == 'run'


def test_post_hook_and_outputs_are_part_of_done(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    parquet = tmp_path / 'cells.parquet'
    parquet.write_bytes(b'p' * 10)
    write_recs(out, rec('preflight', 1, im.qpdata), rec('0', 2, im.qpdata), rec('2', 3, im.qpdata),
               rec('C_cells', 4, im.qpdata), rec('C_geojson', 5, im.qpdata))
    sel = ['C_cells']
    assert actions(plan(im, sel, out)) == {'C_cells': 'run'}   # no post record
    v3ctl.write_json(v3ctl.site_dir(out, 'img', SITE) / 'C_cells.post.json',
                     {'status': 'pass', 'attempt_id': 'a4', 'outputs': [v3ctl.stat_of(parquet)]})
    assert actions(plan(im, sel, out)) == {'C_cells': 'skip'}
    parquet.unlink()
    assert actions(plan(im, sel, out)) == {'C_cells': 'run'}   # output gone


def test_optional_dependency_only_counts_once_it_ran(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata),
               rec('0', 4, im.qpdata), rec('2', 5, im.qpdata), rec('3', 6, im.qpdata))
    assert actions(plan(im, ['3'], out)) == {'3': 'skip'}
    write_recs(out, rec('R_import', 7, im.qpdata))              # reviewed islets imported after 3 ran
    assert actions(plan(im, ['3'], out)) == {'3': 'run'}


def test_drift_record_blocks_unless_allowed(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    d = write_recs(out, rec('preflight', 1, im.qpdata))
    v3ctl.write_json(d / '1.drift.json', {'server_json': {}})
    assert plan(im, ['1'], out)['state'] == 'drift'
    assert plan(im, ['1'], out, allow_drift=True)['state'] == 'ok'


def test_sentinel_without_finish_is_not_done(tmp_path):
    """finding 11: a pass sentinel an orphaned QuPath wrote after the driver died was never validated."""
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata))
    write_recs(out, rec('1', 2, im.qpdata), finish=False)
    p = plan(im, ['1'], out)
    assert actions(p) == {'1': 'run'} and 'never validated by finish' in reasons(p)['1']
    v3ctl.write_json(v3ctl.site_dir(out, 'img', SITE) / '1.finish.json', {'attempt_id': 'a2', 'status': 'DRIFT'})
    assert 'finish said DRIFT' in reasons(plan(im, ['1'], out))['1']
    v3ctl.write_json(v3ctl.site_dir(out, 'img', SITE) / '1.finish.json', {'attempt_id': 'other', 'status': 'PASS'})
    assert actions(plan(im, ['1'], out)) == {'1': 'run'}


def scripts_now():
    return v3ctl.resolve_all_scripts(STEPS, REPO, 'v3', {})


def with_script(r, sid, scripts):
    r['script'] = {'source_path': scripts[sid]['path'], 'source_sha256': scripts[sid]['sha256']}
    return r


def test_stale_when_the_script_changed(tmp_path):
    """finding 6: a record made by another script (path or content) is not done under the current plan."""
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    sc = scripts_now()
    write_recs(out, with_script(rec('preflight', 1, im.qpdata), 'preflight', sc), with_script(rec('1', 2, im.qpdata), '1', sc))
    assert actions(plan(im, ['1'], out, scripts=sc)) == {'1': 'skip'}
    orig = v3ctl.resolve_all_scripts(STEPS, REPO, 'orig', {})
    p = plan(im, ['1'], out, scripts=orig)                        # same records, the plan now uses the original
    assert actions(p) == {'1': 'run'} and 'script changed' in reasons(p)['1']
    sc2 = dict(sc)
    sc2['1'] = dict(sc['1'], sha256='0' * 64)                    # same path, edited content
    assert 'script changed' in reasons(plan(im, ['1'], out, scripts=sc2))['1']
    write_recs(out, with_script(rec('4', 3, im.qpdata), '4', sc))
    sc3 = dict(sc)
    sc3['1'] = {'error': 'gone'}                                 # an upstream ('after') step whose script vanished
    assert 'cannot be resolved' in reasons(plan(im, ['4'], out, scripts=sc3))['4']
    sc4 = dict(sc)
    sc4['preflight'] = {'error': 'gone'}                         # 'needs' only gates starting: 1 stays done
    assert actions(plan(im, ['1'], out, scripts=sc4)) == {'1': 'skip'}


def test_accepted_external_edit_makes_readers_stale(tmp_path):
    """finding 1: after --accept-external every step that reads the edited image re-runs (never-policy
    exports, 3, 8, J_*, R_export), the reviewed steps (1, 4, 0, 2) do not."""
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    names = ['preflight', '1', '4', '0', '2', 'C_cells', 'C_geojson', '3', '8', 'J_context', 'J_annotations',
             'J_islets', 'R_export']
    d = write_recs(out, *[rec(s, i + 1, im.qpdata) for i, s in enumerate(names)])
    for s in ('C_cells', 'J_context'):
        v3ctl.write_json(d / f'{s}.post.json', {'status': 'pass', 'attempt_id': f'a{names.index(s) + 1}', 'outputs': []})
    sel = [s for s in STEPS if s != 'R_import']
    assert set(actions(plan(im, sel, out)).values()) == {'skip'}
    v3ctl.write_json(d / 'external.14.json', {'schema': 'run_v3.provenance/1', 'kind': 'external', 'step': 'external',
                                              'status': 'pass', 'seq': 14, 'note': 'gate 1', 'rerun_from': [],
                                              'qpdata': {'before': v3ctl.stat_of(im.qpdata), 'after': v3ctl.stat_of(im.qpdata)}})
    p = plan(im, sel, out)
    a = actions(p)
    assert {s for s, x in a.items() if x == 'run'} == {'preflight', 'C_cells', 'C_geojson', '3', '8', 'J_context',
                                                         'J_annotations', 'J_islets', 'R_export'}
    assert {s for s, x in a.items() if x == 'skip'} == {'1', '4', '0', '2'}
    assert 'external edit was accepted' in reasons(p)['R_export']


def test_external_rerun_from_is_recorded_and_propagates(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    names = ['preflight', '1', '4', '0', '2']
    d = write_recs(out, *[rec(s, i + 1, im.qpdata) for i, s in enumerate(names)])
    v3ctl.write_json(d / 'external.6.json', {'schema': 'run_v3.provenance/1', 'kind': 'external', 'step': 'external',
                                             'status': 'pass', 'seq': 6, 'note': 'fixed cells', 'rerun_from': ['2'],
                                             'qpdata': {'before': v3ctl.stat_of(im.qpdata), 'after': v3ctl.stat_of(im.qpdata)}})
    p = plan(im, v3ctl.selected_steps(STEPS, ['I', 'C'], []), out)
    a = actions(p)
    assert a['2'] == 'run' and 'rerun-from 2' in reasons(p)['2']
    assert a['0'] == 'skip' and a['1'] == 'skip' and a['preflight'] == 'run'     # preflight: on_external=stale


def test_force_outside_selection_is_an_error():
    """finding 7"""
    class A:
        force, rerun_from = '8', ''
    with pytest.raises(SystemExit, match='not among the selected steps'):
        v3ctl.plan_force(A(), STEPS, ['preflight', '1', '4'])
    A.force, A.rerun_from = '', '4'
    assert v3ctl.plan_force(A(), STEPS, ['preflight', '1', '4']) == {'4'}


# ---- provenance per site, merge, chain (finding 4) ----------------------------------------------
def test_records_of_another_site_block_until_merged(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata), site='ws')
    write_recs(out, rec('0', 4, im.qpdata), site='hpg')
    p = plan(im, ['0'], out, site='ws')
    assert p['state'] == 'unmerged' and 'hpg' in p['message'] and set(actions(p).values()) == {'blocked'}


def run_merge(out, project, site, *extra):
    return subprocess.run([sys.executable, str(HERE.parent / 'v3ctl.py'), '--project', str(project), '--out', str(out),
                           '--site', site, 'merge-provenance', '--images', 'img', '--lockdir', str(out / 'locks'),
                           '--shared-lockdir', '', *extra], capture_output=True, text=True)


def make_project(tmp_path):
    proj = tmp_path / 'proj' / 'project.qpproj'
    proj.write_text(json.dumps({'images': [{'imageName': 'img', 'entryID': 1, 'serverBuilder': {}}]}))
    return proj


def test_merge_provenance_mirrors_and_validates_the_chain(tmp_path):
    im = make_image(tmp_path)
    proj = make_project(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('1', 2, im.qpdata), rec('4', 3, im.qpdata), site='ws')
    before = v3ctl.stat_of(im.qpdata)
    im.qpdata.write_bytes(b'z' * 300)                              # hpg ran step 0 and saved; the file came back
    after = v3ctl.stat_of(im.qpdata)
    write_recs(out, rec('0', 4, im.qpdata, qpdata={'before': before, 'after': after}), site='hpg')
    r = run_merge(out, proj, 'ws', '--check')
    assert r.returncode == 0 and 'would merge' in r.stdout
    assert not (v3ctl.site_dir(out, 'img', 'ws') / 'merged' / 'hpg').exists()
    r = run_merge(out, proj, 'ws')
    assert r.returncode == 0 and 'merged (hpg: +2 new' in r.stdout, r.stdout + r.stderr
    p = plan(im, ['0', '2'], out, site='ws')
    assert p['state'] == 'ok' and 'hpg:0' in p['message'] and actions(p) == {'0': 'skip', '2': 'run'}
    # the other site re-runs step 0 later: its new record must be merged before anything plans
    newer = v3ctl.stat_of(im.qpdata)
    write_recs(out, rec('0', 5, im.qpdata, qpdata={'before': after, 'after': newer}, attempt_id='a5b'), site='hpg')
    assert plan(im, ['0'], out, site='ws')['state'] == 'unmerged'


def test_merge_refuses_when_the_data_did_not_come_along(tmp_path):
    im = make_image(tmp_path)
    proj = make_project(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), site='ws')
    start = v3ctl.stat_of(im.qpdata)
    other = dict(start, size=999)                                 # hpg's file after its save, not synced here
    write_recs(out, rec('0', 2, im.qpdata, qpdata={'before': start, 'after': other}), site='hpg')
    r = run_merge(out, proj, 'ws')
    assert r.returncode == 1 and 'NOT merged' in r.stdout and 'custody_broken' in r.stdout
    assert 'older copy' in r.stdout                               # it says which record the file matches


def test_duplicate_seq_across_sites_is_a_conflict(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('0', 2, im.qpdata), site='ws')
    d = v3ctl.site_dir(out, 'img', 'ws') / 'merged' / 'hpg'
    d.mkdir(parents=True)
    v3ctl.write_json(d / '0.json', rec('0', 2, im.qpdata, attempt_id='other'))       # the same seq on hpg
    recs = v3ctl.load_records(v3ctl.prov_dir(out, 'img'), 'ws')
    state, _, msg = v3ctl.custody(recs, v3ctl.stat_of(im.qpdata))
    assert state == 'conflict' and 'share a sequence number' in msg and 'including the highest' in msg
    p = plan(im, ['0'], out, site='ws')
    assert p['state'] == 'custody_conflict' and set(actions(p).values()) == {'blocked'}


def test_identical_record_reachable_twice_counts_once(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    d = write_recs(out, rec('preflight', 1, im.qpdata), site='ws')
    (d / 'history').mkdir()
    v3ctl.write_json(d / 'history' / 'preflight.seq1.1.json', rec('preflight', 1, im.qpdata))
    assert plan(im, ['preflight'], out, site='ws')['state'] == 'ok'


def test_custody_chain_break_is_a_conflict(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    s = v3ctl.stat_of(im.qpdata)
    elsewhere = dict(s, size=12345)
    write_recs(out, rec('preflight', 1, im.qpdata, qpdata={'before': s, 'after': elsewhere}),
               rec('1', 2, im.qpdata), site='ws')                  # 2 started from a file 1 did not leave
    p = plan(im, ['1'], out, site='ws')
    assert p['state'] == 'custody_conflict' and 'chain broken' in p['message']


def test_failed_record_that_changed_the_file_breaks_the_chain(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    s = v3ctl.stat_of(im.qpdata)
    write_recs(out, rec('1', 1, im.qpdata, status='fail', qpdata={'before': dict(s, size=1), 'after': s}), site='ws')
    assert plan(im, ['1'], out, site='ws')['state'] == 'custody_conflict'


def test_foreign_in_flight_marker_blocks_unless_that_attempt_finished(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), site='ws')
    hd = v3ctl.site_dir(out, 'img', 'hpg')
    hd.mkdir(parents=True)
    v3ctl.write_json(hd / '.running.0.json', {'attempt_id': 'z1', 'started': 'now'})
    assert plan(im, ['0'], out, site='ws')['state'] == 'foreign_in_flight'
    v3ctl.write_json(hd / '0.finish.json', {'attempt_id': 'z1', 'status': 'FAIL'})       # it finished there
    p = plan(im, ['0'], out, site='ws')
    assert p['state'] == 'unmerged'                               # (the finish record itself is not merged yet)


def test_legacy_flat_layout_is_refused(tmp_path):
    im = make_image(tmp_path)
    out = tmp_path / 'out'
    v3ctl.write_json(v3ctl.prov_dir(out, 'img') / 'preflight.json', rec('preflight', 1, im.qpdata))
    assert plan(im, ['1'], out)['state'] == 'legacy_layout'


def test_site_name_from_environment(monkeypatch):
    monkeypatch.setenv('V3_SITE', 'hpg')
    assert v3ctl.resolve_site(None) == 'hpg'
    monkeypatch.setenv('V3_SITE', 'merged')
    with pytest.raises(SystemExit):
        v3ctl.resolve_site(None)
    monkeypatch.delenv('V3_SITE')
    assert v3ctl.resolve_site(None) == os.uname().nodename.split('.')[0]


# ---- accept-external (findings 1, 11) -------------------------------------------------------------
def test_accept_external_records_rerun_from_under_the_image_lock(tmp_path):
    im = make_image(tmp_path)
    proj = make_project(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), rec('0', 2, im.qpdata), site='ws')
    im.qpdata.write_bytes(b'gui' * 50)                            # edited in the GUI
    cmd = [sys.executable, str(HERE.parent / 'v3ctl.py'), '--project', str(proj), '--out', str(out), '--site', 'ws',
           'accept-external', '--images', 'img', '--note', 'edited cells', '--rerun-from', '2',
           '--lockdir', str(tmp_path / 'locks'), '--shared-lockdir', '']
    (tmp_path / 'locks').mkdir()
    with open(tmp_path / 'locks' / 'img.lock', 'a') as held:     # a run holds the image
        fcntl.flock(held, fcntl.LOCK_EX)
        r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode != 0 and 'image lock' in r.stderr
    assert not list(v3ctl.site_dir(out, 'img', 'ws').glob('external.*.json'))
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    ext = json.loads((v3ctl.site_dir(out, 'img', 'ws') / 'external.3.json').read_text())
    assert ext['rerun_from'] == ['2'] and 'C_cells' in ext['invalidates'] and ext['qpdata']['after']['size'] == 150
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0 and 'nothing to accept' in r.stdout
    bad = cmd[:cmd.index('--rerun-from') + 1] + ['nope'] + cmd[cmd.index('--rerun-from') + 2:]
    assert subprocess.run(bad, capture_output=True, text=True).returncode != 0


# ---- prepare / finish (findings 2, 9, 11) ---------------------------------------------------------
def prepared(tmp_path, capsys, deadline='1900000000'):
    im = make_image(tmp_path)
    proj = make_project(tmp_path)
    out = tmp_path / 'out'
    worker = tmp_path / 'w' / 'x'
    v3ctl.main(['--project', str(proj), '--out', str(out), '--site', 'ws', '--run-id', 'r1', 'prepare', '--image', 'img',
                '--step', '0', '--workdir', str(worker), '--tmp', str(tmp_path), '--deadline-epoch', deadline,
                '--save-margin-min', '7'])
    lines = capsys.readouterr().out.strip().splitlines()
    return im, proj, out, worker, dict(ln.split('=', 1) for ln in lines), lines


def test_prepare_prints_only_known_keys_and_writes_the_step_config(tmp_path, capsys):
    im, proj, out, worker, kv, lines = prepared(tmp_path, capsys)
    assert [ln.split('=', 1)[0] for ln in lines] == list(v3ctl.PREPARE_KEYS)
    cfg = json.loads((worker / 'step_config.json').read_text())
    assert cfg['qpdata_before'] == v3ctl.stat_of(im.qpdata)
    assert cfg['deadline_epoch_ms'] == 1900000000000 and cfg['save_margin_min'] == 7.0 and cfg['site'] == 'ws'
    assert len(cfg['repo_git']['commit'] or '') == 40 and isinstance(cfg['repo_git']['dirty'], bool)
    assert cfg['args'] == [str(worker / 'step_config.json')]          # {CFG}
    assert Path(kv['QARGS_FILE']).read_bytes() == b'-a\0' + str(worker / 'step_config.json').encode() + b'\0'
    body = Path(kv['WSCRIPT']).read_text()
    assert "'deadline_epoch_ms': 1900000000000L" in body and 'REFUSED' in body
    assert (v3ctl.site_dir(out, 'img', 'ws') / '.running.0.json').exists()


def test_prepare_refuses_when_custody_broke_after_the_plan(tmp_path, capsys):
    im = make_image(tmp_path)
    proj = make_project(tmp_path)
    out = tmp_path / 'out'
    write_recs(out, rec('preflight', 1, im.qpdata), site='ws')
    im.qpdata.write_bytes(b'q' * 7)
    with pytest.raises(SystemExit) as e:
        v3ctl.main(['--project', str(proj), '--out', str(out), '--site', 'ws', 'prepare', '--image', 'img', '--step', '0',
                    '--workdir', str(tmp_path / 'w'), '--tmp', str(tmp_path)])
    assert e.value.code == 2 and 'custody_broken' in capsys.readouterr().err


def finish(proj, out, worker, attempt, log='', code='0'):
    try:
        v3ctl.main(['--project', str(proj), '--out', str(out), '--site', 'ws', 'finish', '--image', 'img', '--step', '0',
                    '--attempt', attempt, '--workdir', str(worker), '--exit-code', code, '--log', log])
    except SystemExit as e:
        return e.code


def test_finish_writes_its_verdict_and_names_the_bkp(tmp_path, capsys):
    im, proj, out, worker, kv, _ = prepared(tmp_path, capsys)
    sd = v3ctl.site_dir(out, 'img', 'ws')
    sentinel = rec('0', int(kv['SEQ']), im.qpdata, attempt_id=kv['ATTEMPT'])
    sentinel['script'] = {'instantiated_sha256': v3ctl.sha256_file(Path(kv['WSCRIPT']))}
    v3ctl.write_json(sd / '0.json', sentinel)
    v3ctl.bkp_path(im.qpdata).write_bytes(b'b')
    assert finish(proj, out, worker, kv['ATTEMPT']) == 0
    fin = json.loads((sd / '0.finish.json').read_text())
    assert fin['status'] == 'PASS' and fin['attempt_id'] == kv['ATTEMPT'] and 'data.qpdata.bkp exists' in fin['reason']
    assert fin['qpdata_bkp']['exists'] and not (sd / '.running.0.json').exists()


def test_finish_without_sentinel_quotes_the_refusal(tmp_path, capsys):
    im, proj, out, worker, kv, _ = prepared(tmp_path, capsys)
    log = tmp_path / 'step.log'
    log.write_text(f"===== x run r attempt OLD seq 1\n[run_v3] REFUSED: old\n===== y run r attempt {kv['ATTEMPT']} seq 2\n"
                   f"[run_v3] REFUSED: step 0: not saving: data.qpdata changed since the step began\n")
    assert finish(proj, out, worker, kv['ATTEMPT'], log=str(log), code='1') == 1
    out_txt = capsys.readouterr().out
    assert 'no sentinel' in out_txt and 'changed since the step began' in out_txt and 'REFUSED: old' not in out_txt
    assert not (v3ctl.site_dir(out, 'img', 'ws') / '0.finish.json').exists()


def test_run_log_keeps_an_older_header(tmp_path):
    (tmp_path / 'run_log.csv').write_text('run_id,step,status\n')
    v3ctl.append_run_log(tmp_path, {'run_id': 'r', 'step': 's', 'status': 'PASS', 'site': 'ws'})
    rows = list(csv.reader(open(tmp_path / 'run_log.csv')))
    assert rows == [['run_id', 'step', 'status'], ['r', 's', 'PASS']]
    v3ctl.append_run_log(tmp_path / 'new', {'run_id': 'r', 'site': 'ws'})
    assert next(csv.reader(open(tmp_path / 'new' / 'run_log.csv'))) == v3ctl.RUN_LOG_FIELDS


# ---- shared mkdir locks (finding 5) -------------------------------------------------------------
def test_shared_lock_acquire_release_and_stale(tmp_path):
    root = tmp_path / 'locks'
    me = v3ctl.lock_owner(os.getpid(), 'r1')
    ok, _ = v3ctl.shared_lock_acquire(root, 'img', me)
    assert ok
    ok, msg = v3ctl.shared_lock_acquire(root, 'img', v3ctl.lock_owner(os.getppid(), 'r2'))
    assert not ok and 'held by' in msg                            # a live owner keeps it
    ok, msg = v3ctl.shared_lock_release(root, 'img', os.getppid())
    assert not ok and 'not ours' in msg
    assert v3ctl.shared_lock_release(root, 'img', os.getpid())[0]
    # an owner whose process is gone: the lock is broken and taken over
    p = subprocess.Popen(['true'])
    p.wait()
    dead = dict(v3ctl.lock_owner(os.getpid(), 'r3'), pid=p.pid, pid_start=None)
    root.joinpath('img.lock.d').mkdir()
    v3ctl.write_json(root / 'img.lock.d' / 'owner.json', dead)
    ok, msg = v3ctl.shared_lock_acquire(root, 'img', me)
    assert ok and json.loads((root / 'img.lock.d' / 'owner.json').read_text())['pid'] == os.getpid()
    # an owner on another host without a SLURM job cannot be checked from here: held
    other = dict(dead, host='elsewhere', slurm_job=None)
    v3ctl.write_json(root / 'img.lock.d' / 'owner.json', other)
    ok, msg = v3ctl.shared_lock_acquire(root, 'img', me)
    assert not ok and 'another host' in msg


def test_shared_lock_owner_on_ended_slurm_job_is_stale(tmp_path, monkeypatch):
    bin_ = tmp_path / 'bin'
    bin_.mkdir()
    (bin_ / 'squeue').write_text('#!/bin/sh\necho "slurm_load_jobs error: Invalid job id specified" >&2\nexit 1\n')
    (bin_ / 'squeue').chmod(0o755)
    monkeypatch.setenv('PATH', f'{bin_}:{os.environ["PATH"]}')
    alive, why = v3ctl.owner_alive({'host': 'c1', 'pid': 1, 'slurm_job': '123_4', 'time': 't'}, time.time())
    assert not alive and 'ended' in why
    (bin_ / 'squeue').write_text('#!/bin/sh\necho RUNNING\n')
    alive, why = v3ctl.owner_alive({'host': 'c1', 'pid': 1, 'slurm_job': '123_4', 'time': 't'}, time.time())
    assert alive and 'RUNNING' in why


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
    res = v3ctl.tsv_to_parquet(tsv, tmp_path / 'x.parquet', A(tmp_path), ('Cell: DAPI: Mean', 'Missing col'))
    assert (res['rows'], res['object_id_distinct'], res['columns']) == (2, 2, 5)
    import pyarrow.parquet as pq
    t = pq.read_table(res['part']).to_pydict()
    assert t['Centroid X µm'] == [1.5, 2.5]
    assert str(t['Cell: DAPI: Mean'][0]) == 'nan' and t['Cell: DAPI: Mean'][1] is None
    assert res['string_columns'] == ['Image', 'Object ID', 'Parent']
    # finding 2: the post hook compares this with the row count -- NaN and empty both count as missing
    assert res['complete_counts'] == {'Cell: DAPI: Mean': 0, 'Missing col': None}


def test_tsv_with_text_in_numeric_column_fails(tmp_path):
    pytest.importorskip('duckdb')
    tsv = tmp_path / 'bad.tsv'
    tsv.write_text('Object ID\tSomething new\nid-1\thello\n', encoding='utf-8')
    with pytest.raises(Exception):
        v3ctl.tsv_to_parquet(tsv, tmp_path / 'bad.parquet', A(tmp_path))


def test_cells_post_hook_requires_dapi_on_every_row(tmp_path, capsys):
    pytest.importorskip('duckdb')
    im, proj, out, worker, kv, _ = prepared(tmp_path, capsys)
    sd = v3ctl.site_dir(out, 'img', 'ws')
    tsv = tmp_path / 'img.cells.tsv'
    tsv.write_text('Object ID\tCell: DAPI: Mean\nid-1\t5.0\nid-2\tNaN\n', encoding='utf-8')
    r = rec('C_cells', 9, im.qpdata, attempt_id='att', transient_outputs=[v3ctl.stat_of(tsv)], counts={'detections': 2})
    v3ctl.write_json(sd / 'C_cells.json', r)
    with pytest.raises(SystemExit) as e:
        v3ctl.main(['--project', str(proj), '--out', str(out), '--site', 'ws', 'post', '--image', 'img', '--step', 'C_cells',
                    '--attempt', 'att', '--tmp', str(tmp_path), '--duckdb-mem-gb', '1', '--duckdb-threads', '1'])
    assert e.value.code == 1
    post = json.loads((sd / 'C_cells.post.json').read_text())
    assert post['status'] == 'fail' and "'Cell: DAPI: Mean' is null/NaN on 1 of 2 rows" in post['error']
    assert not (out / 'cells' / 'img.parquet').exists() and tsv.exists()
