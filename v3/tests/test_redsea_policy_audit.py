"""Tests for scripts/v3/redsea_policy_audit.py (synthetic data only; reads no project data).

Run with the Senior venv (it has numpy/pandas/pyarrow, and the end-to-end test imports the
central donor-exclusion list from the Senior worktree, read-only):

    /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python -m pytest -q \\
        scripts/v3/tests/test_redsea_policy_audit.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import redsea_policy_audit as rp  # noqa: E402

RNG = np.random.default_rng(20260929)


# --------------------------------------------------------------------------- column mapping
def test_raw_column_uses_the_cells_table_spelling():
    assert rp.raw_column("Pan_Cytokeratin", "Cell") == "Pan_Cytokeratin"
    assert rp.raw_column("Histone_H3_(p_Ser28)", "Nucleus") == "Nucleus__Histone_H3_p_Ser28"
    assert rp.raw_column("FOXP3", "Nucleus") == "Nucleus__FOXP3"


def test_registry_rows_default_measurement_source(tmp_path):
    reg = tmp_path / "reg.json"
    reg.write_text(json.dumps({"markers": [
        {"name": "INS", "kind": "type", "uses": ["broad"], "compartment": "Cell",
         "spillover_policy": "redsea_corrected"},
        {"name": "FOXP3", "kind": "type", "uses": ["deferred_type"], "compartment": "Nucleus",
         "spillover_policy": "uncompensated"}]}))
    rows = rp.load_registry_rows(reg).set_index("marker")
    assert rows.at["INS", "measurement_source"] == "redsea"
    assert rows.at["INS", "redsea_column"] == "Cell__INS"
    assert rows.at["FOXP3", "measurement_source"] == "qupath_raw"


# --------------------------------------------------------------------------- positivity rules
def _mixture(n_bg=20000, n_pos=600, bg_log=3.0, pos_log=7.0):
    bg = RNG.lognormal(bg_log, 0.5, n_bg)
    pos = RNG.lognormal(pos_log, 0.4, n_pos)
    return np.concatenate([bg, pos]), n_bg


def test_typing_rule_is_run_step3_formula():
    v, _ = _mixture()
    fit = rp.typing_rule(v)
    h, e = np.histogram(np.log1p(v[v > 0]), bins=512)
    floor = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))
    top = float(np.percentile(v, 99.9))
    assert fit.floor == pytest.approx(floor)
    assert fit.top == pytest.approx(top)
    assert fit.threshold == pytest.approx(floor + 0.10 * (top - floor))
    # positive iff scaled >= 0.10, as score.meets_requirements
    scaled = np.clip((v - floor) / (top - floor), 0, 1)
    assert np.array_equal(rp.positive(v, fit), scaled >= 0.10 - 1e-12)


def test_otsu_matches_skimage():
    skf = pytest.importorskip("skimage.filters")
    x = np.log1p(_mixture(5000, 3000)[0])
    assert rp.otsu_threshold(x) == pytest.approx(skf.threshold_otsu(x, nbins=512))


def test_rule_falls_back_to_otsu_when_the_floor_changes_population():
    # abundant marker: 30% background, 70% positive, background mode dominant in raw
    bg = RNG.lognormal(3.0, 0.15, 30000)
    pos = RNG.lognormal(7.0, 0.6, 70000)
    raw = np.concatenate([bg, pos])
    rs = raw.copy()
    rs[:25000] = 0.0            # REDSEA's clamp removes most of the background
    rule, fr, fs, note = rp.choose_rule(raw, rs)
    assert rule == "O" and "unstable" in note
    # a marker REDSEA only rescales keeps the typing rule
    rule2, *_ = rp.choose_rule(raw, raw * 1.1)
    assert rule2 == "T"


def test_positives_at_bar_moves_typing_markers_only():
    v, _ = _mixture()
    t = rp.typing_rule(v)
    o = rp.otsu_rule(v)
    fits = {("T1", "raw"): t, ("T1", "rs"): t, ("O1", "raw"): o, ("O1", "rs"): o}
    vals = {k: v for k in fits}
    base = rp.positives_at_bar(vals, fits, {"T1": "T", "O1": "O"}, rp.REQUIRE_MIN, v.size)
    strict = rp.positives_at_bar(vals, fits, {"T1": "T", "O1": "O"}, 0.30, v.size)
    assert np.array_equal(base[("T1", "raw")], rp.positive(v, t))
    assert strict[("T1", "raw")].sum() < base[("T1", "raw")].sum()
    assert np.array_equal(strict[("O1", "raw")], base[("O1", "raw")])


def test_count_matched_takes_the_k_largest_finite_values():
    v = np.array([5.0, np.nan, 1.0, 9.0, 7.0])
    assert rp.count_matched(v, 2).tolist() == [False, False, False, True, True]
    assert rp.count_matched(v, 10).sum() == 4          # NaN is never positive


# --------------------------------------------------------------------------- verdicts
@pytest.mark.parametrize("args,expected", [
    ((10, 500, 50, 0.5, 0.99, 5.0, True), "insufficient"),
    ((500, 500, 5, 0.5, 0.99, 5.0, True), "nothing_to_fix"),
    ((500, 500, 50, 0.5, 0.75, 5.0, True), "harm"),
    ((500, 500, 50, -0.3, 0.95, 5.0, True), "worse"),
    ((500, 500, 50, 0.3, 0.95, 5.0, True), "benefit"),
    ((500, 500, 50, 0.3, 0.85, 5.0, True), "tradeoff"),
    ((500, 500, 50, 0.3, 0.95, 1.2, True), "no_benefit"),       # not selective
    ((500, 500, 50, 0.15, 0.95, 5.0, True), "small_benefit"),
    ((500, 500, 50, 0.05, 0.95, 5.0, True), "no_benefit"),
    ((500, 500, 50, 0.5, 0.99, 5.0, False), "insufficient"),    # no usable rule
    ((500, 500, 50, 0.5, float("nan"), 5.0, True), "insufficient"),   # no single positive
])
def test_verdict_boundaries(args, expected):
    assert rp.verdict(*args) == expected


def _focal(rows):
    base = {"gate": "", "rule": "T", "rule_partner": "T", "copos_drop": 0.3, "retention": 0.95,
            "retention_bright": 0.99, "selectivity": 4.0, "copos_drop_k": 0.3,
            "cosine_raw": 0.2, "cosine_marker_corrected": 0.1}
    return pd.DataFrame([{**base, **r} for r in rows])


def test_marker_summary_keep_judgement_raw_and_donor_dependent():
    keep = rp.summarize_marker(_focal(
        [{"donor_id": d, "pair": "X|Y", "partner": "Y", "tier": 1, "verdict": "benefit",
          "verdict_k": "benefit"} for d in "abc"]))
    assert keep["proposal"] == "keep REDSEA" and keep["n_donors_evaluable"] == 3
    mixed = rp.summarize_marker(_focal(
        [{"donor_id": d, "pair": "X|Y", "partner": "Y", "tier": 1, "verdict": "benefit",
          "verdict_k": "benefit"} for d in "abc"]
        + [{"donor_id": d, "pair": "X|Z", "partner": "Z", "tier": 1, "verdict": "harm",
            "verdict_k": "harm"} for d in "abc"]))
    assert mixed["proposal"] == "judgement"
    raw = rp.summarize_marker(_focal(
        [{"donor_id": d, "pair": "X|Y", "partner": "Y", "tier": 1, "verdict": v,
          "verdict_k": v} for d, v in zip("abc", ["no_benefit", "nothing_to_fix",
                                                  "small_benefit"])]))
    assert raw["proposal"] == "raw"
    split = rp.summarize_marker(_focal(
        [{"donor_id": d, "pair": "X|Y", "partner": "Y", "tier": 1, "verdict": v,
          "verdict_k": v} for d, v in zip("abcd", ["benefit", "benefit", "no_benefit",
                                                   "no_benefit"])]))
    assert split["proposal"] == "judgement" and "some donors" in split["reason"]
    only_t2 = rp.summarize_marker(_focal(
        [{"donor_id": d, "pair": "V|W", "partner": "W", "tier": 2, "verdict": "no_benefit",
          "verdict_k": "no_benefit"} for d in "ab"]))
    assert only_t2["proposal"] == "raw" and only_t2["evidence_tier"] == 2
    few = rp.summarize_marker(_focal(
        [{"donor_id": "a", "pair": "X|Y", "partner": "Y", "tier": 1, "verdict": "benefit",
          "verdict_k": "benefit"}]))
    assert "provisional" in few["flags"]


# --------------------------------------------------------------------------- synthetic donors
MARKERS = ("INS", "GCG", "SST", "CD3e", "CD20")


def _registry(path: Path) -> Path:
    path.write_text(json.dumps({"markers": [
        {"name": m, "kind": "type", "uses": ["broad"], "compartment": "Cell",
         "spillover_policy": "redsea_corrected"} for m in MARKERS]}))
    return path


def _donor(root: Path, donor: str, n: int = 6000, fix_ins: bool = True) -> None:
    """Beta/alpha/delta/T/B/other cells; alpha cells carry INS spillover in the raw values.

    The "REDSEA" variant removes that spillover (fix_ins) and otherwise rescales by 1.05, so
    INS should earn REDSEA and every other marker should not.
    """
    kinds = RNG.choice(6, size=n, p=[0.08, 0.10, 0.05, 0.10, 0.10, 0.57])
    raw = {m: RNG.lognormal(3.0, 0.4, n) for m in MARKERS}
    for k, m in enumerate(MARKERS):
        own = kinds == k
        raw[m][own] = RNG.lognormal(7.5, 0.3, own.sum())
    spill = (kinds == 1) & (RNG.random(n) < 0.4)
    raw["INS"][spill] = RNG.lognormal(6.3, 0.2, spill.sum())
    red = {m: raw[m] * 1.05 for m in MARKERS}
    if fix_ins:
        red["INS"] = red["INS"].copy()
        red["INS"][spill] = RNG.lognormal(3.0, 0.4, spill.sum())
    oid = np.array([f"{donor}-{i}" for i in range(n)])
    perm = RNG.permutation(n)                   # REDSEA rows in a different order
    for sub, frame in (
        ("cells", pd.DataFrame({"object_id": oid, **raw})),
        ("redsea", pd.DataFrame({"object_id": oid[perm], "donor_id": donor,
                                 **{f"Cell__{m}": red[m][perm] for m in MARKERS}})),
        ("geometry_qc", pd.DataFrame({"donor_id": donor, "object_id": oid,
                                      "analysis_eligible": RNG.random(n) > 0.02})),
    ):
        part = root / sub / f"donor_id={donor}"
        part.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(part / ("data_0.parquet" if sub != "geometry_qc" else "data.parquet"),
                         index=False)


def test_process_donor_credits_the_corrected_marker_only(tmp_path):
    _donor(tmp_path, "D1")
    reg = rp.load_registry_rows(_registry(tmp_path / "reg.json"))
    res = rp.process_donor(tmp_path, "D1", reg)
    focal = res["focal"].set_index(["marker", "partner"])
    ins = focal.loc[("INS", "GCG")]
    assert ins["verdict"] == "benefit"
    assert ins["copos_drop"] > 0.8 and ins["retention"] > 0.95
    gcg = focal.loc[("GCG", "INS")]              # a pure rescale changes nothing
    assert gcg["verdict"] in ("no_benefit", "nothing_to_fix")
    assert abs(gcg["copos_drop"]) < 0.05
    for q in rp.SWEEP:                           # the typing-bar sensitivity is recorded
        col = f"verdict_q{int(round(q * 100))}"
        assert set(res["focal"][col]) <= set(rp.VERDICT_GROUP) | {"insufficient"}
    pairs = res["pairs"].set_index("pair")
    assert pairs.at["CD4|CD8@CD3e+", "status"].startswith("not_acquired")
    assert pairs.at["INS|GCG", "copos_union_redsea_raw"] < pairs.at["INS|GCG",
                                                                     "copos_union_raw_raw"]


def test_universes_must_match(tmp_path):
    _donor(tmp_path, "D1")
    qc = tmp_path / "geometry_qc" / "donor_id=D1" / "data.parquet"
    t = pd.read_parquet(qc)
    t.iloc[:-1].to_parquet(qc, index=False)
    with pytest.raises(ValueError, match="universes differ"):
        rp.open_donor(tmp_path, "D1")


SENIOR = rp.DEFAULT_SENIOR_PKG


@pytest.mark.skipif(not (SENIOR / "phenocycler" / "cohort.py").exists(),
                    reason="Senior worktree (central exclusion list) not present")
def test_incremental_run_gates_on_the_log_and_skips_done_donors(tmp_path, capsys):
    data, out = tmp_path / "data", tmp_path / "out"
    for d in ("D1", "D2", "6457"):
        _donor(data, d, n=4000)
    reg = _registry(tmp_path / "reg.json")
    runlog = tmp_path / "run.log"
    runlog.write_text("[donor-done] redsea/D1\n[donor-done] redsea/6457\n"
                      "[donor-run] redsea/D2\n")
    argv = ["--data-dir", str(data), "--log", str(runlog), "--out", str(out),
            "--registry", str(reg), "--senior-pkg", str(SENIOR), "--status-csv",
            str(tmp_path / "none.csv")]
    assert rp.main(argv) == 0
    first = capsys.readouterr().out
    assert "[done] D1" in first
    assert "[wait] D2: REDSEA not reported done" in first
    assert "[wait] 6457: centrally excluded" in first
    assert (out / "per_donor" / "donor_id=D1" / "meta.json").exists()
    assert not (out / "per_donor" / "donor_id=D2").exists()

    with runlog.open("a") as fh:
        fh.write("[donor-done] redsea/D2\n")
    assert rp.main(argv) == 0
    second = capsys.readouterr().out
    assert "[skip] D1: already done" in second and "[done] D2" in second
    pairs = pd.read_csv(out / "pairs_by_donor.tsv", sep="\t", dtype={"donor_id": str})
    assert set(pairs["donor_id"]) == {"D1", "D2"}
    summary = pd.read_csv(out / "markers_summary.tsv", sep="\t").set_index("marker")
    assert summary.at["INS", "n_donors_evaluable"] == 2
    readme = (out / "README.md").read_text()
    assert "| D1 |" in readme and "| D2 |" in readme and "6457" in readme
