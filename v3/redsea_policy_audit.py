#!/usr/bin/env python3
"""Per-marker evidence for the REDSEA spillover-correction policy ("option 3").

The rule the user approved on 2026-09-29:

    Keep REDSEA correction for a marker only where it reduces biologically impossible
    co-positives, without destroying that marker's real single-positive signal.

This script measures both halves of that sentence for every marker that takes part in a
declared mutually-exclusive pair, per donor, raw vs REDSEA, and proposes a per-marker policy.
It writes evidence only.  It never edits ``marker_registry.json``, the Senior repo, or
anything under ``results/v3/senior``; the policy change is a separate PR after review.

Inputs (read-only), per donor partition ``donor_id=<id>`` under ``--data-dir``:

* ``cells/``        raw QuPath means.  ``Cell: <m>: Mean`` is stored as the bare sanitized
                    marker name (``Vimentin``); ``Nucleus: <m>: Mean`` as ``Nucleus__<m>``.
* ``redsea/``       REDSEA output, ``{Nucleus,Cytoplasm,Membrane,Cell}__<registry name>``.
                    The expression stage reads ``<compartment>__<name>`` for every
                    ``redsea_corrected`` / ``redsea_passthrough`` marker, so that column
                    (the registry's ``measurement_column``) is the "REDSEA" variant here.
* ``geometry_qc/``  ``analysis_eligible`` (the expression stage's ``qc_analysis_eligible``).

A donor is read only when its REDSEA partition is complete: the partition holds exactly one
parquet AND the run log has ``[donor-done] redsea/<id>`` (pipelined run) or
``[<id>] DONE in`` (stage-only run).  Donors already processed from the same data directory
are skipped, so the script can be rerun as more partitions land.

Run from ``/scratch/io60/work/senior-pr17/core/phenocycler`` with the Senior venv::

    python /home/smith6jt/IO60panc2nd/scripts/v3/redsea_policy_audit.py            # full run
    python .../redsea_policy_audit.py --data-dir .../pilot_data \\
        --log .../pilot_redsea.log --donors 6533                                    # pilot

Outputs under ``--out`` (default ``results/v3/redsea_policy``): ``pairs_by_donor.tsv``,
``marker_evidence_by_donor.tsv``, ``thresholds_by_donor.tsv``, ``markers_summary.tsv`` and
``README.md``; per-donor caches live in ``per_donor/donor_id=<id>/``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import re
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

#: A REDSEA job shares this host: cap every thread pool before numpy/pyarrow start.
THREADS = int(os.environ.get("REDSEA_AUDIT_THREADS", "4"))
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, str(THREADS))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

pa.set_cpu_count(THREADS)
pa.set_io_thread_count(THREADS)

# ----------------------------------------------------------------------------- locations
V3 = Path("/home/smith6jt/IO60panc2nd/results/v3")
DEFAULT_DATA = V3 / "senior" / "data"
DEFAULT_LOG = V3 / "senior" / "full_run_20260929_pipelined.log"
DEFAULT_OUT = V3 / "redsea_policy"
DEFAULT_SENIOR_PKG = Path("/scratch/io60/work/senior-pr17/core/phenocycler")
DEFAULT_STATUS_CSV = Path(
    "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/donor_acquisition_batch.csv")
SCRIPT_VERSION = "redsea-policy-audit-v1"


# ----------------------------------------------------------------------------- exclusive pairs
@dataclass(frozen=True)
class Pair:
    """Two markers no single cell should carry together.

    ``tier`` 1 pairs decide the proposal; tier 2 pairs have a documented way to co-occur and
    are reported as context only.  ``gate`` restricts the pair to cells positive for the gate
    marker in BOTH variants (so the population is the same on both sides of the comparison).
    """

    a: str
    b: str
    tier: int
    why: str
    caveat: str = ""
    gate: str | None = None

    @property
    def key(self) -> str:
        return f"{self.a}|{self.b}" + (f"@{self.gate}+" if self.gate else "")


#: Declared biology, not a finding.  PR #18 (``uniform_audit.DEFAULT_PAIRS``) contributes
#: INS/GCG, INS/SST, GCG/SST, CD3e/CD20, CD3e/CD79a and CD68/CD3e; its Ker8_18/CD45 pair is
#: not usable because CD45 is not in this panel.  The registry declares no pairs
#: (``mutually_exclusive`` is false for every marker, by its no-verified-typing rule).
PAIRS: tuple[Pair, ...] = (
    # endocrine
    Pair("INS", "GCG", 1, "beta vs alpha cells", "rare bihormonal cells, a little commoner in T1D"),
    Pair("INS", "SST", 1, "beta vs delta cells", "rare polyhormonal cells"),
    Pair("GCG", "SST", 1, "alpha vs delta cells", "rare polyhormonal cells"),
    # lymphoid lineages
    Pair("CD3e", "CD20", 1, "T vs B cells", "a small CD20-dim T-cell subset exists"),
    Pair("CD3e", "CD79a", 1, "T vs B/plasma cells"),
    Pair("CD3e", "CD68", 1, "T cells vs macrophages", "efferocytosed T-cell fragments"),
    Pair("CD4", "CD8", 1, "CD4 vs CD8 T cells, scored inside CD3e+ cells only",
         "CD4+CD8+ T cells are ~1-3% of peripheral T cells", gate="CD3e"),
    Pair("CD4", "CD20", 1, "B cells never express CD4"),
    Pair("CD8", "CD20", 1, "B cells never express CD8"),
    # immune / endothelial / stromal vs epithelial
    Pair("CD3e", "E_cadherin", 1, "T cells vs epithelium (intraepithelial T cells touch, "
         "but do not express, E-cadherin)"),
    Pair("CD20", "E_cadherin", 1, "B cells vs epithelium"),
    Pair("CD31", "E_cadherin", 1, "endothelium vs epithelium"),
    Pair("CD34", "E_cadherin", 1, "endothelium / CD34+ stroma vs epithelium"),
    Pair("CD3e", "Pan_Cytokeratin", 1, "T cells vs keratin+ epithelium"),
    Pair("CD20", "Pan_Cytokeratin", 1, "B cells vs keratin+ epithelium"),
    Pair("Ker8_18", "CD3e", 1, "K8/18 epithelium vs T cells",
         "lymph-node K8/18+ reticular cells are not T cells"),
    Pair("Keratin_5", "CD3e", 1, "K5+ epithelium vs T cells"),
    Pair("EpCAM", "CD3e", 1, "epithelium vs T cells"),
    Pair("EpCAM", "CD31", 1, "epithelium vs endothelium"),
    Pair("CD31", "INS", 1, "islet capillary endothelium vs beta cells"),
    Pair("CD34", "CD3e", 1, "endothelium (incl. HEV) / CD34+ stroma vs mature T cells"),
    Pair("SMA", "CD3e", 1, "smooth muscle / pericytes / myofibroblasts vs T cells"),
    Pair("SMA", "INS", 1, "islet pericytes / arteriolar smooth muscle vs beta cells"),
    Pair("Podoplanin", "CD20", 1, "lymphatic endothelium / FRC / FDC vs B cells"),
    Pair("B3TUBB", "CD3e", 1, "nerve fibres vs T cells"),
    Pair("CD206", "CD3e", 1, "CD206+ macrophages / sinus endothelium vs T cells"),
    Pair("CD163", "CD3e", 1, "CD163+ macrophages vs T cells"),
    Pair("CD209", "CD3e", 1, "DC-SIGN+ DCs / macrophages vs T cells"),
    Pair("MPO", "CD3e", 1, "neutrophils vs T cells"),
    # tier 2: plausible exclusivity with a documented exception -> context only
    Pair("Vimentin", "INS", 2, "mesenchymal / endothelial / immune vs beta cells",
         "vimentin reported in stressed endocrine and ductal cells"),
    Pair("Vimentin", "E_cadherin", 2, "mesenchymal vs epithelial",
         "EMT; vimentin+ ductal cells in pancreatitis and T1D"),
    Pair("Caveolin", "CD3e", 2, "caveolae-rich endothelium / stroma vs T cells",
         "caveolin-1 reported in activated T cells"),
    Pair("Caveolin", "CD20", 2, "caveolae-rich endothelium / stroma vs B cells",
         "caveolin-1 reported in B cells"),
    Pair("Iba1", "CD3e", 2, "macrophages vs T cells", "AIF-1 (Iba1) reported in activated T cells"),
    Pair("CD11c", "E_cadherin", 2, "DCs / macrophages vs epithelium",
         "E-cadherin+ DC subsets; E-cadherin induced on IL-4-activated macrophages"),
    Pair("SMA", "CD31", 2, "vascular smooth muscle vs endothelium",
         "endothelial-to-mesenchymal transition in injury and fibrosis"),
)

#: Candidate pairs considered and NOT used, with the reason (written to the README).
REJECTED_PAIRS: tuple[tuple[str, str], ...] = (
    ("CD45 / E-cadherin", "CD45 is not in this panel"),
    ("CD68, Iba1 / INS", "islet macrophages capture beta-cell granules (real co-localisation)"),
    ("CD68 / CD20", "tingible-body macrophages engulf B cells in germinal centres"),
    ("CD68, CD163, CD206, Iba1 / E-cadherin", "E-cadherin is induced on IL-4-activated macrophages"),
    ("CD56 / CD3e, CD56 / hormones", "NKT cells; islet endocrine cells express NCAM (CD56)"),
    ("CD57 / CD3e", "CD57 marks effector T cells as well as NK cells"),
    ("CD11c / CD3e, CD11c / CD20", "CD11c+ effector T cells and CD11c+ age-associated B cells"),
    ("CD21 / CD3e", "CD21 is expressed by a subset of T cells"),
    ("CD66 / immune", "CEACAM1 (CD66a) on activated T cells, CD66b on granulocytes; clone ambiguous"),
    ("CD141 / immune", "thrombomodulin marks cDC1 and monocytes as well as endothelium"),
    ("M2Gal3 / epithelial", "galectin-3 is expressed by epithelium and macrophages"),
    ("Caveolin / INS", "beta cells express caveolin-1"),
    ("Podoplanin / CD3e", "podoplanin is expressed by Th17 cells"),
)

# ----------------------------------------------------------------------------- positivity rules
#: The typing pipeline's rule (scripts/celltype/run_step3.py + score.py): a marker is "shown"
#: when (v - floor) / (top - floor) >= REQUIRE_MIN, floor = the log1p mode of v > 0 (512
#: bins, bin 0 ignored), top = the p99.9.  Typing fits floor/top cohort-wide on UniFORM-
#: aligned values; UniFORM is one multiplicative scale per (donor, marker), so fitting them
#: per donor here is the same rule on each donor's own scale.
REQUIRE_MIN = 0.10
SCALE_Q = 99.9
HIST_BINS = 512
MIN_VALUES_FOR_FLOOR = 1000
#: A rule is used for a (donor, marker) only if its threshold moves < 2x between raw and
#: REDSEA and, where the floor is a material part of the threshold, the floor does too.
#: A larger jump means the rule's background anchor changed population between variants
#: (e.g. REDSEA's zero clamp removes enough background that the log1p mode lands on the
#: positive population), so the two variants would not be thresholded by the same rule.
RULE_STABILITY = 2.0
OTSU_BINS = 512
#: Sensitivity only: the same decision at stricter typing bars.  Typing sweeps REQUIRE_MIN
#: rather than fixing a cut per marker, so a proposal that flips with it is flagged.  Only
#: typing-rule (T) markers move; an Otsu-scored (donor, marker) keeps its threshold.
SWEEP = (0.20, 0.30)

# ----------------------------------------------------------------------------- decision thresholds
MIN_POSITIVE = 50        # raw positives needed for focal and partner
MIN_COPOS = 20           # raw co-positive cells needed before a reduction can be measured
BENEFIT_DROP = 0.20      # net relative drop of the co-positive rate (focal corrected alone)
SMALL_DROP = 0.10        # 10-20% selective, loss-free drops are reported as small_benefit
SELECTIVITY_MIN = 2.0    # focal calls cleared on co-positive cells vs lost on single positives
RETENTION_OK = 0.90      # single-positive retention needed for a clean benefit
RETENTION_HARM = 0.80    # below this REDSEA destroys single-positive signal
WORSE_RISE = 0.20        # co-positive rate rising by this fraction = REDSEA makes it worse
BRIGHT_RETENTION_FLAG = 0.95
MAJORITY = 2.0 / 3.0     # share of evaluable donors that must agree
MIN_DONORS_FINAL = 3     # fewer evaluable donors -> the verdict is provisional

VERDICT_GROUP = {
    "benefit": "benefit",
    "tradeoff": "tradeoff",
    "harm": "adverse",
    "worse": "adverse",
    "small_benefit": "no_gain",      # below the bar; counted separately and flagged
    "no_benefit": "no_gain",
    "nothing_to_fix": "no_gain",
}


def log(msg: str) -> None:
    print(msg, flush=True)


# ----------------------------------------------------------------------------- registry
def sanitize(name: str) -> str:
    """``phenocycler.cells_parquet.sanitize``: the cells-table spelling of a marker."""
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def raw_column(name: str, compartment: str) -> str:
    """Raw QuPath column in ``data/cells`` for a registry marker and compartment."""
    s = sanitize(name)
    return s if compartment == "Cell" else f"{compartment}__{s}"


def load_registry_rows(registry_path: Path) -> pd.DataFrame:
    """One row per registry marker, with the columns this audit reads for it."""
    doc = json.loads(Path(registry_path).read_text())
    rows = []
    for m in doc["markers"]:
        policy = m["spillover_policy"]
        source = m.get("measurement_source",
                       "qupath_raw" if policy == "uncompensated" else "redsea")
        rows.append({
            "marker": m["name"], "kind": m["kind"], "uses": ",".join(m["uses"]),
            "compartment": m["compartment"], "spillover_policy": policy,
            "measurement_source": source,
            "redsea_column": f"{m['compartment']}__{m['name']}",
            "raw_column": raw_column(m["name"], m["compartment"]),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- donors
DONE_PATTERNS = (
    re.compile(r"\[donor-done\] redsea/(\S+)"),
    re.compile(r"^\[(\S+)\] DONE in \d+s", re.M),
)


def donors_done_in_log(log_path: Path) -> set[str]:
    text = Path(log_path).read_text(errors="replace")
    done: set[str] = set()
    for pat in DONE_PATTERNS:
        done.update(pat.findall(text))
    return done


def partition_file(root: Path, donor: str) -> Path | None:
    """The single complete parquet of a donor partition, or None."""
    part = Path(root) / f"donor_id={donor}"
    if not part.is_dir():
        return None
    files = [f for f in part.glob("*.parquet") if not f.name.startswith(".")]
    return files[0] if len(files) == 1 else None


def ready_donors(data_dir: Path, log_path: Path | None, requested: list[str] | None,
                 exclude_filter) -> tuple[list[str], dict[str, str]]:
    """Donors whose REDSEA is complete, plus why each other candidate was not read."""
    redsea_root = data_dir / "redsea"
    # every donor ingested (cells) is a candidate, so the ones still waiting are listed too
    present = sorted({p.name.split("=", 1)[1] for sub in ("cells", "redsea")
                      for p in (data_dir / sub).glob("donor_id=*") if p.is_dir()})
    candidates = requested if requested else present
    done_log = donors_done_in_log(log_path) if log_path is not None else None
    eligible = set(exclude_filter(candidates))
    ready, skipped = [], {}
    for d in candidates:
        if d not in eligible:
            skipped[d] = "centrally excluded (phenocycler.cohort.DONOR_EXCLUSIONS)"
        elif partition_file(redsea_root, d) is None:
            skipped[d] = "REDSEA partition absent or incomplete"
        elif done_log is not None and d not in done_log:
            skipped[d] = "REDSEA not reported done in the run log"
        elif partition_file(data_dir / "cells", d) is None:
            skipped[d] = "cells partition absent"
        elif partition_file(data_dir / "geometry_qc", d) is None:
            skipped[d] = "geometry_qc partition absent"
        else:
            ready.append(d)
    return ready, skipped


def donor_status_map(path: Path | None) -> dict[str, str]:
    """``donor_id -> disease.status`` for context; unknown donors map to ''."""
    if path is None or not Path(path).exists():
        return {}
    t = pd.read_csv(path, dtype=str)
    if "donor_id" not in t or "disease.status" not in t:
        return {}
    return {_norm(d): s for d, s in zip(t["donor_id"], t["disease.status"])}


def _norm(d: str) -> str:
    s = str(d).strip()
    return str(int(s)) if s.isdecimal() else s


# ----------------------------------------------------------------------------- rules
@dataclass(frozen=True)
class RuleFit:
    rule: str          # "T" (typing) or "O" (Otsu on log1p of v > 0)
    floor: float       # background anchor (T: log1p mode; O: background-class centre)
    top: float         # p99.9
    threshold: float   # positive iff v >= threshold (T) / v > threshold (O)
    ok: bool


def typing_rule(values) -> RuleFit:
    """The typing pipeline's positivity rule, fitted on one variant of one donor."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return RuleFit("T", np.nan, np.nan, np.nan, False)
    pos = v[v > 0]
    floor, ok = 0.0, pos.size >= MIN_VALUES_FOR_FLOOR
    if ok:
        h, e = np.histogram(np.log1p(pos), bins=HIST_BINS)
        floor = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))
    top = float(np.percentile(v, SCALE_Q))
    span = top - floor
    span = span if np.isfinite(span) and span > 0 else 1.0
    return RuleFit("T", floor, top, floor + REQUIRE_MIN * span, ok)


def otsu_threshold(x: np.ndarray, nbins: int = OTSU_BINS) -> float:
    """``skimage.filters.threshold_otsu`` (same histogram and tie convention)."""
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0 or np.all(x == x[0]):
        return float("nan")
    hist, edges = np.histogram(x, bins=nbins, range=(x.min(), x.max()))
    centers = (edges[:-1] + edges[1:]) / 2.0
    hist = hist.astype(np.float64)
    w1 = np.cumsum(hist)
    w2 = np.cumsum(hist[::-1])[::-1]
    m1 = np.cumsum(hist * centers) / np.where(w1 > 0, w1, 1.0)
    m2 = (np.cumsum((hist * centers)[::-1]) / np.where(w2[::-1] > 0, w2[::-1], 1.0))[::-1]
    var12 = w1[:-1] * w2[1:] * (m1[:-1] - m2[1:]) ** 2
    return float(centers[int(np.argmax(var12))])


def otsu_rule(values) -> RuleFit:
    """Otsu on log1p(v) among v > 0; positive iff v > threshold."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    pos = v[v > 0]
    if pos.size < MIN_VALUES_FOR_FLOOR:
        return RuleFit("O", np.nan, np.nan, np.nan, False)
    lx = np.log1p(pos)
    t = otsu_threshold(lx)
    if not np.isfinite(t):
        return RuleFit("O", np.nan, np.nan, np.nan, False)
    below = lx[lx <= t]
    floor = float(np.expm1(below.mean())) if below.size else 0.0
    return RuleFit("O", floor, float(np.percentile(v, SCALE_Q)), float(np.expm1(t)), True)


def _ratio(a: float, b: float) -> float:
    if not (np.isfinite(a) and np.isfinite(b)) or a <= 0 or b <= 0:
        return float("inf")
    return max(a, b) / min(a, b)


def rule_stable(raw: RuleFit, rs: RuleFit) -> bool:
    if not (raw.ok and rs.ok):
        return False
    if _ratio(raw.threshold, rs.threshold) > RULE_STABILITY:
        return False
    floor_matters = max(raw.floor / raw.threshold, rs.floor / rs.threshold) >= 0.5
    return not (floor_matters and _ratio(raw.floor, rs.floor) > RULE_STABILITY)


def choose_rule(raw_values, rs_values) -> tuple[str, RuleFit, RuleFit, str]:
    """Typing rule where it is stable across variants, else Otsu, else none."""
    tr, ts = typing_rule(raw_values), typing_rule(rs_values)
    if rule_stable(tr, ts):
        return "T", tr, ts, ""
    note = (f"typing rule unstable (floor {tr.floor:.4g}->{ts.floor:.4g}, "
            f"threshold {tr.threshold:.4g}->{ts.threshold:.4g})")
    orr, ors = otsu_rule(raw_values), otsu_rule(rs_values)
    if rule_stable(orr, ors):
        return "O", orr, ors, note
    return "none", tr, ts, note + "; Otsu unstable too"


def positive(values: np.ndarray, fit: RuleFit) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        if fit.rule == "T":
            return np.isfinite(v) & (v >= fit.threshold)
        return np.isfinite(v) & (v > fit.threshold)


def count_matched(values: np.ndarray, k: int) -> np.ndarray:
    """The k highest finite values (ties broken by position): REDSEA at raw's positive count."""
    v = np.asarray(values, dtype=np.float64)
    out = np.zeros(v.size, dtype=bool)
    idx = np.flatnonzero(np.isfinite(v))
    k = int(min(max(k, 0), idx.size))
    if k:
        order = np.argsort(-v[idx], kind="stable")[:k]
        out[idx[order]] = True
    return out


def scaled(values: np.ndarray, fit: RuleFit) -> np.ndarray:
    """``score.apply_scale`` on the rule's own floor/top (threshold-light measures)."""
    span = fit.top - fit.floor
    span = span if np.isfinite(span) and span > 0 else 1.0
    return np.clip((np.asarray(values, dtype=np.float64) - fit.floor) / span, 0.0, 1.0)


# ----------------------------------------------------------------------------- metrics
def copos_union(a: np.ndarray, b: np.ndarray) -> float:
    u = int((a | b).sum())
    return float((a & b).sum()) / u if u else float("nan")


def _frac(mask_num: np.ndarray, mask_den: np.ndarray) -> float:
    d = int(mask_den.sum())
    return float((mask_num & mask_den).sum()) / d if d else float("nan")


def svd_ratio(x: np.ndarray, y: np.ndarray) -> float:
    """PR #18 ``exclusivity_ratio``: sigma2/sigma1 of the uncentred [x, y] matrix."""
    gram = np.array([[x @ x, x @ y], [x @ y, y @ y]], dtype=np.float64)
    lam = np.clip(np.linalg.eigvalsh(gram), 0.0, None)
    return float(np.sqrt(lam[0] / lam[1])) if lam[1] > 0 else float("nan")


def cosine(x: np.ndarray, y: np.ndarray) -> float:
    """Uncentred cosine of the pair: 0 for perfectly exclusive, balance-free."""
    d = float(np.sqrt((x @ x) * (y @ y)))
    return float(x @ y) / d if d > 0 else float("nan")


def selectivity(cleared: float, retention: float) -> float:
    lost = 1.0 - retention
    if not (np.isfinite(cleared) and np.isfinite(lost)):
        return float("nan")
    return min(cleared / max(lost, 1e-3), 999.0)


def verdict(n_focal: int, n_partner: int, n_copos: int, drop: float, retention: float,
            select: float, rules_ok: bool) -> str:
    """One (donor, focal marker, partner) verdict; thresholds are the module constants."""
    if not rules_ok or n_focal < MIN_POSITIVE or n_partner < MIN_POSITIVE:
        return "insufficient"
    if n_copos < MIN_COPOS:
        return "nothing_to_fix"
    if not np.isfinite(retention):          # no raw single positive to retain
        return "insufficient"
    if retention < RETENTION_HARM:
        return "harm"
    if np.isfinite(drop) and drop <= -WORSE_RISE:
        return "worse"
    selective = np.isfinite(select) and select >= SELECTIVITY_MIN
    if np.isfinite(drop) and drop >= BENEFIT_DROP and selective:
        return "benefit" if retention >= RETENTION_OK else "tradeoff"
    if np.isfinite(drop) and drop >= SMALL_DROP and selective and retention >= RETENTION_OK:
        return "small_benefit"
    return "no_benefit"


def _median_ratio(v_rs, v_raw, mask) -> float:
    if not mask.any():
        return float("nan")
    x, y = v_raw[mask], v_rs[mask]
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    return float(np.median(y[ok] / x[ok])) if ok.any() else float("nan")


def focal_metrics(prefix: str, f_raw, f_rs, f_rsk, p_raw, p_rs, v_raw, v_rs=None) -> dict:
    """Metrics for correcting the focal marker alone (partner held raw)."""
    cu_raw = copos_union(f_raw, p_raw)
    single = f_raw & ~p_raw
    cop = f_raw & p_raw
    out = {}
    if v_rs is not None:      # how much REDSEA changes the focal value where it should act
        out[f"{prefix}_ratio_in_copos"] = _median_ratio(v_rs, v_raw, cop)
        out[f"{prefix}_ratio_in_single"] = _median_ratio(v_rs, v_raw, single)
        out[f"{prefix}_zeroed_in_copos"] = float((v_rs[cop] == 0).mean()) if cop.any() \
            else float("nan")
    out.update({
        f"{prefix}_n_single_raw": int(single.sum()),
        f"{prefix}_cu_corrected": copos_union(f_rs, p_raw),
        f"{prefix}_cu_corrected_k": copos_union(f_rsk, p_raw),
    })
    for tag, fr in (("", f_rs), ("_k", f_rsk)):
        cu_f = copos_union(fr, p_raw)
        ret = _frac(fr, single)
        cleared = 1.0 - _frac(fr, cop) if cop.any() else float("nan")
        out[f"{prefix}_drop{tag}"] = (1.0 - cu_f / cu_raw) if cu_raw and np.isfinite(cu_raw) \
            and cu_raw > 0 else float("nan")
        out[f"{prefix}_retention{tag}"] = ret
        out[f"{prefix}_cleared{tag}"] = cleared
        out[f"{prefix}_selectivity{tag}"] = selectivity(cleared, ret)
    fin = v_raw[f_raw]
    cut = float(np.median(fin)) if fin.size else float("nan")
    bright = single & (v_raw >= cut)
    out[f"{prefix}_retention_bright"] = _frac(f_rs, bright)
    out[f"{prefix}_n_bright_raw"] = int(bright.sum())
    # the co-positive share of the focal marker's own calls (impossible-call fraction)
    out[f"{prefix}_impossible_share_raw"] = _frac(p_raw, f_raw)
    out[f"{prefix}_impossible_share_corrected"] = _frac(p_raw, f_rs)
    return out


def pair_row(pair: Pair, universe: np.ndarray, vals: dict, fits: dict, pos: dict) -> dict:
    """Every raw/REDSEA combination for one pair on one donor, as one wide row."""
    a, b = pair.a, pair.b
    U = universe
    A = {k: pos[(a, k)][U] for k in ("raw", "rs", "rsk")}
    B = {k: pos[(b, k)][U] for k in ("raw", "rs", "rsk")}
    row = {"n_cells": int(U.sum())}
    for m, P in (("a", A), ("b", B)):
        row[f"n_{m}_raw"] = int(P["raw"].sum())
        row[f"n_{m}_redsea"] = int(P["rs"].sum())
    combos = {"raw_raw": (A["raw"], B["raw"]), "redsea_raw": (A["rs"], B["raw"]),
              "raw_redsea": (A["raw"], B["rs"]), "redsea_redsea": (A["rs"], B["rs"])}
    for name, (x, y) in combos.items():
        row[f"n_copos_{name}"] = int((x & y).sum())
        row[f"copos_union_{name}"] = copos_union(x, y)
    row.update(focal_metrics("a", A["raw"], A["rs"], A["rsk"], B["raw"], B["rs"],
                             vals[(a, "raw")][U], vals[(a, "rs")][U]))
    row.update(focal_metrics("b", B["raw"], B["rs"], B["rsk"], A["raw"], A["rs"],
                             vals[(b, "raw")][U], vals[(b, "rs")][U]))
    # threshold-light: PR #18 SVD ratio and the uncentred cosine on each rule's own scale
    sa = {k: scaled(vals[(a, k)][U], fits[(a, k)]) for k in ("raw", "rs")}
    sb = {k: scaled(vals[(b, k)][U], fits[(b, k)]) for k in ("raw", "rs")}
    for name, (x, y) in {"raw_raw": (sa["raw"], sb["raw"]), "redsea_raw": (sa["rs"], sb["raw"]),
                         "raw_redsea": (sa["raw"], sb["rs"]),
                         "redsea_redsea": (sa["rs"], sb["rs"])}.items():
        row[f"svd_ratio_{name}"] = svd_ratio(x, y)
        row[f"cosine_{name}"] = cosine(x, y)
    return row


def positives_at_bar(vals: dict, fits: dict, rule_used: dict, q: float, n: int) -> dict:
    """Raw/REDSEA positives with the typing bar at ``q`` (T-rule markers only move)."""
    out = {}
    for m, rule in rule_used.items():
        for var in ("raw", "rs"):
            f = fits[(m, var)]
            if rule == "T":
                span = f.top - f.floor
                span = span if np.isfinite(span) and span > 0 else 1.0
                out[(m, var)] = positive(vals[(m, var)], RuleFit("T", f.floor, f.top,
                                                                  f.floor + q * span, True))
            elif rule == "O":
                out[(m, var)] = positive(vals[(m, var)], f)
            else:
                out[(m, var)] = np.zeros(n, bool)
    return out


def sweep_metrics(pair: Pair, pos_q: dict, vals: dict, base_universe: np.ndarray) -> dict:
    """Focal verdict inputs for both markers of a pair at one sweep bar."""
    U = base_universe.copy()
    if pair.gate:
        U &= pos_q[(pair.gate, "raw")] & pos_q[(pair.gate, "rs")]
    out = {}
    for focal, partner in ((pair.a, pair.b), (pair.b, pair.a)):
        f_raw, f_rs = pos_q[(focal, "raw")][U], pos_q[(focal, "rs")][U]
        p_raw = pos_q[(partner, "raw")][U]
        m = focal_metrics("x", f_raw, f_rs, f_rs, p_raw, pos_q[(partner, "rs")][U],
                          vals[(focal, "raw")][U])
        out[focal] = {"n_f": int(f_raw.sum()), "n_p": int(p_raw.sum()),
                      "n_cp": int((f_raw & p_raw).sum()), "drop": m["x_drop"],
                      "retention": m["x_retention"], "selectivity": m["x_selectivity"]}
    return out


# ----------------------------------------------------------------------------- per donor
def pair_markers(pairs=PAIRS) -> list[str]:
    names: list[str] = []
    for p in pairs:
        for m in (p.a, p.b, p.gate):
            if m and m not in names:
                names.append(m)
    return names


@dataclass
class DonorFrame:
    """One donor's QC-eligible cells, aligned on the cells table's ``object_id`` order."""

    files: dict
    cells_cols: set
    rs_cols: set
    order_rs: np.ndarray       # position in the REDSEA file of each cells-table row
    keep: np.ndarray           # QC-eligible, in cells-table order
    n_total: int

    def columns(self, which: str, names: list[str]) -> dict[str, np.ndarray]:
        """Float columns of one file, in cells-table order, QC-eligible rows only."""
        if not names:
            return {}
        t = pq.read_table(self.files[which], columns=list(dict.fromkeys(names)))
        out = {}
        for name in names:
            v = np.asarray(t.column(name).to_numpy(zero_copy_only=False), dtype=np.float64)
            if which == "redsea":
                v = v[self.order_rs]
            out[name] = v[self.keep]
        return out


def open_donor(data_dir: Path, donor: str) -> DonorFrame:
    """Align cells, REDSEA and geometry QC on ``object_id``; the universes must be equal."""
    files = {k: partition_file(data_dir / k, donor) for k in ("cells", "redsea", "geometry_qc")}
    oid_c = pd.Index(pq.read_table(files["cells"], columns=["object_id"]).column(0)
                     .to_pandas().astype(str))
    oid_r = pd.Index(pq.read_table(files["redsea"], columns=["object_id"]).column(0)
                     .to_pandas().astype(str))
    qc = pq.read_table(files["geometry_qc"],
                       columns=["object_id", "analysis_eligible"]).to_pandas()
    if not (oid_c.is_unique and oid_r.is_unique and qc["object_id"].is_unique):
        raise ValueError(f"donor {donor}: duplicate object_id")
    order_rs = oid_r.get_indexer(oid_c)
    elig = pd.Series(qc["analysis_eligible"].fillna(False).astype(bool).to_numpy(),
                     index=qc["object_id"].astype(str)).reindex(oid_c)
    if (order_rs < 0).any() or elig.isna().any() or len(oid_r) != len(oid_c) \
            or len(qc) != len(oid_c):
        raise ValueError(f"donor {donor}: object universes differ (cells {len(oid_c)}, "
                         f"redsea {len(oid_r)}, geometry_qc {len(qc)})")
    return DonorFrame(files=files, cells_cols=set(pq.read_schema(files["cells"]).names),
                      rs_cols=set(pq.read_schema(files["redsea"]).names), order_rs=order_rs,
                      keep=elig.to_numpy(bool), n_total=len(oid_c))


def read_donor(data_dir: Path, donor: str, reg: pd.DataFrame, markers: list[str]
               ) -> tuple[dict, list[str], dict]:
    """Raw and REDSEA arrays of the acquired ``markers`` over QC-eligible cells."""
    dfm = open_donor(data_dir, donor)
    r = reg.set_index("marker")
    cand = [m for m in markers if m in r.index and r.at[m, "raw_column"] in dfm.cells_cols
            and r.at[m, "redsea_column"] in dfm.rs_cols]
    raw = dfm.columns("cells", [r.at[m, "raw_column"] for m in cand])
    red = dfm.columns("redsea", [r.at[m, "redsea_column"] for m in cand])
    vals, have = {}, []
    for m in cand:
        x = raw[r.at[m, "raw_column"]]
        if not np.isfinite(x).any():      # union column of a marker this panel lacks
            continue
        vals[(m, "raw")], vals[(m, "rs")] = x, red[r.at[m, "redsea_column"]]
        have.append(m)
    meta = {"cells_file": str(dfm.files["cells"]), "redsea_file": str(dfm.files["redsea"]),
            "geometry_qc_file": str(dfm.files["geometry_qc"]),
            "n_cells_total": int(dfm.n_total), "n_cells_qc": int(dfm.keep.sum())}
    return vals, have, meta


def clamp_rates(data_dir: Path, donor: str, reg: pd.DataFrame, batch: int = 16) -> pd.DataFrame:
    """For every REDSEA-sourced marker: share of raw>0 QC cells REDSEA clamps to exactly 0."""
    dfm = open_donor(data_dir, donor)
    rows = []
    sel = reg[reg["measurement_source"] == "redsea"]
    todo = [(m, rc, sc) for m, rc, sc in zip(sel["marker"], sel["raw_column"],
                                              sel["redsea_column"])
            if rc in dfm.cells_cols and sc in dfm.rs_cols]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        raw = dfm.columns("cells", [c[1] for c in chunk])
        red = dfm.columns("redsea", [c[2] for c in chunk])
        for m, rc, sc in chunk:
            x, y = raw[rc], red[sc]
            if not np.isfinite(x).any():       # declared column, marker not acquired
                continue
            ok = np.isfinite(x) & np.isfinite(y)
            rawpos = ok & (x > 0)
            rows.append({"marker": m,
                         "frac_zeroed_given_raw_pos": float(((y == 0) & rawpos).sum())
                         / max(int(rawpos.sum()), 1),
                         "median_ratio_redsea_raw": float(np.median(y[rawpos] / x[rawpos]))
                         if rawpos.any() else float("nan")})
        del raw, red
    return pd.DataFrame(rows)


def process_donor(data_dir: Path, donor: str, reg: pd.DataFrame, pairs=PAIRS) -> dict:
    t0 = time.time()
    markers = pair_markers(pairs)
    vals, have, meta = read_donor(data_dir, donor, reg, markers)
    n = int(meta["n_cells_qc"])
    fits, pos, rule_used, notes = {}, {}, {}, {}
    th_rows = []
    for m in have:
        x, y = vals[(m, "raw")], vals[(m, "rs")]
        rule, fr, fs, note = choose_rule(x, y)
        rule_used[m], notes[m] = rule, note
        fits[(m, "raw")], fits[(m, "rs")] = fr, fs
        pos[(m, "raw")] = positive(x, fr) if rule != "none" else np.zeros(n, bool)
        pos[(m, "rs")] = positive(y, fs) if rule != "none" else np.zeros(n, bool)
        pos[(m, "rsk")] = count_matched(y, int(pos[(m, "raw")].sum())) if rule != "none" \
            else np.zeros(n, bool)
        for var, f in (("raw", fr), ("redsea", fs)):
            p = pos[(m, "raw" if var == "raw" else "rs")]
            th_rows.append({"donor_id": donor, "marker": m, "variant": var, "rule": rule,
                            "floor": f.floor, "top": f.top, "threshold": f.threshold,
                            "n_positive": int(p.sum()), "frac_positive": float(p.mean()) if n
                            else float("nan"), "rule_note": note})
    sweep_pos = {q: positives_at_bar(vals, fits, rule_used, q, n) for q in SWEEP}
    pair_rows, focal_rows = [], []
    for pair in pairs:
        need = [pair.a, pair.b] + ([pair.gate] if pair.gate else [])
        missing = [m for m in need if m not in have]
        base = {"donor_id": donor, "pair": pair.key, "a": pair.a, "b": pair.b,
                "gate": pair.gate or "", "tier": pair.tier}
        if missing:
            pair_rows.append({**base, "status": "not_acquired:" + ",".join(missing)})
            continue
        rules_ok = all(rule_used[m] != "none" for m in need)
        base.update({"rule_a": rule_used[pair.a], "rule_b": rule_used[pair.b]})
        universe = np.ones(n, bool)
        for m in (pair.a, pair.b):
            universe &= np.isfinite(vals[(m, "raw")]) & np.isfinite(vals[(m, "rs")])
        base_universe = universe.copy()
        if pair.gate:
            universe &= pos[(pair.gate, "raw")] & pos[(pair.gate, "rs")]
        if not rules_ok:
            pair_rows.append({**base, "status": "no_usable_rule",
                              "n_cells": int(universe.sum())})
            continue
        row = pair_row(pair, universe, vals, fits, pos)
        sweep = {q: sweep_metrics(pair, sweep_pos[q], vals, base_universe) for q in SWEEP}
        for m, focal, partner in (("a", pair.a, pair.b), ("b", pair.b, pair.a)):
            o = "b" if m == "a" else "a"
            n_f, n_p = row[f"n_{m}_raw"], row[f"n_{o}_raw"]
            n_cp = row["n_copos_raw_raw"]
            v = verdict(n_f, n_p, n_cp, row[f"{m}_drop"], row[f"{m}_retention"],
                        row[f"{m}_selectivity"], rules_ok)
            vk = verdict(n_f, n_p, n_cp, row[f"{m}_drop_k"], row[f"{m}_retention_k"],
                         row[f"{m}_selectivity_k"], rules_ok)
            row[f"verdict_{m}"], row[f"verdict_{m}_k"] = v, vk
            corr = "redsea_raw" if m == "a" else "raw_redsea"
            focal_rows.append({
                "donor_id": donor, "marker": focal, "partner": partner, "pair": pair.key,
                "gate": pair.gate or "", "tier": pair.tier, "rule": rule_used[focal],
                "rule_partner": rule_used[partner], "n_cells": row["n_cells"],
                "n_marker_raw": n_f, "n_marker_redsea": row[f"n_{m}_redsea"],
                "n_partner_raw": n_p, "n_copos_raw": n_cp,
                "n_copos_marker_corrected": row[f"n_copos_{corr}"],
                "copos_union_raw": row["copos_union_raw_raw"],
                "copos_union_marker_corrected": row[f"{m}_cu_corrected"],
                "copos_union_both_corrected": row["copos_union_redsea_redsea"],
                "copos_drop": row[f"{m}_drop"], "retention": row[f"{m}_retention"],
                "retention_bright": row[f"{m}_retention_bright"],
                "cleared": row[f"{m}_cleared"], "selectivity": row[f"{m}_selectivity"],
                "redsea_ratio_in_copos": row[f"{m}_ratio_in_copos"],
                "redsea_ratio_in_single": row[f"{m}_ratio_in_single"],
                "redsea_zeroed_in_copos": row[f"{m}_zeroed_in_copos"],
                "impossible_share_raw": row[f"{m}_impossible_share_raw"],
                "impossible_share_corrected": row[f"{m}_impossible_share_corrected"],
                "copos_drop_k": row[f"{m}_drop_k"], "retention_k": row[f"{m}_retention_k"],
                "selectivity_k": row[f"{m}_selectivity_k"],
                "cosine_raw": row["cosine_raw_raw"], "cosine_marker_corrected":
                row[f"cosine_{corr}"],
                "svd_ratio_raw": row["svd_ratio_raw_raw"],
                "svd_ratio_marker_corrected": row[f"svd_ratio_{corr}"],
                "verdict": v, "verdict_k": vk})
            for q in SWEEP:
                sq, tag = sweep[q][focal], f"q{int(round(q * 100))}"
                vq = verdict(sq["n_f"], sq["n_p"], sq["n_cp"], sq["drop"], sq["retention"],
                             sq["selectivity"], rules_ok)
                focal_rows[-1].update({f"n_marker_raw_{tag}": sq["n_f"],
                                       f"n_copos_raw_{tag}": sq["n_cp"],
                                       f"copos_drop_{tag}": sq["drop"],
                                       f"retention_{tag}": sq["retention"],
                                       f"selectivity_{tag}": sq["selectivity"],
                                       f"verdict_{tag}": vq})
                row[f"verdict_{m}_{tag}"] = vq
        pair_rows.append({**base, "status": "ok", **row})
    meta.update({"donor_id": donor, "data_dir": str(data_dir), "markers_read": have,
                 "seconds": round(time.time() - t0, 1)})
    return {"pairs": pd.DataFrame(pair_rows), "focal": pd.DataFrame(focal_rows),
            "thresholds": pd.DataFrame(th_rows), "meta": meta}


def write_donor_cache(out: Path, donor: str, result: dict, clamps: pd.DataFrame) -> None:
    """Atomic per-donor cache: files land in a temp dir that is renamed into place."""
    final = out / "per_donor" / f"donor_id={donor}"
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=f".donor_id={donor}.", dir=final.parent))
    try:
        result["pairs"].to_csv(tmp / "pairs.tsv", sep="\t", index=False)
        result["focal"].to_csv(tmp / "focal.tsv", sep="\t", index=False)
        th = result["thresholds"]
        if len(clamps):
            th = th.merge(clamps, on="marker", how="left")
        th.to_csv(tmp / "thresholds.tsv", sep="\t", index=False)
        clamps.assign(donor_id=donor).to_csv(tmp / "clamp_rates.tsv", sep="\t", index=False)
        meta = dict(result["meta"], script_version=SCRIPT_VERSION,
                    written=_dt.datetime.now().isoformat(timespec="seconds"))
        (tmp / "meta.json").write_text(json.dumps(meta, indent=1))
        if final.exists():
            shutil.rmtree(final)
        os.replace(tmp, final)
    finally:
        if tmp.exists():
            shutil.rmtree(tmp)


def cached_source(out: Path, donor: str) -> str | None:
    meta = out / "per_donor" / f"donor_id={donor}" / "meta.json"
    if not meta.exists():
        return None
    try:
        return json.loads(meta.read_text()).get("data_dir")
    except json.JSONDecodeError:
        return None


# ----------------------------------------------------------------------------- aggregation
TEXT_COLUMNS = ("gate", "rule", "rule_partner", "rule_a", "rule_b", "rule_note", "status",
                "pair", "a", "b", "marker", "partner", "verdict", "verdict_k", "verdict_a",
                "verdict_b", "verdict_a_k", "verdict_b_k", "variant",
                *(f"verdict_q{int(round(q * 100))}" for q in SWEEP),
                *(f"verdict_{m}_q{int(round(q * 100))}" for m in "ab" for q in SWEEP))


def _read_all(out: Path, name: str) -> pd.DataFrame:
    frames = []
    for p in sorted((out / "per_donor").glob(f"donor_id=*/{name}")):
        try:
            f = pd.read_csv(p, sep="\t", dtype={"donor_id": str})
        except pd.errors.EmptyDataError:
            continue
        for c in TEXT_COLUMNS:
            if c in f:
                f[c] = f[c].fillna("").astype(str)
        if len(f):
            frames.append(f)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def majority(verdicts: list[str]) -> tuple[str, int, dict]:
    """Majority verdict group over evaluable donors, or 'donor-dependent'."""
    groups = [VERDICT_GROUP[v] for v in verdicts if v in VERDICT_GROUP]
    n = len(groups)
    counts = {g: groups.count(g) for g in dict.fromkeys(groups)}
    if not n:
        return "insufficient", 0, counts
    best = max(counts, key=counts.get)
    return (best if counts[best] / n >= MAJORITY else "donor-dependent"), n, counts


def _label(g: pd.DataFrame) -> str:
    gate = g["gate"].iloc[0]
    return g["partner"].iloc[0] + (f"@{gate}+" if gate else "")


def decide(t1: pd.DataFrame, t2: pd.DataFrame, vcol: str = "verdict",
           kcol: str | None = "verdict_k") -> dict:
    """The marker-level decision from per-(donor, partner) verdicts in column ``vcol``."""
    major, major_k, share, parts = {}, {}, {}, []
    for _, g in t1.groupby("pair", sort=True):
        name = _label(g)
        maj, n, counts = majority(list(g[vcol]))
        major[name] = maj
        if kcol:
            major_k[name] = majority(list(g[kcol]))[0]
        share[name] = {k: v / n for k, v in counts.items()} if n else {}
        parts.append(f"{name}:{maj}(" + "/".join(f"{k}={v}" for k, v in counts.items()) + ")")
    t2major = {_label(g): majority(list(g[vcol]))[0] for _, g in t2.groupby("pair", sort=True)}
    flags: list[str] = []
    majors = set(major.values())
    donor_dep = sorted(p for p, m in major.items() if m == "donor-dependent")
    dep_benefit = [p for p in donor_dep if share[p].get("benefit", 0) >= 1 / 3]
    if not len(t1):
        proposal, reason = "raw", "no tier-1 exclusive pair (cannot earn REDSEA)"
        if "benefit" in t2major.values():
            proposal, reason = "judgement", "only tier-2 pairs, and they show a benefit"
    elif majors <= {"insufficient"}:
        proposal, reason = "raw", "pairs present but never evaluable (too few positives)"
    elif "benefit" in majors:
        conflict = majors & {"adverse", "tradeoff", "donor-dependent"}
        if conflict:
            proposal = "judgement"
            reason = "benefit on some partners, " + "/".join(sorted(conflict)) + " on others"
        else:
            proposal, reason = "keep REDSEA", "fewer impossible co-positives, single positives kept"
        mism = [p for p, m in major.items()
                if kcol and (m == "benefit") != (major_k.get(p) == "benefit")]
        if mism:
            flags.append("benefit not reproduced at matched positive count: " + ",".join(mism))
            if proposal == "keep REDSEA":
                proposal, reason = "judgement", reason + "; not reproduced at matched count"
    elif dep_benefit:
        proposal, reason = "judgement", "benefit in some donors only (" + ",".join(dep_benefit) + ")"
    elif "tradeoff" in majors:
        proposal, reason = "judgement", "co-positives fall but >10% of single positives lost"
    elif "adverse" in majors:
        proposal, reason = "raw", "REDSEA destroys single positives or adds co-positives"
    elif donor_dep:
        proposal, reason = "raw", "donors disagree and none shows a benefit majority"
    else:
        proposal, reason = "raw", "no measurable reduction of impossible co-positives"
    if donor_dep:
        flags.append("donor-dependent on " + ",".join(donor_dep))
    return {"proposal": proposal, "reason": reason, "flags": flags, "major": major,
            "share": share, "parts": parts, "majors": majors,
            "benefit_partners": [p for p, m in major.items() if m == "benefit"]}


def summarize_marker(ev: pd.DataFrame) -> dict:
    """Per-marker proposal from its tier-1 focal rows (all donors); tier 2 as context."""
    out: dict = {}
    t1 = ev[ev["tier"] == 1]
    t2 = ev[ev["tier"] == 2]
    d = decide(t1, t2)
    t2parts = []
    for _, g in t2.groupby("pair", sort=True):
        maj, _, counts = majority(list(g["verdict"]))
        t2parts.append(f"{_label(g)}:{maj}(" + "/".join(f"{k}={v}" for k, v in counts.items())
                       + ")")
    # medians come from the tier-1 rows; a marker with only tier-2 pairs shows those instead
    basis = t1 if len(t1) else t2
    evaluable = basis[basis["verdict"] != "insufficient"]
    out["evidence_tier"] = (1 if len(t1) else 2) if len(basis) else ""
    out["n_donors_evaluable"] = int(evaluable["donor_id"].nunique())
    out["tier1_evidence"] = "; ".join(d["parts"])
    out["tier2_evidence"] = "; ".join(t2parts)
    for col, src in (("median_copos_drop", "copos_drop"), ("median_retention", "retention"),
                     ("median_retention_bright", "retention_bright"),
                     ("median_selectivity", "selectivity"),
                     ("median_copos_drop_k", "copos_drop_k"),
                     ("median_cosine_raw", "cosine_raw"),
                     ("median_cosine_corrected", "cosine_marker_corrected")):
        v = evaluable[src].replace([np.inf, -np.inf], np.nan).dropna() \
            if len(evaluable) and src in evaluable else pd.Series(dtype=float)
        out[col] = float(v.median()) if len(v) else float("nan")
    proposal, reason, flags = d["proposal"], d["reason"], list(d["flags"])
    small = [p for p in d["major"] if d["share"][p].get("no_gain", 0)
             and _small_share(t1, p) >= 0.5]
    if small and proposal == "raw":
        flags.append("borderline: selective loss-free 10-20% drop on " + ",".join(small))
    out["borderline_partners"] = ",".join(small) if proposal == "raw" else ""
    # the same decision at stricter typing bars (sensitivity; T-rule markers only move)
    for q in SWEEP:
        col = f"verdict_q{int(round(q * 100))}"
        if col in ev and len(t1):
            dq = decide(t1, t2, col, None)
            out[f"proposal_q{int(round(q * 100))}"] = dq["proposal"]
            if dq["proposal"] != proposal:
                flags.append(f"at typing bar {q:g}: {dq['proposal']}"
                             + (f" (benefit on {','.join(dq['benefit_partners'])})"
                                if dq["benefit_partners"] else ""))
    if out["n_donors_evaluable"] and out["n_donors_evaluable"] < MIN_DONORS_FINAL:
        flags.append(f"provisional: {out['n_donors_evaluable']} evaluable donor(s)")
    rb = out["median_retention_bright"]
    if np.isfinite(rb) and rb < BRIGHT_RETENTION_FLAG:
        flags.append(f"loses bright single positives (median {rb:.2f})")
    elif "adverse" in d["majors"] and np.isfinite(rb):
        flags.append(f"losses are dim calls (bright retention {rb:.2f})")
    unstable = sorted(ev.loc[ev["rule"] != "T", "donor_id"].unique())
    out["typing_rule_unstable_donors"] = len(unstable)
    if unstable:
        flags.append(f"typing rule's floor jumps under REDSEA in {len(unstable)} donor(s): "
                     + ",".join(unstable) + " (scored with Otsu there)")
    out.update(proposal=proposal, reason=reason, flags="; ".join(flags))
    return out


def _small_share(t1: pd.DataFrame, partner_label: str) -> float:
    """Share of a partner's evaluable donors whose verdict is small_benefit."""
    lab = t1["partner"] + np.where(t1["gate"] != "", "@" + t1["gate"] + "+", "")
    g = t1[(lab == partner_label) & (t1["verdict"] != "insufficient")]
    return float((g["verdict"] == "small_benefit").mean()) if len(g) else 0.0


def build_summary(reg: pd.DataFrame, focal: pd.DataFrame, clamps: pd.DataFrame,
                  pairs=PAIRS) -> pd.DataFrame:
    rows = []
    t1 = {m for p in pairs if p.tier == 1 for m in (p.a, p.b)}
    t2 = {m for p in pairs if p.tier == 2 for m in (p.a, p.b)}
    for r in reg.itertuples(index=False):
        m = r.marker
        partners1 = sorted({(p.b if p.a == m else p.a) + (f"@{p.gate}+" if p.gate else "")
                            for p in pairs if p.tier == 1 and m in (p.a, p.b)})
        partners2 = sorted({p.b if p.a == m else p.a for p in pairs
                            if p.tier == 2 and m in (p.a, p.b)})
        row = {"marker": m, "kind": r.kind, "uses": r.uses, "compartment": r.compartment,
               "registry_policy": r.spillover_policy, "in_tier1_pair": m in t1,
               "in_tier2_pair": m in t2, "tier1_partners": ",".join(partners1),
               "tier2_partners": ",".join(partners2)}
        c = clamps[clamps["marker"] == m] if len(clamps) else pd.DataFrame()
        row["n_donors_acquired"] = int(c["donor_id"].nunique()) if len(c) else 0
        row["median_frac_zeroed_by_redsea"] = float(c["frac_zeroed_given_raw_pos"].median()) \
            if len(c) else float("nan")
        row["median_ratio_redsea_raw"] = float(c["median_ratio_redsea_raw"].median()) \
            if len(c) else float("nan")
        if r.spillover_policy != "redsea_corrected":
            row.update(proposal=f"unchanged ({r.spillover_policy})",
                       reason="REDSEA does not compensate this marker", flags="")
        else:
            ev = focal[focal["marker"] == m] if len(focal) else pd.DataFrame()
            if len(ev):
                row.update(summarize_marker(ev))
            else:
                reason = ("no exclusive pair (cannot earn REDSEA under option 3)"
                          if not partners1 and not partners2 else
                          "pairs declared but no donor measured them yet")
                row.update(proposal="raw" if not partners1 else "pending",
                           reason=reason, flags="", n_donors_evaluable=0)
                if not partners1 and partners2:
                    row["reason"] = "only tier-2 pairs; no donor evaluated yet"
        rows.append(row)
    cols = ["marker", "kind", "uses", "compartment", "registry_policy", "proposal", "reason",
            "flags", "n_donors_acquired", "n_donors_evaluable", "evidence_tier",
            "median_frac_zeroed_by_redsea", "median_ratio_redsea_raw",
            "median_copos_drop", "median_retention", "median_retention_bright",
            "median_selectivity", "median_copos_drop_k", "median_cosine_raw",
            "median_cosine_corrected", "borderline_partners", "typing_rule_unstable_donors",
            "proposal_q20", "proposal_q30",
            "tier1_partners", "tier1_evidence", "tier2_partners", "tier2_evidence",
            "in_tier1_pair", "in_tier2_pair"]
    s = pd.DataFrame(rows)
    return s[[c for c in cols if c in s.columns]]


# ----------------------------------------------------------------------------- README
def pair_label(key: str) -> str:
    """Markdown-safe pair name: ``CD4|CD8@CD3e+`` -> ``CD4 / CD8 (in CD3e+)``."""
    pair, _, gate = key.partition("@")
    return pair.replace("|", " / ") + (f" (in {gate})" if gate else "")


def _s(x) -> str:
    return x if isinstance(x, str) else ""


def _fmt(x, nd=2):
    return "" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


HEADLINE_PAIRS = ("INS|GCG", "INS|SST", "GCG|SST", "CD3e|CD20", "CD3e|CD79a", "CD3e|CD68",
                  "CD4|CD8@CD3e+", "CD3e|E_cadherin", "CD31|E_cadherin")


def headline(summary: pd.DataFrame, focal: pd.DataFrame, pairs: pd.DataFrame | None) -> list[str]:
    """Data-derived headline numbers (regenerated on every run)."""
    lines = ["## Headline (computed from the donors above)", ""]
    ev = focal[(focal["tier"] == 1) & (focal["verdict"] != "insufficient")] if len(focal) \
        else focal
    if len(ev):
        q = lambda s: f"{s.median():.2f} (IQR {s.quantile(.25):.2f}-{s.quantile(.75):.2f})"  # noqa: E731
        lines += [
            f"- Over {len(ev)} tier-1 (donor, marker, partner) rows, correcting the marker alone "
            f"changes the impossible co-positive rate by a median {q(ev['copos_drop'])} "
            f"(positive = fewer), keeps {q(ev['retention'])} of its raw single positives, with "
            f"selectivity {q(ev['selectivity'].clip(upper=999))} and a matched-count drop of "
            f"{q(ev['copos_drop_k'])}.",
            "- Verdict counts: " + ", ".join(f"{k} {v}" for k, v in
                                           ev["verdict"].value_counts().items()) + ".",
        ]
        if "redsea_ratio_in_copos" in ev:
            rc, zc = ev["redsea_ratio_in_copos"].median(), ev["redsea_zeroed_in_copos"].median()
            lines += [
                "- Inside the impossible co-positive cells, REDSEA/raw of the focal value is "
                f"{q(ev['redsea_ratio_in_copos'])}, clamped to 0 in "
                f"{q(ev['redsea_zeroed_in_copos'])} of them; in its single positives "
                f"{q(ev['redsea_ratio_in_single'])}."
                + (" REDSEA barely treats these cells as spillover recipients, so most "
                   "co-positives come from something its 1-px lateral model does not see "
                   "(segmentation, overlap in z, or a permissive positivity bar)."
                   if 0.9 <= rc <= 1.1 and zc < 0.05 else "")]
    if pairs is not None and len(pairs) and "copos_union_raw_raw" in pairs:
        ok = pairs[pairs["status"] == "ok"]
        lines += ["", "Raw co-positive rate among (A+ | B+), and after correcting both markers "
                  "(median over donors):", "", "| pair | raw | both REDSEA | joint drop |",
                  "|---|---:|---:|---:|"]
        for key in HEADLINE_PAIRS:
            g = ok[ok["pair"] == key]
            if len(g):
                a, b = g["copos_union_raw_raw"].median(), g["copos_union_redsea_redsea"].median()
                jd = (1 - g["copos_union_redsea_redsea"] / g["copos_union_raw_raw"]).median()
                lines.append(f"| {pair_label(key)} | {a:.3f} | {b:.3f} | {jd:+.2f} |")
    z = summary.loc[summary["registry_policy"] == "redsea_corrected",
                    "median_frac_zeroed_by_redsea"].dropna()
    if len(z):
        top = summary.loc[z.sort_values(ascending=False).index[:6]]
        lines += ["", "- REDSEA zeroes a median " + ", ".join(
            f"{r.marker} {r.median_frac_zeroed_by_redsea * 100:.1f}%"
            for r in top.itertuples(index=False)) + " of raw>0 QC cells (highest six)."]
        rr = summary.loc[summary["registry_policy"] == "redsea_corrected",
                         "median_ratio_redsea_raw"].dropna()
        if len(rr):
            lines += [f"- Elsewhere REDSEA moves values little: the per-cell REDSEA/raw ratio "
                      f"has a median of {rr.median():.3f} across corrected markers (range of "
                      f"marker medians {rr.min():.3f}-{rr.max():.3f})."]
    return lines + [""]


def write_readme(out: Path, summary: pd.DataFrame, donors: list[dict], args, focal: pd.DataFrame,
                 skipped: dict, pairs_df: pd.DataFrame | None = None) -> None:
    corrected = summary[summary["registry_policy"] == "redsea_corrected"]
    lines = [
        "# REDSEA spillover-correction policy audit (option 3)",
        "",
        f"Generated {_dt.datetime.now().isoformat(timespec='minutes')} by "
        "`scripts/v3/redsea_policy_audit.py` (" + SCRIPT_VERSION + "). Evidence only: the "
        "registry is untouched, and the policy change is a separate PR after review.",
        "",
        "**Rule under test** (user-approved \"option 3\"): keep REDSEA correction for a marker "
        "only where it reduces biologically impossible co-positives without destroying that "
        "marker's real single-positive signal. Same-slide Xenium RNA agreement (donors 6533 and "
        "HDL112) is a secondary check run separately and is not used here.",
        "",
        "## Donors included",
        "",
        "| donor | status | source | QC cells |",
        "|---|---|---|---:|",
    ]
    for d in donors:
        lines.append(f"| {d['donor_id']} | {d.get('status', '')} | "
                     f"`{Path(d['data_dir']).name}` | {d.get('n_cells_qc', 0):,} |")
    have = {d["donor_id"]: Path(d["data_dir"]).name for d in donors}
    unknown = [d["donor_id"] for d in donors if not d.get("status")]
    if unknown:
        lines += ["", "Status is context only, from `donor_acquisition_batch.csv`, which does "
                  "not list " + ", ".join(unknown) + "."]
    if skipped:
        lines += ["", "Not read on the last run: " + "; ".join(
            f"{k} ({v}" + (f"; included from `{have[k]}`" if k in have else "") + ")"
            for k, v in sorted(skipped.items())) + "."]
    lines += [
        "",
    ]
    lines += headline(summary, focal, pairs_df)
    lines += [
        "## Proposal (REDSEA-corrected markers)",
        "",
        "| marker | proposal | reason | donors | zeroed | copos drop | retention | bright ret. "
        "| selectivity | drop @k | at bar 0.2 / 0.3 | evidence | flags |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    order = {"keep REDSEA": 0, "judgement": 1, "pending": 2, "raw": 3}
    cs = corrected.assign(_o=corrected["proposal"].map(order).fillna(9),
                          _t=~corrected["in_tier1_pair"]).sort_values(["_t", "_o", "marker"])
    paired = cs[cs["in_tier1_pair"] | cs["in_tier2_pair"]]
    for r in paired.itertuples(index=False):
        n_ev = getattr(r, "n_donors_evaluable", 0)
        n_ev = int(n_ev) if isinstance(n_ev, (int, float)) and np.isfinite(n_ev) else 0
        z = r.median_frac_zeroed_by_redsea
        ev1, ev2 = _s(getattr(r, "tier1_evidence", "")), _s(getattr(r, "tier2_evidence", ""))
        evidence = ev1 if ev1 else ("tier 2 only: " + ev2 if ev2 else "")
        lines.append(
            f"| {r.marker} | {r.proposal} | {r.reason} | {n_ev} | "
            f"{_fmt(z * 100, 1) + '%' if np.isfinite(z) else ''} | "
            f"{_fmt(getattr(r, 'median_copos_drop', np.nan))} | "
            f"{_fmt(getattr(r, 'median_retention', np.nan))} | "
            f"{_fmt(getattr(r, 'median_retention_bright', np.nan))} | "
            f"{_fmt(getattr(r, 'median_selectivity', np.nan), 1)} | "
            f"{_fmt(getattr(r, 'median_copos_drop_k', np.nan))} | "
            f"{_s(getattr(r, 'proposal_q20', ''))} / {_s(getattr(r, 'proposal_q30', ''))} | "
            f"{evidence} | {_s(r.flags)} |")
    nopair = cs[~(cs["in_tier1_pair"] | cs["in_tier2_pair"])]
    lines += ["", "**Raw by rule (no exclusive pair, so they cannot earn REDSEA under option 3)**, "
              "with the median share of raw>0 cells REDSEA zeroes: "
              + ", ".join(f"{r.marker} ({_fmt(r.median_frac_zeroed_by_redsea * 100, 1)}%)"
                          if np.isfinite(r.median_frac_zeroed_by_redsea) else r.marker
                          for r in nopair.itertuples(index=False)) + "."]
    judg = cs[cs["proposal"] == "judgement"]["marker"].tolist()
    dep = cs[(cs["proposal"] != "judgement")
             & cs["flags"].fillna("").astype(str).str.contains("donor-dependent")]
    dep = [f"{r.marker} ({r.proposal}; {r.reason})" for r in dep.itertuples(index=False)]
    border = cs[(cs["proposal"] == "raw")
                & (cs["borderline_partners"].fillna("").astype(str) != "")]["marker"].tolist() \
        if "borderline_partners" in cs else []
    lines += [
        "",
        "Markers needing the user's judgement: " + (", ".join(judg) if judg else "none") + ".",
        "",
        "Donor-dependent evidence behind a raw proposal (flagged): "
        + ("; ".join(dep) if dep else "none") + ".",
        "",
        "Borderline (raw under the 20% bar, but REDSEA gives a selective, loss-free 10-20% "
        "drop on at least one partner in most donors; a 10% bar would keep REDSEA): "
        + (", ".join(border) if border else "none") + ".",
        "",
        "Unchanged: " + ", ".join(summary.loc[summary["registry_policy"] != "redsea_corrected",
                                              "marker"]) + " (REDSEA does not compensate them).",
        "",
    ]
    lines += [
        "## Where the values live, and what a policy change needs",
        "",
        "- **Raw**: `data/cells/donor_id=<id>/data_0.parquet`, QuPath means. "
        "`Cell: <m>: Mean` is the bare sanitized name (`Vimentin`, `Histone_H3_p_Ser28`); "
        "`Nucleus: <m>: Mean` is `Nucleus__<m>`.",
        "- **REDSEA**: `data/redsea/donor_id=<id>/data_0.parquet`, "
        "`{Nucleus,Cytoplasm,Membrane,Cell}__<registry name>` (means, clamped at 0; a stored 0 "
        "is a clamp). The raw re-rasterised mean is not stored for corrected markers.",
        "- **Expression** (`phenocycler/expression.py::_resolve_measurement`): per registry "
        "marker, `measurement_source` redsea -> REDSEA `<compartment>__<name>`; qupath_raw -> "
        "cells `<compartment>__<name>`, then `Raw__...`, then the bare name for Cell. "
        "`measurement_source` defaults to qupath_raw only for `uncompensated`, else redsea. "
        "Typing reads this table. REDSEA itself compensates exactly the `redsea_corrected` "
        "markers (`marker_taxonomy.spillover_corrected_mask`).",
        "- **Switching a Cell marker to raw needs only an expression rerun** if the registry "
        "sets `spillover_policy: uncompensated` (source qupath_raw, reads the cells table), "
        "**but** `validate_registry` currently rejects `uncompensated` unless the compartment "
        "is Nucleus, so that PR must relax that check. `panel.NUCLEUS_MEASURED` "
        "in the typing lane is a literal, so it does not follow the policy. The alternative, "
        "`redsea_passthrough`, "
        "keeps the source on REDSEA and needs a REDSEA rerun (qptiff re-read; intermediates are "
        "not saved). `Histone_H3_(p_Ser28)` cannot use the qupath_raw route as is: expression "
        "looks for `Nucleus__Histone_H3_(p_Ser28)` but the cells table spells it "
        "`Nucleus__Histone_H3_p_Ser28`.",
        "- The pipeline skips a stage whose output exists, and `write_donor_partition` refuses "
        "to overwrite, so move `data/expression` (and its availability audit) aside, then run "
        "`python -m phenocycler.pipeline run --only expression --no-export`.",
        "",
        "## Exclusive pairs",
        "",
        "Tier 1 decides the proposal; tier 2 has a documented way to co-occur and is context "
        "only. CD4/CD8 is scored inside cells that are CD3e+ in both variants.",
        "",
        "| pair | tier | why exclusive | caveat |",
        "|---|---|---|---|",
    ]
    for p in PAIRS:
        lines.append(f"| {pair_label(p.key)} | {p.tier} | {p.why} | {p.caveat} |")
    lines += ["", "Considered and rejected:", ""]
    lines += [f"- {a}: {b}." for a, b in REJECTED_PAIRS]
    lines += [
        "- Ubiquitous or multi-lineage markers (metabolic enzymes, Beta_actin, HLA_A, CD99, "
        "CD44, Bcl_2, b_Catenin1) and activation/state markers (CD38, CD39, CD107a, "
        "Granzyme_B, HLA_DR, ICOS, IDO1, iNOS, LAG3, PD_1, PD_L1, VISTA, Histone H3 pSer28) have "
        "no conservative exclusive partner. CD11b is assay QC only; IAPP failed QC.",
        "",
        "## Method",
        "",
        "- Cells: `geometry_qc.analysis_eligible`, finite in both variants for the pair.",
        "- **Positivity, same rule on both variants, per donor.** Default is the typing "
        "pipeline's own rule (`run_step3.build_scale` + `score.meets_requirements`): positive "
        f"iff (v - floor)/(p{SCALE_Q} - floor) >= {REQUIRE_MIN}, floor = log1p mode of v>0 "
        "(512 bins). Typing fits it cohort-wide on UniFORM-aligned values; UniFORM is one "
        "scale per (donor, marker), so the per-donor fit is the same rule on each donor's own "
        "scale. The rule is used only where its threshold (and, when material, its floor) "
        f"moves < {RULE_STABILITY:g}x between raw and REDSEA; otherwise the log1p mode has "
        "changed population (REDSEA's zero clamp removes background, so the mode lands on the "
        "positive peak: E-cadherin and Vimentin on 6533) and Otsu on log1p(v>0) is used for "
        "that (donor, marker), with the same stability check. `thresholds_by_donor.tsv` "
        "records the rule, floor, top and threshold of every (donor, marker, variant).",
        "- **Attribution.** Each marker is scored as the *focal* marker with its partner held "
        "raw, so its own correction is isolated; both-corrected rates are reported too.",
        "- `copos_drop` = 1 - copos(focal REDSEA, partner raw) / copos(raw, raw), with copos "
        "= |A+ & B+| / |A+ | B+|.",
        "- `retention` = share of raw focal+ partner- cells still focal+ after REDSEA; "
        "`retention_bright` = the same for those above the median raw value of raw focal+ cells.",
        "- `cleared` = share of raw co-positive cells that lose the focal call; `selectivity` = "
        "cleared / (1 - retention). A uniform dimming of the channel gives ~1.",
        "- Sensitivity (added by this audit, flagged): `*_k` repeats the scoring with REDSEA "
        "positives = the top-k REDSEA cells, k = raw positive count, i.e. at matched positive "
        "count, so a benefit that is only a threshold shift is exposed.",
        "- Sensitivity (added by this audit, flagged): the whole decision is repeated with the "
        "typing bar at " + " and ".join(f"{q:g}" for q in SWEEP) + " instead of "
        f"{REQUIRE_MIN:g} (typing sweeps this bar rather than cutting per marker); Otsu-scored "
        "(donor, marker) cells keep their threshold. A proposal that changes is flagged.",
        "- Threshold-light: PR #18's SVD ratio (sigma2/sigma1 of the uncentred scaled pair) "
        "and the uncentred cosine (0 = perfectly exclusive, balance-free), each on the rule's "
        "own floor/p99.9 scale.",
        "",
        "## Decision thresholds",
        "",
        f"Per (donor, focal marker, partner): `insufficient` if either marker has < "
        f"{MIN_POSITIVE} raw positives or no usable rule; `nothing_to_fix` if < {MIN_COPOS} raw "
        f"co-positives; `harm` if retention < {RETENTION_HARM:.2f}; `worse` if copos_drop <= "
        f"-{WORSE_RISE:.2f}; `benefit` if copos_drop >= {BENEFIT_DROP:.2f} and selectivity >= "
        f"{SELECTIVITY_MIN:g} and retention >= {RETENTION_OK:.2f} (`tradeoff` if retention is "
        f"{RETENTION_HARM:.2f}-{RETENTION_OK:.2f}); `small_benefit` if the drop is "
        f"{SMALL_DROP:.2f}-{BENEFIT_DROP:.2f} with the same selectivity and retention >= "
        f"{RETENTION_OK:.2f} (counts as no gain, flagged `borderline`); otherwise `no_benefit`.",
        "",
        f"- {BENEFIT_DROP:.0%} net drop: a material share of impossible calls. REDSEA's own-rim "
        "reinforcement alone raises the p99.9, and with it the typing threshold, by ~5-20% "
        "(`thresholds_by_donor.tsv`), which trims co-positives and single positives alike; a "
        "smaller drop is not distinguishable from that.",
        f"- Selectivity >= {SELECTIVITY_MIN:g}: REDSEA must remove the marker from impossible "
        "cells at least twice as often as from plausible ones; otherwise it is dimming the "
        "channel, not removing spillover.",
        f"- Retention >= {RETENTION_OK:.2f} for a clean keep, < {RETENTION_HARM:.2f} = destroys "
        "signal: some single positives are themselves spillover onto non-partner neighbours, "
        "so 1.0 is not expected, but losing > 10% of a marker's calls is material for typing "
        "and > 20% is destructive.",
        f"- Across donors, a partner's verdict is the group shared by >= {MAJORITY:.0%} of "
        f"evaluable donors (benefit / tradeoff / adverse = harm+worse / no_gain = no_benefit+"
        f"nothing_to_fix), else `donor-dependent`. Fewer than {MIN_DONORS_FINAL} evaluable "
        "donors marks the verdict provisional.",
        "- Marker: **keep REDSEA** if >= 1 tier-1 partner is `benefit` and none is adverse, "
        "tradeoff or donor-dependent, and the benefit survives the matched-count check; "
        "**raw** if no partner is `benefit` (or the marker has no tier-1 pair: under option 3 "
        "it cannot earn REDSEA); **judgement** otherwise.",
        "",
    ]
    lines += [
        "## Files",
        "",
        "- `pairs_by_donor.tsv`: one row per donor x pair, every raw/REDSEA combination "
        "(`copos_union_<A variant>_<B variant>`), per-marker retention/cleared/selectivity "
        "(prefix `a_`/`b_`), matched-count (`_k`), SVD ratio and cosine, verdicts.",
        "- `marker_evidence_by_donor.tsv`: the same, one row per donor x focal marker x partner.",
        "- `thresholds_by_donor.tsv`: rule, floor, top, threshold and positives per donor x "
        "marker x variant, plus REDSEA's clamp rate.",
        "- `markers_summary.tsv`: every registry marker, its pairs, evidence and proposal.",
        "- `per_donor/donor_id=<id>/`: the per-donor cache the tables are rebuilt from.",
        "",
        "## Rerun",
        "",
        "```bash",
        "cd /scratch/io60/work/senior-pr17/core/phenocycler",
        "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python \\",
        "    /home/smith6jt/IO60panc2nd/scripts/v3/redsea_policy_audit.py",
        "```",
        "",
        "It reads only donors whose REDSEA partition is present and logged `[donor-done] "
        "redsea/<id>`, skips donors already done from the same data directory (a donor done "
        "from `pilot_data` is recomputed once the full run has it), and rebuilds every table. "
        "`--force` recomputes all; `--summarize-only` rebuilds tables from the caches.",
    ]
    (out / "README.md").write_text("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------- main
def _atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(f".{path.name}.tmp")
    frame.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)


def rebuild_tables(out: Path, reg: pd.DataFrame, args, skipped: dict,
                   status: dict[str, str]) -> pd.DataFrame:
    pairs = _read_all(out, "pairs.tsv")
    focal = _read_all(out, "focal.tsv")
    th = _read_all(out, "thresholds.tsv")
    clamps = _read_all(out, "clamp_rates.tsv")
    for t in (pairs, focal, th):
        if len(t):      # disease status is context only (donor_acquisition_batch.csv)
            t.insert(1, "donor_status", t["donor_id"].map(lambda d: status.get(_norm(d), "")))
    _atomic_tsv(pairs, out / "pairs_by_donor.tsv")
    _atomic_tsv(focal, out / "marker_evidence_by_donor.tsv")
    _atomic_tsv(th, out / "thresholds_by_donor.tsv")
    summary = build_summary(reg, focal, clamps)
    _atomic_tsv(summary, out / "markers_summary.tsv")
    donors = []
    for p in sorted((out / "per_donor").glob("donor_id=*/meta.json")):
        m = json.loads(p.read_text())
        m["status"] = status.get(_norm(m["donor_id"]), "")
        donors.append(m)
    write_readme(out, summary, donors, args, focal, skipped, pairs)
    return summary


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--log", type=Path, default=None,
                    help=f"run log gating completion (default {DEFAULT_LOG} for the default "
                         "data dir; required otherwise unless --no-log-gate)")
    ap.add_argument("--no-log-gate", action="store_true",
                    help="trust complete partitions without a log line (tests only)")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--senior-pkg", type=Path, default=DEFAULT_SENIOR_PKG)
    ap.add_argument("--registry", type=Path, default=None,
                    help="default: <senior-pkg>/phenocycler/marker_registry.json")
    ap.add_argument("--status-csv", type=Path, default=DEFAULT_STATUS_CSV)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--summarize-only", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    sys.path.insert(0, str(args.senior_pkg))
    from phenocycler.cohort import filter_eligible_donors   # the one exclusion authority

    registry = args.registry or args.senior_pkg / "phenocycler" / "marker_registry.json"
    reg = load_registry_rows(registry)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    status = donor_status_map(args.status_csv)
    data_dir = args.data_dir.resolve()
    log_path = args.log
    if log_path is None and not args.no_log_gate:
        if data_dir == DEFAULT_DATA.resolve():
            log_path = DEFAULT_LOG
        else:
            raise SystemExit("--log is required for a non-default --data-dir "
                             "(or --no-log-gate for synthetic tests)")
    skipped: dict[str, str] = {}
    if not args.summarize_only:
        ready, skipped = ready_donors(data_dir, None if args.no_log_gate else log_path,
                                      args.donors, filter_eligible_donors)
        todo = []
        for d in ready:
            src = cached_source(out, d)
            if src == str(data_dir) and not args.force:
                log(f"[skip] {d}: already done from {data_dir.name}")
            else:
                todo.append(d)
        for d, why in sorted(skipped.items()):
            log(f"[wait] {d}: {why}")
        for d in todo:
            t = time.time()
            log(f"[run] {d} from {data_dir}")
            res = process_donor(data_dir, d, reg)
            clamps = clamp_rates(data_dir, d, reg)
            write_donor_cache(out, d, res, clamps)
            log(f"[done] {d}: {res['meta']['n_cells_qc']:,} QC cells, "
                f"{(res['pairs']['status'] == 'ok').sum()} pairs scored, {time.time() - t:.0f}s")
    summary = rebuild_tables(out, reg, args, skipped, status)
    counts = summary["proposal"].value_counts().to_dict()
    log(f"[summary] {out}/markers_summary.tsv: {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
