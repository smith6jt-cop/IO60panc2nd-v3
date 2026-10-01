"""Tests for scripts/v3/eval_6457_cycles.py (synthetic images only; never touches the share).

Run:  <Senior>/.venv/bin/python -m pytest -q scripts/v3/tests/test_eval_6457_cycles.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import eval_6457_cycles as ev  # noqa: E402

RNG = np.random.default_rng(7)


def _nuclei(shape=(1024, 1024), sigma=3):
    img = gaussian_filter(RNG.random(shape), sigma)
    return ((img - img.min()) / (img.max() - img.min()) * 4000 + 100).astype(np.uint16)


@pytest.mark.parametrize("sigma", [1.5, 3, 6])
def test_measure_pair_reports_the_shift_that_aligns_the_moving_tile_to_its_reference(sigma):
    # sigma 6 is a hazy tile: the phase-normalised estimate without a window returned 0 there
    ref = _nuclei(sigma=sigma)
    mov = np.roll(ref, (3, -2), axis=(0, 1)).astype(np.float32) + RNG.normal(0, 30, ref.shape).astype(np.float32)
    tissue = np.ones((2, 2), dtype=bool)
    d = ev.measure_pair(ref, mov, tissue, 512)
    assert len(d) == 4
    assert d["dy_um"].to_numpy() == pytest.approx(-3 * ev.PX_UM, abs=0.06)   # content sits 3 px low: move it up
    assert d["dx_um"].to_numpy() == pytest.approx(2 * ev.PX_UM, abs=0.06)
    assert d["shift_um"].to_numpy() == pytest.approx(np.hypot(3, 2) * ev.PX_UM, abs=0.06)
    same = ev.measure_pair(ref, ref, tissue, 512)
    assert same["ncc"].to_numpy() == pytest.approx(1.0)
    assert same["shift_um"].to_numpy() == pytest.approx(0.0, abs=1e-6)


def test_measure_pair_skips_tiles_without_tissue():
    ref = _nuclei()
    tissue = np.array([[True, False], [False, False]])
    assert len(ev.measure_pair(ref, ref, tissue, 512)) == 1


def test_relative_sharpness_falls_when_the_same_tile_is_blurred_and_ignores_brightness():
    t = _nuclei((512, 512))
    assert ev.rel_sharpness(gaussian_filter(t.astype(np.float32), 2)) < 0.5 * ev.rel_sharpness(t)
    assert ev.rel_sharpness(t.astype(np.float32) * 3) == pytest.approx(ev.rel_sharpness(t), rel=1e-4)


def _field(tile, tx, ty, a1, a2, b1, b2, noise=0.0, n=40):
    rows, cols = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    x = (cols.ravel() + 0.5) * tile * ev.PX_UM
    y = (rows.ravel() + 0.5) * tile * ev.PX_UM
    x0, y0 = x.mean(), y.mean()
    dx = tx + a1 * (x - x0) + a2 * (y - y0) + RNG.normal(0, noise, x.size)
    dy = ty + b1 * (x - x0) + b2 * (y - y0) + RNG.normal(0, noise, x.size)
    return pd.DataFrame({"row": rows.ravel(), "col": cols.ravel(), "dx_um": dx, "dy_um": dy})


def test_affine_fit_recovers_translation_rotation_and_scale():
    theta = np.radians(0.010)                           # 10 millidegrees
    d = _field(512, 1.5, -0.5, 200e-6, -theta, theta, 200e-6)
    f = ev.affine_fit(d, 512)
    assert f["tx_um"] == pytest.approx(1.5, abs=1e-6) and f["ty_um"] == pytest.approx(-0.5, abs=1e-6)
    assert f["rotation_mdeg"] == pytest.approx(10.0, rel=1e-4)
    assert f["scale_ppm"] == pytest.approx(200.0, rel=1e-4)
    assert f["stretch_ppm"] == pytest.approx(0.0, abs=1e-3)
    assert f["residual_um_p95"] == pytest.approx(0.0, abs=1e-6)


def test_affine_fit_leaves_a_non_affine_field_in_the_residual_and_survives_bad_tiles():
    d = _field(512, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, noise=0.05)
    d.loc[::50, "dx_um"] = 40.0                         # 2 % of tiles: a failed phase correlation
    f = ev.affine_fit(d, 512)
    assert abs(f["tx_um"]) < 0.02                       # the 40 um tiles did not drag the fit
    assert 0.02 <= f["tiles_outlier_share"] < 0.05      # they, and a little of the noise tail, were set aside
    bend = _field(512, 0, 0, 0, 0, 0, 0)
    x = (bend["col"] + 0.5) * 512 * ev.PX_UM
    bend["dx_um"] = 3.0 * np.sin(2 * np.pi * x / x.max())   # a wave: no affine model fits it
    assert ev.affine_fit(bend, 512)["residual_um_p50"] > 0.5


def test_affine_fit_with_too_few_tiles_fits_nothing():
    d = _field(512, 1, 1, 0, 0, 0, 0, n=5)
    assert ev.affine_fit(d.iloc[:10], 512) == {"tiles_fit": 10}
