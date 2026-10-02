"""uniform_bright_end.py: a gain moves both ends of a donor's distribution, an offset only the background."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "senior"))
import uniform_bright_end as B  # noqa: E402

RNG = np.random.default_rng(7)


def donor(background: float, positives: float, n: int = 60000, n_pos: int = 600) -> np.ndarray:
    """A background peak at `background` and a few hundred bright cells around `positives`."""
    return np.concatenate([background * np.exp(RNG.normal(0, 0.25, n)), positives * np.exp(RNG.normal(0, 0.15, n_pos))])


def per_donor(kind: str) -> pd.DataFrame:
    backgrounds = {"a": 40.0, "b": 80.0, "c": 160.0, "d": 20.0, "0112": 10.0, "6347": 12.0, "6533": 9.0, "6580": 11.0}
    vals, scale = {}, {}
    for d, bg in backgrounds.items():
        # gain: the bright cells sit 50x over the donor's own background; offset: at 4,000 whatever the background
        vals[d] = donor(bg, 50.0 * bg if kind == "gain" else 4000.0)
        scale[d] = 80.0 / bg                      # the factor that lines the backgrounds up on donor b
    return pd.DataFrame(B.donor_rows("M", vals, scale, "automatic_registration"))


def test_the_factor_lines_the_backgrounds_up_in_both_cases():
    for kind in ("gain", "offset"):
        per = per_donor(kind)
        assert np.allclose(per["mode_aligned"], 80.0, rtol=0.08), kind
        assert set(per.loc[per["four"], "donor_id"]) == set(B.FOUR)


def test_a_gain_is_corrected_at_the_bright_end_too():
    s = B.spread_table(per_donor("gain"), min_donors=4).iloc[0]
    assert s["slope_top_on_background"] > 0.9
    assert s["top_sd_after_over_before"] < 0.2
    g = B.group_table(per_donor("gain")).iloc[0]
    assert 0.8 < g["top_four_over_others_after"] < 1.25
    assert not g["four_tops_all_above_others"]


def test_an_offset_is_not_and_the_factor_moves_the_bright_ends_apart():
    per = per_donor("offset")
    s = B.spread_table(per, min_donors=4).iloc[0]
    assert abs(s["slope_top_on_background"]) < 0.1
    assert s["top_sd_after_over_before"] > 5
    g = B.group_table(per).iloc[0]
    assert 0.8 < g["top_four_over_others_before"] < 1.25       # equally bright in the stored counts
    assert g["top_four_over_others_after"] > 3                 # and several times apart once the backgrounds are lined up
    assert g["four_tops_all_above_others"] and not g["four_tops_before_all_above"]


def test_the_level_in_stored_counts_follows_the_factor():
    per = per_donor("offset")
    lv = B.level_table(per, "M", floor=80.0, span=3920.0, level=0.10).set_index("donor_id")
    # floor + 0.10 x span = 472 aligned; donor b is the reference (factor 1), donor c has factor 0.5, 0112 has factor 8
    assert np.isclose(lv.loc["b", "count_at_level"], 472.0, rtol=1e-6)
    assert np.isclose(lv.loc["c", "count_at_level"], 944.0, rtol=1e-6)
    assert np.isclose(lv.loc["0112", "count_at_level"], 59.0, rtol=1e-6)
    assert np.allclose(lv["unshifted_over_background"] + lv["mode"], 472.0)
