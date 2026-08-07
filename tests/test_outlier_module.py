"""
tests/test_outlier_module.py
----------------------------
Unit tests for wavy/outlier_module.py.

All tests use purely synthetic data; no real model files or network access
is required.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for CI

from wavy.outlier_module import (
    find_outliers,
    _detect_outliers,
    outlier_class,
    PATTERN_REGISTRY,
    register_pattern,
    _detect_checkerboard,
    _detect_garden_sprinkler,
    _detect_source_term_ringing,
    _detect_hs_collapse,
    _detect_spinup_insufficient,
    required_spinup_hours,
)

# ---------------------------------------------------------------------------
# Shared fixture — a minimal cco-like object backed by a synthetic dataset
# ---------------------------------------------------------------------------


class _MockCCO:
    """
    Minimal stand-in for a collocation_class instance.
    Exposes only the attributes that outlier_module actually reads:
        .vars   – xr.Dataset
        .model  – str
        .nID    – str
        .varalias – str
        .sd, .ed, .region, .name – str / None
    """

    def __init__(self, ds: xr.Dataset):
        self.vars = ds
        self.model = "mock_model"
        self.nID = "mock_nID"
        self.varalias = "Hs"
        self.sd = "2023-01-01"
        self.ed = "2023-01-02"
        self.region = None
        self.name = None


def _make_times(n: int, start: str = "2023-01-01", freq: str = "1min"):
    """Return a DatetimeIndex of length n."""
    return pd.date_range(start=start, periods=n, freq=freq)


def _build_cco(obs_hs, model_hs, times=None, extra_vars: dict | None = None):
    """
    Build a _MockCCO with obs_Hs, model_Hs, obs_lons, obs_lats, time.
    obs_hs and model_hs must be array-like of equal length.
    """
    n = len(obs_hs)
    if times is None:
        times = _make_times(n)
    ds = xr.Dataset(
        {
            "obs_Hs": ("time", np.asarray(obs_hs, dtype=float)),
            "model_Hs": ("time", np.asarray(model_hs, dtype=float)),
            "obs_lons": ("time", np.linspace(0.0, 10.0, n)),
            "obs_lats": ("time", np.linspace(50.0, 60.0, n)),
        },
        coords={"time": times},
    )
    if extra_vars:
        for k, v in extra_vars.items():
            ds[k] = ("time", np.asarray(v))
    return _MockCCO(ds)


@pytest.fixture
def mock_cco():
    """
    Synthetic cco with 100 normal points and 5 obvious positive outliers
    (bias ≈ +15 m each, injected at indices 10, 25, 50, 75, 90).
    The background bias is ~0 m with std ~0.1 m, so a z-score threshold
    of 3 reliably flags exactly these 5 points.
    """
    rng = np.random.default_rng(seed=42)
    n = 100
    obs = rng.normal(loc=2.0, scale=0.5, size=n)
    mod = obs + rng.normal(loc=0.0, scale=0.1, size=n)  # bias ≈ 0

    outlier_indices = [10, 25, 50, 75, 90]
    for i in outlier_indices:
        mod[i] = obs[i] + 15.0  # large positive bias → outlier

    return _build_cco(obs, mod)


# ---------------------------------------------------------------------------
# 1. find_outliers – z-score
# ---------------------------------------------------------------------------


def test_find_outliers_zscore(mock_cco):
    """find_outliers with method='zscore' must flag exactly the 5 injected points."""
    df, stats = find_outliers(mock_cco, n_std=3, method="zscore")

    assert df is not None, "Expected outliers to be found"
    assert stats["method"] == "zscore"
    assert stats["n_outliers"] == 5

    # The flagged times must correspond to indices 10, 25, 50, 75, 90
    expected_times = mock_cco.vars["time"].values[[10, 25, 50, 75, 90]]
    found_times = pd.DatetimeIndex(df["time"]).values
    assert set(expected_times) == set(found_times)


# ---------------------------------------------------------------------------
# 2. find_outliers – mad
# ---------------------------------------------------------------------------


def test_find_outliers_mad(mock_cco):
    """find_outliers with method='mad' must flag the same 5 injected points."""
    df, stats = find_outliers(mock_cco, n_std=3, method="mad")

    assert df is not None
    assert stats["method"] == "mad"
    # MAD-based z-score of bias=+15 m on a near-zero distribution is enormous
    assert stats["n_outliers"] == 5

    expected_times = mock_cco.vars["time"].values[[10, 25, 50, 75, 90]]
    found_times = pd.DatetimeIndex(df["time"]).values
    assert set(expected_times) == set(found_times)


# ---------------------------------------------------------------------------
# 3. outlier_class.populate()
# ---------------------------------------------------------------------------


def test_outlier_class_populate(mock_cco):
    """populate() must produce an xr.Dataset with n_outliers == 5."""
    outo = outlier_class(mock_cco, n_std=3, method="zscore").populate()

    assert isinstance(
        outo.vars, xr.Dataset
    ), "outo.vars must be an xr.Dataset after populate()"
    assert outo.stats["n_outliers"] == 5
    assert "obs_Hs" in outo.vars
    assert "model_Hs" in outo.vars
    assert "bias" in outo.vars
    assert len(outo.vars["time"]) == 5


# ---------------------------------------------------------------------------
# 4. detect_numerical_patterns – checkerboard
# ---------------------------------------------------------------------------


def test_detect_checkerboard():
    """
    Build a synthetic outo where colidx_y alternates between even (LOW Hs)
    and odd (HIGH Hs).  The checkerboard criterion must be detected.
    """
    n = 20
    # Alternating grid indices: 800, 801, 802, 803, …  (even, odd, even, odd …)
    colidx_y = np.arange(800, 800 + n)  # 800,801,802,...
    # model_Hs alternates: even cells → HIGH, odd cells → LOW
    model_hs = np.where(colidx_y % 2 == 0, 13.0, 1.5)  # even=HIGH, odd=LOW
    obs_hs = np.full(n, 5.0)  # constant obs

    times = _make_times(n)
    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", model_hs - obs_hs),
            "colidx_y": ("time", colidx_y),
        },
        coords={"time": times},
    )

    # Build the outo directly without going through populate()
    cco = _build_cco(obs_hs, model_hs, times=times, extra_vars={"colidx_y": colidx_y})
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    report = outo.detect_numerical_patterns(checkerboard_threshold=0.70)

    assert report["checkerboard"] is not None
    assert (
        report["checkerboard"]["detected"] is True
    ), f"Checkerboard should be detected; got {report['checkerboard']}"
    assert (
        report["checkerboard"]["amplitude_m"] > 5.0
    ), f"Amplitude should be > 5 m; got {report['checkerboard']['amplitude_m']:.3f}"
    assert "checkerboard" in report["summary"].lower()


# ---------------------------------------------------------------------------
# 5. detect_numerical_patterns – garden-sprinkler
# ---------------------------------------------------------------------------


def test_detect_garden_sprinkler():
    """
    Build two tight time clusters 30 min apart, each with clearly bimodal
    model_Hs (4 LOW @ 1.5 m + 5 HIGH @ 13 m = 9 points).

    For a perfectly bimodal 50/50 split the bimodal score is exactly 2.0;
    using 4 LOW + 5 HIGH (odd-sized window, unequal split) yields ~2.012 so
    the strict ``> 2.0`` threshold is satisfied.
    """
    # Cluster A: 9 points, t=0..8 min
    t_A = pd.date_range("2023-01-01 00:00", periods=9, freq="1min")
    # Cluster B: 9 points, t=30..38 min (22-min gap > window_min=10)
    t_B = pd.date_range("2023-01-01 00:30", periods=9, freq="1min")
    times = t_A.append(t_B)

    n = len(times)
    # Each cluster: 4 LOW (1.5 m) + 5 HIGH (13 m)  →  bimodal score ≈ 2.012
    model_hs = np.array([1.5] * 4 + [13.0] * 5 + [1.5] * 4 + [13.0] * 5)
    obs_hs = np.full(n, 5.0)

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", model_hs - obs_hs),
        },
        coords={"time": times},
    )

    cco = _build_cco(obs_hs, model_hs, times=times)
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    report = outo.detect_numerical_patterns(
        window_min=10,
        bimodal_score_threshold=2.0,
    )

    assert report["garden_sprinkler"] is not None
    assert (
        report["garden_sprinkler"]["detected"] is True
    ), f"Garden-sprinkler should be detected; got {report['garden_sprinkler']}"
    assert len(report["garden_sprinkler"]["windows"]) >= 1
    assert "garden" in report["summary"].lower()


# ---------------------------------------------------------------------------
# 6. detect_numerical_patterns – no patterns
# ---------------------------------------------------------------------------


def test_no_patterns():
    """
    Smooth, unimodal outliers (all bias = +4 m, no grid-index alternation,
    no bimodal clusters) must return detected=False for both patterns.
    """
    n = 20
    obs_hs = np.full(n, 2.0)
    model_hs = np.full(n, 6.0)  # constant +4 m bias — no pattern
    colidx_y = np.arange(100, 100 + n)  # monotonically increasing, no alternation

    times = _make_times(n)
    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", model_hs - obs_hs),
            "colidx_y": ("time", colidx_y),
        },
        coords={"time": times},
    )

    cco = _build_cco(obs_hs, model_hs, times=times, extra_vars={"colidx_y": colidx_y})
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    report = outo.detect_numerical_patterns(
        window_min=30,
        checkerboard_threshold=0.70,
        bimodal_score_threshold=2.0,
    )

    # Checkerboard: all bias values have the same sign (+4), so no alternation
    assert (
        report["checkerboard"]["detected"] is False
    ), f"Checkerboard should NOT be detected; got {report['checkerboard']}"

    # Garden-sprinkler: all model_Hs identical → std=0 → no flagging
    assert (
        report["garden_sprinkler"]["detected"] is False
    ), f"Garden-sprinkler should NOT be detected; got {report['garden_sprinkler']}"


# ---------------------------------------------------------------------------
# 7. detect_numerical_patterns – raises when vars is None
# ---------------------------------------------------------------------------


def test_detect_patterns_raises_when_no_vars(mock_cco):
    """detect_numerical_patterns() must raise ValueError when vars is None."""
    outo = outlier_class.__new__(outlier_class)
    outo.cco = mock_cco
    outo.vars = None  # no outliers detected
    outo.stats = {}
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = None
    outo.nID = None
    outo.varalias = None
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    with pytest.raises(ValueError, match="outo.vars is None"):
        outo.detect_numerical_patterns()


# ---------------------------------------------------------------------------
# 8. pattern_report stored on self
# ---------------------------------------------------------------------------


def test_pattern_report_stored(mock_cco):
    """After detect_numerical_patterns(), self.pattern_report must be set."""
    outo = outlier_class(mock_cco, n_std=3, method="zscore").populate()
    assert outo.vars is not None

    # Inject colidx_y so checkerboard can run (monotone → not detected)
    n = len(outo.vars["time"])
    outo.vars["colidx_y"] = ("time", np.arange(100, 100 + n))

    report = outo.detect_numerical_patterns()

    assert outo.pattern_report is report
    assert "checkerboard" in report
    assert "garden_sprinkler" in report
    assert "summary" in report


# ---------------------------------------------------------------------------
# 9. Pattern registry — built-in patterns present
# ---------------------------------------------------------------------------


def test_registry_has_builtin_patterns():
    """PATTERN_REGISTRY must contain all four built-in patterns."""
    for key in (
        "checkerboard",
        "garden_sprinkler",
        "source_term_ringing",
        "hs_collapse",
    ):
        assert key in PATTERN_REGISTRY, f"'{key}' missing from PATTERN_REGISTRY"
        entry = PATTERN_REGISTRY[key]
        assert "detect_fn" in entry
        assert "solutions" in entry
        assert "references" in entry
        assert len(entry["solutions"]) > 0
        assert len(entry["references"]) > 0


# ---------------------------------------------------------------------------
# 10. Register a custom pattern and verify it runs
# ---------------------------------------------------------------------------


def test_register_custom_pattern(mock_cco):
    """A user-registered pattern must appear in detect_numerical_patterns output."""

    def _my_detector(outo, *, my_threshold=0.5):
        """Custom stub detector — never fires."""
        return {
            "detected": False,
            "contribution_pct": 0.0,
            "affected_idx": [],
            "my_metric": my_threshold,
        }

    register_pattern(
        "test_custom",
        name="Test custom pattern",
        instability_id="test",
        description="Always returns not-detected.",
        detect_fn=_my_detector,
        solutions=["No fix needed."],
        references=["Test reference."],
        default_params={"my_threshold": 0.42},
    )

    outo = outlier_class(mock_cco, n_std=3, method="zscore").populate()
    n = len(outo.vars["time"])
    outo.vars["colidx_y"] = ("time", np.arange(100, 100 + n))

    report = outo.detect_numerical_patterns()

    assert "test_custom" in report, "Custom pattern not in report"
    assert report["test_custom"]["detected"] is False
    assert report["test_custom"]["my_metric"] == 0.42

    # Clean up so other tests are not affected
    del PATTERN_REGISTRY["test_custom"]


# ---------------------------------------------------------------------------
# 11. detect_numerical_patterns – source-term ringing
# ---------------------------------------------------------------------------


def _make_ringing_outo(n=20, period=1, amplitude=4.0):
    """
    Build an outo whose model_Hs oscillates every `period` index with
    `amplitude` peak-to-peak, all within a 5-min window.
    """
    times = _make_times(n, freq="30s")  # 30 s cadence → 10 min for 20 pts
    obs_hs = np.full(n, 2.0)
    # Ringing: alternating HIGH/LOW every `period` point
    model_hs = np.where(
        np.arange(n) % (2 * period) < period, 2.0 + amplitude / 2, 2.0 - amplitude / 2
    )

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 1, n)),
            "obs_lats": ("time", np.linspace(50, 51, n)),
            "bias": ("time", model_hs - obs_hs),
        },
        coords={"time": times},
    )
    cco = _build_cco(obs_hs, model_hs, times=times)
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None
    return outo


def test_detect_source_term_ringing():
    """
    A rapidly oscillating model_Hs series (alternating every 30 s by ±2 m)
    must trigger the source-term-ringing detector.
    """
    outo = _make_ringing_outo(n=20, amplitude=4.0)
    result = _detect_source_term_ringing(
        outo,
        ringing_window_min=10,
        ringing_threshold=0.60,
        ringing_amplitude_m=1.0,
    )
    assert (
        result["detected"] is True
    ), f"Source-term ringing should be detected; got {result}"
    assert result["mean_amplitude_m"] > 1.0


# ---------------------------------------------------------------------------
# 12. detect_numerical_patterns – Hs collapse
# ---------------------------------------------------------------------------


def _make_collapse_outo(n=20, n_collapsed=5):
    """Build an outo with n_collapsed points where model_Hs ≈ 0 but obs is large."""
    times = _make_times(n)
    obs_hs = np.full(n, 3.0)
    model_hs = np.full(n, 3.0)
    # inject collapsed points
    model_hs[:n_collapsed] = 0.01

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", model_hs - obs_hs),
        },
        coords={"time": times},
    )
    cco = _build_cco(obs_hs, model_hs, times=times)
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None
    return outo, n_collapsed


def test_detect_hs_collapse():
    """Points with model_Hs ≈ 0 and obs_Hs > 0.5 m must be flagged."""
    outo, n_collapsed = _make_collapse_outo(n=20, n_collapsed=5)
    result = _detect_hs_collapse(
        outo,
        hs_collapse_max_model=0.05,
        hs_collapse_min_obs=0.5,
    )
    assert result["detected"] is True, f"Hs collapse should be detected; got {result}"
    assert result["n_collapsed"] == n_collapsed
    assert result["mean_obs_hs_m"] == pytest.approx(3.0)


def test_hs_collapse_not_detected_for_calm():
    """Near-zero model_Hs that matches calm obs must NOT be flagged."""
    n = 10
    times = _make_times(n)
    obs_hs = np.full(n, 0.1)  # genuinely calm
    model_hs = np.full(n, 0.02)

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.zeros(n)),
            "obs_lats": ("time", np.zeros(n)),
            "bias": ("time", model_hs - obs_hs),
        },
        coords={"time": times},
    )
    cco = _build_cco(obs_hs, model_hs, times=times)
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {}
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = None
    outo.nID = None
    outo.varalias = None
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    result = _detect_hs_collapse(
        outo, hs_collapse_max_model=0.05, hs_collapse_min_obs=0.5
    )
    assert result["detected"] is False


# ---------------------------------------------------------------------------
# 13. suggest_fixes — output contains solutions and references
# ---------------------------------------------------------------------------


def test_suggest_fixes_checkerboard():
    """suggest_fixes() on a detected checkerboard must mention key solutions."""
    n = 20
    colidx_y = np.arange(800, 800 + n)
    model_hs = np.where(colidx_y % 2 == 0, 13.0, 1.5)
    obs_hs = np.full(n, 5.0)
    times = _make_times(n)

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", model_hs - obs_hs),
            "colidx_y": ("time", colidx_y),
        },
        coords={"time": times},
    )
    cco = _build_cco(obs_hs, model_hs, times=times, extra_vars={"colidx_y": colidx_y})
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    report = outo.detect_numerical_patterns(checkerboard_threshold=0.70)
    text = outo.suggest_fixes()

    assert isinstance(text, str)
    assert len(text) > 0
    # Should mention the CFL / DTXY fix
    assert "DTXY" in text or "DTMAX" in text or "CFL" in text
    # Should cite Tolman (1992)
    assert "Tolman" in text and "1992" in text


def test_suggest_fixes_no_report(mock_cco):
    """suggest_fixes() with no pattern_report must print an informational message."""
    outo = outlier_class.__new__(outlier_class)
    outo.cco = mock_cco
    outo.vars = None
    outo.stats = {}
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = None
    outo.nID = None
    outo.varalias = None
    outo.sd = None
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    text = outo.suggest_fixes()
    assert "No pattern report" in text


# ---------------------------------------------------------------------------
# 14. patterns kwarg filters the registry
# ---------------------------------------------------------------------------


def test_detect_patterns_filter(mock_cco):
    """passing patterns=['checkerboard'] must skip all other detectors."""
    outo = outlier_class(mock_cco, n_std=3, method="zscore").populate()
    n = len(outo.vars["time"])
    outo.vars["colidx_y"] = ("time", np.arange(100, 100 + n))

    report = outo.detect_numerical_patterns(patterns=["checkerboard"])

    assert "checkerboard" in report
    # Other patterns must NOT appear
    for key in (
        "garden_sprinkler",
        "source_term_ringing",
        "hs_collapse",
        "spinup_insufficient",
    ):
        assert key not in report, f"'{key}' should not be in report when filtered"


# ---------------------------------------------------------------------------
# 17. required_spinup_hours — physics formula
# ---------------------------------------------------------------------------


def test_required_spinup_hours_formula():
    """
    Verify deep-water group velocity formula: c_g = g/(4π·f).
    For f = 0.1 Hz and 1 000 km domain: T ≈ 35.6 h.
    """
    g = 9.81
    f = 0.1
    c_g = g / (4 * np.pi * f)  # ≈ 7.806 m/s
    expected = 1000e3 / c_g / 3600.0  # ≈ 35.6 h

    h = required_spinup_hours(1000.0, f_low_hz=f)
    assert h == pytest.approx(expected, rel=1e-6)


def test_required_spinup_hours_scaling():
    """Larger domain → longer time; higher frequency → longer time."""
    h_base = required_spinup_hours(1000.0, f_low_hz=0.1)
    h_2x = required_spinup_hours(2000.0, f_low_hz=0.1)
    h_highf = required_spinup_hours(1000.0, f_low_hz=0.2)  # higher f → slower c_g

    assert h_2x == pytest.approx(2 * h_base, rel=1e-6)
    assert h_highf > h_base  # slower group velocity → more time needed


# ---------------------------------------------------------------------------
# 18. _detect_spinup_insufficient — detected
# ---------------------------------------------------------------------------


def _make_spinup_outo(n_early=15, n_late=15, run_start_str="2023-01-01 00:00"):
    """
    Synthetic outo where the early period shows:
    - model_Hs growing from 0.5 m to 2.5 m (energy slope > thresh)
    - mean early bias ≈ −1.5 m vs mean late bias ≈ 0 m (ratio >> 2)
    """
    run_start = pd.Timestamp(run_start_str)
    t_early = [run_start + pd.Timedelta(hours=h) for h in np.linspace(1, 47, n_early)]
    t_late = [run_start + pd.Timedelta(hours=h) for h in np.linspace(49, 96, n_late)]
    times = pd.DatetimeIndex(t_early + t_late)
    n = n_early + n_late

    obs_hs = np.full(n, 3.0)
    # Growing model Hs in early period; stable at obs level in late period
    mod_hs_early = np.linspace(0.5, 2.5, n_early)
    mod_hs_late = np.full(n_late, 3.0)
    model_hs = np.concatenate([mod_hs_early, mod_hs_late])
    bias_arr = model_hs - obs_hs

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", bias_arr),
        },
        coords={"time": times},
    )
    cco = _build_cco(obs_hs, model_hs, times=times)
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {
        "lo": -100,
        "hi": 100,
        "center": 0.0,
        "n_outliers": n,
        "n_total": n,
        "pct_outliers": 100.0,
        "method": "zscore",
    }
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = "mock"
    outo.nID = None
    outo.varalias = "Hs"
    outo.sd = run_start_str
    outo.ed = None
    outo.region = None
    outo.pattern_report = None
    return outo, run_start_str


def test_detect_spinup_insufficient_detected():
    """
    Growing model_Hs and strongly negative early bias must trigger detection.
    """
    outo, run_start_str = _make_spinup_outo(n_early=15, n_late=15)
    result = _detect_spinup_insufficient(
        outo,
        spinup_end_hours=48.0,
        rel_slope_thresh=0.02,
        early_bias_factor=2.0,
        min_points=4,
        run_start=run_start_str,
    )
    assert result["detected"] is True, f"Spin-up should be detected; got {result}"
    # At least one sub-test must fire
    assert result["energy_still_growing"] or result["bias_drift_detected"]
    # All early points (15) are affected
    assert len(result["affected_idx"]) == 15
    # early bias more negative than late
    assert result["early_mean_bias_m"] < result["late_mean_bias_m"]


# ---------------------------------------------------------------------------
# 19. _detect_spinup_insufficient — not detected (flat bias)
# ---------------------------------------------------------------------------


def test_detect_spinup_not_detected_flat():
    """
    Flat model_Hs and flat bias over the full period → not detected.
    """
    n = 30
    run_start = pd.Timestamp("2023-01-01 00:00")
    times = pd.DatetimeIndex(
        [run_start + pd.Timedelta(hours=h) for h in np.linspace(1, 96, n)]
    )
    obs_hs = np.full(n, 2.0)
    model_hs = np.full(n, 2.5)  # constant +0.5 m bias — no trend
    bias_arr = model_hs - obs_hs

    ds = xr.Dataset(
        {
            "obs_Hs": ("time", obs_hs),
            "model_Hs": ("time", model_hs),
            "obs_lons": ("time", np.linspace(0, 5, n)),
            "obs_lats": ("time", np.linspace(50, 55, n)),
            "bias": ("time", bias_arr),
        },
        coords={"time": times},
    )
    cco = _build_cco(obs_hs, model_hs, times=times)
    outo = outlier_class.__new__(outlier_class)
    outo.cco = cco
    outo.vars = ds
    outo.stats = {}
    outo.obs_var = "obs_Hs"
    outo.mod_var = "model_Hs"
    outo.n_std = 3
    outo.method = "zscore"
    outo.model = None
    outo.nID = None
    outo.varalias = None
    outo.sd = str(run_start)
    outo.ed = None
    outo.region = None
    outo.pattern_report = None

    result = _detect_spinup_insufficient(
        outo,
        spinup_end_hours=48.0,
        rel_slope_thresh=0.02,
        early_bias_factor=2.0,
        min_points=4,
        run_start=str(run_start),
    )
    assert (
        result["detected"] is False
    ), f"Should NOT detect spin-up on flat data; got {result}"
    assert result["affected_idx"] == []


# ---------------------------------------------------------------------------
# 20. detect_numerical_patterns includes spinup_insufficient in registry
# ---------------------------------------------------------------------------


def test_registry_includes_spinup():
    """PATTERN_REGISTRY must now include spinup_insufficient with 5 built-ins."""
    assert "spinup_insufficient" in PATTERN_REGISTRY
    entry = PATTERN_REGISTRY["spinup_insufficient"]
    assert "run_start" in entry["default_params"]
    assert len(entry["solutions"]) >= 3
    assert len(entry["references"]) >= 2


def test_detect_spinup_via_registry(mock_cco):
    """
    detect_numerical_patterns() must include spinup_insufficient in output
    (not-detected for mock_cco which has no temporal trend).
    """
    outo = outlier_class(mock_cco, n_std=3, method="zscore").populate()
    n = len(outo.vars["time"])
    outo.vars["colidx_y"] = ("time", np.arange(100, 100 + n))

    report = outo.detect_numerical_patterns(
        patterns=["spinup_insufficient"],
        run_start="2023-01-01",
    )
    assert "spinup_insufficient" in report
    # mock_cco has no temporal structure → should not fire
    assert report["spinup_insufficient"]["detected"] is False
