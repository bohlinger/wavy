#!/usr/bin/env python3
"""
wavy/outlier_patterns.py
------------------------
Built-in numerical-artifact pattern detectors for the wavy outlier module.

Each detector is a plain function with the contract::

    fn(outo: outlier_class, **params) -> dict

where the returned dict must contain at minimum:

    detected         : bool
    contribution_pct : float
    affected_idx     : list[int]

Detectors are registered via ``register_pattern()`` and are automatically
picked up by ``outlier_class.detect_numerical_patterns()`` and
``outlier_class.suggest_fixes()``.

Adding a new pattern
--------------------
1. Write a detection function (follow the contract above).
2. Call ``register_pattern(key, ...)`` at module level here.
3. Optionally extend the summary-builder labels in
   ``outlier_class.detect_numerical_patterns()`` if you want a
   letter-labelled section in the text report.

Users can also register patterns from their own code *after* import:

    from wavy.outlier_patterns import register_pattern
    register_pattern("my_pattern", ..., detect_fn=my_fn, ...)
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pattern Registry — open for user contributions
# ---------------------------------------------------------------------------

PATTERN_REGISTRY: dict = {}


def register_pattern(
    key: str,
    *,
    name: str,
    instability_id: str,
    description: str,
    detect_fn,
    solutions: list,
    references: list,
    default_params: dict | None = None,
) -> None:
    """
    Register a numerical-pattern detector with the global registry.

    Registered patterns are automatically picked up by
    ``outlier_class.detect_numerical_patterns()`` and
    ``outlier_class.suggest_fixes()``.

    Parameters
    ----------
    key            : str      – unique registry key (e.g. ``"checkerboard"``).
    name           : str      – human-readable pattern name.
    instability_id : str      – doc section reference (e.g. ``"1"``, ``"7"``).
    description    : str      – one-line description.
    detect_fn      : callable – ``fn(outo, **params) → dict`` with at minimum
                                keys ``detected`` (bool),
                                ``contribution_pct`` (float),
                                ``affected_idx`` (list[int]).
    solutions      : list[str] – known WW3 / model configuration fixes.
    references     : list[str] – key literature references.
    default_params : dict     – default keyword arguments forwarded to
                                ``detect_fn``.

    Example
    -------
    >>> def my_detector(outo, *, my_threshold=0.5):
    ...     affected = []      # fill with matching indices
    ...     return {"detected": bool(affected),
    ...             "contribution_pct": 0.0, "affected_idx": affected}
    >>> register_pattern(
    ...     "my_pattern",
    ...     name="My custom pattern",
    ...     instability_id="custom",
    ...     description="Detects my custom artifact.",
    ...     detect_fn=my_detector,
    ...     solutions=["Reduce time step."],
    ...     references=["Smith et al. (2000)."],
    ...     default_params={"my_threshold": 0.5},
    ... )
    """
    if key in PATTERN_REGISTRY:
        logger.warning("Pattern key '%s' is already registered; overwriting.", key)
    PATTERN_REGISTRY[key] = {
        "name": name,
        "instability_id": instability_id,
        "description": description,
        "detect_fn": detect_fn,
        "solutions": list(solutions),
        "references": list(references),
        "default_params": dict(default_params or {}),
    }


# ---------------------------------------------------------------------------
# Pattern detection functions  (operate on an outlier_class instance)
# ---------------------------------------------------------------------------


def _detect_checkerboard(outo, *, checkerboard_threshold=0.70):
    """
    Detect checkerboard / 2Δx spatial instability.

    A checkerboard arises when model Hs at even and odd ``colidx_y``
    (or ``colidx_x``) colocation-grid indices systematically differ —
    the spatial signature of a geographic CFL violation (instability #1).

    Parameters
    ----------
    outo                   : outlier_class
    checkerboard_threshold : float – minimum alternation fraction to flag
                                     (default 0.70).

    Returns
    -------
    dict : detected, parity_fraction, amplitude_m,
           contribution_pct, affected_idx.
    """
    n_total = int(outo.vars.sizes["time"])
    mod_hs = np.asarray(outo.vars[outo.mod_var])
    obs_hs = np.asarray(outo.vars[outo.obs_var])

    _empty = {
        "detected": False,
        "parity_fraction": 0.0,
        "amplitude_m": 0.0,
        "contribution_pct": 0.0,
        "affected_idx": [],
    }

    # Resolve grid index from outo.vars (preferred) or warn
    grid_idx = None
    if "colidx_y" in outo.vars:
        grid_idx = np.asarray(outo.vars["colidx_y"], dtype=int)
    elif "colidx_x" in outo.vars:
        grid_idx = np.asarray(outo.vars["colidx_x"], dtype=int)
    elif outo.cco is not None and (
        "colidx_y" in outo.cco.vars or "colidx_x" in outo.cco.vars
    ):
        # outo was saved before colidx_* was carried through populate().
        # Fall back to time-matching outo to cco to extract grid indices.
        logger.info(
            "colidx_* not in outo.vars; extracting from cco.vars by "
            "time-matching (legacy fallback)."
        )
        try:
            outo_ts = outo.vars.time.values.astype("datetime64[s]")
            cco_ts = pd.DatetimeIndex(outo.cco.vars.time.values).values.astype(
                "datetime64[s]"
            )
            cco_map = {t: i for i, t in enumerate(cco_ts)}
            cco_idx = np.array([cco_map.get(t, -1) for t in outo_ts])
            n_miss = int((cco_idx < 0).sum())
            if n_miss:
                logger.warning(
                    "%d outo timestamps unmatched in cco; "
                    "those will use index 0 (may affect checkerboard score).",
                    n_miss,
                )
            col_key = "colidx_y" if "colidx_y" in outo.cco.vars else "colidx_x"
            all_cols = np.asarray(outo.cco.vars[col_key])
            grid_idx = all_cols[np.clip(cco_idx, 0, len(all_cols) - 1)].astype(int)
        except Exception as exc:
            logger.warning(
                "colidx fallback from cco.vars failed (%s); "
                "skipping checkerboard detection.",
                exc,
            )
            return _empty
    else:
        logger.warning(
            "colidx_x / colidx_y not found in outo.vars or cco.vars; "
            "skipping checkerboard detection.",
        )
        return _empty

    if len(grid_idx) < 2:
        return _empty

    bias_vals = mod_hs - obs_hs
    parity = grid_idx % 2
    diffs = np.diff(grid_idx)
    sign_changes = np.diff(np.sign(bias_vals))
    consecutive_step1 = np.abs(diffs) == 1
    alternating = (sign_changes != 0) & consecutive_step1

    n_transitions = int(consecutive_step1.sum())
    n_alternating = int(alternating.sum())
    parity_fraction = n_alternating / n_transitions if n_transitions > 0 else 0.0
    detected = parity_fraction >= checkerboard_threshold

    even_mask = parity == 0
    odd_mask = parity == 1
    mean_even = float(np.mean(mod_hs[even_mask])) if even_mask.any() else float("nan")
    mean_odd = float(np.mean(mod_hs[odd_mask])) if odd_mask.any() else float("nan")
    amplitude = (
        abs(mean_even - mean_odd)
        if np.isfinite(mean_even) and np.isfinite(mean_odd)
        else 0.0
    )

    if detected:
        affected_set = set()
        for k in range(len(diffs)):
            if consecutive_step1[k] and alternating[k]:
                affected_set.add(k)
                affected_set.add(k + 1)
        affected_idx = sorted(affected_set)
    else:
        affected_idx = []

    contribution_pct = 100.0 * len(affected_idx) / n_total if n_total > 0 else 0.0

    return {
        "detected": detected,
        "parity_fraction": float(parity_fraction),
        "amplitude_m": float(amplitude),
        "contribution_pct": float(contribution_pct),
        "affected_idx": affected_idx,
    }


def _detect_garden_sprinkler(outo, *, window_min=30, bimodal_score_threshold=2.0):
    """
    Detect the Garden Sprinkler Effect (instability #7).

    Within a ±window_min-minute time window that contains ≥ 4 outlier
    points, compute the bimodal score on model_Hs:
    ``bimodal_score = (mean_high − mean_low) / std_all``.
    Flag windows where ``bimodal_score > bimodal_score_threshold``.

    Note: for a perfectly equal 50/50 bimodal split the score equals
    exactly 2.0.  Windows with an odd point count (e.g. 4 LOW + 5 HIGH)
    yield ~2.01 and satisfy the strict ``>`` criterion.

    Parameters
    ----------
    outo                    : outlier_class
    window_min              : int   – half-window in minutes (default 30).
    bimodal_score_threshold : float – minimum score to flag (default 2.0).

    Returns
    -------
    dict : detected, windows, contribution_pct, affected_idx.
    """
    n_total = int(outo.vars.sizes["time"])
    mod_hs = np.asarray(outo.vars[outo.mod_var])
    times = pd.DatetimeIndex(outo.vars["time"].values)
    times_sec = np.array([t.timestamp() for t in times])
    window_sec = window_min * 60.0

    flagged_windows = []
    affected_set = set()

    for i in range(n_total):
        dt = np.abs(times_sec - times_sec[i])
        in_win = np.where(dt <= window_sec)[0]
        if len(in_win) < 4:
            continue
        if i != int(in_win[0]):
            continue

        hs_win = mod_hs[in_win]
        std_all = float(np.std(hs_win))
        if std_all < 1e-12:
            continue

        sorted_hs = np.sort(hs_win)
        half = len(sorted_hs) // 2
        mean_low = float(np.mean(sorted_hs[:half]))
        mean_high = float(np.mean(sorted_hs[half:]))
        bimodal_score = (mean_high - mean_low) / std_all

        if bimodal_score > bimodal_score_threshold:
            center_time = times[in_win[len(in_win) // 2]]
            flagged_windows.append(
                {
                    "center_time": center_time,
                    "bimodal_score": float(bimodal_score),
                    "level_low_m": float(mean_low),
                    "level_high_m": float(mean_high),
                    "n_points": int(len(in_win)),
                    "indices": in_win.tolist(),
                }
            )
            affected_set.update(in_win.tolist())

    affected_idx = sorted(affected_set)
    detected = len(flagged_windows) > 0
    contribution_pct = 100.0 * len(affected_idx) / n_total if n_total > 0 else 0.0

    return {
        "detected": detected,
        "windows": flagged_windows,
        "contribution_pct": float(contribution_pct),
        "affected_idx": affected_idx,
    }


def _detect_source_term_ringing(
    outo,
    *,
    ringing_window_min: int = 5,
    ringing_threshold: float = 0.60,
    ringing_amplitude_m: float = 0.5,
):
    """
    Detect source-term stiffness / temporal ringing (instability #4).

    In short time windows, high-frequency oscillation of model_Hs
    (alternating sign of first differences) with non-trivial amplitude is
    the signature of stiff source-term integration under strong wind forcing
    or in very shallow water.

    Parameters
    ----------
    outo                 : outlier_class
    ringing_window_min   : int   – window length in minutes (default 5).
    ringing_threshold    : float – minimum fraction of alternating first
                                   differences to flag (default 0.60).
    ringing_amplitude_m  : float – minimum model_Hs range in the window [m]
                                   (default 0.5 m).

    Returns
    -------
    dict : detected, n_flagged_windows, mean_amplitude_m,
           contribution_pct, affected_idx.
    """
    n_total = int(outo.vars.sizes["time"])
    mod_hs = np.asarray(outo.vars[outo.mod_var])
    times = pd.DatetimeIndex(outo.vars["time"].values)
    times_sec = np.array([t.timestamp() for t in times])
    window_sec = ringing_window_min * 60.0

    flagged_windows = []
    affected_set = set()

    for i in range(n_total):
        dt = np.abs(times_sec - times_sec[i])
        in_win = np.where(dt <= window_sec)[0]
        if len(in_win) < 4:
            continue
        if i != int(in_win[0]):
            continue

        hs_win = mod_hs[in_win]
        amplitude = float(np.max(hs_win) - np.min(hs_win))
        if amplitude < ringing_amplitude_m:
            continue

        diffs = np.diff(hs_win)
        signs_nz = np.sign(diffs)
        signs_nz = signs_nz[signs_nz != 0]
        if len(signs_nz) < 2:
            continue

        n_alternations = int(np.sum(np.diff(signs_nz) != 0))
        ringing_fraction = n_alternations / (len(signs_nz) - 1)

        if ringing_fraction >= ringing_threshold:
            center_time = times[in_win[len(in_win) // 2]]
            flagged_windows.append(
                {
                    "center_time": center_time,
                    "ringing_fraction": float(ringing_fraction),
                    "amplitude_m": float(amplitude),
                    "n_points": int(len(in_win)),
                    "indices": in_win.tolist(),
                }
            )
            affected_set.update(in_win.tolist())

    affected_idx = sorted(affected_set)
    detected = len(flagged_windows) > 0
    amplitudes = [w["amplitude_m"] for w in flagged_windows]
    mean_amplitude = float(np.mean(amplitudes)) if amplitudes else 0.0
    contribution_pct = 100.0 * len(affected_idx) / n_total if n_total > 0 else 0.0

    return {
        "detected": detected,
        "n_flagged_windows": len(flagged_windows),
        "mean_amplitude_m": mean_amplitude,
        "contribution_pct": float(contribution_pct),
        "affected_idx": affected_idx,
    }


def _detect_hs_collapse(
    outo,
    *,
    hs_collapse_max_model: float = 0.05,
    hs_collapse_min_obs: float = 0.5,
):
    """
    Detect near-zero / collapsed model Hs (instability #10).

    Negative spectral energy densities are clipped to zero by WW3's
    non-negativity limiter, producing unrealistically low (near-zero)
    model Hs values where the satellite observes non-trivial wave heights.

    Parameters
    ----------
    outo                  : outlier_class
    hs_collapse_max_model : float – model_Hs threshold below which the field
                                    is considered collapsed [m] (default 0.05).
    hs_collapse_min_obs   : float – minimum obs_Hs to rule out genuinely calm
                                    conditions [m] (default 0.5).

    Returns
    -------
    dict : detected, n_collapsed, mean_obs_hs_m, contribution_pct, affected_idx.
    """
    n_total = int(outo.vars.sizes["time"])
    mod_hs = np.asarray(outo.vars[outo.mod_var])
    obs_hs = np.asarray(outo.vars[outo.obs_var])

    mask = (mod_hs < hs_collapse_max_model) & (obs_hs > hs_collapse_min_obs)
    affected_idx = np.where(mask)[0].tolist()
    detected = len(affected_idx) > 0

    contribution_pct = 100.0 * len(affected_idx) / n_total if n_total > 0 else 0.0
    mean_obs = float(np.mean(obs_hs[mask])) if detected else 0.0

    return {
        "detected": detected,
        "n_collapsed": len(affected_idx),
        "mean_obs_hs_m": mean_obs,
        "contribution_pct": float(contribution_pct),
        "affected_idx": affected_idx,
    }


def required_spinup_hours(
    domain_size_km: float,
    f_low_hz: float = 1.0 / 20.0,
) -> float:
    """
    Estimate the minimum model spin-up duration from domain geometry.

    The longest-period (lowest-frequency) swell travels at the highest
    group velocity and therefore determines how quickly energy from a
    distant boundary or storm can reach an interior validation point.
    This gives a *per-domain* lower bound on spin-up rather than the
    classic fixed 36–48 h rule of thumb.

    Deep-water linear wave theory:
        c_g = g / (4π · f_low)     [group velocity, m/s]
        T_required = domain_size / c_g   [seconds → hours]

    Parameters
    ----------
    domain_size_km : float – characteristic domain size, or distance from
                             the open-ocean boundary to the validation point
                             [km].
    f_low_hz       : float – lowest swell frequency of interest [Hz].
                             Default 1/20 Hz (20 s period swell).

    Returns
    -------
    hours : float – minimum spin-up duration [h].

    Notes
    -----
    For T = 20 s swell (f = 0.05 Hz): c_g ≈ 15.6 m/s → 56 km/h.
    A 1 500 km domain therefore needs at least ~26.7 h.
    For T = 10 s swell (f = 0.10 Hz): c_g ≈ 7.8 m/s → 28 km/h.
    A 1 500 km domain then needs ~53.4 h — already beyond the 48 h default.

    Examples
    --------
    >>> required_spinup_hours(1500, f_low_hz=1/10)   # 10-s swell, 1 500 km
    53.4
    >>> required_spinup_hours(500, f_low_hz=1/20)    # 20-s swell, 500 km
    8.9
    """
    g = 9.81  # m s⁻²
    c_g = g / (4.0 * np.pi * f_low_hz)  # deep-water group velocity [m/s]
    return (domain_size_km * 1e3) / c_g / 3600.0


def _detect_spinup_insufficient(
    outo,
    *,
    spinup_end_hours: float = 48.0,
    trend_window_h: float = 6.0,
    rel_slope_thresh: float = 0.02,
    early_bias_factor: float = 2.0,
    min_points: int = 4,
    run_start: str | None = None,
):
    """
    Detect outliers caused by insufficient model spin-up (instability #0).

    Two complementary diagnostics are applied from the collocated dataset:

    * **Test A — energy growth rate**: if model Hs² is still growing at
      more than ``rel_slope_thresh`` per hour in the early period (measured
      as ``d(ln Hs²)/dt``), the model has not equilibrated — local wind-sea
      is still developing.

    * **Test B — temporal bias drift**: if the mean bias in the early period
      is at least ``early_bias_factor`` times more negative than in the late
      period, and it improves over time, the model field was systematically
      underdeveloped at spin-up end.  This proxy catches both local wind-sea
      disequilibrium *and* remote swell in transit (which would otherwise
      require frequency-resolved spectral output to detect directly).

    Use ``required_spinup_hours()`` for an a-priori domain-size-based
    estimate before deciding on the nominal spin-up duration.

    Parameters
    ----------
    outo              : outlier_class
    spinup_end_hours  : float – nominal spin-up duration from model run
                                start [h] (default 48).
    trend_window_h    : float – trailing window length [h]; reserved for
                                future spectral extension (default 6).
    rel_slope_thresh  : float – relative energy slope threshold [/h].
                                If ``d(ln Hs²)/dt > rel_slope_thresh`` in
                                the early period, energy is still growing
                                (default 0.02 = 2 % per hour).
    early_bias_factor : float – |early_mean_bias| / |late_mean_bias| ratio
                                above which bias drift is flagged
                                (default 2.0).
    min_points        : int   – minimum outlier points per period for each
                                test to run (default 4).
    run_start         : str, optional – ISO-8601 model run start time.
                                If None, resolution order is:
                                ``outo.sd`` → ``cco.sd`` → first
                                observation time (warning issued).

    Returns
    -------
    dict : detected, energy_still_growing, bias_drift_detected,
           early_mean_bias_m, late_mean_bias_m, energy_slope_per_hr,
           bias_ratio, n_early, n_late, affected_idx,
           contribution_pct, run_start_used.
    """
    n_total = int(outo.vars.sizes["time"])
    mod_hs = np.asarray(outo.vars[outo.mod_var])
    bias = np.asarray(outo.vars["bias"])
    times = pd.DatetimeIndex(outo.vars["time"].values)

    _empty: dict = {
        "detected": False,
        "energy_still_growing": False,
        "bias_drift_detected": False,
        "early_mean_bias_m": float("nan"),
        "late_mean_bias_m": float("nan"),
        "energy_slope_per_hr": float("nan"),
        "bias_ratio": float("nan"),
        "n_early": 0,
        "n_late": 0,
        "affected_idx": [],
        "contribution_pct": 0.0,
        "run_start_used": "",
    }

    # ---- Resolve model run start ----------------------------------------
    t_ref: pd.Timestamp | None = None

    if run_start is not None:
        try:
            t_ref = pd.Timestamp(run_start)
        except Exception:
            logger.warning(
                "Could not parse run_start='%s'; falling back to outo.sd.", run_start
            )

    if t_ref is None:
        for attr_src in (outo, getattr(outo, "cco", None)):
            sd_val = getattr(attr_src, "sd", None) if attr_src is not None else None
            if sd_val is not None:
                try:
                    t_ref = pd.Timestamp(sd_val)
                    break
                except Exception:
                    pass

    if t_ref is None:
        t_ref = times[0]
        logger.warning(
            "run_start not supplied and not found on outo/cco.sd; "
            "using first observation time (%s) as t = 0 for spin-up "
            "detection. Pass run_start='YYYY-MM-DD HH' for accuracy.",
            t_ref,
        )

    run_start_str = str(t_ref)

    # ---- Elapsed time from run start [h] ---------------------------------
    elapsed_h = np.array([(t - t_ref).total_seconds() / 3600.0 for t in times])

    valid_mask = elapsed_h >= 0.0
    if valid_mask.sum() < 2:
        return {**_empty, "run_start_used": run_start_str}

    early_mask = valid_mask & (elapsed_h < spinup_end_hours)
    late_mask = valid_mask & (elapsed_h >= spinup_end_hours)
    n_early = int(early_mask.sum())
    n_late = int(late_mask.sum())

    # ---- Test A: relative energy growth rate in the early period --------
    energy_still_growing = False
    energy_slope_per_hr = float("nan")

    if n_early >= min_points:
        t_e = elapsed_h[early_mask]
        hs_e = np.maximum(mod_hs[early_mask], 1e-3)  # guard log(0)
        log_e = np.log(hs_e**2)
        try:
            coeffs = np.polyfit(t_e, log_e, 1)
            energy_slope_per_hr = float(coeffs[0])
            energy_still_growing = energy_slope_per_hr > rel_slope_thresh
        except np.linalg.LinAlgError:
            pass

    # ---- Test B: temporal bias drift ------------------------------------
    bias_drift_detected = False
    early_mean_bias = float("nan")
    late_mean_bias = float("nan")
    bias_ratio = float("nan")

    if n_early >= min_points:
        early_mean_bias = float(np.mean(bias[early_mask]))
    if n_late >= min_points:
        late_mean_bias = float(np.mean(bias[late_mask]))

    if n_early >= min_points and n_late >= min_points:
        late_ref = max(abs(late_mean_bias), 0.1)  # prevent division by ~0
        bias_ratio = abs(early_mean_bias) / late_ref
        bias_drift_detected = (
            early_mean_bias < 0  # spin-up → model underestimates
            and early_mean_bias < late_mean_bias  # bias improves over time
            and bias_ratio >= early_bias_factor
        )

    # ---- Combine --------------------------------------------------------
    detected = energy_still_growing or bias_drift_detected
    affected_idx = sorted(np.where(early_mask)[0].tolist()) if detected else []
    contribution_pct = 100.0 * len(affected_idx) / n_total if n_total > 0 else 0.0

    return {
        "detected": detected,
        "energy_still_growing": energy_still_growing,
        "bias_drift_detected": bias_drift_detected,
        "early_mean_bias_m": early_mean_bias,
        "late_mean_bias_m": late_mean_bias,
        "energy_slope_per_hr": energy_slope_per_hr,
        "bias_ratio": float(bias_ratio) if np.isfinite(bias_ratio) else bias_ratio,
        "n_early": n_early,
        "n_late": n_late,
        "affected_idx": affected_idx,
        "contribution_pct": float(contribution_pct),
        "run_start_used": run_start_str,
    }


# ---------------------------------------------------------------------------
# Built-in pattern registrations
# ---------------------------------------------------------------------------

register_pattern(
    "checkerboard",
    name="Checkerboard / 2Δx spatial instability",
    instability_id="1",
    description=(
        "Alternating high/low model Hs on consecutive colocation-grid cells — "
        "signature of a geographic CFL violation."
    ),
    detect_fn=_detect_checkerboard,
    solutions=[
        "Reduce the propagation time step DTXY (or DTMAX) in ww3_grid.inp so "
        "the Courant number C_g * Δt / Δx < 1 everywhere.",
        "Use the implicit PDLIB solver for unstructured grids — removes the "
        "CFL restriction on Δx at the cost of a sparse linear solve.",
        "Smooth very small grid cells near the poles or in narrow straits, or "
        "increase Δx locally to relax the CFL constraint.",
    ],
    references=[
        "Tolman, H.L. (1992). Effects of numerics on the physics in a "
        "third-generation wind-wave model. J. Phys. Oceanogr. 22, 1770–1786.",
        "WW3 Development Group. WAVEWATCH III User Manual (NOAA/NCEP), "
        "section on time stepping and CFL.",
    ],
    default_params={"checkerboard_threshold": 0.70},
)

register_pattern(
    "garden_sprinkler",
    name="Garden Sprinkler Effect (GSE)",
    instability_id="7",
    description=(
        "Discretely banded Hs patterns downstream of distant storms — caused "
        "by finite directional/frequency resolution spreading swell as "
        "discrete 'jets' rather than a continuous field."
    ),
    detect_fn=_detect_garden_sprinkler,
    solutions=[
        "Increase directional spectral resolution (NDIR in ww3_grid.inp); "
        "doubling NDIR is the most robust fix but roughly doubles cost.",
        "Activate the GSE divergence alleviation scheme (Tolman 2002b): "
        "set FLAGTR > 0 in the propagation input.",
        "Apply a diffusion tensor in spectral space (Booij & Holthuijsen 1987): "
        "non-zero SDMAX / DDMAX in the propagation input.",
        "Use directional averaging (Lavrenov & Onvlee 1995) as a cheaper "
        "alternative on structured grids.",
        "Note: GSE alleviation is currently unavailable for unstructured grids.",
    ],
    references=[
        "Booij, N. and Holthuijsen, L.H. (1987). Propagation of ocean waves in "
        "discrete spectral wave models. J. Comput. Phys. 68, 307–326.",
        "Lavrenov, I.V. and Onvlee, J. (1995). On the directional spreading in "
        "discrete spectral wave models. J. Phys. Oceanogr. 25, 62–71.",
        "Tolman, H.L. (2002b). Alleviating the Garden Sprinkler Effect in wind "
        "wave models. Ocean Modelling 4, 269–289.",
    ],
    default_params={"window_min": 30, "bimodal_score_threshold": 2.0},
)

register_pattern(
    "source_term_ringing",
    name="Source-term stiffness / temporal ringing",
    instability_id="4",
    description=(
        "High-frequency oscillation of model Hs in time at isolated points — "
        "signature of stiff source-term integration under strong wind forcing "
        "or in very shallow water."
    ),
    detect_fn=_detect_source_term_ringing,
    solutions=[
        "Reduce DTMIN (minimum source-term sub-step) so the adaptive algorithm "
        "can shrink the step more aggressively in high-forcing cells.",
        "Ensure the fully implicit source-term integration scheme "
        "(Hargreaves & Annan 2000) is active — this is the WW3 default but "
        "can be disabled by the IMPLCT switch.",
        "Reduce the overall model time step in domains with hurricane-force "
        "winds or very shallow nested grids.",
        "Apply the spectral change limiter (Tolman 2002a) to cap growth per "
        "time step — this can introduce bias in rapidly evolving spectra.",
    ],
    references=[
        "Hargreaves, J.C. and Annan, J.D. (2000). Comments on 'Improvement of "
        "the short-fetch behaviour in the WAM model'. "
        "J. Atmos. Ocean. Technol. 17, 498–503.",
        "Tolman, H.L. (2002a). Testing of WAVEWATCH III version 2.22 in "
        "NCEP's NWW3 ocean wave forecasting system (NOAA/NCEP Tech Note 214).",
        "WW3 Development Group. WAVEWATCH III User Manual (NOAA/NCEP), "
        "section on source term time stepping and DTMIN.",
    ],
    default_params={
        "ringing_window_min": 5,
        "ringing_threshold": 0.60,
        "ringing_amplitude_m": 0.5,
    },
)

register_pattern(
    "hs_collapse",
    name="Near-zero Hs collapse (negative-energy clipping)",
    instability_id="10",
    description=(
        "Model Hs near zero where the satellite observes non-trivial waves — "
        "caused by negative spectral energy densities clipped to zero "
        "after an over-dissipating source-term update."
    ),
    detect_fn=_detect_hs_collapse,
    solutions=[
        "Reduce the source-term time step (DTMIN) so individual spectral bins "
        "are not over-dissipated in a single update.",
        "Verify the non-negativity clipping switch is active (it is the "
        "default; check that no custom build flag disabled it).",
        "Inspect the whitecapping / dissipation parameterisation coefficients "
        "(BETAMAX, SDSBR, etc.) — over-tuned dissipation often drives this.",
        "Increase source-term time step resolution in very fine coastal nests.",
    ],
    references=[
        "WW3 Development Group. WAVEWATCH III User Manual (NOAA/NCEP), "
        "section on non-negativity of spectral densities.",
        "Tolman, H.L. (1992). Effects of numerics on the physics in a "
        "third-generation wind-wave model. J. Phys. Oceanogr. 22, 1770–1786.",
    ],
    default_params={
        "hs_collapse_max_model": 0.05,
        "hs_collapse_min_obs": 0.5,
    },
)

register_pattern(
    "spinup_insufficient",
    name="Insufficient model spin-up",
    instability_id="0",
    description=(
        "Model Hs systematically underestimated early in the run — "
        "wind-sea not yet equilibrated or remote swell still in transit "
        "at the end of the nominal spin-up period."
    ),
    detect_fn=_detect_spinup_insufficient,
    solutions=[
        "Run required_spinup_hours(domain_size_km) for a physics-based "
        "minimum estimate before deciding on spin-up length. For domains "
        "> 500 km with 10-s swell the required time often exceeds 48 h.",
        "Use a warm start from a long prior run (e.g. a multi-year "
        "reanalysis or global model) rather than a calm-sea IC, to "
        "pre-populate the swell field and avoid zero-energy at t=0.",
        "Extend the spin-up period: 72–96 h is safer for swell-dominated "
        "or large open-ocean domains.",
        "Archive 1-D frequency spectra at validation points in addition to "
        "bulk Hs, so the low-frequency swell-arrival diagnostic can be "
        "applied directly to spectral output (Test B in the full framework).",
        "If your run straddles a storm onset, verify that both the storm "
        "growth and the swell propagation delay are fully covered by the "
        "spin-up window.",
    ],
    references=[
        "Tolman, H.L. (1992). Effects of numerics on the physics in a "
        "third-generation wind-wave model. J. Phys. Oceanogr. 22, 1770–1786.",
        "WW3 Development Group. WAVEWATCH III User Manual (NOAA/NCEP), "
        "section on initial conditions and spin-up.",
        "required_spinup_hours() in this module implements the a-priori "
        "geometric estimate: T_req = domain_size / c_g(f_low), "
        "c_g = g / (4π f_low) [deep-water linear wave theory].",
    ],
    default_params={
        "spinup_end_hours": 48.0,
        "trend_window_h": 6.0,
        "rel_slope_thresh": 0.02,
        "early_bias_factor": 2.0,
        "min_points": 4,
        "run_start": None,
    },
)
