#!/usr/bin/env python3
"""
wavy/outlier_module.py
----------------------
Identify and visualise outliers in a wavy collocation_class object.

bias = model_Hs - obs_Hs  (point-wise). Three interchangeable detection
methods are supported via the `method` argument (see `_detect_outliers`):

    "zscore"  (default)  |bias_i - mean_bias|        > n_std * std_bias
    "mad"                |modified z-score (median/MAD)| > n_std
    "iqr"                bias_i outside [Q1 - n_std*IQR, Q3 + n_std*IQR]

"mad" and "iqr" are robust to the outliers themselves skewing the
threshold (the classic masking problem with plain Z-scores) and are
recommended when bias distributions are heavy-tailed (e.g. during storms).


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import timedelta
import os
import logging


import pickle
import dill
import xarray as xr
from copy import deepcopy

import cmocean
from wavy.model_module import model_class as mc
from wavy.utils import hour_rounder
from wavy.gridder_module import gridder_class as gc
from wavy.grid_stats import apply_metric
from wavy.outlier_patterns import (
    PATTERN_REGISTRY,
    register_pattern,
    required_spinup_hours,
    _detect_checkerboard,
    _detect_garden_sprinkler,
    _detect_source_term_ringing,
    _detect_hs_collapse,
    _detect_spinup_insufficient,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _round_to_hour(t):
    """
    Round a timestamp to the nearest hour using wavy's own hour_rounder,
    ensuring consistency with the model file loader logic.
    """
    t = pd.Timestamp(t)
    # hour_rounder expects a datetime-like with .minute / .hour / .replace
    from datetime import datetime
    dt = datetime(t.year, t.month, t.day, t.hour, t.minute,
                  t.second, t.microsecond)
    return pd.Timestamp(hour_rounder(dt))


def _draw_colored_track(ax, lons, lats, values, cmap, norm, projection,
                         linewidth=2.5, zorder=10):
    """
    Draw the satellite track as a continuous line whose colour varies with
    `values`.  Coordinates are projected to the map's native system so that
    a plain LineCollection (no cartopy transform) is placed correctly for
    any projection.

    Returns the LineCollection (attach a colorbar to it).
    """
    pts = projection.transform_points(ccrs.PlateCarree(),
                                       np.asarray(lons),
                                       np.asarray(lats))
    x, y = pts[:, 0], pts[:, 1]
    values = np.asarray(values)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(values)
    if valid.sum() < 2:
        return None

    x, y, values = x[valid], y[valid], values[valid]

    # build [N-1, 2, 2] segment array
    points   = np.stack([x, y], axis=1).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    seg_vals = (values[:-1] + values[1:]) / 2.0   # midpoint colour

    lc = LineCollection(segments, cmap=cmap, norm=norm,
                        linewidth=linewidth, zorder=zorder,
                        capstyle="round", joinstyle="round")
    lc.set_array(seg_vals)
    ax.add_collection(lc)
    return lc


def _detect_outliers(bias_all, n_std=3, method="zscore"):
    """
    Core outlier-detection logic, shared by `find_outliers`,
    `plot_bias_distribution`, and `plot_bias_timeseries`, so that all three
    always agree on exactly which points are flagged.

    Parameters
    ----------
    bias_all : np.ndarray   – point-wise bias (model - obs), may contain NaN.
    n_std    : float        – threshold multiplier. Meaning depends on
                              `method` (sigma for "zscore", modified
                              z-score cutoff for "mad", IQR multiplier
                              "k" for "iqr"). Default 3.
    method   : {"zscore", "mad", "iqr"}
        "zscore" – classic mean/std threshold. Simple, but large outliers
                   inflate std and can mask other outliers ("masking").
        "mad"    – median + MAD-based modified z-score (0.6745 factor
                   makes it ~comparable to a standard Z-score under
                   normality). Robust to the outliers themselves.
        "iqr"    – Tukey-style fence: outside [Q1 - n_std*IQR,
                   Q3 + n_std*IQR]. Also robust; a common default for
                   `n_std` here is 1.5 (mild) to 3 (extreme), not 3σ.

    Returns
    -------
    outlier_mask : np.ndarray[bool]  – same shape as bias_all.
    stats : dict – mean_bias, std_bias, median_bias, mad_bias, q1, q3, iqr,
                   center, threshold, lo, hi, method, n_outliers, n_total,
                   pct_outliers.
    """
    bias_all   = np.asarray(bias_all, dtype=float)
    mask_valid = np.isfinite(bias_all)
    valid      = bias_all[mask_valid]

    if valid.size == 0:
        raise ValueError("No finite bias values to run outlier detection on.")

    mean_bias   = float(np.mean(valid))
    std_bias    = float(np.std(valid))
    median_bias = float(np.median(valid))
    mad_bias    = float(np.median(np.abs(valid - median_bias)))
    q1, q3      = (float(v) for v in np.percentile(valid, [25, 75]))
    iqr         = q3 - q1

    if method == "zscore":
        center    = mean_bias
        threshold = n_std * std_bias
        lo, hi    = center - threshold, center + threshold
        outlier_mask = mask_valid & (np.abs(bias_all - center) > threshold)

    elif method == "mad":
        center = median_bias
        scale  = mad_bias if mad_bias > 1e-12 else 1e-12
        modified_z = 0.6745 * (bias_all - center) / scale
        threshold  = n_std  # in modified-z units, not metres
        lo = center - (n_std * scale / 0.6745)
        hi = center + (n_std * scale / 0.6745)
        outlier_mask = mask_valid & (np.abs(modified_z) > n_std)

    elif method == "iqr":
        center    = median_bias
        threshold = n_std * iqr  # n_std acts as the IQR multiplier "k"
        lo, hi    = q1 - threshold, q3 + threshold
        outlier_mask = mask_valid & ((bias_all < lo) | (bias_all > hi))

    else:
        raise ValueError(
            f"Unknown method '{method}'; choose 'zscore', 'mad', or 'iqr'."
        )

    n_outliers = int(outlier_mask.sum())
    n_total    = int(mask_valid.sum())

    stats = {
        "method":       method,
        "mean_bias":    mean_bias,
        "std_bias":     std_bias,
        "median_bias":  median_bias,
        "mad_bias":     mad_bias,
        "q1":           q1,
        "q3":           q3,
        "iqr":          iqr,
        "center":       center,
        "threshold":    threshold,
        "lo":           lo,
        "hi":           hi,
        "n_outliers":   n_outliers,
        "n_total":      n_total,
        "pct_outliers": 100.0 * n_outliers / max(1, n_total),
    }
    return outlier_mask, stats


def _format_threshold_message(stats, n_std):
    """Human-readable one-liner describing the detection threshold used."""
    method = stats["method"]
    if method == "zscore":
        rule = (f"|bias − {stats['mean_bias']:+.3f}| > "
                f"{n_std} × {stats['std_bias']:.3f} = {stats['threshold']:.3f} m")
    elif method == "mad":
        rule = (f"|modified z-score (median={stats['median_bias']:+.3f}, "
                f"MAD={stats['mad_bias']:.3f})| > {n_std}  "
                f"[≈ outside {stats['lo']:.3f} … {stats['hi']:.3f} m]")
    else:  # iqr
        rule = (f"bias outside [Q1−{n_std}×IQR, Q3+{n_std}×IQR] = "
                f"[{stats['lo']:.3f}, {stats['hi']:.3f}] m "
                f"(Q1={stats['q1']:.3f}, Q3={stats['q3']:.3f})")
    return (
        f"[{method}] threshold = {rule}\n"
        f"  → {stats['n_outliers']} / {stats['n_total']} points flagged "
        f"({stats['pct_outliers']:.2f} %)"
    )


def _lon_extent(lons, margin_deg=5.0, global_span_threshold=180.0):
    """
    Compute a longitude (lonmin, lonmax) extent that correctly handles
    antimeridian (dateline) crossing, e.g. a track running
    [..., 179, -179, -178, ...]. Naively taking min/max of such values
    gives a near-global bounding box instead of the small span the track
    actually covers.

    Works by normalising to [0, 360), finding the single largest angular
    gap between sorted points (assumed to be "outside" the track), and
    unwrapping around that gap. cartopy's PlateCarree correctly renders
    extents that stray outside the canonical [-180, 180] range, so the
    result does not need to be re-wrapped.

    This unwrap logic assumes a single contiguous cluster of points (a
    track segment, one storm event). For data scattered across most/all
    of the globe (e.g. a global collocation dataset), there is no
    meaningful gap to unwrap around, and forcing one produces an extent
    that is *technically* not global (e.g. [176.5, 530]) but close
    enough to trigger a known cartopy/shapely gridliner bug: a degenerate
    map-boundary polygon (`GEOSException: Points of LinearRing do not
    form a closed linestring`). If the data's angular coverage exceeds
    `global_span_threshold` degrees, fall back to the plain global extent
    instead, since a "zoom" spanning more than that isn't a useful
    zoom anyway.
    """
    lons = np.mod(np.asarray(lons, dtype=float), 360.0)
    lons = lons[np.isfinite(lons)]
    if lons.size == 0:
        return -margin_deg, margin_deg
    if lons.size == 1:
        return float(lons[0]) - margin_deg, float(lons[0]) + margin_deg

    lons_sorted = np.sort(lons)
    gaps  = np.diff(np.concatenate([lons_sorted, [lons_sorted[0] + 360.0]]))
    split = int(np.argmax(gaps))
    largest_gap = float(gaps[split])

    if (360.0 - largest_gap) > global_span_threshold:
        return -180.0, 180.0

    lons_unwrapped = np.concatenate([
        lons_sorted[split + 1:],
        lons_sorted[:split + 1] + 360.0,
    ])
    lonmin = float(lons_unwrapped[0]) - margin_deg
    lonmax = float(lons_unwrapped[-1]) + margin_deg

    # if the span comfortably fits within a normal -180..180 window,
    # shift it back for a more readable extent
    if lonmax > 180.0 and lonmin - 360.0 >= -180.0:
        lonmin -= 360.0
        lonmax -= 360.0

    return lonmin, lonmax


def _save_figure(fig, save_path, dpi=200):
    """
    Save a figure to disk, creating parent directories as needed.

    Defensively retries once with all gridlines removed if the initial
    save fails. This guards against a known cartopy/shapely gridliner bug
    (`GEOSException: Points of LinearRing do not form a closed
    linestring`) triggered when a map's boundary polygon becomes
    degenerate for certain projection/extent combinations — so that one
    bad gridliner doesn't take down an entire batch run.
    """
    os.makedirs(os.path.dirname(os.path.abspath(str(save_path))), exist_ok=True)
    try:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    except Exception as exc:
        logger.warning(
            f"savefig failed ({exc.__class__.__name__}: {exc}); "
            f"retrying '{save_path}' with gridlines removed"
        )
        for ax in fig.axes:
            for gl in list(getattr(ax, "_gridliners", [])):
                try:
                    gl.remove()
                except Exception:
                    pass
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    logger.info("Saved → %s", save_path)

def load_from_dill(path):
    """Load any dill-serialized object (e.g. a cco) from disk."""
    with open(path, "rb") as fh:
        return dill.load(fh)


# ---------------------------------------------------------------------------
# 1.  Outlier detection
# ---------------------------------------------------------------------------

def find_outliers(cco, n_std=3, method="zscore", obs_var=None, mod_var=None):
    """
    Identify outlier collocated pairs from the point-wise bias.

    Parameters
    ----------
    cco      : collocation_class  – populated collocation object.
    n_std    : float              – threshold multiplier (default 3).
                                    Meaning depends on `method` — see
                                    `_detect_outliers`.
    method   : {"zscore", "mad", "iqr"} – detection method (default
                                    "zscore"). "mad" or "iqr" are
                                    recommended when the bias distribution
                                    is heavy-tailed, since a few large
                                    outliers can otherwise inflate the
                                    plain std and mask other outliers.
    obs_var  : str, optional      – obs variable name; auto-detected if None.
    mod_var  : str, optional      – model variable name; auto-detected if None.

    Returns
    -------
    df_outliers : pd.DataFrame or None
        Columns: time, obs_lons, obs_lats, <obs_var>, <mod_var>,
        bias, [model_time].
    stats : dict
        method, mean_bias, std_bias, median_bias, mad_bias, q1, q3, iqr,
        center, threshold, lo, hi, n_outliers, n_total, pct_outliers,
        obs_var, mod_var.
    """
    if obs_var is None:
        obs_var = "obs_Hs" if "obs_Hs" in cco.vars else "obs_values"
    if mod_var is None:
        mod_var = "model_Hs" if "model_Hs" in cco.vars else "model_values"
    if obs_var not in cco.vars or mod_var not in cco.vars:
        raise KeyError(
            f"obs_var='{obs_var}' / mod_var='{mod_var}' not found in "
            f"cco.vars; available: {list(cco.vars.data_vars)}"
        )

    obs   = np.asarray(cco.vars[obs_var])
    mod   = np.asarray(cco.vars[mod_var])
    times = pd.DatetimeIndex(cco.vars["time"].values)
    lons  = np.asarray(cco.vars["obs_lons"])
    lats  = np.asarray(cco.vars["obs_lats"])

    bias_all = mod - obs
    outlier_mask, stats = _detect_outliers(bias_all, n_std=n_std, method=method)
    stats["obs_var"] = obs_var
    stats["mod_var"] = mod_var

    logger.info("%s", _format_threshold_message(stats, n_std))

    if stats["n_outliers"] == 0:
        return None, stats

    idx = np.where(outlier_mask)[0]
    df = pd.DataFrame({
        "time":     times[idx],
        "obs_lons": lons[idx],
        "obs_lats": lats[idx],
        obs_var:    obs[idx],
        mod_var:    mod[idx],
        "bias":     bias_all[idx],
    })
    if "model_time" in cco.vars:
        df["model_time"] = pd.DatetimeIndex(
            np.asarray(cco.vars["model_time"])[idx]
        )

    return df.sort_values("time").reset_index(drop=True), stats


# ---------------------------------------------------------------------------
# 2.  Single-event diagnostic map
# ---------------------------------------------------------------------------

def plot_outlier_map(
    cco,
    outlier_time,
    model_nID=None,
    model_name=None,
    track_window_min=30,
    projection=None,
    bb=None,
    margin_deg=5.0,
    n_std_label=3,
    obs_var=None,
    mod_var=None,
    track_color_by="bias",
    vmin_hs=0,
    vmax_hs=14,
    levels_incr=0.5,
    save_path=None,
):
    """
    Single-panel diagnostic map for one outlier event.

    Draws the 2D model Hs field (contourf) with a colour-mapped satellite
    track (±track_window_min minutes) overlaid, and the outlier location
    marked with a star.

    Parameters
    ----------
    cco              : collocation_class
    outlier_time     : datetime-like
    model_nID        : str, optional     – model nID from model_cfg.yaml;
                                           falls back to cco.model.
    track_window_min : int               – half-window in minutes (def 30).
    projection       : cartopy.crs       – e.g.
                         ccrs.Stereographic(central_latitude=90,
                                            central_longitude=-30,
                                            true_scale_latitude=90)
                         Defaults to PlateCarree.
    bb               : (lonmin,lonmax,latmin,latmax) – explicit map extent;
                                           auto-computed if None.
    margin_deg       : float             – auto-extent margin (def 5 °).
    n_std_label      : float             – for annotation only.
    obs_var/mod_var  : str, optional     – override variable names.
    track_color_by   : "bias" or "obs"   – what to colour the track with.
    vmin_hs/vmax_hs  : float             – model Hs colour scale [m].
    levels_incr      : float             – Hs contour step [m].
    save_path        : str/Path, optional – saved at 200 dpi if provided.

    Returns
    -------
    fig : matplotlib.Figure or None
    """
    if obs_var is None:
        obs_var = "obs_Hs" if "obs_Hs" in cco.vars else "obs_values"
    if mod_var is None:
        mod_var = "model_Hs" if "model_Hs" in cco.vars else "model_values"
    if model_nID is None:
        model_nID = getattr(cco, "model", None)
    if model_nID is None:
        raise ValueError("model_nID must be supplied or cco.model must be set.")

    outlier_time = pd.Timestamp(outlier_time)

    # ---- 1. satellite track window ----------------------------------------
    t_lo = outlier_time - timedelta(minutes=track_window_min)
    t_hi = outlier_time + timedelta(minutes=track_window_min)
    ds_track = cco.vars.sel(time=slice(t_lo, t_hi))

    if len(ds_track.time) == 0:
        logger.info(
            "[plot_outlier_map] No collocated points in "
            "±%d min around %s", track_window_min, outlier_time
        )
        return None

    track_lons = np.asarray(ds_track["obs_lons"])
    track_lats = np.asarray(ds_track["obs_lats"])
    track_obs  = np.asarray(ds_track[obs_var])
    track_mod  = np.asarray(ds_track[mod_var])
    track_bias = track_mod - track_obs

    # ---- 2. nearest outlier point -----------------------------------------
    times_track = pd.DatetimeIndex(ds_track.time.values)
    nearest_idx = int(np.argmin(np.abs(times_track - outlier_time)))
    out_lon  = float(track_lons[nearest_idx])
    out_lat  = float(track_lats[nearest_idx])
    out_obs  = float(track_obs[nearest_idx])
    out_mod  = float(track_mod[nearest_idx])
    out_bias = out_mod - out_obs

    # ---- 3. map extent ----------------------------------------------------
    if bb is not None:
        lonmin, lonmax, latmin, latmax = bb
    else:
        lonmin, lonmax = _lon_extent(track_lons, margin_deg=margin_deg)
        latmin = max(float(np.min(track_lats)) - margin_deg,  -90)
        latmax = min(float(np.max(track_lats)) + margin_deg,   90)

    # ---- 4. load 2D model field -------------------------------------------
    # IMPORTANT: round to the nearest hour so that generate_bestguess_leadtime
    # and _get_model_filedate produce a valid file path.  Sub-hour timestamps
    # cause the internal hour_rounder to overshoot the init_time and return None.
    mco = None
    model_time_h = _round_to_hour(outlier_time)
    # 'name' is required when the model path template uses a name-based
    # substitution (e.g. ensemble member / config variant).
    # Priority: explicit model_name arg > cco.name > omit.
    _mc_name = model_name if model_name is not None else getattr(cco, "name", None)
    _mc_kwargs = {
        "nID": model_nID,
        "sd":  str(model_time_h),
        "ed":  str(model_time_h),
    }
    if _mc_name is not None:
        _mc_kwargs["name"] = _mc_name

    logger.info(
        "[plot_outlier_map] Loading model '%s'%s at %s "
        "(rounded from %s) …",
        model_nID,
        f" name='{_mc_name}'" if _mc_name else "",
        model_time_h,
        outlier_time.strftime("%H:%M:%S"),
    )
    try:
        mco = mc(**_mc_kwargs).populate()
    except Exception as exc:
        logger.warning("Could not load model field: %s", exc)

    # ---- 5. single-panel figure -------------------------------------------
    if projection is None:
        projection = ccrs.PlateCarree()

    fig = plt.figure(figsize=(10, 9))
    ax  = fig.add_subplot(1, 1, 1, projection=projection)

    levels_hs = np.arange(vmin_hs, vmax_hs + levels_incr, levels_incr)
    cmap_hs   = cmocean.cm.amp
    norm_hs   = mpl.colors.BoundaryNorm(levels_hs, cmap_hs.N)

    land = cfeature.GSHHSFeature(
        scale="i", levels=[1], facecolor=cfeature.COLORS["land"]
    )
    ax.add_feature(land, facecolor="burlywood", alpha=0.5, zorder=5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, zorder=6)
    try:
        ax.set_extent([lonmin, lonmax, latmin, latmax],
                       crs=ccrs.PlateCarree())
    except Exception as exc:
        logger.warning(f"set_extent failed: {exc}")

        # draw_labels=True causes a shapely GEOSException with polar / stereographic
    # projections when cartopy tries to convert the circular map boundary into a
    # LinearRing.  Only safe for PlateCarree where the boundary is rectangular.
    _draw_labels = isinstance(projection, ccrs.PlateCarree)
    gl = ax.gridlines(draw_labels=_draw_labels, linewidth=0.4,
                       color="gray", alpha=0.5, linestyle="--")
    if _draw_labels:
        gl.top_labels   = False
        gl.right_labels = False

    # ---- model 2D Hs contourf ---------------------------------------------
    if mco is not None:
        try:
            Mlons = mco.vars["lons"].values
            Mlats = mco.vars["lats"].values
            if Mlons.ndim == 1:
                Mlons, Mlats = np.meshgrid(Mlons, Mlats)
            Mhs_da = mco.vars["Hs"]
            if "time" in Mhs_da.dims:
                Mhs = Mhs_da.sel(time=model_time_h, method="nearest")
            else:
                Mhs = Mhs_da
            cf = ax.contourf(
                Mlons, Mlats, Mhs.values,
                levels=levels_hs, cmap=cmap_hs, norm=norm_hs,
                transform=ccrs.PlateCarree(), extend="max",
            )
            fig.colorbar(cf, ax=ax, label="Hs model [m]",
                         fraction=0.033, pad=0.04,
                         ticks=levels_hs[::4])
        except Exception as exc:
            logger.warning("Could not plot model field: %s", exc)

    # ---- satellite track (continuous coloured line) -----------------------
    if track_color_by == "obs":
        color_vals  = track_obs
        track_cmap  = cmocean.cm.amp
        track_norm  = norm_hs
        track_label = obs_var.replace("obs_", "Sat obs ") + " [m]"
    else:  # "bias"
        color_vals   = track_bias
        track_cmap   = mpl.cm.RdBu_r
        bias_abs_max = max(float(np.abs(track_bias).max()), 0.5)
        track_norm   = mpl.colors.TwoSlopeNorm(
            vmin=-bias_abs_max, vcenter=0.0, vmax=bias_abs_max
        )
        track_label = "Bias model−obs [m]"

    lc = _draw_colored_track(
        ax, track_lons, track_lats, color_vals,
        cmap=track_cmap, norm=track_norm,
        projection=projection,
        linewidth=2.5, zorder=12,
    )
    if lc is not None:
        fig.colorbar(lc, ax=ax, label=track_label,
                     fraction=0.033, pad=0.09, shrink=0.8)

    # ---- outlier highlight ------------------------------------------------
    ax.scatter(
        out_lon, out_lat,
        c="red", s=200, marker="*",
        edgecolors="black", linewidths=0.8,
        transform=ccrs.PlateCarree(), zorder=15,
        label=f"Outlier (|bias|>{n_std_label}σ)",
    )
    ax.legend(loc="lower right", fontsize=9, framealpha=0.85)

    ax.set_title(
        f"Outlier — {model_nID}\n"
        f"{outlier_time.strftime('%Y-%m-%d %H:%M UTC')}  "
        f"({out_lat:.2f}°N, {out_lon:.2f}°E)  "
        f"bias={out_bias:+.3f} m  obs={out_obs:.2f} m  mod={out_mod:.2f} m",
        fontsize=9,
    )
    fig.tight_layout()

    if save_path is not None:
        _save_figure(fig, save_path)

    return fig


# ---------------------------------------------------------------------------
# 4.  Bias distribution — histogram with threshold(s) marked
# ---------------------------------------------------------------------------

def plot_bias_distribution(
    cco,
    n_std=3,
    method="zscore",
    obs_var=None,
    mod_var=None,
    bins=50,
    save_path=None,
):
    """
    Histogram of the point-wise bias with the detection threshold(s)
    overlaid, so the effect of `n_std`/`method` can be sanity-checked
    before generating a batch of per-event maps.

    Parameters
    ----------
    cco             : collocation_class
    n_std, method   : same meaning as in `find_outliers`.
    obs_var/mod_var : str, optional – override variable names.
    bins            : int – histogram bin count (default 50).
    save_path       : str/Path, optional – saved at 200 dpi if provided.

    Returns
    -------
    fig : matplotlib.Figure
    """
    if obs_var is None:
        obs_var = "obs_Hs" if "obs_Hs" in cco.vars else "obs_values"
    if mod_var is None:
        mod_var = "model_Hs" if "model_Hs" in cco.vars else "model_values"

    obs = np.asarray(cco.vars[obs_var])
    mod = np.asarray(cco.vars[mod_var])
    bias_all = mod - obs

    outlier_mask, stats = _detect_outliers(bias_all, n_std=n_std, method=method)
    logger.info("%s", _format_threshold_message(stats, n_std))

    finite_bias = bias_all[np.isfinite(bias_all)]

    fig, ax = plt.subplots(figsize=(8, 5))

    # Exactly 50 bins between -5 m and +5 m
    bin_edges = np.linspace(-2.5, 2.5, bins + 1)

    ax.hist(
        finite_bias,
        bins=bin_edges,
        color="steelblue",
        alpha=0.75,
        edgecolor="white",
        linewidth=0.4,
        zorder=2,
    )

    ax.axvline(stats["lo"], color="crimson", linewidth=1.2, linestyle="--",
               label=f"threshold lo = {stats['lo']:.3f} m")
    ax.axvline(stats["hi"], color="crimson", linewidth=1.2, linestyle="--",
               label=f"threshold hi = {stats['hi']:.3f} m")
    ax.axvline(stats["center"], color="black", linewidth=1.0, linestyle="-",
               label=f"center = {stats['center']:.3f} m")
    
    # Shade the flagged tails
    left_tail = finite_bias[finite_bias < stats["lo"]]
    right_tail = finite_bias[finite_bias > stats["hi"]]

    if left_tail.size:
        ax.hist(
            left_tail,
            bins=bin_edges,
            color="crimson",
            alpha=0.6,
            zorder=2,
        )

    if right_tail.size:
        ax.hist(
            right_tail,
            bins=bin_edges,
            color="crimson",
            alpha=0.6,
            zorder=2,
        )

    ax.set_xlim(-2.5, 2.5)
    ax.set_xlabel("Bias, model − obs [m]")
    ax.set_ylabel("Count")
    ax.set_title(
        f"Bias distribution — method='{method}', n_std={n_std}\n"
        f"{stats['n_outliers']} / {stats['n_total']} points flagged "
        f"({stats['pct_outliers']:.2f} %)",
        fontsize=10,
    )
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    fig.tight_layout()

    if save_path is not None:
        _save_figure(fig, save_path)

    return fig


# ---------------------------------------------------------------------------
# 5.  Time series — obs vs. model, and bias, with outliers marked
# ---------------------------------------------------------------------------

def plot_bias_timeseries(
    cco,
    n_std=3,
    method="zscore",
    obs_var=None,
    mod_var=None,
    save_path=None,
):
    """
    Two-panel time series: (top) obs vs. model Hs, (bottom) bias with the
    detection threshold band shaded. Outlier points are highlighted in
    both panels. Useful for seeing whether outliers are isolated spikes
    (bad collocations, spurious retrievals) or part of a systematic
    drift / storm event that the model missed.

    Parameters
    ----------
    cco             : collocation_class
    n_std, method   : same meaning as in `find_outliers`.
    obs_var/mod_var : str, optional – override variable names.
    save_path       : str/Path, optional – saved at 200 dpi if provided.

    Returns
    -------
    fig : matplotlib.Figure
    """
    if obs_var is None:
        obs_var = "obs_Hs" if "obs_Hs" in cco.vars else "obs_values"
    if mod_var is None:
        mod_var = "model_Hs" if "model_Hs" in cco.vars else "model_values"

    obs   = np.asarray(cco.vars[obs_var])
    mod   = np.asarray(cco.vars[mod_var])
    times = pd.DatetimeIndex(cco.vars["time"].values)
    bias_all = mod - obs

    outlier_mask, stats = _detect_outliers(bias_all, n_std=n_std, method=method)
    logger.info("%s", _format_threshold_message(stats, n_std))

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 7), sharex=True,
        gridspec_kw={"height_ratios": [1.3, 1]},
    )

    # ---- top: obs vs model ------------------------------------------------
    ax1.plot(times, obs, color="tab:blue", linewidth=1.0, label=obs_var, zorder=2)
    ax1.plot(times, mod, color="tab:orange", linewidth=1.0, label=mod_var, zorder=2)
    ax1.scatter(times[outlier_mask], obs[outlier_mask], color="red", s=25,
                zorder=5, label="outlier (obs)")
    ax1.scatter(times[outlier_mask], mod[outlier_mask], facecolors="none",
                edgecolors="red", s=25, zorder=5, label="outlier (model)")
    ax1.set_ylabel("Hs [m]")
    ax1.legend(loc="upper right", fontsize=8, ncol=2, framealpha=0.85)
    ax1.set_title(
        f"Obs vs. model, and bias — method='{method}', n_std={n_std}  "
        f"({stats['n_outliers']} / {stats['n_total']} flagged, "
        f"{stats['pct_outliers']:.2f} %)",
        fontsize=10,
    )

    # ---- bottom: bias with threshold band ---------------------------------
    ax2.axhspan(stats["lo"], stats["hi"], color="seagreen", alpha=0.12, zorder=1,
                label=f"threshold band ({stats['lo']:.2f} … {stats['hi']:.2f} m)")
    ax2.axhline(stats["center"], color="black", linewidth=1.0, zorder=2)
    ax2.plot(times, bias_all, color="grey", linewidth=0.8, zorder=3)
    ax2.scatter(times[outlier_mask], bias_all[outlier_mask], color="red", s=25,
                zorder=5, label="outlier")
    ax2.set_ylabel("Bias, model − obs [m]")
    ax2.set_xlabel("Time")
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.85)

    fig.autofmt_xdate()
    fig.tight_layout()

    if save_path is not None:
        _save_figure(fig, save_path)

    return fig


# ---------------------------------------------------------------------------
# 6.  Batch wrapper — top N events ranked by severity
# ---------------------------------------------------------------------------

def plot_all_outliers(
    cco,
    key,
    model_nID=None,
    model_name=None,
    n_std=3,
    method="zscore",
    track_window_min=30,
    n_top=10,
    projection=None,
    bb=None,
    margin_deg=5.0,
    track_color_by="bias",
    make_overview=True,
    make_distribution=True,
    make_timeseries=True,
    save_dir=None,
):
    """
    Detect all outliers, group them into independent events, rank by severity
    (max |bias| per event), and plot the top `n_top`. Optionally also
    produces the three summary diagnostics (overview map, bias histogram,
    time series) so a single call gives a complete picture.

    Parameters
    ----------
    cco              : collocation_class
    key              : str             – label used in saved filenames.
    model_nID        : str, optional   – falls back to cco.model.
    n_std            : float           – threshold multiplier (default 3);
                                          meaning depends on `method`.
    method           : {"zscore", "mad", "iqr"} – detection method
                                          (default "zscore").
    track_window_min : int             – grouping + track half-window (def 30).
    n_top            : int             – per-event maps to plot (default 10).
    projection       : cartopy.crs, optional
    bb               : (lonmin,lonmax,latmin,latmax), optional
    margin_deg       : float           – auto-extent margin (default 5 °).
    track_color_by   : "bias" or "obs"
    make_overview    : bool  – also produce `plot_outlier_overview` (def True).
    make_distribution: bool  – also produce `plot_bias_distribution` (def True).
    make_timeseries  : bool  – also produce `plot_bias_timeseries` (def True).
    save_dir         : str/Path, optional – PNGs saved here.

    Returns
    -------
    figs        : list[matplotlib.Figure]
        Summary diagnostics (overview/distribution/timeseries, if enabled)
        come first, followed by the per-event maps in rank order.
    df_outliers : pd.DataFrame or None
    stats       : dict
    """
    df_outliers, stats = find_outliers(cco, n_std=n_std, method=method)

    figs = []
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)

    if make_distribution:
        sp = os.path.join(save_dir, f"outlier_distribution_{key}.png") if save_dir else None
        fig_d = plot_bias_distribution(cco, n_std=n_std, method=method, save_path=sp)
        figs.append(fig_d)
        plt.close(fig_d)

    if make_timeseries:
        sp = os.path.join(save_dir, f"outlier_timeseries_{key}.png") if save_dir else None
        fig_t = plot_bias_timeseries(cco, n_std=n_std, method=method, save_path=sp)
        figs.append(fig_t)
        plt.close(fig_t)

    if make_overview and df_outliers is not None:
        sp = os.path.join(save_dir, f"outlier_overview_{key}.png") if save_dir else None
        _tmp = outlier_class.__new__(outlier_class)
        _tmp.cco     = cco
        _tmp.stats   = stats
        _tmp.n_std   = n_std
        _tmp.method  = method
        _tmp.obs_var = stats.get("obs_var", "obs_Hs")
        _tmp.mod_var = stats.get("mod_var", "model_Hs")
        _tmp.vars = xr.Dataset(
            {
                "obs_lons": (["time"], df_outliers["obs_lons"].values),
                "obs_lats": (["time"], df_outliers["obs_lats"].values),
                "bias":     (["time"], df_outliers["bias"].values),
            },
            coords={"time": df_outliers["time"].values},
        )
        fig_o = _tmp.plot_overview(
            projection=projection, bb=bb, margin_deg=margin_deg,
            save_path=sp,
        )
        if fig_o is not None:
            figs.append(fig_o)
            plt.close(fig_o)

    if df_outliers is None:
        logger.info("[plot_all_outliers] No outliers found.")
        return figs, None, stats

    times = pd.DatetimeIndex(df_outliers["time"])

    # ---- group into independent events ------------------------------------
    used        = np.zeros(len(times), dtype=bool)
    event_times = []
    event_scores = []

    for i, t in enumerate(times):
        if not used[i]:
            window_mask = np.array([
                abs((t2 - t).total_seconds()) <= track_window_min * 60
                for t2 in times
            ])
            event_times.append(t)
            event_scores.append(
                float(np.abs(df_outliers.loc[window_mask, "bias"].values).max())
            )
            used[window_mask] = True

    # ---- rank descending by max |bias|, keep top n_top -------------------
    order        = np.argsort(event_scores)[::-1]
    top_times    = [event_times[i]  for i in order[:n_top]]
    top_scores   = [event_scores[i] for i in order[:n_top]]

    logger.info(
        "[plot_all_outliers] %d event(s) detected — "
        "plotting top %d by max |bias|:",
        len(event_times), len(top_times),
    )
    for rank, (t, s) in enumerate(zip(top_times, top_scores), start=1):
        logger.info(
            "  #%2d  %s  max|bias|=%.3f m",
            rank, pd.Timestamp(t).strftime("%Y-%m-%d %H:%M"), s,
        )

    # ---- produce per-event maps (appended after the summary diagnostics) --
    for rank, t in enumerate(top_times, start=1):
        t_str = pd.Timestamp(t).strftime("%Y%m%dT%H%M")
        sp = (
            os.path.join(save_dir,
                         f"outlier_top{rank:02d}_{key}_{t_str}.png")
            if save_dir is not None
            else None
        )
        logger.info("  --- #%d/%d at %s ---", rank, len(top_times), pd.Timestamp(t))
        fig = plot_outlier_map(
            cco=cco,
            outlier_time=t,
            model_nID=model_nID,
            model_name=model_name,
            track_window_min=track_window_min,
            projection=projection,
            bb=bb,
            margin_deg=margin_deg,
            n_std_label=n_std,
            track_color_by=track_color_by,
            save_path=sp,
        )
        if fig is not None:
            figs.append(fig)
            plt.close(fig)

    return figs, df_outliers, stats


# ---------------------------------------------------------------------------
# Outlier class  (outo)
# ---------------------------------------------------------------------------


class outlier_class:
    """
    An outlier object (``outo``) is a structured subset of a collocation
    object: it stores only the flagged observation/model pairs together with
    their coordinates, time, bias, and detection metadata.

    The xarray Dataset at ``outo.vars`` mirrors ``cco.vars`` column names so
    the existing plotting helpers can be reused directly.

    Parameters
    ----------
    cco     : collocation_class   – populated collocation object.
    n_std   : float               – threshold multiplier (default 3).
    method  : {"zscore","mad","iqr"} – detection method (default "zscore").
    obs_var : str, optional       – obs variable; auto-detected if None.
    mod_var : str, optional       – model variable; auto-detected if None.

    Examples
    --------
    >>> outo = outlier_class(cco, n_std=3, method="mad").populate()
    >>> outo.write_to_nc("outliers_2023.nc")
    """

    def __init__(self, cco, n_std=3, method="zscore",
                 obs_var=None, mod_var=None):
        logger.debug("Initializing outlier_class")
        self.cco    = cco
        self.stats  = None
        self.vars   = None
        self.n_std  = n_std
        self.method = method
        self.obs_var = obs_var or (
            "obs_Hs" if "obs_Hs" in cco.vars else "obs_values"
        )
        self.mod_var = mod_var or (
            "model_Hs" if "model_Hs" in cco.vars else "model_values"
        )
        # mirror key metadata from cco
        self.model    = getattr(cco, "model",    None)
        self.nID      = getattr(cco, "nID",      None)
        self.varalias = getattr(cco, "varalias", None)
        self.sd       = getattr(cco, "sd",       None)
        self.ed       = getattr(cco, "ed",       None)
        self.region   = getattr(cco, "region",   None)
        self.pattern_report = None
        logger.debug(
            "outlier_class initialised: nID=%s model=%s method=%s n_std=%s",
            self.nID, self.model, method, n_std,
        )

    # ------------------------------------------------------------------
    # Core method
    # ------------------------------------------------------------------

    def populate(self, **kwargs):
        """
        Run outlier detection and build the xarray Dataset of flagged points.

        All collocated variables that exist in ``cco.vars`` are preserved in
        ``outo.vars``; ``bias`` (model − obs) is added.  Detection statistics
        are stored both in ``outo.stats`` and as ``outo.vars`` global
        attributes so they survive a round-trip through NetCDF.

        Returns
        -------
        outo : outlier_class   – new instance with .vars and .stats set.
        """
        new = deepcopy(self)
        logger.debug("populate(): running outlier detection")

        cco     = new.cco
        obs_var = new.obs_var
        mod_var = new.mod_var

        obs      = np.asarray(cco.vars[obs_var])
        mod      = np.asarray(cco.vars[mod_var])
        times    = pd.DatetimeIndex(cco.vars["time"].values)
        lons     = np.asarray(cco.vars["obs_lons"])
        lats     = np.asarray(cco.vars["obs_lats"])
        bias_all = mod - obs

        outlier_mask, stats = _detect_outliers(
            bias_all, n_std=new.n_std, method=new.method
        )
        stats["obs_var"] = obs_var
        stats["mod_var"] = mod_var
        new.stats = stats

        logger.info("%s", _format_threshold_message(stats, new.n_std))

        if stats["n_outliers"] == 0:
            logger.info("No outliers found – outo.vars is None.")
            new.vars = None
            return new

        idx = np.where(outlier_mask)[0]

        # ---- core variables (always present) --------------------------
        ds_dict = {
            "obs_lons": (["time"], lons[idx]),
            "obs_lats": (["time"], lats[idx]),
            obs_var:    (["time"], obs[idx]),
            mod_var:    (["time"], mod[idx]),
            "bias":     (["time"], bias_all[idx]),
        }
        coords = {"time": times[idx]}

        # ---- optional variables carried from cco.vars -----------------
        for extra in ("dist", "model_time", "model_lons", "model_lats",
                      "colidx_x", "colidx_y"):
            if extra in cco.vars:
                ds_dict[extra] = (["time"], np.asarray(cco.vars[extra])[idx])

        # ---- detection stats as global attrs (survive NetCDF save) ----
        ds_attrs = {k: v for k, v in stats.items()
                    if isinstance(v, (str, int, float))}
        ds_attrs["title"] = "outlier_class dataset"

        new.vars = xr.Dataset(ds_dict, coords=coords, attrs=ds_attrs)
        new.sd = str(pd.Timestamp(new.vars.time.values.min()))
        new.ed = str(pd.Timestamp(new.vars.time.values.max()))

        logger.info(
            "%d outlier(s) detected (%.2f %% of %d valid points).",
            stats["n_outliers"], stats["pct_outliers"], stats["n_total"],
        )
        return new

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def write_to_pickle(self, path):
        """Serialize the whole outlier_class object (incl. parent cco)."""
        os.makedirs(os.path.dirname(os.path.abspath(str(path))),
                    exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump(self, fh)
        logger.info("Saved outlier object (pickle) → %s", path)

    def write_to_nc(self, path):
        """
        Write ``outo.vars`` to a NetCDF file.
        Detection stats are preserved as global attributes.
        """
        if self.vars is None:
            logger.info("Nothing to save – vars is None (no outliers).")
            return
        os.makedirs(os.path.dirname(os.path.abspath(str(path))),
                    exist_ok=True)
        self.vars.to_netcdf(path)
        logger.info("Saved outlier vars (NetCDF) → %s", path)

    @classmethod
    def read_from_pickle(cls, path):
        """Load a previously serialized outlier_class object."""
        with open(path, "rb") as fh:
            return pickle.load(fh)

    @classmethod
    def read_from_dill(cls, path):
        """Load a previously serialized outlier_class object from a dill file."""
        with open(path, "rb") as fh:
            return dill.load(fh)

    @classmethod
    def read_from_nc(cls, path):
        """
        Reconstruct a minimal outlier_class from a saved NetCDF.
        ``outo.vars`` and ``outo.stats`` are populated from the file;
        ``outo.cco`` will be None (the raw cco is not stored in the NC).
        """
        obj = cls.__new__(cls)
        obj.vars   = xr.open_dataset(path)
        obj.stats  = dict(obj.vars.attrs)
        obj.cco    = None
        obj.n_std  = obj.stats.get("n_std",  3)
        obj.method = obj.stats.get("method", "zscore")
        return obj

    @classmethod
    def from_dill(cls, path, n_std=3, method="zscore",
                obs_var=None, mod_var=None):
        """
        Load a collocation object from a dill file and construct an
        outlier_class ready to populate.

        Parameters
        ----------
        path    : str/Path  – path to the .dill file.
        n_std, method, obs_var, mod_var : same as __init__.

        Returns
        -------
        outo : outlier_class  – not yet populated; call .populate() next.

        Example
        -------
        >>> outo = outlier_class.from_dill(
        ...     "BETAMAX_165_IC1__global__spinup_36h.dill",
        ...     n_std=3, method="mad",
        ... ).populate()
        """

        with open(path, "rb") as fh:
            cco = dill.load(fh)
        return cls(cco, n_std=n_std, method=method,
                obs_var=obs_var, mod_var=mod_var)

    # ------------------------------------------------------------------
    # Analysis and plotting methods
    # ------------------------------------------------------------------

    def plot_distribution(self, bins=50, save_path=None):
        """
        Histogram of the full bias (from cco) with detection thresholds.
        Requires self.cco to be set.
        """
        if self.cco is None:
            raise ValueError("cco is None; load the outo with read_from_dill.")
        return plot_bias_distribution(
            self.cco, n_std=self.n_std, method=self.method,
            obs_var=self.obs_var, mod_var=self.mod_var,
            bins=bins, save_path=save_path,
        )

    def plot_overview(self, projection=None, bb=None, margin_deg=5.0,
                  point_size=40, show_full_track=False,
                  mode="scatter",
                  vmin=None, vmax=None,
                  density_res=(1.0, 1.0),
                  save_path=None):
        """
        Map of all outlier locations.

        Parameters
        ----------
        projection   : cartopy.crs, optional – default PlateCarree.
        bb           : (lonmin,lonmax,latmin,latmax), optional
        mode         : "scatter" (default) or "density"
            "scatter" – one dot per outlier, coloured by bias.
            "density" – 2-D histogram of outlier counts per grid cell.
        vmin, vmax   : float, optional – bias colour limits for scatter mode.
                        Auto-computed (symmetric around 0) if not provided.
        density_res  : tuple(float, float) – grid resolution for density mode (default (1.0, 1.0)).
        point_size   : float – scatter marker size (default 40).
        show_full_track : bool – draw the full cco track in grey (default False).
        save_path    : str/Path, optional
        """
        if self.vars is None:
            logger.info("No outliers to plot.")
            return None
        if projection is None:
            projection = ccrs.PlateCarree()

        lons = self.vars["obs_lons"].values
        lats = self.vars["obs_lats"].values
        bias = self.vars["bias"].values

        # ---- map extent -------------------------------------------------------
        if bb is not None:
            lonmin, lonmax, latmin, latmax = bb
        else:
            lonmin, lonmax = _lon_extent(lons, margin_deg=margin_deg)
            latmin = max(float(lats.min()) - margin_deg, -90)
            latmax = min(float(lats.max()) + margin_deg,  90)

        # ---- figure -----------------------------------------------------------
        fig = plt.figure(figsize=(10, 8))
        ax  = fig.add_subplot(1, 1, 1, projection=projection)

        land = cfeature.GSHHSFeature(scale="i", levels=[1],
                                    facecolor=cfeature.COLORS["land"])
        ax.add_feature(land, facecolor="burlywood", alpha=0.5, zorder=2)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.6, zorder=3)

        if show_full_track and self.cco is not None:
            try:
                all_lons = np.asarray(self.cco.vars["obs_lons"])
                all_lats = np.asarray(self.cco.vars["obs_lats"])
                ax.plot(all_lons, all_lats, color="grey", linewidth=0.6,
                        alpha=0.5, transform=ccrs.PlateCarree(), zorder=4)
            except Exception as exc:
                logger.warning(f"could not draw full track: {exc}")

        try:
            ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
        except Exception as exc:
            logger.warning(f"set_extent failed: {exc}")

        extent_span = lonmax - lonmin
        if extent_span < 350.0:
            _draw_labels = isinstance(projection, ccrs.PlateCarree)
            gl = ax.gridlines(draw_labels=_draw_labels, linewidth=0.4,
                            color="gray", alpha=0.5, linestyle="--")
            if _draw_labels:
                gl.top_labels = False
                gl.right_labels = False

        # ---- scatter or density -----------------------------------------------
        if mode == "density":
            gco = gc(
                lons=lons, lats=lats, values=bias,
                bb=(lonmin, lonmax, latmin, latmax),
                res=density_res,
                varalias=self.obs_var,
            )
            gridvar, lon_grid, lat_grid = apply_metric(gco=gco)

            count_grid = np.ma.masked_equal(gridvar["nov"], 0)

            cf = ax.pcolormesh(
                lon_grid + density_res[0] / 2,   # shift to cell centres
                lat_grid + density_res[1] / 2,
                count_grid,
                cmap="YlOrRd",
                transform=ccrs.PlateCarree(),
                zorder=5,
            )
            fig.colorbar(cf, ax=ax, label="Number of outliers per cell",
                        fraction=0.033, pad=0.06, shrink=0.85)

        else:  # scatter
            if vmin is not None and vmax is not None:
                # diverging norm only when 0 lies strictly between vmin and vmax
                if vmin < 0 < vmax:
                    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
                else:
                    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            else:
                bias_abs_max = max(float(np.abs(bias).max()), 0.5)
                norm = mpl.colors.TwoSlopeNorm(
                    vmin=-bias_abs_max, vcenter=0.0, vmax=bias_abs_max
                )
            sc = ax.scatter(
                lons, lats, c=bias, cmap=mpl.cm.RdBu_r, norm=norm,
                s=point_size, edgecolors="black", linewidths=0.4,
                transform=ccrs.PlateCarree(), zorder=10,
            )
            fig.colorbar(sc, ax=ax, label="Bias model−obs [m]",
                        fraction=0.033, pad=0.06, shrink=0.85)

        # ---- title ------------------------------------------------------------
        method_used = (self.stats.get("method", self.method)
                    if self.stats else self.method)
        n_tot   = self.stats["n_total"] if self.stats else None
        subtitle = f"{len(lons)} outlier(s)"
        if n_tot:
            subtitle += f" / {n_tot} points ({100.0 * len(lons) / n_tot:.2f} %)"
        ax.set_title(
            f"Outlier overview — method='{method_used}', n_std={self.n_std}"
            f"  [{mode}]\n{subtitle}",
            fontsize=10,
        )
        fig.tight_layout()

        if save_path is not None:
            _save_figure(fig, save_path)

        return fig

    def group_events(self, track_window_min=30):
        """
        Group outlier points into independent events by time proximity.

        Two points belong to the same event when their timestamps are
        within ``track_window_min`` minutes.  Events are seeded greedily
        in chronological order and ranked descending by max |bias|.

        Returns
        -------
        events : pd.DataFrame
            Columns: event_id, peak_time, max_bias, n_points, indices.
        """
        if self.vars is None:
            return pd.DataFrame(
                columns=["event_id", "peak_time", "max_bias",
                        "n_points", "indices"]
            )
        times   = pd.DatetimeIndex(self.vars.time.values)
        bias    = self.vars["bias"].values
        used    = np.zeros(len(times), dtype=bool)
        records = []
        event_id = 0

        for i, t in enumerate(times):
            if not used[i]:
                dt_sec = np.array([
                    abs((t2 - t).total_seconds()) for t2 in times
                ])
                in_window  = dt_sec <= track_window_min * 60
                idx_window = np.where(in_window)[0]
                abs_bias   = np.abs(bias[in_window])
                peak_local = int(np.argmax(abs_bias))
                records.append({
                    "event_id":  event_id,
                    "peak_time": times[idx_window[peak_local]],
                    "max_bias":  float(abs_bias.max()),
                    "n_points":  int(in_window.sum()),
                    "indices":   idx_window.tolist(),
                })
                used[in_window] = True
                event_id += 1

        return (
            pd.DataFrame(records)
            .sort_values("max_bias", ascending=False)
            .reset_index(drop=True)
        )

    # ------------------------------------------------------------------
    # Numerical pattern detection
    # ------------------------------------------------------------------

    def detect_numerical_patterns(
        self,
        window_min: int = 30,
        checkerboard_threshold: float = 0.70,
        bimodal_score_threshold: float = 2.0,
        patterns: list | None = None,
        **pattern_kwargs,
    ) -> dict:
        """
        Scan the outlier set for registered numerical model artifacts.

        Iterates over every entry in the module-level ``PATTERN_REGISTRY``
        (or a user-selected subset via ``patterns``), calls each detector,
        collects results, computes RMSE/bias impact, and prints a summary.

        New patterns can be added without touching this method: call
        ``register_pattern()`` at module level and they will automatically
        appear here.

        Parameters
        ----------
        window_min               : half-window for garden-sprinkler (and
                                   source-term-ringing) grouping [min].
        checkerboard_threshold   : minimum alternation fraction for checkerboard.
        bimodal_score_threshold  : minimum bimodal score for garden-sprinkler.
        patterns                 : list of registry keys to run, or None for all.
        **pattern_kwargs         : any additional keyword args forwarded to the
                                   relevant detector (override default_params).

        Returns
        -------
        report : dict
            One key per registered pattern (e.g. ``"checkerboard"``,
            ``"garden_sprinkler"``, ``"source_term_ringing"``,
            ``"hs_collapse"``), plus ``"summary"`` (str).
            Each pattern value is the dict returned by its ``detect_fn``.
            ``self.pattern_report`` is set to the returned dict.
        """
        if self.vars is None:
            raise ValueError(
                "outo.vars is None – run .populate() before "
                "detect_numerical_patterns()."
            )

        n_total = int(self.vars.sizes["time"])
        times   = pd.DatetimeIndex(self.vars["time"].values)

        # Consolidate all caller-supplied kwargs (explicit + **pattern_kwargs)
        all_kwargs: dict = {
            "window_min":               window_min,
            "checkerboard_threshold":   checkerboard_threshold,
            "bimodal_score_threshold":  bimodal_score_threshold,
            **pattern_kwargs,
        }

        # Determine which patterns to run
        keys_to_run = (
            [k for k in PATTERN_REGISTRY if k in patterns]
            if patterns is not None
            else list(PATTERN_REGISTRY)
        )

        pattern_results: dict = {}
        for key in keys_to_run:
            entry = PATTERN_REGISTRY[key]
            # Merge: start from registered defaults, override with caller kwargs
            params = dict(entry["default_params"])
            for p in list(params):
                if p in all_kwargs:
                    params[p] = all_kwargs[p]
            try:
                pattern_results[key] = entry["detect_fn"](self, **params)
            except Exception as exc:
                logger.warning(
                    "Pattern detector '%s' raised %s: %s — skipping.",
                    key, exc.__class__.__name__, exc,
                )
                pattern_results[key] = {
                    "detected": False,
                    "contribution_pct": 0.0,
                    "affected_idx": [],
                    "error": str(exc),
                }

        # ------------------------------------------------------------------
        # Impact on model statistics (RMSE / bias, full collocated dataset)
        # ------------------------------------------------------------------
        impact_lines: list = []
        if self.cco is not None:
            try:
                full_obs      = np.asarray(self.cco.vars[self.obs_var])
                full_mod      = np.asarray(self.cco.vars[self.mod_var])
                full_bias_arr = full_mod - full_obs
                valid         = np.isfinite(full_bias_arr)

                cco_times_s  = np.array(
                    pd.DatetimeIndex(self.cco.vars["time"].values)
                    .astype("int64") // 10**9
                )
                outo_times_s = np.array(times.astype("int64") // 10**9)

                def _rmse(arr):
                    a = arr[np.isfinite(arr)]
                    return float(np.sqrt(np.mean(a ** 2))) if len(a) else float("nan")

                def _bias_stat(arr):
                    a = arr[np.isfinite(arr)]
                    return float(np.mean(a)) if len(a) else float("nan")

                rmse_all     = _rmse(full_bias_arr[valid])
                bias_all_val = _bias_stat(full_bias_arr[valid])

                impact_lines.append("\n## Impact on model statistics")
                impact_lines.append(
                    f"RMSE (all outliers)       = {rmse_all:.3f} m"
                )
                impact_lines.append(
                    f"Bias (all outliers)       = {bias_all_val:+.3f} m"
                )

                def _map_to_cco(local_idx):
                    target_ts = outo_times_s[np.asarray(local_idx)]
                    return np.isin(cco_times_s, target_ts)

                for key, res in pattern_results.items():
                    if not res.get("detected") or not res.get("affected_idx"):
                        continue
                    excl  = _map_to_cco(res["affected_idx"])
                    keep  = valid & ~excl
                    label = PATTERN_REGISTRY[key]["name"][:18].rstrip()
                    impact_lines.append(
                        f"RMSE (excl. {label:<18}) = {_rmse(full_bias_arr[keep]):.3f} m  "
                        f"(Δ = {_rmse(full_bias_arr[keep]) - rmse_all:+.3f} m)"
                    )
                    impact_lines.append(
                        f"Bias (excl. {label:<18}) = {_bias_stat(full_bias_arr[keep]):+.3f} m  "
                        f"(Δ = {_bias_stat(full_bias_arr[keep]) - bias_all_val:+.3f} m)"
                    )
            except Exception as exc:
                logger.warning("Impact statistics computation failed: %s", exc)

        # ------------------------------------------------------------------
        # Build human-readable summary
        # ------------------------------------------------------------------
        lines = [
            "## Numerical Pattern Detection Report",
            f"   Total outliers analysed : {n_total}",
            f"   Patterns checked        : {', '.join(keys_to_run)}",
        ]

        _section_labels = {
            "checkerboard":        "A",
            "garden_sprinkler":    "B",
            "source_term_ringing": "C",
            "hs_collapse":         "D",
            "spinup_insufficient": "E",
        }
        # Assign letters: registered order, falling back to sequential index
        letter_idx = 0
        used_letters = set()

        for key in keys_to_run:
            if key in _section_labels:
                letter = _section_labels[key]
            else:
                while chr(ord("A") + letter_idx) in used_letters:
                    letter_idx += 1
                letter = chr(ord("A") + letter_idx)
                letter_idx += 1
            used_letters.add(letter)

            entry = PATTERN_REGISTRY[key]
            res   = pattern_results[key]
            lines.append(f"\n### {letter}. {entry['name']} (instability #{entry['instability_id']})")

            if "error" in res:
                lines.append(f"   ERROR: {res['error']}")
                continue

            if res["detected"]:
                lines.append("   DETECTED")
                lines.append(
                    f"   contribution_pct : {res['contribution_pct']:.1f} %  "
                    f"({len(res['affected_idx'])} affected point(s))"
                )
                # Pattern-specific extra lines
                if key == "checkerboard":
                    lines.append(
                        f"   parity_fraction  : {res['parity_fraction']:.3f}"
                    )
                    lines.append(
                        f"   amplitude_m      : {res['amplitude_m']:.3f} m"
                    )
                elif key == "garden_sprinkler":
                    lines.append(
                        f"   flagged windows  : {len(res['windows'])}"
                    )
                    for w in res["windows"]:
                        lines.append(
                            f"     [{pd.Timestamp(w['center_time']).strftime('%Y-%m-%d %H:%M')}]"
                            f"  score={w['bimodal_score']:.2f}"
                            f"  low={w['level_low_m']:.2f} m"
                            f"  high={w['level_high_m']:.2f} m"
                            f"  n={w['n_points']}"
                        )
                elif key == "source_term_ringing":
                    lines.append(
                        f"   flagged windows  : {res['n_flagged_windows']}"
                    )
                    lines.append(
                        f"   mean amplitude   : {res['mean_amplitude_m']:.3f} m"
                    )
                elif key == "hs_collapse":
                    lines.append(
                        f"   collapsed points : {res['n_collapsed']}"
                    )
                    lines.append(
                        f"   mean obs_Hs      : {res['mean_obs_hs_m']:.3f} m"
                    )
                elif key == "spinup_insufficient":
                    lines.append(
                        f"   energy_still_growing : {res['energy_still_growing']}"
                    )
                    lines.append(
                        f"   bias_drift_detected  : {res['bias_drift_detected']}"
                    )
                    lines.append(
                        f"   n_early / n_late     : {res['n_early']} / {res['n_late']}"
                    )
                    if np.isfinite(res["early_mean_bias_m"]):
                        lines.append(
                            f"   early_mean_bias      : {res['early_mean_bias_m']:+.3f} m"
                        )
                    if np.isfinite(res["late_mean_bias_m"]):
                        lines.append(
                            f"   late_mean_bias       : {res['late_mean_bias_m']:+.3f} m"
                        )
                    if np.isfinite(res["energy_slope_per_hr"]):
                        lines.append(
                            f"   energy_slope_per_hr  : {res['energy_slope_per_hr']:+.4f} /h"
                        )
                    lines.append(
                        f"   run_start_used       : {res['run_start_used']}"
                    )
            else:
                lines.append("   Not detected.")

        lines.extend(impact_lines)
        summary = "\n".join(lines)
        logger.info("%s", summary)

        report = {**pattern_results, "summary": summary}
        self.pattern_report = report
        return report

    # ------------------------------------------------------------------
    # Fix suggestions
    # ------------------------------------------------------------------

    def suggest_fixes(self, pattern_report: dict | None = None) -> str:
        """
        Print known solutions and key references for every detected pattern.

        Parameters
        ----------
        pattern_report : dict, optional – output of ``detect_numerical_patterns()``.
                          Defaults to ``self.pattern_report``.

        Returns
        -------
        text : str  – the full formatted message (also printed to stdout).
        """
        report = pattern_report if pattern_report is not None else self.pattern_report
        if report is None:
            msg = (
                "No pattern report available. "
                "Run detect_numerical_patterns() first."
            )
            logger.info("%s", msg)
            return msg

        detected_keys = [
            k for k in report
            if k != "summary"
            and isinstance(report[k], dict)
            and report[k].get("detected")
        ]

        if not detected_keys:
            msg = "No numerical patterns were detected — no fixes to suggest."
            logger.info("%s", msg)
            return msg

        lines = ["=" * 70, "  Suggested fixes for detected numerical instabilities",
                 "=" * 70]

        for key in detected_keys:
            if key not in PATTERN_REGISTRY:
                continue
            entry   = PATTERN_REGISTRY[key]
            res     = report[key]
            pct     = res.get("contribution_pct", 0.0)
            n_aff   = len(res.get("affected_idx", []))

            lines.append(f"\n{'─' * 70}")
            lines.append(
                f"  Pattern : {entry['name']}"
                f"  (instability #{entry['instability_id']})"
            )
            lines.append(f"  Impact  : {pct:.1f} % of outliers ({n_aff} point(s))")
            lines.append(f"  What    : {entry['description']}")
            lines.append("")
            lines.append("  SOLUTIONS")
            for i, sol in enumerate(entry["solutions"], start=1):
                wrapped = f"    {i}. {sol}"
                lines.append(wrapped)
            lines.append("")
            lines.append("  REFERENCES")
            for ref in entry["references"]:
                lines.append(f"    • {ref}")

        lines.append(f"\n{'=' * 70}")
        text = "\n".join(lines)
        logger.info("%s", text)
        return text

    def sel_event(self, event_id, events_df):
        """
        Return a copy of this outo with vars restricted to one event.

        Parameters
        ----------
        event_id  : int            – value from events_df["event_id"].
        events_df : pd.DataFrame   – output of group_events().
        """
        row = events_df.set_index("event_id").loc[event_id]
        new = deepcopy(self)
        new.vars = self.vars.isel(time=row["indices"])
        return new

    def plot_event_timeseries(self, event_id, events_df,
                            window_h=3, save_path=None):
        """
        Two-panel time series for one grouped event:
        (top) obs vs. model Hs, (bottom) bias with the threshold band.
        Outlier points belonging to the event are highlighted in red.

        Parameters
        ----------
        event_id  : int            – from events_df["event_id"].
        events_df : pd.DataFrame   – output of group_events().
        window_h  : float          – hours of context on each side of the
                                    event peak (default 3).
        save_path : str/Path, optional
        """
        if self.cco is None:
            raise ValueError("cco is None; full time series requires the parent cco.")

        row  = events_df.set_index("event_id").loc[event_id]
        peak = pd.Timestamp(row["peak_time"])
        t_lo = peak - pd.Timedelta(hours=window_h)
        t_hi = peak + pd.Timedelta(hours=window_h)

        cco_times = pd.DatetimeIndex(self.cco.vars.time.values)
        mask = (cco_times >= t_lo) & (cco_times <= t_hi)
        ds = self.cco.vars.isel(time=np.where(mask)[0])

        obs      = np.asarray(ds[self.obs_var])
        mod      = np.asarray(ds[self.mod_var])
        times    = pd.DatetimeIndex(ds.time.values)
        bias_all = mod - obs

        # match outlier timestamps at second precision to avoid ns drift
        outo_event   = self.sel_event(event_id, events_df)
        out_times_s  = outo_event.vars.time.values.astype("datetime64[s]")
        times_s      = times.values.astype("datetime64[s]")
        out_mask     = np.isin(times_s, out_times_s)

        lo, hi, center = self.stats["lo"], self.stats["hi"], self.stats["center"]

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(12, 6), sharex=True,
            gridspec_kw={"height_ratios": [1.3, 1]},
        )

        ax1.plot(times, obs, color="tab:blue",   linewidth=1.0, label=self.obs_var)
        ax1.plot(times, mod, color="tab:orange", linewidth=1.0, label=self.mod_var)
        if out_mask.any():
            ax1.scatter(times[out_mask], obs[out_mask], color="red",
                        s=35, zorder=5, label="outlier (obs)")
            ax1.scatter(times[out_mask], mod[out_mask], facecolors="none",
                        edgecolors="red", s=35, zorder=5, label="outlier (model)")
        ax1.set_ylabel("Hs [m]")
        ax1.legend(loc="upper right", fontsize=8, ncol=2, framealpha=0.85)
        ax1.set_title(
            f"Event #{event_id} — peak {peak.strftime('%Y-%m-%d %H:%M UTC')}  "
            f"max|bias|={row['max_bias']:.3f} m  ({row['n_points']} outlier point(s))",
            fontsize=10,
        )

        ax2.axhspan(lo, hi, color="seagreen", alpha=0.12, zorder=1,
                    label=f"threshold [{lo:.2f}, {hi:.2f}] m")
        ax2.axhline(center, color="black", linewidth=0.8, zorder=2)
        ax2.plot(times, bias_all, color="grey", linewidth=0.8, zorder=3)
        if out_mask.any():
            ax2.scatter(times[out_mask], bias_all[out_mask], color="red",
                        s=35, zorder=5, label="outlier")
        ax2.set_ylabel("Bias model−obs [m]")
        ax2.set_xlabel("Time")
        ax2.legend(loc="upper right", fontsize=8, framealpha=0.85)

        fig.autofmt_xdate()
        fig.tight_layout()

        if save_path is not None:
            _save_figure(fig, save_path)

        return fig

    def plot_track_over_model(self,
                          event_id=None, events_df=None,
                          idx=None,
                          model_nID=None, model_name=None,
                          track_window_min=30,
                          projection=None, bb=None, margin_deg=5.0,
                          track_color_by="obs",
                          vmin_hs=0, vmax_hs=14, levels_incr=0.5,
                          save_path=None):
        """
        Plot the satellite track over the 2D model Hs field with all outlier
        points in the window highlighted.

        The center time (model snapshot + track window) is determined from
        either a grouped event or a single point index in outo.vars.

        Parameters
        ----------
        event_id         : int, optional   – event ID from group_events();
                                            mutually exclusive with idx.
        events_df        : pd.DataFrame    – required when event_id is given.
        idx              : int, optional   – integer index into outo.vars;
                                            uses that point's time as center.
        model_nID        : str, optional   – falls back to self.model.
        model_name       : str, optional   – falls back to self.cco.name.
        track_window_min : int             – half-window in minutes (default 30).
        projection       : cartopy.crs, optional – default PlateCarree.
        bb               : (lonmin,lonmax,latmin,latmax), optional
        margin_deg       : float           – auto-extent margin (default 5°).
        track_color_by   : "obs" or "bias" – variable used to colour the track
                                            (default "obs").
        vmin_hs/vmax_hs  : float           – model Hs colour scale [m].
        levels_incr      : float           – Hs contour step [m].
        save_path        : str/Path, optional

        Returns
        -------
        fig : matplotlib.Figure or None
        """
        if self.cco is None:
            raise ValueError(
                "cco is None; plot_track_over_model requires the parent cco."
            )
        if self.vars is None:
            logger.info("No outliers — nothing to plot.")
            return None

        # ---- resolve center time ------------------------------------------
        if event_id is not None:
            if events_df is None:
                raise ValueError("events_df is required when event_id is given.")
            row = events_df.set_index("event_id").loc[event_id]
            center_time = pd.Timestamp(row["peak_time"])
        elif idx is not None:
            center_time = pd.Timestamp(self.vars.isel(time=idx)["time"].values)
        else:
            raise ValueError("Provide either event_id (+ events_df) or idx.")

        # ---- model identifiers --------------------------------------------
        if model_nID is None:
            model_nID = self.model
        if model_nID is None:
            raise ValueError("model_nID must be supplied or outo.model must be set.")
        _mc_name = (model_name if model_name is not None
                    else getattr(self.cco, "name", None))

        # ---- satellite track from cco -------------------------------------
        t_lo = center_time - timedelta(minutes=track_window_min)
        t_hi = center_time + timedelta(minutes=track_window_min)
        ds_track = self.cco.vars.sel(time=slice(t_lo, t_hi))

        if len(ds_track.time) == 0:
            logger.info(
                "[plot_track_over_model] No collocated points in "
                "±%d min around %s", track_window_min, center_time
            )
            return None

        track_lons = np.asarray(ds_track["obs_lons"])
        track_lats = np.asarray(ds_track["obs_lats"])
        track_obs  = np.asarray(ds_track[self.obs_var])
        track_mod  = np.asarray(ds_track[self.mod_var])
        track_bias = track_mod - track_obs

        # ---- locate outlier points in this window (from outo.vars) --------
        out_times_s   = self.vars.time.values.astype("datetime64[s]")
        track_times_s = ds_track.time.values.astype("datetime64[s]")
        out_in_window = np.isin(track_times_s, out_times_s)
        out_lons      = track_lons[out_in_window]
        out_lats      = track_lats[out_in_window]

        # ---- map extent ---------------------------------------------------
        if bb is not None:
            lonmin, lonmax, latmin, latmax = bb
        else:
            lonmin, lonmax = _lon_extent(track_lons, margin_deg=margin_deg)
            latmin = max(float(np.min(track_lats)) - margin_deg, -90)
            latmax = min(float(np.max(track_lats)) + margin_deg,  90)

        # ---- load model field ---------------------------------------------
        model_time_h = _round_to_hour(center_time)
        _mc_kwargs = {"nID": model_nID,
                    "sd":  str(model_time_h),
                    "ed":  str(model_time_h)}
        if _mc_name is not None:
            _mc_kwargs["name"] = _mc_name

        logger.info(
            "[plot_track_over_model] Loading '%s'%s at %s …",
            model_nID,
            f" name='{_mc_name}'" if _mc_name else "",
            model_time_h,
        )
        mco = None
        try:
            mco = mc(**_mc_kwargs).populate()
        except Exception as exc:
            logger.warning("Could not load model field: %s", exc)

        # ---- build figure with the requested projection -----------------
        if projection is None:
            projection = ccrs.PlateCarree()

        levels_hs = np.arange(vmin_hs, vmax_hs + levels_incr, levels_incr)
        cmap_hs   = cmocean.cm.amp
        norm_hs   = mpl.colors.BoundaryNorm(levels_hs, cmap_hs.N)

        fig = plt.figure(figsize=(10, 9))
        ax  = fig.add_subplot(1, 1, 1, projection=projection)

        land = cfeature.GSHHSFeature(
            scale="i", levels=[1], facecolor=cfeature.COLORS["land"]
        )
        ax.add_feature(land, facecolor="burlywood", alpha=0.5, zorder=5)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.6, zorder=6)
        try:
            ax.set_extent([lonmin, lonmax, latmin, latmax],
                           crs=ccrs.PlateCarree())
        except Exception as exc:
            logger.warning(f"set_extent failed: {exc}")

        _draw_labels = isinstance(projection, ccrs.PlateCarree)
        gl = ax.gridlines(draw_labels=_draw_labels, linewidth=0.4,
                           color="gray", alpha=0.5, linestyle="--")
        if _draw_labels:
            gl.top_labels   = False
            gl.right_labels = False

        # ---- model 2D Hs contourf -----------------------------------------
        if mco is not None:
            try:
                Mlons = mco.vars["lons"].values
                Mlats = mco.vars["lats"].values
                if Mlons.ndim == 1:
                    Mlons, Mlats = np.meshgrid(Mlons, Mlats)
                Mhs_da = mco.vars["Hs"]
                if "time" in Mhs_da.dims:
                    Mhs = Mhs_da.sel(time=model_time_h, method="nearest")
                else:
                    Mhs = Mhs_da
                cf = ax.contourf(
                    Mlons, Mlats, Mhs.values,
                    levels=levels_hs, cmap=cmap_hs, norm=norm_hs,
                    transform=ccrs.PlateCarree(), extend="max",
                )
                fig.colorbar(cf, ax=ax, label="Hs model [m]",
                             fraction=0.033, pad=0.04,
                             ticks=levels_hs[::4])
            except Exception as exc:
                logger.warning("Could not plot model field: %s", exc)

        # ---- satellite track (coloured by track_color_by) -----------------
        if track_color_by == "bias":
            color_vals   = track_bias
            track_cmap   = mpl.cm.RdBu_r
            bias_abs_max = max(float(np.abs(track_bias).max()), 0.5)
            track_norm   = mpl.colors.TwoSlopeNorm(
                vmin=-bias_abs_max, vcenter=0.0, vmax=bias_abs_max
            )
            track_label = "Bias model\u2212obs [m]"
        else:  # "obs"
            color_vals  = track_obs
            track_cmap  = cmocean.cm.amp
            track_norm  = norm_hs
            track_label = self.obs_var.replace("obs_", "Sat obs ") + " [m]"

        lc = _draw_colored_track(
            ax, track_lons, track_lats, color_vals,
            cmap=track_cmap, norm=track_norm,
            projection=projection,
            linewidth=2.0, zorder=12,
        )
        if lc is not None:
            fig.colorbar(lc, ax=ax, label=track_label,
                         fraction=0.033, pad=0.09, shrink=0.8)

        # ---- outlier highlights -------------------------------------------
        if out_lons.size > 0:
            ax.scatter(
                out_lons, out_lats, c="red", s=120, marker="*",
                edgecolors="black", linewidths=0.8,
                transform=ccrs.PlateCarree(), zorder=15,
                label=f"outlier ({out_lons.size})",
            )
            ax.legend(loc="lower right", fontsize=9, framealpha=0.85)

        label = (f"event #{event_id}" if event_id is not None
                 else f"point #{idx}")
        ax.set_title(
            f"Track over model \u2014 {model_nID}  [{label}]\n"
            f"{center_time.strftime('%Y-%m-%d %H:%M UTC')}  "
            f"\u00b1{track_window_min} min  |  {out_lons.size} outlier(s) in window",
            fontsize=9,
        )
        fig.tight_layout()

        if save_path is not None:
            _save_figure(fig, save_path)

        return fig
