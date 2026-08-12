"""
Tests for wavy.quicklookmod.

Coverage:
  - Pure utility functions: _check_projection, _set_lonlat_minmax,
    _set_polar_extent
  - quicklook_class_sat._check_varalias (via a minimal stub)
  - quicklook_class_sat._PLOT_REGISTRY structure validation
  - quicklook_class_sat._build_plot_context (via populated satellite object)
  - needs_colloc guard in quicklook()
  - Integration: plot_timeseries and plot_map via satellite_class.quicklook()
    using the local L3 test data
"""

import matplotlib

matplotlib.use("Agg")  # non-interactive backend; must come before pyplot import
import matplotlib.pyplot as plt
import numpy as np
import pytest
import cartopy.crs as ccrs
from unittest.mock import MagicMock, patch

from wavy.quicklookmod import (
    PlotContext,
    _check_projection,
    _set_lonlat_minmax,
    _set_polar_extent,
    quicklook_class_sat,
)
from wavy.satellite_module import satellite_class as sc

# ------------------------------------------------------------------ #
# Helpers
# ------------------------------------------------------------------ #


def _make_quicklook_stub(varalias, units):
    """Return a bare quicklook_class_sat instance with only the
    attributes needed by _check_varalias."""
    obj = quicklook_class_sat.__new__(quicklook_class_sat)
    obj.varalias = varalias
    obj.units = units
    return obj


# ------------------------------------------------------------------ #
# _check_projection
# ------------------------------------------------------------------ #


def test_check_projection_none_returns_platecarree():
    proj = _check_projection(None)
    assert isinstance(proj, ccrs.PlateCarree)


def test_check_projection_specified_projection_returned():
    mercator = ccrs.Mercator()
    result = _check_projection(mercator)
    assert result is mercator


def test_check_projection_north_polar_stereo_returned():
    nps = ccrs.NorthPolarStereo()
    result = _check_projection(nps)
    assert result is nps


# ------------------------------------------------------------------ #
# _set_lonlat_minmax
# ------------------------------------------------------------------ #


def test_set_lonlat_minmax_defaults_to_array_extremes():
    lons = np.array([10.0, 20.0, 30.0])
    lats = np.array([50.0, 60.0, 70.0])
    lonmax, lonmin, latmax, latmin = _set_lonlat_minmax(lons, lats)
    assert lonmax == pytest.approx(30.0)
    assert lonmin == pytest.approx(10.0)
    assert latmax == pytest.approx(70.0)
    assert latmin == pytest.approx(50.0)


def test_set_lonlat_minmax_kwargs_override_all():
    lons = np.array([10.0, 20.0, 30.0])
    lats = np.array([50.0, 60.0, 70.0])
    lonmax, lonmin, latmax, latmin = _set_lonlat_minmax(
        lons, lats, lonmax=40.0, lonmin=0.0, latmax=80.0, latmin=40.0
    )
    assert lonmax == pytest.approx(40.0)
    assert lonmin == pytest.approx(0.0)
    assert latmax == pytest.approx(80.0)
    assert latmin == pytest.approx(40.0)


def test_set_lonlat_minmax_partial_kwarg_override():
    lons = np.array([10.0, 20.0, 30.0])
    lats = np.array([50.0, 60.0, 70.0])
    lonmax, lonmin, latmax, latmin = _set_lonlat_minmax(lons, lats, lonmax=40.0)
    assert lonmax == pytest.approx(40.0)
    assert lonmin == pytest.approx(10.0)  # still from array


def test_set_lonlat_minmax_2d_arrays():
    lons = np.array([[0.0, 10.0], [20.0, 30.0]])
    lats = np.array([[40.0, 50.0], [60.0, 70.0]])
    lonmax, lonmin, latmax, latmin = _set_lonlat_minmax(lons, lats)
    assert lonmax == pytest.approx(30.0)
    assert lonmin == pytest.approx(0.0)
    assert latmax == pytest.approx(70.0)
    assert latmin == pytest.approx(40.0)


# ------------------------------------------------------------------ #
# _set_polar_extent
# ------------------------------------------------------------------ #


def test_set_polar_extent_northern_hemisphere_latmax_is_90():
    proj = ccrs.NorthPolarStereo()
    extent = _set_polar_extent(None, 30.0, -30.0, 80.0, 60.0, proj)
    assert extent[3] == 90


def test_set_polar_extent_northern_hemisphere_latmin_clamped():
    proj = ccrs.NorthPolarStereo()
    extent = _set_polar_extent(None, 30.0, -30.0, 80.0, 60.0, proj)
    assert extent[2] >= 30  # clamped to max(-180 is lon, lat >= 30)


def test_set_polar_extent_southern_hemisphere_latmin_is_minus_90():
    proj = ccrs.SouthPolarStereo()
    extent = _set_polar_extent(None, 30.0, -30.0, -60.0, -80.0, proj)
    assert extent[2] == -90


def test_set_polar_extent_southern_hemisphere_latmax_clamped():
    proj = ccrs.SouthPolarStereo()
    extent = _set_polar_extent(None, 30.0, -30.0, -60.0, -80.0, proj)
    assert extent[3] <= -30  # clamped


def test_set_polar_extent_returns_four_elements():
    proj = ccrs.NorthPolarStereo()
    extent = _set_polar_extent(None, 30.0, -30.0, 80.0, 60.0, proj)
    assert len(extent) == 4


def test_set_polar_extent_lon_clamped_to_pm180():
    proj = ccrs.NorthPolarStereo()
    extent = _set_polar_extent(None, 175.0, -175.0, 80.0, 60.0, proj)
    assert extent[0] >= -180
    assert extent[1] <= 180


def test_set_polar_extent_mid_latitude_stereographic():
    """Stereographic with mid-latitude centre falls through to else branch."""
    proj = ccrs.Stereographic(central_latitude=45.0)
    extent = _set_polar_extent(None, 30.0, -30.0, 60.0, 40.0, proj)
    assert len(extent) == 4
    # no clamping to ±90 in the else branch
    assert extent[2] < extent[3]


# ------------------------------------------------------------------ #
# quicklook_class_sat._check_varalias
# ------------------------------------------------------------------ #


def test_check_varalias_string_returned_directly():
    stub = _make_quicklook_stub("Hs", "m")
    va, units = stub._check_varalias()
    assert va == "Hs"
    assert units == "m"


def test_check_varalias_list_defaults_to_first_element():
    stub = _make_quicklook_stub(["Hs", "Tm01"], ["m", "s"])
    va, units = stub._check_varalias()
    assert va == "Hs"
    assert units == "m"


def test_check_varalias_list_kwarg_selects_second_element():
    stub = _make_quicklook_stub(["Hs", "Tm01"], ["m", "s"])
    va, units = stub._check_varalias(varalias="Tm01")
    assert va == "Tm01"
    assert units == "s"


def test_check_varalias_list_invalid_kwarg_raises():
    stub = _make_quicklook_stub(["Hs", "Tm01"], ["m", "s"])
    with pytest.raises(AssertionError):
        stub._check_varalias(varalias="U10")


def test_check_varalias_list_non_string_kwarg_raises():
    stub = _make_quicklook_stub(["Hs", "Tm01"], ["m", "s"])
    with pytest.raises(AssertionError):
        stub._check_varalias(varalias=["Hs"])


# ------------------------------------------------------------------ #
# _PLOT_REGISTRY structure
# ------------------------------------------------------------------ #

_REQUIRED_SPEC_KEYS = {"method", "default", "needs_colloc"}


def test_plot_registry_is_nonempty():
    assert len(quicklook_class_sat._PLOT_REGISTRY) > 0


def test_plot_registry_all_entries_have_required_keys():
    for kwarg, spec in quicklook_class_sat._PLOT_REGISTRY.items():
        assert _REQUIRED_SPEC_KEYS <= set(spec.keys()), (
            f"Registry entry {kwarg!r} is missing keys: "
            f"{_REQUIRED_SPEC_KEYS - set(spec.keys())}"
        )


def test_plot_registry_methods_exist_on_class():
    for kwarg, spec in quicklook_class_sat._PLOT_REGISTRY.items():
        assert hasattr(
            quicklook_class_sat, spec["method"]
        ), f"Registry entry {kwarg!r} references missing method {spec['method']!r}"


def test_plot_registry_default_and_needs_colloc_are_bool():
    for kwarg, spec in quicklook_class_sat._PLOT_REGISTRY.items():
        assert isinstance(
            spec["default"], bool
        ), f"Registry entry {kwarg!r}: 'default' must be bool"
        assert isinstance(
            spec["needs_colloc"], bool
        ), f"Registry entry {kwarg!r}: 'needs_colloc' must be bool"


# ------------------------------------------------------------------ #
# PlotContext dataclass
# ------------------------------------------------------------------ #


def test_plot_context_instantiation():
    ctx = PlotContext(
        varalias="Hs",
        units="m",
        plot_var=np.array([1.0, 2.0]),
        plot_lons=np.array([0.0, 1.0]),
        plot_lats=np.array([60.0, 61.0]),
        plot_var_obs=None,
        plot_var_model=None,
        fs=12,
        cmap=None,
        projection=None,
    )
    assert ctx.varalias == "Hs"
    assert ctx.mode == "comb"  # default


def test_plot_context_mode_override():
    ctx = PlotContext(
        varalias="Hs",
        units="m",
        plot_var=None,
        plot_lons=None,
        plot_lats=None,
        plot_var_obs=None,
        plot_var_model=None,
        fs=12,
        cmap=None,
        projection=None,
        mode="indiv",
    )
    assert ctx.mode == "indiv"


# ------------------------------------------------------------------ #
# Integration tests – satellite_class.quicklook()
# ------------------------------------------------------------------ #


@pytest.fixture
def populated_sco(test_data):
    """satellite_class object populated from local L3 test data."""
    sco = sc(
        sd="2022-2-1 12",
        ed="2022-2-1 12",
        nID="cmems_L3_NRT",
        name="s3a",
        varalias="Hs",
        twin=30,
    )
    return sco.populate(path=str(test_data / "L3/s3a"))


def test_quicklook_timeseries_returns_figure_and_axes(populated_sco):
    fig, ax = populated_sco.quicklook(ts=True, show=False)
    assert fig is not None
    assert ax is not None
    plt.close(fig)


def test_quicklook_timeseries_ylabel_contains_varalias(populated_sco):
    fig, ax = populated_sco.quicklook(ts=True, show=False)
    ylabel = ax.get_ylabel()
    assert "Hs" in ylabel
    plt.close(fig)


def test_build_plot_context_returns_plot_context(populated_sco):
    ctx = populated_sco._build_plot_context()
    assert isinstance(ctx, PlotContext)
    assert ctx.varalias == "Hs"
    assert ctx.plot_var is not None
    assert ctx.plot_lons is not None
    assert ctx.plot_lats is not None


def test_build_plot_context_no_colloc_data_has_none_obs_model(populated_sco):
    ctx = populated_sco._build_plot_context()
    # Satellite-only data: no obs_/model_ prefix → obs and model are None
    assert ctx.plot_var_obs is None
    assert ctx.plot_var_model is None


def test_quicklook_needs_colloc_guard_skips_without_data(populated_sco):
    """sc=True on non-collocated data should log a warning, not raise."""
    fig, ax = populated_sco.quicklook(sc=True, show=False)
    # No collocated data → guard skips; quicklook returns None when no plot made
    # (or returns the last valid (fig, ax) from a non-colloc plot if also requested)
    # Here only sc=True is set so nothing is plotted.
    assert fig is None and ax is None


def _mock_land():
    """Return a MagicMock that stands in for cfeature.GSHHSFeature.

    GSHHSFeature downloads data from the internet; the mock keeps tests
    fully offline by returning an empty geometry iterable.
    """
    mock = MagicMock()
    mock.intersecting_geometries.return_value = []
    return mock


def test_quicklook_map_returns_figure_and_axes(populated_sco):
    with patch("wavy.quicklookmod.cfeature.GSHHSFeature", return_value=_mock_land()):
        fig, ax = populated_sco.quicklook(m=True, show=False)
    assert fig is not None
    assert ax is not None
    plt.close(fig)


def test_quicklook_map_custom_projection(populated_sco):
    proj = ccrs.Stereographic(central_latitude=90, central_longitude=0)
    with patch("wavy.quicklookmod.cfeature.GSHHSFeature", return_value=_mock_land()):
        fig, ax = populated_sco.quicklook(m=True, show=False, projection=proj)
    assert fig is not None
    plt.close(fig)
