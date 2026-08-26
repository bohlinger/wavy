import sys
import os
from datetime import datetime
import pytest
import numpy as np
import xarray as xr
from types import SimpleNamespace

from wavy.wconfig import load_or_default
from wavy.satellite_module import satellite_class as sc
from wavy import satellite_module
from wavy.insitu_module import poi_class as pc
from wavy.model_module import model_class as mc
from wavy.errors import (
    SatelliteFileNotFoundError,
    SatelliteReadError,
    SatelliteProcessingError,
    SatelliteVariableError,
    SatellitePathTemplateError,
    RegionNotDefinedError,
    GridRetrievalError,
)


@pytest.mark.need_credentials
def test_collectors_cmems_L3(tmpdir):
    sco = sc(sd="2023-2-1 12", ed="2023-2-1 12", nID="cmems_L3_NRT", name="s3a")
    sco.download(path=tmpdir, nproc=4)
    # check if files were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist)) if ".nc" in filelist[i]]
    assert len(nclist) >= 1


# @pytest.mark.need_credentials
# def test_collectors_cci_v3_20Hz(tmpdir):
#    sco = sc(sd='2020-2-1 12', ed='2020-2-1 12',
#             nID='L2_20Hz_s3a', name='s3a')
#    sco.download(path=tmpdir, nproc=8)
#    # check if files were download to tmp directory
#    filelist = os.listdir(tmpdir)
#    nclist = [i for i in range(len(filelist))
#              if '.nc' in filelist[i]]
#    assert len(nclist) >= 1

# @pytest.mark.need_credentials
# def test_collectors_cci_v1_01Hz(tmpdir):
#    sco = sc(sd='2018-1-1', ed='2018-1-1',
#             nID='CCIv1_L3', name='multi')
#    sco.download(path=tmpdir, nproc=8)
#    # check if files were download to tmp directory
#    filelist = os.listdir(tmpdir)
#    nclist = [i for i in range(len(filelist))
#              if '.nc' in filelist[i]]
#    assert len(nclist) >= 1

# @pytest.mark.need_credentials
# def test_collectors_aviso(tmpdir):
#    sco = sc(sd='2020-2-1 12', ed='2020-2-1 12',
#             nID='L2_20Hz_s3a', name='s3a')
#    sco.download(path=tmpdir, nproc=8)
#    # check if files were download to tmp directory
#    filelist = os.listdir(tmpdir)
#    nclist = [i for i in range(len(filelist))
#              if '.nc' in filelist[i]]
#    assert len(nclist) >= 1


def test_manually_specified_reader(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = "s3a"
    varalias = "Hs"
    twin = 30
    nID = "cmems_L3_NRT"
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(reader="read_local_ncfiles", path=str(test_data / "L3/s3a"))
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    print(vlst)
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    print(flst)
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()


def test_default_reader(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = "s3a"
    varalias = "Hs"
    twin = 30
    nID = "cmems_L3_NRT"
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(path=str(test_data / "L3/s3a"))
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()


def test_sco_multivar(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = "s3a"
    varalias = ["Hs", "U"]
    twin = 30
    nID = "cmems_L3_NRT"
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(path=str(test_data / "L3/s3a"))
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()
    assert len(sco.vars["time"]) > 0
    assert len(sco.vars.keys()) == 4
    assert not all(np.isnan(v) for v in sco.vars["Hs"])
    assert not all(np.isnan(v) for v in sco.vars["U"])


def test_polygon_region(test_data):
    sd = "2022-2-01 01"
    ed = "2022-2-03 23"
    name = "s3a"
    varalias = "Hs"
    twin = 30
    nID = "cmems_L3_NRT"
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(path=str(test_data / "L3/s3a"))
    sco = sco.crop_to_region("NordicSeas")
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()


def test_rectangular_region(test_data):
    sd = "2022-2-01 01"
    ed = "2022-2-03 23"
    name = "s3a"
    varalias = "Hs"
    nID = "cmems_L3_NRT"
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias)
    # read data
    sco = sco.populate(path=str(test_data / "L3/s3a"))
    sco = sco.crop_to_region("Sulafj")
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()


def test_direct_input_custom_region(test_data):
    sd = "2022-2-01 01"
    ed = "2022-2-03 23"
    name = "s3a"
    varalias = "Hs"
    nID = "cmems_L3_NRT"
    region_dict = {
        "name": "custom",
        "region": {
            "llcrnrlon": -180.0,
            "llcrnrlat": -90.0,
            "urcrnrlon": 180.0,
            "urcrnrlat": 90.0,
        },
    }

    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias)
    # read data
    sco = sco.populate(path=str(test_data / "L3/s3a"), region=region_dict)
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()
    assert sco.region == "custom"


def test_poi_storm_track(test_data):
    import pandas as pd
    from wavy.utils import parse_date

    # read track
    f = pd.read_csv(test_data / "track/Katrina_track.csv")
    # convert dates to datetime
    dt = [parse_date(d) for d in f.date]
    lons = f.lon.values
    lats = f.lat.values

    # define poi dictionary for track
    poi_dict = {"time": dt, "lons": lons, "lats": lats}

    # initialize poi class object
    pco = pc(poi_dict, nID="Katrina", name="Katrina", varalias="Hs")

    # retrievals
    sco = sc(
        twin=180,
        distlim=200,
        name="multi",
        nID="CCIv1_L3",
        varalias="Hs",  # default
        poi=pco,
    )

    sco = sco.populate(path=str(test_data / "L3/multi"))

    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    print(vlst)
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()


# make test for get closest only


# make test for reading 20Hz
def test_manually_specified_reader_CCIv3_20Hz(test_data):
    sd = "2019-3-24 15"
    ed = "2019-3-24 16"
    name = "s3a"
    varalias = "Hs"
    nID = "L2_20Hz_s3a"

    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias)
    # populate
    sco = sco.populate(path=str(test_data / "CCIv3_20Hz"))
    # adjustments
    sco = sco.crop_to_region("BarentsSea")
    assert sco.__class__.__name__ == "satellite_class"
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if "__" not in n]
    assert len(flst) >= 46
    assert type(sco.vars) == xr.core.dataset.Dataset
    assert not "error" in vars(sco).keys()
    # check if some data was imported
    assert len(sco.vars["time"]) > 0
    # check that not all data is nan
    assert not all(np.isnan(v) for v in sco.vars["time"])
    assert not all(np.isnan(v) for v in sco.vars["Hs"])
    assert not all(np.isnan(v) for v in sco.vars["lons"])
    assert not all(np.isnan(v) for v in sco.vars["lats"])


# def test_write_to_nc(test_data):
# write to nc
# sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
# check if created -> assert
# read nc
# check if varalias assert


def test_populate_no_files_found(tmp_path):
    """
    populate() should raise SatelliteFileNotFoundError - not crash on
    a bare IndexError from self.pathlst[0] - when no candidate files
    exist for the requested period/path.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)

    # tmp_path is guaranteed empty - no .nc files inside
    with pytest.raises(
        SatelliteFileNotFoundError, match="No accessible satellite files"
    ):
        sco.populate(path=str(tmp_path))

    assert not hasattr(sco, "vars")


def test_populate_reader_error(monkeypatch, test_data):
    """
    populate() should wrap a failure inside _get_sat_ts() (the reader
    stage) as SatelliteReadError, chaining the original exception.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)

    def broken_get_sat_ts(self, **kwargs):
        raise ValueError("simulated reader failure")

    monkeypatch.setattr(
        satellite_module.satellite_class, "_get_sat_ts", broken_get_sat_ts
    )

    with pytest.raises(SatelliteReadError) as excinfo:
        sco.populate(path=str(test_data / "L3/s3a"))

    assert isinstance(excinfo.value.__cause__, ValueError)
    assert "simulated reader failure" in str(excinfo.value.__cause__)
    assert not hasattr(sco, "vars")


def test_populate_processing_error(monkeypatch, test_data):
    """
    populate() should wrap a post-processing failure as
    SatelliteProcessingError, distinct from a reader failure.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)

    def broken_change_stdvarname(self, **kwargs):
        raise KeyError("simulated CF standard-name failure")

    monkeypatch.setattr(
        satellite_module.satellite_class,
        "_change_stdvarname_to_cfname",
        broken_change_stdvarname,
    )

    with pytest.raises(SatelliteProcessingError) as excinfo:
        sco.populate(path=str(test_data / "L3/s3a"))

    assert isinstance(excinfo.value.__cause__, KeyError)
    assert not hasattr(sco, "vars")


def test_change_varname_missing_coordinate_raises(monkeypatch):
    """
    _change_varname_to_aliases() should raise SatelliteVariableError
    (not a masked NameError) when a required coordinate ('time',
    'lons', 'lats') cannot be resolved to a source variable name.

    Uses a minimal fake object (not a real populated satellite_class)
    so the failure point is isolated from any real data/config.
    """

    def fake_get_filevarname(item, variable_def, cfgdict, meta, **kwargs):
        if item == "time":
            raise KeyError("no mapping for 'time' in this config")
        return item

    monkeypatch.setattr(satellite_module, "get_filevarname", fake_get_filevarname)

    # minimal stand-in for a satellite_class instance: a dataset with
    # only a data variable ('Hs') and no 'time'/'lons'/'lats' present,
    # so the coordinate-renaming loop is forced to look them up
    fake_self = SimpleNamespace(
        nID="cmems_L3_NRT",  # a real nID so satellite_dict[nID] resolves
        varalias=["Hs"],
        vars=xr.Dataset({"Hs": ("time", [1.0, 2.0])}),
        meta={},
    )

    with pytest.raises(SatelliteVariableError, match="coordinate 'time'"):
        satellite_module.satellite_class._change_varname_to_aliases(fake_self)


def test_change_varname_already_normalized_coordinate_is_skipped(monkeypatch):
    """
    Regression test for a bug found in test_manually_specified_reader_CCIv3_20Hz:
    when a reader (e.g. read_local_20Hz_files, via build_xr_ds_multivar)
    already returns a dataset with 'time' as a coordinate (not a data
    variable), the old check `c in list(new.vars.keys())` failed to
    detect it - Dataset.keys() only reflects data variables, not
    coordinates - so the code went on to resolve a stale source name
    from self.meta (e.g. 'time_echo_sar_ku') and tried to rename a key
    that no longer existed, raising an opaque xarray ValueError.

    With the fix (`c in new.vars.variables`), an already-present
    coordinate must be detected and skipped, even though
    get_filevarname() would resolve to a completely different,
    no-longer-present name.
    """

    def fake_get_filevarname(item, variable_def, cfgdict, meta, **kwargs):
        # simulates stale metadata describing the *original* raw
        # NetCDF name, which the reader has already normalized away
        if item == "time":
            return "time_echo_sar_ku"
        return item

    monkeypatch.setattr(satellite_module, "get_filevarname", fake_get_filevarname)

    ds = xr.Dataset(
        data_vars={"Hs": ("time", [1.0, 2.0])},
        coords={
            "time": (
                "time",
                np.array(["2019-01-01", "2019-01-02"], dtype="datetime64[ns]"),
            ),
            "lats": ("time", [60.0, 61.0]),
            "lons": ("time", [5.0, 6.0]),
        },
    )

    fake_self = SimpleNamespace(
        nID="cmems_L3_NRT",
        varalias=["Hs"],
        vars=ds,
        meta={},
    )

    # should NOT raise - 'time' is already present and must be
    # detected via new.vars.variables, not new.vars.keys()
    result = satellite_module.satellite_class._change_varname_to_aliases(fake_self)
    assert "time" in result.vars.variables


def test_change_varname_stale_source_name_raises(monkeypatch):
    """
    Complementary case: if a required coordinate is genuinely absent
    under both its alias name AND its resolved source name, that's a
    real mismatch and should raise SatelliteVariableError with a
    clear message - not silently pass through to xarray's rename()
    and fail with an opaque ValueError.
    """

    def fake_get_filevarname(item, variable_def, cfgdict, meta, **kwargs):
        if item == "time":
            return "time_echo_sar_ku"  # does not exist anywhere below
        return item

    monkeypatch.setattr(satellite_module, "get_filevarname", fake_get_filevarname)

    # 'time' is neither present as an alias nor under the resolved
    # source name - genuinely unresolvable
    ds = xr.Dataset(
        data_vars={"Hs": ("t", [1.0, 2.0])},
        coords={"lats": ("t", [60.0, 61.0])},
    )

    fake_self = SimpleNamespace(
        nID="cmems_L3_NRT",
        varalias=["Hs"],
        vars=ds,
        meta={},
    )

    with pytest.raises(SatelliteVariableError, match="does not exist in the dataset"):
        satellite_module.satellite_class._change_varname_to_aliases(fake_self)


def test_get_files_bad_path_template(monkeypatch):
    """
    _get_files() should raise SatellitePathTemplateError immediately
    on a deterministic config problem (malformed src_tmplt/strsub),
    rather than repeating the identical failure for every date in
    sd..ed and silently returning no files.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-2 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)

    def broken_make_subdict(*args, **kwargs):
        raise TypeError("simulated malformed strsub config")

    monkeypatch.setattr(satellite_module, "make_subdict", broken_make_subdict)

    with pytest.raises(
        SatellitePathTemplateError, match="Could not build a local path"
    ):
        sco._get_files(dict_for_sub=vars(sco), path=None, wavy_path=None)


def test_region_not_defined(test_data):
    """
    crop_to_region() should raise RegionNotDefinedError for an
    unknown region name, not sys.exit() the whole process.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)
    sco = sco.populate(path=str(test_data / "L3/s3a"))

    with pytest.raises(RegionNotDefinedError, match="not defined"):
        sco.crop_to_region("ThisRegionDoesNotExist")


def test_grid_retrieval_error(monkeypatch, test_data):
    """
    crop_to_region() with a model-grid region should raise
    GridRetrievalError (not an unhandled/cryptic exception) if the
    model grid cannot be retrieved via either the requested date or
    the model's fallback grid_date.
    """

    def broken_wrapper(self, fc_date, leadtime, **kwargs):
        raise ValueError("simulated: no accessible model file for grid")

    monkeypatch.setattr(mc, "_make_model_filename_wrapper", broken_wrapper)

    sd = "2022-2-01 01"
    ed = "2022-2-03 23"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)
    sco = sco.populate(path=str(test_data / "L3/s3a"))

    # 'ww3_4km' is a model-grid region (defined in model_cfg.yaml,
    # not in region_cfg.yaml's poly/geojson sections)
    with pytest.raises(GridRetrievalError, match="Could not retrieve the model grid"):
        sco.crop_to_region("ww3_4km")
