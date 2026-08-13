from types import SimpleNamespace
import pytest
import xarray as xr

from wavy.satellite_module import satellite_class as sc
from wavy.model_module import model_class as mc
from wavy.collocation_module import collocation_class as cc
from wavy import collocation_module
from wavy.insitu_module import insitu_class as ic
from wavy.insitu_module import poi_class as pc

from wavy.errors import (
    CollocationInputError,
    CollocationBuildError,
    CollocationRunError,
    ModelFileSearchError,
)

# include possibility for collocating different variable
# varalias = 'Hs', 'U', aso...


def test_sat_collocation_and_validation(test_data, tmpdir):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = "s3a"
    varalias = "Hs"
    twin = 30
    nID = "cmems_L3_NRT"
    model = "ww3_4km"
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(reader="read_local_ncfiles", path=str(test_data / "L3/s3a"))
    # crop to region
    sco = sco.crop_to_region(model)

    # collocate
    cco = cc(oco=sco, model=model, leadtime="best", distlim=6).populate()
    assert len(vars(cco).keys()) == 22
    assert len(cco.vars.keys()) == 10

    # validate


def test_cco_multivar(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = "s3a"
    varalias = "Hs"
    twin = 30
    nID = "cmems_L3_NRT"
    model = "ww3_4km"
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, name=name, varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(reader="read_local_ncfiles", path=str(test_data / "L3/s3a"))
    # crop to region
    sco = sco.crop_to_region(model)

    # collocate
    cco = cc(
        oco=sco, model=model, leadtime="best", distlim=6, varalias=["Hs", "Tm01"]
    ).populate()
    assert len(vars(cco).keys()) == 22
    assert len(cco.vars.keys()) == 11


def test_insitu_collocation_and_validation(test_data, tmpdir):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    varalias = "Hs"
    twin = 30
    model = "ww3_4km"
    nID = "D_Breisundet_wave"
    name = "wavescan"

    # init insitu_object and check for polygon region
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name, twin=twin)

    # read data
    ico = ico.populate()

    # collocate
    cco = cc(oco=ico, model=model, leadtime="best", distlim=6).populate()
    assert len(vars(cco).keys()) == 22
    assert len(cco.vars.keys()) == 10

    # validate


def test_insitu_collocation_leadtime(test_data, tmpdir):
    sd = "2024-01-01 10"
    ed = "2024-01-01 19"
    varalias = "Hs"
    twin = 30
    model = "ww3_4km"
    nID = "D_Breisundet_wave"
    name = "wavescan"

    # init insitu_object and check for polygon region
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name, twin=twin)

    # read data
    ico = ico.populate()

    # collocate
    cco = cc(oco=ico, model=model, leadtime=10, twin=9).populate()
    assert len(vars(cco).keys()) == 22
    assert len(cco.vars.keys()) == 10
    assert len(cco.vars.time) == 2


def test_poi_collocation():
    # define poi dictionary for track
    dt = ["2023-7-1", "2023-7-2", "2023-7-3"]
    lats = [56.5, 59.3, 64.3]
    lons = [3.5, 1.8, 4.2]
    poi_dict = {"time": dt, "lons": lons, "lats": lats}

    # init poi_class
    pco = pc(poi_dict)

    # collocate
    cco = cc(oco=pco, model="ww3_4km", leadtime="best").populate()
    assert len(vars(cco).keys()) == 22
    assert len(cco.vars.keys()) == 10


#    # write to nc
#    cco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
#    # test validation
#    cco.validate_collocated_values()
#
# def test_insitu_collocation_and_validation():
#    sd = "2021-8-2 01"
#    ed = "2021-8-2 03"
#    nID = 'D_Breisundet_wave'
#    sensor = 'wavescan'
#    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,sensor=sensor)
#    # collocate
#    cco = cc(model='mwam4',obs_obj_in=ico,distlim=6,
#             leadtime='best',date_incr=1)
#    # test validation
#    cco.validate_collocated_values()


def test_collocate_observations(test_data):
    from wavy.collocation_module import collocate_observations

    sd = "2023-07-04"
    ed = "2023-07-05"
    ico = ic(sd=sd, ed=ed, nID="MO_Draugen_monthly", name="Draugen").populate(
        path=str(test_data / "insitu/monthly/Draugen/")
    )
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a").populate(
        path=str(test_data / "L3/s3a")
    )

    ico_colloc, sco_colloc = collocate_observations(ico, sco)
    print(ico_colloc)
    print(len(ico_colloc.vars.keys()))
    assert len(ico_colloc.vars.keys()) == 3
    assert len(ico_colloc.vars.time.values) > 0
    print(sco_colloc)
    print(len(sco_colloc.vars.keys()))
    assert len(sco_colloc.vars.keys()) == 4
    assert len(sco_colloc.vars.time.values) > 0


# ----------------------------------------------------------------------- #
# --- new error-handling tests (wavy.errors) ----------------------------- #
# ----------------------------------------------------------------------- #


def test_collocate_no_observations_raises():
    """
    collocate() should raise CollocationInputError when oco is None.

    Regression test: the original check used `and` where it needed
    `or` - since `self.oco is None` is True, short-circuit evaluation
    of `and` still evaluated `self.oco.vars[...]`, crashing with
    AttributeError on the exact case this check exists to catch.
    Calling collocate() directly on a minimal fake object isolates
    this precondition logic from the rest of the pipeline.
    """
    fake_self = SimpleNamespace(oco=None, model="ww3_4km", method="closest")

    with pytest.raises(CollocationInputError, match="no observation values"):
        cc.collocate(fake_self)


def test_collocate_no_model_raises():
    """
    collocate() should raise CollocationInputError when no model is
    specified, distinct from the 'no observations' case.
    """
    fake_oco = SimpleNamespace(
        varalias=["Hs"],
        stdvarname=["Hs"],
        vars=xr.Dataset({"Hs": ("time", [1.0, 2.0])}),
    )
    fake_self = SimpleNamespace(
        oco=fake_oco,
        varalias_obs=["Hs"],
        model=None,
        method="closest",
    )

    with pytest.raises(CollocationInputError, match="no model specified"):
        cc.collocate(fake_self)


def test_collocate_unknown_method_raises():
    """
    collocate() should raise CollocationInputError for an unrecognized
    'method', rather than leaving results_dict undefined and raising
    an opaque UnboundLocalError at `return results_dict`.
    """
    fake_oco = SimpleNamespace(
        varalias=["Hs"],
        stdvarname=["Hs"],
        vars=xr.Dataset({"Hs": ("time", [1.0, 2.0])}),
    )
    fake_self = SimpleNamespace(
        oco=fake_oco,
        varalias_obs=["Hs"],
        model="ww3_4km",
        method="not_a_real_method",
    )

    with pytest.raises(CollocationInputError, match="Unknown collocation method"):
        cc.collocate(fake_self)


def test_populate_propagates_collocation_input_error(test_data):
    """
    populate() should propagate CollocationInputError unchanged
    (not re-wrap it as CollocationRunError) when a precondition fails
    - here, model=None via the real constructor/populate() path.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)
    sco = sco.populate(reader="read_local_ncfiles", path=str(test_data / "L3/s3a"))
    sco = sco.crop_to_region("ww3_4km")

    cco = cc(oco=sco, model=None, leadtime="best", distlim=6)

    with pytest.raises(CollocationInputError, match="no model specified"):
        cco.populate()


def test_build_xr_dataset_missing_key_raises():
    """
    _build_xr_dataset() should raise CollocationBuildError (not an
    opaque NameError/UnboundLocalError from referencing an
    unassigned 'ds') when results_dict is missing a required key.
    """
    fake_self = SimpleNamespace(
        varalias_obs=["Hs"], varalias_mod=["Hs"], nID="fake_nID", model="ww3_4km"
    )

    # deliberately incomplete - missing 'dist', 'obs_lons', etc.
    results_dict = {"obs_time": [1, 2]}

    with pytest.raises(
        CollocationBuildError, match="Could not assemble the collocated"
    ):
        cc._build_xr_dataset(fake_self, results_dict)


def test_populate_wraps_unexpected_error(monkeypatch, test_data):
    """
    populate() should wrap a genuinely unexpected failure as
    CollocationRunError, chaining the original exception, rather than
    silently setting new.vars = None / new.error = e and returning a
    seemingly-successful object.
    """
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    sco = sc(sd=sd, ed=ed, nID="cmems_L3_NRT", name="s3a", varalias="Hs", twin=30)
    sco = sco.populate(reader="read_local_ncfiles", path=str(test_data / "L3/s3a"))
    sco = sco.crop_to_region("ww3_4km")

    cco = cc(oco=sco, model="ww3_4km", leadtime="best", distlim=6)

    def broken_collocate(self, **kwargs):
        raise ValueError("simulated unexpected collocation failure")

    monkeypatch.setattr(cc, "collocate", broken_collocate)

    with pytest.raises(CollocationRunError) as excinfo:
        cco.populate()

    assert isinstance(excinfo.value.__cause__, ValueError)
    assert not hasattr(cco, "error")


def test_get_model_filename_returns_none_on_search_exhausted(monkeypatch):
    """
    Regression test for the interaction with the model_module.py
    fix: _make_model_filename_wrapper() now raises
    ModelFileSearchError instead of looping forever once its search
    is exhausted. get_model_filename() must translate that back to
    None, since find_valid_fc_dates_for_model_and_leadtime() relies
    on a None return to filter out dates with no available file.
    """

    def broken_wrapper(self, fc_date, leadtime, **kwargs):
        raise ModelFileSearchError("simulated: search exhausted")

    monkeypatch.setattr(mc, "_make_model_filename_wrapper", broken_wrapper)

    result = collocation_module.get_model_filename("ww3_4km", "2022-2-1 12", "best")
    assert result is None


def test_find_valid_fc_dates_raises_when_none_valid(monkeypatch):
    """
    find_valid_fc_dates_for_model_and_leadtime() filters out
    individual dates with no available file (via get_model_filename()
    translating ModelFileSearchError to None - see the previous
    test), but if that leaves zero valid dates at all, it raises
    ModelFileSearchError itself rather than silently returning an
    empty list - consistent with the "fail loud on total failure"
    pattern used elsewhere (e.g. _collocate_track()'s all-dates-failed
    check).
    """

    def broken_wrapper(self, fc_date, leadtime, **kwargs):
        raise ModelFileSearchError("simulated: search exhausted")

    monkeypatch.setattr(mc, "_make_model_filename_wrapper", broken_wrapper)

    fc_dates = ["2022-2-1 12", "2022-2-1 13"]
    with pytest.raises(ModelFileSearchError, match="No model files found"):
        collocation_module.find_valid_fc_dates_for_model_and_leadtime(
            fc_dates, "ww3_4km", "best", "nearest"
        )


def test_collocate_track_all_dates_failed_raises(monkeypatch):
    """
    _collocate_track() (via populate()) should raise
    CollocationRunError when every candidate forecast date fails to
    collocate, rather than silently returning a successful-looking
    but entirely empty dataset.
    """
    sd = "2024-01-01 10"
    ed = "2024-01-01 19"
    varalias = "Hs"
    twin = 30
    model = "ww3_4km"
    nID = "D_Breisundet_wave"
    name = "wavescan"

    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name, twin=twin)
    ico = ico.populate()

    cco = cc(oco=ico, model=model, leadtime=10, twin=9)

    def broken_collocate_field(self, mco, tmp_dict, **kwargs):
        raise ValueError("simulated: no collocation possible")

    monkeypatch.setattr(cc, "_collocate_field", broken_collocate_field)

    with pytest.raises(CollocationRunError, match="could be collocated"):
        cco.populate()
