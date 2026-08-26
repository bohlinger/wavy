from wavy import model_module
from wavy.model_module import model_class as mc
from wavy.errors import (
    ModelFileNotFoundError,
    ModelFileSearchError,
    ModelPathTemplateError,
    ModelProcessingError,
    ModelReadError,
)
import pytest


def test_model_class_init():
    # get_model
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-1 01")
    assert mco.__class__.__name__ == "model_class"


def test_ww3_4km_reader():
    # get_model
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-1 01")
    assert mco.__class__.__name__ == "model_class"
    mco.populate()
    # print(mco.vars)
    assert len(vars(mco).keys()) == 19
    assert len(mco.vars.keys()) == 3


def test_dummy_reader():
    """
    Test that the dummy model_class raises a ModelFileSearchError when populate is called.
    Treats the case where the model class is configured or initialized with the wrong src_tmplt/fl_tmplt, which will cause the model_class to fail to find a model file.
    """
    mco = mc(nID="dummy_model", sd="2023-6-1", ed="2023-6-1 01")
    assert mco.__class__.__name__ == "model_class"

    with pytest.raises(
        ModelFileSearchError, match="Reached maximum number of attempts \\(2\\)"
    ):
        mco.populate(max_iter=2)

    assert not hasattr(mco, "vars")


def test_populate_no_files_found(monkeypatch):
    """
    Test that populate() raises ModelFileNotFoundError when the file
    search completes (no exception from _make_model_filename_wrapper)
    but simply finds nothing - e.g. leadtime resolves to None for
    every fc_date, so list_input_files() legitimately returns an
    empty list.

    We bypass the leadtime search machinery entirely by monkeypatching
    list_input_files() directly, so this test is isolated from the
    ModelFileSearchError path and only exercises the "empty pathlst"
    branch of populate().
    """
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-1 01")

    monkeypatch.setattr(mco, "list_input_files", lambda **kwargs: [])

    with pytest.raises(ModelFileNotFoundError, match="No accessible model files"):
        mco.populate()

    assert not hasattr(mco, "vars")


def test_populate_reader_error(monkeypatch):
    """
    Test that populate() wraps a reader failure in ModelReadError,
    preserving the original exception via `raise ... from e`.

    We let the file search succeed normally (real ww3_4km data is
    available), but monkeypatch _get_model() - the method that
    actually invokes the reader - to simulate the reader itself
    blowing up (e.g. a corrupt file, unexpected variable layout).
    """
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-1 01")

    def broken_get_model(self, **kwargs):
        raise ValueError("simulated reader failure")

    monkeypatch.setattr(model_module.model_class, "_get_model", broken_get_model)

    with pytest.raises(ModelReadError) as excinfo:
        mco.populate()

    # original exception should be chained, not lost
    assert isinstance(excinfo.value.__cause__, ValueError)
    assert "simulated reader failure" in str(excinfo.value.__cause__)
    assert not hasattr(mco, "vars")


def test_populate_processing_error(monkeypatch):
    """
    Test that populate() wraps a post-processing failure (variable
    renaming / CF standard names / convention / longitude formatting)
    in ModelProcessingError, distinct from a reader failure.

    The reader itself is left untouched (real data is read
    successfully); we monkeypatch the first post-processing step to
    fail instead.
    """
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-1 01")

    def broken_change_varname(self, **kwargs):
        raise KeyError("simulated rename failure")

    monkeypatch.setattr(
        model_module.model_class, "_change_varname_to_aliases", broken_change_varname
    )

    with pytest.raises(ModelProcessingError) as excinfo:
        mco.populate()

    assert isinstance(excinfo.value.__cause__, KeyError)
    # assert not hasattr(mco, "vars") ask if this should return the incomplete vars dict or not. I think it should not, but we can discuss.


def test_get_files_bad_path_template(monkeypatch):
    """
    Test that a deterministic config problem (e.g. malformed
    src_tmplt/strsub) raises ModelPathTemplateError immediately,
    rather than retrying identically for every date in sd..ed and
    silently returning an empty file list.
    """
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-2")

    def broken_make_subdict(*args, **kwargs):
        raise TypeError("simulated malformed strsub config")

    monkeypatch.setattr(model_module, "make_subdict", broken_make_subdict)

    with pytest.raises(ModelPathTemplateError, match="Could not build a local path"):
        mco._get_files(dict_for_sub=vars(mco), path=None, wavy_path=None)


@pytest.mark.need_credentials
def test_ww3_unstr_reader():
    # get_model

    bb = (5.8, 6.61, 62.3, 63.1)
    res = (0.01, 0.01)  # lon/lat

    mco = mc(nID="ww3_unstr", sd="2019-3-24 10", ed="2019-3-24 10")
    assert mco.__class__.__name__ == "model_class"
    mco.populate(res=res, bb=bb, interp="nearest")
    # print(mco.vars)

    assert len(vars(mco).keys()) == 19
    assert len(mco.vars.keys()) == 3


def test_NORA3_hc_waves():
    # get_model
    mco = mc(nID="NORA3_hc_waves", sd="2019-1-1", ed="2019-1-1")
    assert mco.__class__.__name__ == "model_class"
    mco.populate()
    # print(mco.vars)
    assert len(vars(mco).keys()) == 19
    assert len(mco.vars.keys()) == 3


def test_mco_multivar():
    # get_model
    mco = mc(nID="ww3_4km", sd="2023-6-1", ed="2023-6-1 00", varalias=["Hs", "U"])
    assert mco.__class__.__name__ == "model_class"
    mco.populate()
    # print(mco.vars)
    assert len(vars(mco).keys()) == 19
    assert len(mco.vars.keys()) == 4


# Fails
# def test_MY_L4_thredds():
#    """
#    Just to check when thredds service stops
#    Test for aggregated reader
#    """
#    #get_model
#    mco = mc(nID='cmems_MY_L4', sd="2021-11-16", ed="2021-11-16")
#    assert mco.__class__.__name__ == 'model_class'
#    mco.populate()
#    print(mco.vars)
#    assert len(vars(mco).keys()) == 18
#    assert len(mco.vars.keys()) == 3
