from wavy.model_module import model_class as mc
import pytest

def test_model_class_init():
    #get_model
    mco = mc(nID='ww3_4km', sd="2023-6-1", ed="2023-6-1 01")
    assert mco.__class__.__name__ == 'model_class'

def test_ww3_4km_reader():
    #get_model
    mco = mc(nID='ww3_4km', sd="2023-6-1", ed="2023-6-1 01")
    assert mco.__class__.__name__ == 'model_class'
    mco.populate()
    print(mco.vars)
    assert len(vars(mco).keys()) == 18
    assert len(mco.vars.keys()) == 3

@pytest.mark.need_credentials
def test_ww3_unstr_reader():
    #get_model

    bb = (5.8, 6.61, 62.3, 63.1)
    res = (0.01, 0.01)  # lon/lat

    mco = mc(nID='ww3_unstr', sd="2019-3-24 10", ed="2019-3-24 10")
    assert mco.__class__.__name__ == 'model_class'
    mco.populate(res=res, bb=bb, interp='nearest')
    print(mco.vars)

    assert len(vars(mco).keys()) == 18
    assert len(mco.vars.keys()) == 3

# Fails
#def test_ARCMFC3_model():
#    #get_model
#    mco = mc(nID='ARCMFC3', sd="2024-1-10", ed="2024-1-10")
#    assert mco.__class__.__name__ == 'model_class'
#    mco.populate()
#    print(mco.vars)
#    assert len(vars(mco).keys()) == 18
#    assert len(mco.vars.keys()) == 3

# not testable
#def test_mwam3():
#    #get_model
#    mco = mc(nID='mwam3', sd="2024-1-10", ed="2024-1-10")
#    assert mco.__class__.__name__ == 'model_class'
#    mco.populate()
#    print(mco.vars)
#    assert len(vars(mco).keys()) == 18
#    assert len(mco.vars.keys()) == 3

def test_NORA3_hc_waves():
    #get_model
    mco = mc(nID='NORA3_hc_waves', sd="2019-1-1", ed="2019-1-1")
    assert mco.__class__.__name__ == 'model_class'
    mco.populate()
    print(mco.vars)
    assert len(vars(mco).keys()) == 18
    assert len(mco.vars.keys()) == 3

# Fails
#def test_MY_L4_thredds():
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

