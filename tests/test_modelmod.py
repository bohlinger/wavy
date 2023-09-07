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
    assert len(vars(mco).keys()) == 17
    assert len(mco.vars.keys()) == 3

def test_ww3_unstr_reader():
    #get_model

    bb = (5.8, 6.61, 62.3, 63.1)
    res = (0.01, 0.01)  # lon/lat

    mco = mc(nID='ww3_unstr', sd="2019-3-24 10", ed="2019-3-24 10")
    assert mco.__class__.__name__ == 'model_class'
    mco.populate(res=res, bb=bb, interp='nearest')
    print(mco.vars)

    assert len(vars(mco).keys()) == 17
    assert len(mco.vars.keys()) == 3
