import pytest
from datetime import datetime, timedelta
import yaml
import numpy as np
import os
from copy import deepcopy

from wavy.insitumod import insitu_class as ic
from wavy.filtermod import apply_land_mask

sd = "2021-8-2 01"
ed = "2021-8-3 12"

#ico = ic(nID,sensor,sd,ed)
#test_dict = deepcopy(ico.vars())

#@pytest.fixture
#def test_data():
#    return os.path.abspath(os.path.join(\
#           os.path.dirname( __file__ ),'data'))

def test_landmask():
    vardict = { 'latitude':[60.12,62.24, 64.08,65.08, 67.65,68.95],
                'longitude':[-23.47,-21.54, -19.32,-17.8, -13.97,-10.99]}
    d,m = apply_land_mask(vardict)
    assert len(m[m==False]) == int(2)

def test_cleaners():
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID,sensor,sd,ed,priorOp='square',cleaner='linearGAM',postOp='root',date_incr=1,filterData=True)
    assert len(vars(ico).keys()) == 15
    assert 'filter' in vars(ico).keys()
    assert 'filterSpecs' in vars(ico).keys()

def test_smoothers():
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID,sensor,sd,ed,smoother='blockMean',date_incr=1,filterData=True)
    assert len(vars(ico).keys()) == 15
    assert 'filter' in vars(ico).keys()
    assert 'filterSpecs' in vars(ico).keys()

