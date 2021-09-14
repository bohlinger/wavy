import pytest
from datetime import datetime, timedelta
import yaml
import numpy as np

from wavy.stationmod import station_class as sc

varalias = 'Hs' # default
sd = datetime(2021,8,2,1)
ed = datetime(2021,8,3,0)

#@pytest.fixture
#def test_data():
#    return os.path.abspath(os.path.join(\
#           os.path.dirname( __file__ ),'data'))

#def test_from_d22():
#    nID = 'draugen'
#    sensor = 'MKIIIradar_1'
#    st_obj = sc(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1,path_local=test_data)
#    assert st_obj.__class__.__name__ == 'station_class'
#    assert len(vars(st_obj).keys()) >= 10
#    assert len(st_obj.vars.keys()) >= 6
#    assert not 'error' in vars(st_obj).keys()

def test_from_nc():
    nID = 'D_Breisundet'
    sensor = 'wavescan'
    st_obj = sc(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1)
    assert st_obj.__class__.__name__ == 'station_class'
    assert len(vars(st_obj).keys()) >= 10
    assert len(st_obj.vars.keys()) >= 6
    assert not 'error' in vars(st_obj).keys()
