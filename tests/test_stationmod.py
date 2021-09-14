from datetime import datetime, timedelta
import numpy as np

from wavy.stationmod import station_class as sc

# challenge is to write a test that does not fail when VPN
# is missing but only gives an informative warning

import yaml
import numpy as np
from datetime import datetime,timedelta
from wavy.stationmod import station_class
varalias = 'Hs' # default
sd = datetime(2021,8,1,1)
ed = datetime(2021,8,2,0)


#def test_from_d22():
#    nID = 'draugen'
#    sensor = 'MKIIIradar_1'
#    st_obj = station_class(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1)
#    assert st_obj.__class__.__name__ == 'station_class'
#    assert len(vars(st_obj).keys()) >= 10
#    assert len(st_obj.vars.keys()) >= 6
#    assert not 'error' in vars(st_obj).keys()

def test_from_nc():
    nID = 'D_Breisundet'
    sensor = 'wavescan'
    st_obj = station_class(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1)
    assert st_obj.__class__.__name__ == 'station_class'
    assert len(vars(st_obj).keys()) >= 10
    assert len(st_obj.vars.keys()) >= 6
    assert not 'error' in vars(st_obj).keys()
