from datetime import datetime

from wavy.insitumod import insitu_class as ic

varalias = 'Hs' # default
sd = datetime(2021,8,2,1)
ed = datetime(2021,8,3,0)

def test_from_d22(test_data):
    nID = 'draugen'
    sensor = 'MKIIIradar_1'
    st_obj = ic(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1,path_local=str(test_data))
    assert st_obj.__class__.__name__ == 'insitu_class'
    assert len(vars(st_obj).keys()) >= 10
    assert len(st_obj.vars.keys()) >= 6
    assert not 'error' in vars(st_obj).keys()

def test_from_nc():
    nID = 'D_Breisundet'
    sensor = 'wavescan'
    st_obj = ic(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1)
    assert st_obj.__class__.__name__ == 'insitu_class'
    assert len(vars(st_obj).keys()) >= 10
    assert len(st_obj.vars.keys()) >= 6
    assert not 'error' in vars(st_obj).keys()
