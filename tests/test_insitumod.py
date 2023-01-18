from wavy.insitumod import insitu_class as ic

varalias = 'Hs' # default
sd = "2021-8-2 01"
ed = "2021-8-3 00"

def test_from_d22(test_data):
    nID = 'draugen'
    sensor = 'MKIIIradar_1'
    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,path_local=str(test_data/"d22"),sensor=sensor)
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) >= 10
    assert len(ico.vars.keys()) >= 6
    assert not 'error' in vars(ico).keys()

def test_from_nc():
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,sensor=sensor)
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) >= 10
    assert len(ico.vars.keys()) >= 6
    assert not 'error' in vars(ico).keys()

def test_from_thredds():
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,fifo='thredds',sensor=sensor)
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) >= 10
    assert len(ico.vars.keys()) >= 6
    assert not 'error' in vars(ico).keys()

def test_to_nc(tmpdir):
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,sensor=sensor)
    ico.write_to_nc(pathtofile=tmpdir.join('test.nc'))
