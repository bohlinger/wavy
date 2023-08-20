from wavy.insitumod import insitu_class as ic

varalias = 'Hs' # default
sd = "2021-8-2 01"
ed = "2021-8-3 00"

#def test_from_nc():
#    nID = 'D_Breisundet_wave'
#    sensor = 'wavescan'
#    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,sensor=sensor)
#    assert ico.__class__.__name__ == 'insitu_class'
#    assert len(vars(ico).keys()) >= 10
#    assert len(ico.vars.keys()) >= 6
#    assert not 'error' in vars(ico).keys()
#
#def test_from_thredds():
#    nID = 'D_Breisundet_wave'
#    sensor = 'wavescan'
#    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,fifo='thredds',sensor=sensor)
#    assert ico.__class__.__name__ == 'insitu_class'
#    assert len(vars(ico).keys()) >= 10
#    assert len(ico.vars.keys()) >= 6
#    assert not 'error' in vars(ico).keys()
#
#def test_to_nc(tmpdir):
#
def test_from_frost_v1():
    nID = 'draugen'
    sensor = 'MKIIIradar_1'
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, sensor=sensor)
    print(ico)
    print(vars(ico).keys())
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) >= 12
    new = ico.populate()
    print(new.vars.keys())
    #assert len(ico.vars.keys()) >= 6
    #assert not 'error' in vars(ico).keys()

#def test_cmems_insitu(test_data):
