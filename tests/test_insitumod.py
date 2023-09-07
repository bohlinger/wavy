from wavy.insitu_module import insitu_class as ic

varalias = 'Hs'  # default
sd = "2021-8-2 01"
ed = "2021-8-3 00"

def test_from_thredds():
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=sensor)
    print(ico)
    print(vars(ico).keys())
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) == 12
    new = ico.populate()
    print(new.vars.keys())
    print(len(new.vars.keys()))
    assert len(new.vars.keys()) == 3

def test_from_thredds_twinID():
    nID = 'A_Sulafjorden_wave'
    twinID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID=nID, twinID=twinID, sd=sd, ed=ed,
             varalias=varalias, name=sensor)
    print(ico)
    print(vars(ico).keys())
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) == 12
    new = ico.populate()
    print(new.vars.keys())
    print(len(new.vars.keys()))
    assert len(new.vars.keys()) == 3


def test_from_frost_v1():
    nID = 'draugen'
    sensor = 'MKIIIradar_1'
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=sensor)
    print(ico)
    print(vars(ico).keys())
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) == 12
    new = ico.populate()
    print(new.vars.keys())
    print(len(new.vars.keys()))
    assert len(new.vars.keys()) == 3

#def test_cmems_insitu(test_data):
#    varalias = 'Hs'  # default
#    sd = "2023-7-2 00"
#    ed = "2023-7-3 00"
#
#    nID = 'MO_Ekofisk_monthly'
#    name = 'NA'
#    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
#    print(ico)
#    print(vars(ico).keys())
#    assert ico.__class__.__name__ == 'insitu_class'
#    assert len(vars(ico).keys()) == 12
#    ico.list_input_files(show=True)
#    new = ico.populate()
#    print(new.vars.keys())
#    print(len(new.vars.keys()))
#    assert len(new.vars.keys()) == 3


#def test_to_nc(tmpdir):
