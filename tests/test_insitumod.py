from wavy.insitu_module import insitu_class as ic
from wavy.insitu_module import poi_class as pc
import pytest
import os
from datetime import datetime
import numpy as np


def test_from_thredds():
    varalias = 'Hs'
    sd = "2021-8-2 01"
    ed = "2021-8-3 00"
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
    varalias = 'Hs'  # default
    sd = "2021-8-2 01"
    ed = "2021-8-3 00"
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

@pytest.mark.need_credentials
def test_from_frost_v1():
    varalias = 'Hs'  # default
    sd = "2021-8-2 01"
    ed = "2021-8-3 00"
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


def test_cmems_insitu_monthly(test_data):
    varalias = 'Hs'  # default
    sd = "2023-7-2 00"
    ed = "2023-7-3 00"
    nID = 'MO_Draugen_monthly'
    name = 'Draugen'
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
    print(ico)
    print(vars(ico).keys())
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) == 12
    ico.list_input_files(show=True)
    new = ico.populate(path=str(test_data/"insitu/monthly/Draugen"))
    new.list_input_files(show=True)
    print(new.vars.keys())
    print(len(new.vars.keys()))
    assert len(new.vars.keys()) == 3
    # check if some data was imported
    assert len(new.vars['time']) > 0
    # check that not all data is nan
    assert not all(np.isnan(v) for v in new.vars['time'])
    assert not all(np.isnan(v) for v in new.vars['Hs'])
    assert not all(np.isnan(v) for v in new.vars['lons'])
    assert not all(np.isnan(v) for v in new.vars['lats'])


def test_cmems_insitu_daily(test_data):
    varalias = 'Hs'  # default
    sd = "2023-8-20 00"
    ed = "2023-8-21 00"
    nID = 'MO_Draugen_daily'
    name = 'Draugen'
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
    print(ico)
    print(vars(ico).keys())
    assert ico.__class__.__name__ == 'insitu_class'
    assert len(vars(ico).keys()) == 12
    ico.list_input_files(show=True)
    new = ico.populate(path=str(test_data/"insitu/daily/Draugen"))
    new.list_input_files(show=True)
    print(new.vars.keys())
    print(len(new.vars.keys()))
    assert len(new.vars.keys()) == 3
    # check if some data was imported
    assert len(new.vars['time']) > 0
    # check that not all data is nan
    assert not all(np.isnan(v) for v in new.vars['time'])
    assert not all(np.isnan(v) for v in new.vars['Hs'])
    assert not all(np.isnan(v) for v in new.vars['lons'])
    assert not all(np.isnan(v) for v in new.vars['lats'])


#@pytest.mark.need_credentials
#def test_insitu_collectors(tmpdir):
#    varalias = 'Hs'  # default
#    sd = "2024-3-10 00"
#    ed = "2024-3-11 00"
#    nID = 'AR_TS_MO_monthly'
#    name = 'Draugen'
#    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
#    ico.download(path=tmpdir, nproc=2)
#    # check if files were download to tmp directory
#    filelist = os.listdir(tmpdir)
#    nclist = [i for i in range(len(filelist))
#              if '.nc' in filelist[i]]
#    assert len(nclist) >= 1


def test_insitu_poi(tmpdir):
    # define poi dictionary for track
    dt = [datetime(2023, 7, 1), datetime(2023, 7, 2), datetime(2023, 7, 3)]
    lats = [56.5, 59.3, 64.3]
    lons = [3.5, 1.8, 4.2]
    poi_dict = {'time': dt, 'lons': lons, 'lats': lats}
    pco = pc(poi_dict)
    assert len(vars(pco)) == 12
    assert len(list(pco.vars.keys())) == 3
