import sys
import os
from datetime import datetime
import pytest

from wavy.wconfig import load_or_default
from wavy.satellite_module import satellite_class as sc

@pytest.mark.need_credentials
def test_collectors_cmems_L3(tmpdir):
    sco = sc(sd='2023-2-1 12', ed='2023-2-1 12',
             nID='cmems_L3_NRT', name='s3a')
    sco.download(path=tmpdir, nproc=4)
    # check if files were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))
              if '.nc' in filelist[i]]
    assert len(nclist) >= 1

@pytest.mark.need_credentials
def test_collectors_cci_v3_20Hz(tmpdir):
    sco = sc(sd='2020-2-1 12', ed='2020-2-1 12',
             nID='L2_20Hz_s3a', name='s3a')
    sco.download(path=tmpdir, nproc=8)
    # check if files were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))
              if '.nc' in filelist[i]]
    assert len(nclist) >= 1

def test_manually_specified_reader(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name,
             varalias=varalias,
             twin=twin)
    # read data
    sco = sco.populate(reader='read_local_ncfiles',
                       path=str(test_data/"L3/s3a"))
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    print(vlst)
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    print(flst)
    assert len(flst) >= 47
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()

def test_default_reader(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name,
             varalias=varalias,
             twin=twin)
    # read data
    sco = sco.populate(path=str(test_data/"L3/s3a"))
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) >= 47
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()

def test_polygon_region(test_data):
    sd = "2022-2-01 01"
    ed = "2022-2-03 23"
    name = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, name=name,
             varalias=varalias,
             twin=twin)
    # read data
    sco = sco.populate(path=str(test_data/"L3/s3a"))
    sco = sco.crop_to_region('NordicSeas')
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) >= 47
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()

def test_rectangular_region(test_data):
    sd = "2022-2-01 01"
    ed = "2022-2-03 23"
    name = 's3a'
    varalias = 'Hs'
    nID = 'cmems_L3_NRT'
    # init satellite_object
    sco = sc(sd=sd, ed=ed, nID=nID, name=name,
             varalias=varalias)
    # read data
    sco = sco.populate(path=str(test_data/"L3/s3a"))
    sco = sco.crop_to_region('Sulafj')
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) >= 47
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()


def test_poi_storm_track(test_data):
    import pandas as pd
    from wavy.utils import parse_date

    # read track
    f = pd.read_csv(test_data/"track/Katrina_track.csv")
    # convert dates to datetime
    dt = [parse_date(d) for d in f.date]
    lons = f.lon.values
    lats = f.lat.values

    # define poi dictionary for track
    poi_dict = {'time': dt, 'lons': lons, 'lats': lats}

    # retrievals
    sco = sc(twin=180, distlim=200, name='multi',
         nID='CCIv1_L3', varalias='Hs',  # default
         poi=poi_dict)

    sco = sco.populate(path=str(test_data/"L3/multi"))

    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    print(vlst)
    assert len(vlst) == 19
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) >= 47
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()

    # get closest only



# def test_write_to_nc(test_data):
    # write to nc
    #sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
    # check if created -> assert
    # read nc
    # check if varalias assert