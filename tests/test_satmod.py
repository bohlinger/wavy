import sys
import os
from datetime import datetime
import pytest

from wavy.wconfig import load_or_default
from wavy.satmod import satellite_class as sc

satellite_dict = load_or_default('satellite_cfg.yaml')

#def test_ftp_files_and_satellite_class_features(tmpdir):
@pytest.mark.need_credentials
def test_collectors_cmems_L3(tmpdir):
    sco = sc(sd='2023-2-1 12', ed='2023-2-1 12',
             nID='cmems_L3_NRT', mission='s3a')
    sco.download(path=tmpdir, nproc=2)
    # check if files were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))
              if '.nc' in filelist[i]]
    assert len(nclist) >= 1

def test_manually_specified_reader(tmpdir, test_data):
    # evoke fct get_remote_files
    sd = "2020-11-1 12"
    ed = "2020-11-1 12"
    mission = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    print(test_data/"L3")
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, mission=mission,
             varalias=varalias,
             twin=twin)
    # read data
    sco.populate(reader='read_local_ncfiles', path=str(test_data/"L3"))
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 20
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) == 32
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()

def test_default_reader(tmpdir, test_data):
    # evoke fct get_remote_files
    sd = "2020-11-1 12"
    ed = "2020-11-1 12"
    mission = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    print(test_data/"L3")
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, mission=mission,
             varalias=varalias,
             twin=twin)
    # read data
    sco.populate(path=str(test_data/"L3"))
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 20
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) == 32
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()

def test_polygon_region(tmpdir, test_data):
    # evoke fct get_remote_files
    sd = "2020-11-01 01"
    ed = "2020-11-03 23"
    mission = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    print(test_data/"L3")
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, mission=mission,
             varalias=varalias,
             twin=twin)
    # read data
    sco.populate(reader='read_local_ncfiles', path=str(test_data/"L3"))
    sco = sco.crop_to_region('NordicSeas')
    print(len(sco.vars['time']))
    assert sco.__class__.__name__ == 'satellite_class'
    # compare number of available variables
    vlst = list(vars(sco).keys())
    assert len(vlst) == 20
    # compare number of available functions
    dlst = dir(sco)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) == 32
    assert type(sco.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(sco).keys()


    # write to nc
    #sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
    # check if created -> assert
    # read nc
    # check if varalias assert
