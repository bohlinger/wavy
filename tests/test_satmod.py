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

def test_readers_and_write_to_nc_cmems_NRT(tmpdir, test_data):
    # evoke fct get_remote_files
    sdate = "2020-11-1 12"
    edate = "2020-11-1 12"
    region = 'NordicSeas'
    mission = 's3a'
    varalias = 'Hs'
    twin = 30
    product = 'cmems_L3_NRT'
    # init satellite_object and check for polygon region
    sco = sc(sdate=sdate,edate=edate,region=region,
             mission=mission,twin=twin,varalias=varalias,
             product=product,path_local=str(test_data/"L3"))
    assert sco.__class__.__name__ == 'satellite_class'
    assert len(vars(sco).keys()) >= 11
    assert len(sco.vars.keys()) >= 6
    assert not 'error' in vars(sco).keys()
    # write to nc
    sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
    # check if created -> assert
    # read nc
    # check if varalias assert
