import sys
import os
from datetime import datetime
import pytest

from wavy.wconfig import load_or_default
import wavy.sat_collectors
from wavy.satmod import satellite_class as sc

sdate = "2020-11-1 12"
sdate_dt = datetime(2020,11,1,12)
edate = "2020-11-1 12"
edate_dt = datetime(2020,11,1,12)
region = 'NordicSeas'
mission = 's3a'
varalias = 'Hs'
twin = 30
nproc = 1
product = 'cmems_L3_NRT'

satellite_dict = load_or_default('satellite_specs.yaml')


@pytest.mark.need_credentials
def test_ftp_files_and_satellite_class_features(tmpdir):
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt,edate=edate_dt,
                            twin=twin,nproc=nproc,
                            product=product,api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1
    # init satellite_object and check for polygon region
    sco = sc(sdate=sdate,edate=edate,region=region,
             mission=mission,twin=twin,varalias=varalias,
             product=product,path_local=tmpdir)
    assert sco.__class__.__name__ == 'satellite_class'
    assert len(vars(sco).keys()) >= 11
    assert len(sco.vars.keys()) >= 6
    assert not 'error' in vars(sco).keys()
    # write to nc
    sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
    # check if created -> assert
    # read nc
    # check if varalias assert
