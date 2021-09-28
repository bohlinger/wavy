import sys
import os
from datetime import datetime, timedelta
import pytest

from wavy.wconfig import load_or_default
import wavy.satmod
from wavy.satmod import satellite_class

sdate = datetime(2020,11,1,12)
edate = datetime(2020,11,1,12)
region = 'NordicSeas'
sat = 's3a'
varalias = 'Hs'
twin = 30
nproc = 1
provider = 'cmems'

satellite_dict = load_or_default('satellite_specs.yaml')


@pytest.mark.need_credentials
def test_ftp_files_and_satellite_class_features(tmpdir):
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':sat}
    wavy.satmod.get_remote_files(tmpdir,
                            sdate, edate, twin, nproc,
                            provider,api_url,sat,dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1
    # init satellite_object and check for polygon region
    sco = satellite_class(sdate=sdate,
                          edate=edate,
                          region=region,
                          sat=sat,
                          twin=twin,
                          varalias=varalias,
                          provider=provider,
                          path_local=tmpdir)
    assert sco.__class__.__name__ == 'satellite_class'
    assert len(vars(sco).keys()) >= 10
    assert len(sco.vars.keys()) >= 6
    assert not 'error' in vars(sco).keys()
    # write to nc
    sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
    # check if created -> assert
    # read nc
    # check if varalias assert
