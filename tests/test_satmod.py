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
instr = 'altimeter'
provider = 'cmems'

satellite_dict = load_or_default('satellite_specs.yaml')


@pytest.mark.need_credentials
def test_ftp_files_and_init_satellite_class(tmpdir):
    # evoke fct get_remote_files
    api_url = None
    sat = 's3a'
    wavy.satmod.get_remote_files(tmpdir,
                            sdate, edate, twin, nproc,
                            instr, provider,api_url,sat)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1
    # init satellite_object and check for polygon region
    sa_obj = satellite_class(sdate=sdate,
                             edate=edate,
                             region=region,
                             sat=sat,
                             twin=twin,
                             varalias=varalias,
                             provider=provider,
                             path_local=tmpdir)
    assert sa_obj.__class__.__name__ == 'satellite_class'
    assert len(vars(sa_obj).keys()) >= 10
    assert len(sa_obj.vars.keys()) >= 6
    assert not 'error' in vars(sa_obj).keys()
    # another one checking for error

