import sys
import os
from datetime import datetime
import pytest

from wavy.wconfig import load_or_default

now = datetime.now()
date = datetime(now.year, now.month, now.day)
region = 'NordicSeas'
sat = 's3a'
varalias = 'Hs'
twin = 60
nproc = 1
instr = 'altimeter'
provider = 'cmems'

satellite_dict = load_or_default('/satellite_specs.yaml')


@pytest.mark.need_credentials
def test_ftp_files_and_init_satellite_class(tmpdir):
    from wavy import satmod

    # make satpath
    path_remote = satellite_dict[instr]\
                    [provider]['remote']['path']\
                    + sat + '/'
    # evoke fct get_remote_files
    satmod.get_remote_files(path_remote, tmpdir,
                            date, date, twin, nproc,
                            instr, provider)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    print(filelist)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1
    # init satellite_object and check for polygon region

    from wavy.satmod import satellite_class

    sa_obj = satellite_class(sdate=date,
                             region=region,
                             sat=sat,
                             varalias=varalias,
                             download_path=tmpdir)
    assert sa_obj.__class__.__name__ == 'satellite_class'
    assert len(vars(sa_obj).keys()) >= 11
    assert len(sa_obj.vars.keys()) >= 6
    assert not 'error' in vars(sa_obj).keys()
    # another one checking for error

