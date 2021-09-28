import sys
import os
from datetime import datetime, timedelta
import pytest

from wavy.wconfig import load_or_default
import wavy.satmod
from wavy.satmod import satellite_class as sc
from wavy.collocmod import collocation_class as cc

sdate = datetime(2020,11,1,12)
edate = datetime(2020,11,1,12)
region = 'NordicSeas'
sat = 's3a'
varalias = 'Hs'
twin = 30
nproc = 1
provider = 'cmems'

satellite_dict = load_or_default('satellite_specs.yaml')


def test_sat_collocation_and_validation(test_data,tmpdir):
    # read sat data
    sco = sc(sdate=sdate,edate=edate,
             region=region,sat=sat,
             twin=twin,varalias=varalias,
             provider=provider,
             path_local=str(test_data))
    # collocate
    cco = cc(model='mwam4',obs_obj_in=sco,distlim=6,
             leadtime='best',date_incr=1)
    assert len(vars(cco).keys()) >= 12
    assert len(cco.vars.keys()) >= 13
    assert not 'error' in vars(sco).keys()
    # write to nc
    cco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
