import sys
import os
from datetime import datetime, timedelta
import pytest

from wavy.wconfig import load_or_default
import wavy.satmod
from wavy.satmod import satellite_class as sc
from wavy.collocmod import collocation_class as cc
from wavy.insitumod import insitu_class as ic

varalias = 'Hs'

satellite_dict = load_or_default('satellite_specs.yaml')


def test_sat_collocation_and_validation(test_data,tmpdir):
    sd = datetime(2020,11,1,12)
    ed = datetime(2020,11,1,12)
    region = 'NordicSeas'
    sat = 's3a'
    provider = 'cmems'
    twin = 30
    # read sat data
    sco = sc(sdate=sd,edate=ed,
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
    # test validation
    cco.validate_collocated_values()

def test_insitu_collocation_and_validation():
    sd = datetime(2021,8,2,1)
    ed = datetime(2021,8,2,3)
    nID = 'D_Breisundet_wave'
    sensor = 'wavescan'
    ico = ic(nID,sensor,sd,ed,varalias=varalias,stwin=1,date_incr=1)
    # collocate
    cco = cc(model='mwam4',obs_obj_in=ico,distlim=6,
             leadtime='best',date_incr=1)
    # test validation
    cco.validate_collocated_values()


