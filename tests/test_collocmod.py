import sys
import os
from datetime import datetime, timedelta
import pytest

from wavy.satmod import satellite_class as sc
from wavy.collocmod import collocation_class as cc
from wavy.insitumod import insitu_class as ic

# include possibility for collocating different variable
# varalias = 'Hs', 'U', aso...

def test_sat_collocation_and_validation(test_data,tmpdir):
    # evoke fct get_remote_files
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    model = 'ww3_4km'
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, name=name,
             varalias=varalias, twin=twin)
    # read data
    sco.populate(reader='read_local_ncfiles', path=str(test_data/"L3/s3a"))
    # crop to region
    sco = sco.crop_to_region(model)

    # collocate
    cco = cc(oco=sco, model=model, leadtime='best', distlim=6)
    assert len(vars(cco).keys()) == 15
    assert len(cco.vars.keys()) == 9
#    # write to nc
#    cco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
#    # test validation
#    cco.validate_collocated_values()
#
#def test_insitu_collocation_and_validation():
#    sd = "2021-8-2 01"
#    ed = "2021-8-2 03"
#    nID = 'D_Breisundet_wave'
#    sensor = 'wavescan'
#    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,sensor=sensor)
#    # collocate
#    cco = cc(model='mwam4',obs_obj_in=ico,distlim=6,
#             leadtime='best',date_incr=1)
#    # test validation
#    cco.validate_collocated_values()
