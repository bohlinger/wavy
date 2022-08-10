import sys
import os
from datetime import datetime
import pytest

from wavy.satmod import satellite_class as sc
from wavy.gridder import gridder_class as gc

def test_gridder_highres(test_data, benchmark):
    sco = sc(sdate="2020-11-1",edate="2020-11-3",region="global",
             path_local=str(test_data/"L3"))
    bb = (-170,170,-75,75)
    res = (.5,.5)
    print("initialize")
    gco = gc(oco=sco,bb=bb,res=res)
    print("assign obs")
    ovals,Midx = gco.get_obs_grid_idx()
    print("compute metric on grid")
    var_gridded = benchmark(gco.apply_metric, Midx,ovals,metric="mean")

def test_gridder_lowres(test_data, benchmark):
    sco = sc(sdate="2020-11-1",edate="2020-11-3",region="global",
             path_local=str(test_data/"L3"))
    bb = (-170,170,-75,75)
    res = (4, 4)
    print("initialize")
    gco = gc(oco=sco,bb=bb,res=res)
    print("assign obs")
    ovals,Midx = gco.get_obs_grid_idx()


    print("compute metric on grid")
    var_gridded = benchmark(gco.apply_metric, Midx,ovals,metric="mean")
