#import sys
#import os
#import numpy as np
#from datetime import datetime
#import pytest
#
#from wavy.satellite_module import satellite_class as sc
#from wavy.gridder import gridder_class as gc
#from wavy.grid_stats import apply_metric

#def test_gridder_lowres(test_data, benchmark):
#    sco = sc(sdate="2020-11-1",edate="2020-11-3",region="global",
#             path_local=str(test_data/"L3"))
#    bb = (-170,170,-75,75)
#    res = (4, 4)
#    print("initialize")
#    gco = gc(oco=sco,bb=bb,res=res)
#    print("assign obs")
#    ovals,mvals,Midx = gco.get_obs_grid_idx()
#
#    ov, olo, ola = apply_metric(gco,metric="mean")
#    gv, glo, gla = apply_metric(gco,metric="mean_group")
#
#    # np.testing.assert_array_equal(olo, glo)
#    # np.testing.assert_array_equal(ola, gla)
#    np.testing.assert_array_almost_equal(ov, gv)
#
#    print("compute metric on grid")
#    var_gridded = benchmark(apply_metric, gco,metric="mean_group")
#
#def test_gridder_highres(test_data, benchmark):
#    sco = sc(sdate="2020-11-1",edate="2020-11-3",region="global",
#             path_local=str(test_data/"L3"))
#    bb = (-170,170,-75,75)
#    res = (.5, .5)
#    print("initialize")
#    gco = gc(oco=sco,bb=bb,res=res)
#    print("assign obs")
#    ovals,mvals,Midx = gco.get_obs_grid_idx()
#
#    ov, olo, ola = apply_metric(gco,metric="mean")
#    gv, glo, gla = apply_metric(gco,metric="mean_group")
#
#    # np.testing.assert_array_equal(olo, glo)
#    # np.testing.assert_array_equal(ola, gla)
#    np.testing.assert_array_almost_equal(ov, gv)
#
#    print("compute metric on grid")
#    var_gridded = benchmark(apply_metric, gco,metric="mean_group")
