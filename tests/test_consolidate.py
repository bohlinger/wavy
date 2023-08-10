import sys
import os
import numpy as np
from datetime import datetime
import pytest

from wavy.consolidate import consolidate_class as cs
from wavy.satmod import satellite_class as sc
from wavy.collocmod import collocation_class as cc

def test_consolidate_scos(test_data, benchmark):
    sco1 = sc(sdate="2020-11-1 12", region="NordicSeas",
              mission='s3a', path_local=str(test_data/"L3"))
    sco2 = sc(sdate="2020-11-1 13", region="NordicSeas",
              mission='s3b', path_local=str(test_data/"L3"))
    cso = cs([sco1, sco2])
    assert len(list(vars(cso).keys())) >= 8
    assert len(list(cso.vars.keys())) == 6

def test_consolidate_ccos(test_data, benchmark):
    sco1 = sc(sdate="2020-11-1 12", region="NordicSeas",
              mission='s3a', path_local=str(test_data/"L3"))
    sco2 = sc(sdate="2020-11-1 13", region="NordicSeas",
              mission='s3b', path_local=str(test_data/"L3"))
    cco1 = cc(model='mwam4', obs_obj_in=sco1, distlim=6,
              leadtime='best', date_incr=1)
    cco2 = cc(model='mwam4', obs_obj_in=sco2, distlim=6,
              leadtime='best', date_incr=1)
    cso = cs([cco1, cco2], obj_type='collocation')
    print(cso.vars.keys())
    assert len(list(vars(cso).keys())) >= 8
    assert len(list(cso.vars.keys())) == 16
