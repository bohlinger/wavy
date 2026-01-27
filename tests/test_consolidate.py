import sys
import os
import numpy as np
from datetime import datetime
import pytest
from wavy.satellite_module import satellite_class as sc
from wavy.consolidate import consolidate_class as cs

def test_consolidate_satellite(test_data):
    # satellite consolidate
    sco1 = sc(sd="2022-2-1",ed ="2022-2-3",region="NordicSeas", nID="cmems_L3_NRT",
              name='s3a').populate(path=str(test_data/"L3/s3a"))
    sco2 = sc(sd="2022-2-1",ed ="2022-2-3",region="NordicSeas", nID="cmems_L3_NRT",
              name='s3b').populate(path=str(test_data/"L3/s3b"))
    
    cso = cs([sco1,sco2])
    
    len(list(vars(cso).keys())) >= 8
    len(list(cso.vars.keys())) == 3

def test_consolidate_satellite_multivar(test_data):

    varalias = ['Hs','U']
    
    # satellite consolidate
    sco1 = sc(sd="2022-2-1",ed ="2022-2-3",region="NordicSeas", nID="cmems_L3_NRT",
              name='s3a', varalias=varalias).populate(path=str(test_data/"L3/s3a"))
    sco2 = sc(sd="2022-2-1",ed ="2022-2-3",region="NordicSeas", nID="cmems_L3_NRT",
              name='s3b', varalias=varalias).populate(path=str(test_data/"L3/s3b"))
    
    cso = cs([sco1,sco2])
    
    len(list(vars(cso).keys())) >= 8
    len(list(cso.vars.keys())) == 3