#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to read model files
files for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library igports
import numpy as np
from datetime import timedelta
import xarray as xr

# own imports
from wavy.ncmod import read_netcdfs
from wavy.wconfig import load_or_default
from wavy.utils import build_xr_ds

# ---------------------------------------------------------------------#
# read yaml config files:
model_dict = load_or_default('model_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')
# ---------------------------------------------------------------------#

def read_ww3_4km(**kwargs):
    pathlst = kwargs.get('pathlst')
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    twin = kwargs.get('twin')

    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)
    # retrieve sliced data
    ds = read_netcdfs(pathlst)
    ds_sort = ds.sortby('time')

    ds_sliced = ds_sort.sel(time=slice(sd, ed))

    return ds_sliced

def read_ww3_unstructured(**kwargs):
    return ds

def read_ww3_unstructured_grid(**kwargs):
    return ds
