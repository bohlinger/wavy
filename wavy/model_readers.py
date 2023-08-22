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
import netCDF4

# own imports
from wavy.wconfig import load_or_default
from wavy.utils import build_xr_ds
import xarray as xr

# ---------------------------------------------------------------------#
# read yaml config files:
model_dict = load_or_default('model_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')
# ---------------------------------------------------------------------#

def read_ww3_4km(**kwargs):
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    fc_dates = kwargs.get('fc_dates')
    varname = kwargs.get('varname')
    ds_lst = []
    # retrieve sliced data
    for i in range(len(fc_dates)):
        d = fc_dates[i]
        p = pathlst[i]
        ds = xr.open_dataset(p)
        ds_sliced = ds.sel({model_dict[nID]['vardef']['time']: d})
        ds_sliced = ds_sliced[[varname,
                               model_dict[nID]['vardef']['lons'],
                               model_dict[nID]['vardef']['lats']]]

        ds_lst.append(ds_sliced)

    print("Concatenate ...")
    combined = xr.concat(ds_lst, model_dict[nID]['vardef']['time'],
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print("... done concatenating")

    return combined

def read_ww3_unstructured(**kwargs):
    return ds

def read_ww3_unstructured_grid(**kwargs):
    return ds
