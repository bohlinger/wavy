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
import netCDF4
from functools import lru_cache

# own imports
from wavy.wconfig import load_or_default
from wavy.utils import build_xr_ds
from wavy.grid_readers import get_gridded_dataset
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.utils import parse_date

# ---------------------------------------------------------------------#
# read yaml config files:
model_dict = load_or_default('model_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')
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

    print(" Concatenate ...")
    combined = xr.concat(ds_lst, model_dict[nID]['vardef']['time'],
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print(" ... done concatenating")

    return combined

@lru_cache(maxsize=32)
def read_single_field_lru(
filestr, varname, lonsname, latsname, timename):
    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    filestr = filestr.replace('\\', '')
    f = netCDF4.Dataset(filestr, 'r')
    # get coordinates and time
    var = f.variables[varname][:]
    lons = f.variables[lonsname][:]
    lats = f.variables[latsname][:]
    time = f.variables[timename]
    timedt = list(
        netCDF4.num2date(time[:], units=time.units))
    f.close()
    return var, lons, lats, timedt

def read_field(**kwargs):
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    fc_dates = kwargs.get('fc_dates')
    varname = kwargs.get('varname')
    timename = model_dict[nID]['vardef']['time']
    lonsname = model_dict[nID]['vardef']['lons']
    latsname = model_dict[nID]['vardef']['lats']

    ds_lst = []
    for i in range(len(pathlst)):
        var_tuple = read_single_field_lru(pathlst[i],
                                          varname,
                                          lonsname, latsname, timename)
        varnames = (varname, 'lons', 'lats', 'time')
        ds_lst.append(build_xr_ds(var_tuple, varnames)\
                      .sel({timename: fc_dates[i]}))

    print(" Concatenate ...")
    combined = xr.concat(ds_lst, timename,
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print(" ... done concatenating")

    return combined

def read_ww3_unstructured_to_grid(**kwargs):
    from wavy.grid_readers import read_ww3_unstructured_to_grid
    ds = read_ww3_unstructured_to_grid(**kwargs)
    return ds
