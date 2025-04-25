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
from wavy.grid_readers import build_xr_ds_grid, build_xr_ds_grid_2D
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.ncmod import read_netcdfs_with_credentials_aggregated
from wavy.utils import parse_date
from wavy.utils import collocate_times
from wavy.credentials import get_credentials

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

    print(' Build dataset')
    print(' dataset ready!')

    return combined

def read_meps(**kwargs):
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

    print(' Build dataset')
    print(' dataset ready!')

    return combined


def read_noresm_making_waves(**kwargs):
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    fc_dates = kwargs.get('fc_dates')
    varname = kwargs.get('varname')
    varalias = kwargs.get('varalias')
    timename = model_dict[nID]['vardef']['time']
    lonsname = model_dict[nID]['vardef']['lons']
    latsname = model_dict[nID]['vardef']['lats']

    ds_lst = []
    for i in range(len(pathlst)):
        ds = xr.open_dataset(pathlst[i])
        ds2 = ds.convert_calendar("all_leap")
        datetimeindex = ds2.indexes[timename].to_datetimeindex()
        ds[timename] = datetimeindex
        del ds2
        var = ds[varname].values
        lons = ds[lonsname].values
        lats = ds[latsname].values
        Mlons, Mlats = np.meshgrid(lons, lats)
        timedt = ds[timename].values
        tlst = [parse_date(str(d)) for d in timedt]
        idx = collocate_times(tlst, target_t=[fc_dates[i]])
        time = np.array([np.array(tlst)[idx[0]]]).reshape((1,))
        varin = var[idx[0], :, :].reshape((1, len(lats), len(lons)))
        ds = build_xr_ds_grid_2D(varin, Mlons, Mlats, time,
                                 lon_grid_coord=lons,
                                 lat_grid_coord=lats,
                                 varstr=varalias)
        ds_lst.append(ds)

    print(" Concatenate ...")
    combined = xr.concat(ds_lst, timename,
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print(" ... done concatenating")

    return combined

def read_remote_ncfiles_aggregated_credentials(**kwargs):
    """
    Wrapping function to read remote opendap satellite netcdf files
    that use credentials
    """
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    twin = kwargs.get('twin')
    varalias = kwargs.get('varalias')
    nID = kwargs.get('nID')
    remoteHostName = kwargs.get('remoteHostName')
    path = kwargs.get('pathlst')[0]

    # varnames
    varname = model_dict[nID]['vardef'][varalias]
    lonsname = model_dict[nID]['vardef']['lons']
    latsname = model_dict[nID]['vardef']['lats']
    timename = model_dict[nID]['vardef']['time']

    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)

    # get credentials
    usr, pw = get_credentials(remoteHostName)

    # retrieve sliced data
    ds = read_netcdfs_with_credentials_aggregated(
            path, remoteHostName, usr, pw)

    ds_sliced = ds.sel(time=slice(sd, ed))
    var_sliced = ds_sliced[[varname, lonsname, latsname]]

    # forge into correct format varalias, lons, lats with dim time
    ds = build_xr_ds_grid(var_sliced[varname],
                          var_sliced[lonsname],
                          var_sliced[latsname],
                          var_sliced[timename],
                          varstr=varalias)

    return ds

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
        ds_lst.append(build_xr_ds(var_tuple, varnames)
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

@lru_cache(maxsize=32)
def read_single_field_lru_ecwam(
filestr, varname, lonsname, latsname, timename):
    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    filestr = filestr.replace('\\', '')
    f = netCDF4.Dataset(filestr, 'r')
    # get coordinates and time
    var = f.variables[varname][:, 0, :, :].squeeze()
    lons = f.variables[lonsname][:]
    lats = f.variables[latsname][:]
    time = f.variables[timename]
    timedt = list(
        netCDF4.num2date(time[:], units=time.units))
    f.close()
    return var, lons, lats, timedt

def read_ecwam(**kwargs):
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    fc_dates = kwargs.get('fc_dates')
    varname = kwargs.get('varname')
    varalias = kwargs.get('varalias')
    timename = model_dict[nID]['vardef']['time']
    lonsname = model_dict[nID]['vardef']['lons']
    latsname = model_dict[nID]['vardef']['lats']

    ds_lst = []
    for i in range(len(pathlst)):
        var, lons, lats, timedt = read_single_field_lru_ecwam(
                                          pathlst[i],
                                          varname,
                                          lonsname, latsname,
                                          timename)
        tlst = [parse_date(str(d)) for d in timedt]
        idx = collocate_times(tlst, target_t=[fc_dates[i]])
        time = np.array([np.array(tlst)[idx[0]]]).reshape((1,))
        varin = var[idx[0], :, :].reshape((1, len(lats), len(lons)))
        ds = build_xr_ds_grid(varin, lons, lats, time,
                              varstr=varalias)
        ds_lst.append(ds)

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


def read_era(**kwargs):
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    fc_dates = kwargs.get('fc_dates')
    varname = kwargs.get('varname')
    ds_lst = []
    # retrieve sliced data
    for i in range(len(fc_dates)):
        d = parse_date(fc_dates[i])
        p = pathlst[i]
        ds = xr.open_dataset(p)
        ds_sliced = ds.sel({model_dict[nID]['vardef']['time']: d},
                            method='nearest')
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

    print(' Build dataset')
    print(' dataset ready!')

    return combined


def read_NORA3_wind(**kwargs):
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    fc_dates = kwargs.get('fc_dates')
    varname = kwargs.get('varname')
    hlevel = kwargs.get('heightlevel', 10)
    ds_lst = []
    # retrieve sliced data
    for i in range(len(fc_dates)):
        d = fc_dates[i]
        p = pathlst[i]
        ds = xr.open_dataset(p)
        ds_sliced = ds.sel({model_dict[nID]['vardef']['time']: d})
        ds_sliced = ds_sliced.sel({'height': hlevel})
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

    print(' Build dataset')
    print(' dataset ready!')

    return combined


#import pyproj
#import xarray as xr
#
#n = xr.open_dataset('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#
#cf_proj = n.projection_stere.attrs
#
#crs_from = pyproj.CRS.from_cf(cf_proj)
#crs_to = pyproj.Proj("EPSG:4326").crs
#
#T = pyproj.Transformer.from_crs(crs_from=crs_from, crs_to=crs_to)
#x=1000
#y=1000
#lon, lat = T.transform(x, y)
#print(lon, lat)
