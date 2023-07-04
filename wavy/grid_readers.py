#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to read satellite altimetry
files for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library igports
import numpy as np
from datetime import timedelta
import netCDF4
import pandas as pd
import xarray as xr

# own imports
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.ncmod import read_netcdfs
#from wavy.ncmod import read_netcdfs_hidefix
from wavy.ncmod import tpe_hidefix
#from wavy.ncmod import read_mf_netcdfs
from wavy.ncmod import read_swim_netcdfs
from wavy.wconfig import load_or_default
from wavy.utils import parse_date, calc_deep_water_T
from wavy.utils import find_included_times, find_included_times_pd
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_cfg.yaml')
model_dict = load_or_default('model_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')

def read_ww3_unstructured(**kwargs):
    """
    Wrapping function to read satellite netcdf files.

    param:
        pathlst - list of paths to be parsed
        product - product as specified in satellite_cfg.yaml
        varalias
        sd - start date (datetime object)
        ed - start date (datetime object)

    return:
        dictionary of variables for the satellite_class object
    """
    print("Reading unstructured grid...")
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    varalias = kwargs.get('varalias')
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')

    # establish coords if defined in config file
    timestr = model_dict[nID]['vardef']['time']
    lonstr = model_dict[nID]['vardef']['lons']
    latstr = model_dict[nID]['vardef']['lats']

    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    varname = get_filevarname(varalias, variable_def,
                              model_dict[nID], ncmeta)
    # retrieve sliced data
    ds_lst = []
    for p in pathlst:
        ds = xr.open_dataset(p)
        # subsample times accroding to temporal range
        if sd == ed:
            ds = ds.sel({timestr: sd})
            times = [ds[timestr].data]
        else:
            ds = ds.sel({timestr: slice(sd, ed)})
            times = ds[timestr].data
        for t in times:
            print(t)
            if (t >= sd and t <= ed):
                # varnames = (varname, lonstr, latstr, timestr)
                if len(times) < 2:
                    var = (ds[varname].data,
                           ds[lonstr].data,
                           ds[latstr].data,
                           t.reshape((1,)))
                else:
                    var = (ds[varname].sel({timestr: t}).data,
                           ds[lonstr].data,
                           ds[latstr].data,
                           t.reshape((1,)))

                new = get_gridded_dataset(var, t, **kwargs)
                ds_lst.append(new)

    print("Concatenate ...")
    combined = xr.concat(ds_lst, 'time',
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print("... done concatenating")

    return combined

def get_gridded_dataset(var, t, **kwargs):
    varalias = kwargs.get('varalias')

    if kwargs.get('interp') is None:
        print("Apply gridding, no interpolation")
        # grid data
        gridvar, lon_grid, lat_grid = \
            grid_point_cloud_ds(var[0], var[1], var[2], t, **kwargs)

        var_means = gridvar['mor']
        # transpose dims
        var_means = np.transpose(var_means)
        field_shape = list(var_means.shape) + [1]
        var_means = var_means.reshape(field_shape)
    else:
        print("Apply gridding with interpolation")
        var_means, lon_grid, lat_grid = \
                grid_point_cloud_interp_ds(
                    var[0], var[1], var[2], **kwargs)
        # transpose dims
        field_shape = list(var_means.shape) + [1]
        var_means = var_means.reshape(field_shape)
    # create xr.dataset
    ds = build_xr_ds_grid(
            var_means,
            np.unique(lon_grid), np.unique(lat_grid),
            t.reshape((1,)), varalias)
    return ds


def build_xr_ds_grid(var_means, lon_grid, lat_grid, t, varalias):
    import xarray as xr

    ds = xr.Dataset({
            varalias: xr.DataArray(
                    data=var_means,
                    dims=['lats', 'lons', 'time'],
                    coords={'lats': lat_grid,
                            'lons': lon_grid,
                            'time': t},
                    )
                },
            attrs={'title': 'wavy dataset'}
        )
    return ds


def grid_point_cloud_ds(values, lons, lats, t, **kwargs):
    from wavy.gridder import gridder_class as gc
    from wavy.grid_stats import apply_metric

    bb = kwargs.get('bb')
    res = kwargs.get('res')

    gco = gc(lons=lons, lats=lats, values=values,
             bb=bb, res=res,
             sdate=t, edate=t)
    gridvar, lon_grid, lat_grid = apply_metric(gco=gco)

    return gridvar, lon_grid, lat_grid

def grid_point_cloud_interp_ds(values, lons, lats, **kwargs):
    from scipy.interpolate import griddata

    res = kwargs.get('res')
    bb = kwargs.get('bb')

    x, y = np.arange(bb[0], bb[1] + res[0], res[0]),\
           np.arange(bb[2], bb[3] + res[1], res[1])
    lon_grid, lat_grid = np.meshgrid(x, y)

    points = [[lons[i], lats[i]] for i in range(len(lons))]
    points = np.array(points)

    gridvar = griddata(points, values, (lon_grid, lat_grid),
                       method=kwargs.get('method', 'linear'))

    return gridvar, lon_grid, lat_grid

