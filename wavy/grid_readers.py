#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to read data directly into grids.
'''
# --- import libraries ------------------------------------------------#
# standard library igports
import numpy as np
import netCDF4
import pandas as pd
import xarray as xr

# own imports
from wavy.ncmod import get_filevarname
from wavy.wconfig import load_or_default
from wavy.utils import parse_date
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_cfg.yaml')
model_dict = load_or_default('model_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')

def read_ww3_unstructured_to_grid(**kwargs):
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
    print(" Reading unstructured grid...")
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    varalias = kwargs.get('varalias')
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    meta = kwargs.get('meta')

    # get meta data
    varstr = get_filevarname(varalias, variable_def,
                              model_dict[nID], meta)
    lonstr = get_filevarname('lons', variable_def,
                              model_dict[nID], meta)
    latstr = get_filevarname('lats', variable_def,
                              model_dict[nID], meta)
    timestr = get_filevarname('time', variable_def,
                              model_dict[nID], meta)


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
            # convert to datetime for comparison
            dt = parse_date(str(t))
            if (dt >= sd and dt <= ed):
                # varnames = (varname, lonstr, latstr, timestr)
                if len(times) < 2:
                    var = (ds[varstr].data,
                           ds[lonstr].data,
                           ds[latstr].data,
                           t.reshape((1,)))
                else:
                    var = (ds[varstr].sel({timestr: t}).data,
                           ds[lonstr].data,
                           ds[latstr].data,
                           t.reshape((1,)))

                new = get_gridded_dataset(var, t, **kwargs)
                ds_lst.append(new)

    print(" Concatenate ...")
    combined = xr.concat(ds_lst, 'time',
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print(" ... done concatenating")

    return combined


def get_gridded_dataset(var, t, **kwargs):
    varstr = kwargs.get('varalias')

    if kwargs.get('interp') is None:
        print(" Apply gridding, no interpolation")
        # grid data
        gridvar, lon_grid, lat_grid = \
            grid_point_cloud_ds(var[0], var[1], var[2], t, **kwargs)

        var_means = gridvar['mor']
        # transpose dims
        var_means = np.transpose(var_means)
        field_shape = (list(var_means.shape[::-1]) + [1])[::-1]
        var_means = var_means.reshape(field_shape)
    else:
        print(" Apply gridding with interpolation")
        var_means, lon_grid, lat_grid = \
                grid_point_cloud_interp_ds(
                    var[0], var[1], var[2], **kwargs)
        # transpose dims
        field_shape = (list(var_means.shape[::-1]) + [1])[::-1]
        var_means = var_means.reshape(field_shape)
    # create xr.dataset
    ds = build_xr_ds_grid(
            var_means,
            np.unique(lon_grid), np.unique(lat_grid),
            t.reshape((1,)),
            varstr=varstr)
    return ds


def build_xr_ds_grid(var_means, lon_grid, lat_grid, t, **kwargs):
    print(" building xarray dataset from grid")
    varstr = kwargs.get('varstr')

    ds = xr.Dataset({
            varstr: xr.DataArray(
                     data=var_means,
                     dims=['time', 'latitude', 'longitude'],
                     coords={'latitude': lat_grid,
                             'longitude': lon_grid,
                             'time': t},
                     attrs=variable_def[varstr],
                     ),
            'lons': xr.DataArray(
                     data=lon_grid,
                     dims=['longitude'],
                     coords={'longitude': lon_grid},
                     attrs=variable_def['lons'],
                     ),
            'lats': xr.DataArray(
                     data=lat_grid,
                     dims=['latitude'],
                     coords={'latitude': lat_grid},
                     attrs=variable_def['lats'],
                     ),
                 },
             attrs={'title': 'wavy dataset'}
         )
    return ds

def build_xr_ds_grid_2D(var_means, lon_grid, lat_grid, t, **kwargs):
    print(" building xarray dataset from grid")
    varstr = kwargs.get('varstr')
    lon_grid_coord = kwargs.get('lon_grid_coord')
    lat_grid_coord = kwargs.get('lat_grid_coord')

    ds = xr.Dataset({
            varstr: xr.DataArray(
                     data=var_means,
                     dims=['time', 'latitude', 'longitude'],
                     coords={'latitude': lat_grid_coord,
                             'longitude': lon_grid_coord,
                             'time': t},
                     attrs=variable_def[varstr],
                     ),
            'lons': xr.DataArray(
                     data=lon_grid,
                     dims=['latitude', 'longitude'],
                     coords={'longitude': lon_grid_coord},
                     attrs=variable_def['lons'],
                     ),
            'lats': xr.DataArray(
                     data=lat_grid,
                     dims=['latitude', 'longitude'],
                     coords={'latitude': lat_grid_coord},
                     attrs=variable_def['lats'],
                     ),
                 },
             attrs={'title': 'wavy dataset'}
         )
    return ds


def grid_point_cloud_ds(values, lons, lats, t, **kwargs):
    print(' gridding point cloud')
    from wavy.gridder_module import gridder_class as gc
    from wavy.grid_stats import apply_metric

    bb = kwargs.get('bb')
    res = kwargs.get('res')

    gco = gc(lons=lons, lats=lats, values=values,
             bb=bb, res=res,
             sdate=t, edate=t)
    gridvar, lon_grid, lat_grid = apply_metric(gco=gco)

    return gridvar, lon_grid, lat_grid

def grid_point_cloud_interp_ds(values, lons, lats, **kwargs):
    print(' gridding point cloud with interpolation')
    from scipy.interpolate import griddata

    res = kwargs.get('res')
    bb = kwargs.get('bb')

    x, y = np.arange(bb[0], bb[1] + res[0]/2, res[0]),\
           np.arange(bb[2], bb[3] + res[1]/2, res[1])
    lon_grid, lat_grid = np.meshgrid(x, y)

    points = [[lons[i], lats[i]] for i in range(len(lons))]
    points = np.array(points)

    gridvar = griddata(points, values, (lon_grid, lat_grid),
                       method=kwargs.get('method', 'linear'))

    return gridvar, lon_grid, lat_grid

