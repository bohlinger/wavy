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
from wavy.ncmod import read_netcdfs_with_credentials_aggregated
from wavy.ncmod import read_swim_netcdfs
from wavy.wconfig import load_or_default
from wavy.utils import parse_date, calc_deep_water_T
from wavy.utils import find_included_times, find_included_times_pd
from wavy.ncmod import build_xr_ds
from wavy.credentials import get_credentials
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')

def read_wavy_ncfiles(**kwargs):
    """
    Wrapping function to read wavy netcdf files.

    param:
        pathlst - list of paths to be parsed
        sd - start date (datetime object)
        ed - start date (datetime object)
        twin - time window (temporal constraint) in minutes

    return:
        dictionary of variables for the satellite_class object
    """
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


def read_local_ncfiles(**kwargs):
    """
    Wrapping function to read satellite netcdf files.
    """
    pathlst = kwargs.get('pathlst')
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    twin = kwargs.get('twin')
    varalias = kwargs.get('varalias')
    nID = kwargs.get('nID')
    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    # varnames
    varname = get_filevarname(varalias, variable_info,
                              satellite_dict[nID], ncmeta)
    lonsname = get_filevarname('lons', variable_info,
                               satellite_dict[nID], ncmeta)
    latsname = get_filevarname('lats', variable_info,
                               satellite_dict[nID], ncmeta)
    timename = get_filevarname('time', variable_info,
                               satellite_dict[nID], ncmeta)
    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)

    # retrieve sliced data
    ds = read_netcdfs(pathlst)

    ds_sort = ds.sortby(timename)
    ds_sliced = ds_sort.sel(time=slice(sd, ed))
    var_sliced = ds_sliced[[varname, lonsname, latsname]]

    # build xr dataset with time as only dim/coords
    varnames = {varalias: varname, 'lons': lonsname,
                'lats': latsname, 'time': timename}
    var = (var_sliced[varname], var_sliced[lonsname],
           var_sliced[latsname], var_sliced[timename])

    ds_new = build_xr_ds(var, varnames, varalias)

    return ds_new


def read_remote_ncfiles_aggregated(**kwargs):
    """
    Wrapping function to read remote opendap satellite netcdf files
    that use credentials
    """
    pathlst = kwargs.get('pathlst')
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    twin = kwargs.get('twin')
    varalias = kwargs.get('varalias')
    nID = kwargs.get('nID')
    remoteHostName = kwargs.get('remoteHostName')
    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    # varnames
    varname = get_filevarname(varalias, variable_info,
                              satellite_dict[nID], ncmeta)
    lonsname = get_filevarname('lons', variable_info,
                               satellite_dict[nID], ncmeta)
    latsname = get_filevarname('lats', variable_info,
                               satellite_dict[nID], ncmeta)
    timename = get_filevarname('time', variable_info,
                               satellite_dict[nID], ncmeta)
    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)

    # get credentials
    usr, pw = get_credentials(remoteHostName)

    # retrieve sliced data
    ds = read_netcdfs_with_credentials_aggregated(
            pathlst, remoteHostName, usr, pw, dim=timename)

    ds_sort = ds.sortby(timename)
    ds_sliced = ds_sort.sel(time=slice(sd, ed))
    var_sliced = ds_sliced[[varname, lonsname, latsname]]

    # build xr dataset with time as only dim/coords
    varnames = {varalias: varname, 'lons': lonsname,
                'lats': latsname, 'time': timename}
    var = (var_sliced[varname], var_sliced[lonsname],
           var_sliced[latsname], var_sliced[timename])

    ds_new = build_xr_ds(var, varnames, varalias)

    return ds_new


def read_local_20Hz_files(**kwargs):
    """
    Wrapping function to read satellite netcdf files.

    param:
        pathlst - list of paths to be parsed
        product - product as specified in satellite_cfg.yaml
        varalias
        sd - start date (datetime object)
        ed - start date (datetime object)
        twin - time window (temporal constraint) in minutes

    return:
        dictionary of variables for the satellite_class object
    """
    pathlst = kwargs.get('pathlst')
    nID = kwargs.get('nID')
    varalias = kwargs.get('varalias')
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    twin = kwargs.get('twin')

    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)
    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    varname = get_filevarname(varalias, variable_info,
                              satellite_dict[nID], ncmeta)
    lonstr = get_filevarname('lons', variable_info,
                              satellite_dict[nID], ncmeta)
    latstr = get_filevarname('lats', variable_info,
                              satellite_dict[nID], ncmeta)
    timestr = get_filevarname('time', variable_info,
                              satellite_dict[nID], ncmeta)


    # retrieve sliced data
    ds = read_netcdfs(pathlst)
    ds_sort = ds.sortby(timestr)

    # get indices for included time period
    nptime = ds_sort[timestr].data
    idx = find_included_times_pd(nptime, sdate=sd, edate=ed)
    var_sliced = ds_sort[varname].values[idx]
    lons = ds_sort[lonstr].values[idx]
    lats = ds_sort[latstr].values[idx]
    nptime = nptime[idx]
    varnames = {varalias: varname, 'lons': lonstr,
                'lats': latstr, 'time': timestr}
    var = (var_sliced, lons, lats, nptime)
    # build xarray ds
    ds_new = build_xr_ds(var, varnames, varalias)
    return ds_new

def read_local_ncfiles_swim(**kwargs):
    """
    Wrapping function to read swim netcdf files.

    param:
        pathlst - list of paths to be parsed
        product - product as specified in satellite_cfg.yaml
        varalias
        sd - start date (datetime object)
        ed - start date (datetime object)
        twin - time window (temporal constraint) in minutes

    return:
        dictionary of variables for the satellite_class object
    """
    pathlst = kwargs.get('pathlst')
    varalias = kwargs.get('varalias')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    twin = kwargs.get('twin')

    # adjust start and end
    sdate = sdate - timedelta(minutes=twin)
    edate = edate + timedelta(minutes=twin)
    # retrieve data
    vardict = read_swim_netcdfs(pathlst, varalias)
    # rm NaN from 'time'
    tmpt = np.array(vardict['time'])
    tmpt = tmpt[~np.isnan(tmpt)]
    # parse time and add to dict
    dtime = [parse_date(d) for d in np.array(tmpt).astype(str)]
    vardict['datetime'] = dtime
    vardict['time_unit'] = variable_info['time']['units']
    vardict['time'] = netCDF4.date2num(vardict['datetime'],
                                       vardict['time_unit'])
    # lon tranformation
    tlons = list(((np.array(vardict['longitude']) - 180) % 360) - 180)
    vardict['longitude'] = tlons
    # find idx for time period
    tidx = find_included_times(dtime, sdate=sdate, edate=edate)
    # adjust dict
    for key in vardict.keys():
        if key != 'meta' and key != 'time_unit':
            vardict[key] = list(np.array(vardict[key])[tidx])
    # if peak wave length transform to peak period
    if kwargs.get('return_var', varalias) == 'Tp':
        Tp = calc_deep_water_T(np.array(vardict[varalias]))
        vardict[variable_info['Tp']['standard_name']] = Tp
        # change varalias to stdvarname from variable_info
    else:
        # change varalias to stdvarname from variable_info
        vardict[variable_info[varalias]['standard_name']] = vardict[varalias]
    # delete varalias key from dict
    vardict.pop(varalias)
    return vardict
