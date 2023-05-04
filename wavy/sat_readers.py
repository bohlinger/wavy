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

# own imports
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.ncmod import read_netcdfs
#from wavy.ncmod import read_netcdfs_hidefix
#from wavy.ncmod import read_mf_netcdfs
from wavy.ncmod import read_swim_netcdfs
from wavy.wconfig import load_or_default
from wavy.utils import parse_date, calc_deep_water_T
from wavy.utils import find_included_times, find_included_times_pd
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')

def read_local_ncfiles(**kwargs):
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
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    twin = kwargs.get('twin')
    ncvar = kwargs.get('ncvar')

    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)
    # retrieve sliced data
    #
    # ds = read_netcdfs_hidefix(pathlst)
    ds = read_netcdfs(pathlst)
    # ds = read_mf_netcdfs(pathlst)
    #
    ds_sort = ds.sortby('time')
    ds_sliced = ds_sort.sel(time=slice(sd, ed))
    var_sliced = ds_sliced[[ncvar]]
    # rename ncvars to respective wavy aliases
    # todo:
    #   attribute cf-standard names to all variables
    #   (only unique time steps?)
    return var_sliced


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

    # establish coords if defined in config file
    timestr = satellite_dict[nID]['vardef']['time']
    lonstr = satellite_dict[nID]['vardef']['lons']
    latstr = satellite_dict[nID]['vardef']['lats']

    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)
    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias, variable_info,
                            satellite_dict[nID], ncmeta)
    # retrieve sliced data
    ds = read_netcdfs(pathlst)
    ds_sort = ds.sortby(timestr)

    # get indices for included time period
    nptime = ds_sort[timestr].data
    idx = find_included_times_pd(nptime, sdate=sd, edate=ed)
    #ds_sliced = ds_sort.sel(time=slice(sd, ed))
    #dtime = [parse_date(str(nptime[idx][i])) for i in range(len(nptime[idx]))]
    dtime = nptime[idx]
    lons = list(((ds_sort[lonstr].data[idx] - 180) % 360) - 180)
    lats = list(ds_sort[latstr].data[idx])

    unxt = (nptime[idx].astype(int) / 10**9)

    # make dict and start with stdvarname for varalias
    stdvarname = variable_info[varalias]['standard_name']
    vardict = {}
    vardict[stdvarname] = list(ds_sort[ncvar].data[idx])
    vardict['longitude'] = lons
    vardict['latitude'] = lats
    vardict['time'] = unxt
    vardict['datetime'] = dtime
    vardict['time_unit'] = variable_info['time']['units']
    return vardict


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
    product = kwargs.get('product')
    varalias = kwargs.get('varalias')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    twin = kwargs.get('twin')

    # adjust start and end
    sdate = sdate - timedelta(minutes=twin)
    edate = edate + timedelta(minutes=twin)
    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias,variable_info,
                            satellite_dict[product],ncmeta)
    # retrieve data
    vardict = read_swim_netcdfs(pathlst,varalias)
    # rm NaN from 'time'
    tmpt = np.array(vardict['time'])
    tmpt = tmpt[~np.isnan(tmpt)]
    # parse time and add to dict
    dtime = [parse_date(d) for d in np.array(tmpt).astype(str)]
    vardict['datetime'] = dtime
    vardict['time_unit'] = variable_info['time']['units']
    vardict['time'] = netCDF4.date2num(vardict['datetime'],vardict['time_unit'])
    # lon tranformation
    tlons = list(((np.array(vardict['longitude']) - 180) % 360) - 180)
    vardict['longitude'] = tlons
    # find idx for time period
    tidx = find_included_times(dtime,sdate=sdate,edate=edate)
    # adjust dict
    for key in vardict.keys():
        if key != 'meta' and key != 'time_unit':
            vardict[key]=list(np.array(vardict[key])[tidx])
    # if peak wave length transform to peak period
    if kwargs.get('return_var',varalias) == 'Tp':
        Tp = calc_deep_water_T(np.array(vardict[varalias]))
        vardict[variable_info['Tp']['standard_name']] = Tp
        # change varalias to stdvarname from variable_info
    else:
        # change varalias to stdvarname from variable_info
        vardict[variable_info[varalias]['standard_name']] = vardict[varalias]
    # delete varalias key from dict
    vardict.pop(varalias)
    return vardict

def read_local_files(**kwargs) -> dict:
    '''
    wrapping function to read altimetry files

    return:
        vardict - dictionary of variables for altimeter data
    '''
    dispatch_reader = {
                #'cmems_L3_NRT': read_local_ncfiles,
                'cmems_L3_NRT': read_local_ncfiles,
                'cmems_L3_s6a': read_local_ncfiles,
                'cmems_L3_MY': read_local_ncfiles,
                'cci_L2P': read_local_ncfiles,
                'cci_L3': read_local_ncfiles,
                'cfo_swim_L2P': read_local_ncfiles_swim,
                'L2_20Hz_s3a': read_local_20Hz_files
                }
    product = kwargs.get('nID')
    # check if product available in dispatcher
    if product in dispatch_reader.keys():
        pass
    else:
        product = 'cmems_L3_NRT'

    vardict = dispatch_reader[product](**kwargs)
    return vardict
