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
import os
import pandas as pd
import zipfile
import tempfile
from tqdm import tqdm
import xarray as xr
from datetime import timedelta
import netCDF4

# own imports
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.ncmod import read_netcdfs, read_netcdfs_zipped_lru
from wavy.ncmod import read_swim_netcdfs
from wavy.wconfig import load_or_default
from wavy.utils import parse_date, find_included_times, calc_deep_water_T
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_specs.yaml')
variable_info = load_or_default('variable_info.yaml')

def unzip_eumetsat(pathlst: list, tmpdir: str):
    """
    Function to unzip eumetsat files prior to reading

    param:
        pathlst - list of paths to zipped files
        tmpdir - temporary folder to unzipped files

    return:
        pathlst_new - new list of paths to unzipped files
    """
    for count, f in enumerate(pathlst):
        zipped = zipfile.ZipFile(f)
        enhanced_measurement = zipped.namelist()[-1]
        extracted = zipped.extract(enhanced_measurement,
                                   path=tmpdir.name)
        fname = extracted.split('/')[-1]
        dignumstr = '_{num:0' + str(len(str(len(pathlst)))) + 'd}.nc'
        # cp extracted file to parent tmp
        cmdstr = ('cp ' + extracted + ' ' + tmpdir.name
                + '/' + fname[0:-3]
                + dignumstr.format(num=count))
                #+ '_{num:04d}.nc'.format(num=count))
        os.system(cmdstr)
        # delete subfolder
        cmdstr = ('rm -r '
                 + os.path.dirname(os.path.realpath(extracted)))
        os.system(cmdstr)
    flst = os.listdir(tmpdir.name)
    pathlst_new = []
    for f in flst:
        pathlst_new.append(os.path.join(tmpdir.name,f))
    return pathlst_new

def read_local_files_eumetsat(**kwargs):
    '''
    Read and concatenate all data to one timeseries for each variable.
    Fct is tailored to EUMETSAT files.
    '''
    pathlst = kwargs.get('pathlst')
    product = kwargs.get('product')
    varalias = kwargs.get('varalias')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    twin = kwargs.get('twin')
    # adjust start and end
    sdate = sdate - timedelta(minutes=twin)
    edate = edate + timedelta(minutes=twin)
    # --- find variable cf names --- #
    print ("Processing " + str(int(len(pathlst))) + " files")
    print (pathlst[0])
    print (pathlst[-1])
    # --- find ncvar cf names --- #
    tmpdir = tempfile.TemporaryDirectory()
    zipped = zipfile.ZipFile(pathlst[0])
    enhanced_measurement = zipped.namelist()[-1]
    extracted = zipped.extract(enhanced_measurement, path=tmpdir.name)
    stdvarname = variable_info[varalias]['standard_name']
    ncmeta = ncdumpMeta(extracted)
    ncvar = get_filevarname(varalias,variable_info,
                            satellite_dict[product],ncmeta)
    latname = get_filevarname('lats',variable_info,
                            satellite_dict[product],ncmeta)
    lonname = get_filevarname('lons',variable_info,
                            satellite_dict[product],ncmeta)
    timename = get_filevarname('time',variable_info,
                            satellite_dict[product],ncmeta)
    tmpdir.cleanup()
    # --- create vardict --- #
    vardict = {}
    ds = read_netcdfs_zipped_lru(pathlst,ncvar,dim=timename)
    ds_sort = ds.sortby(timename)
    ds_sliced = ds_sort.sel({timename:slice(sdate, edate)})
    # make dict and start with stdvarname for varalias
    var_sliced = ds_sliced[ncvar]
    vardict = {}
    vardict[stdvarname] = list(var_sliced.values)
    # add coords to vardict
    # 1. retrieve list of coordinates
    coords_lst = list(var_sliced.coords.keys())
    # 2. iterate over coords_lst
    for varname in coords_lst:
        stdcoordname = ds_sliced[varname].attrs['standard_name']
        if stdcoordname == 'longitude':
            vardict[stdcoordname] = \
                list(((ds_sliced[varname].values - 180) % 360) - 180)
        elif stdcoordname == 'time':
            # convert to unixtime
            df_time = ds_sliced[varname].to_dataframe()
            unxt = (pd.to_datetime(df_time[varname]).view(int) / 10**9)
            vardict[stdcoordname] = unxt.values
            vardict['time_unit'] = variable_info[stdcoordname]['units']
        else:
            vardict[stdcoordname] = list(ds_sliced[varname].values)
    return vardict

def read_local_ncfiles(**kwargs):
    """
    Wrapping function to read satellite netcdf files.

    param:
        pathlst - list of paths to be parsed
        product - product as specified in satellite_specs.yaml
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
    # retrieve sliced data
    ds = read_netcdfs(pathlst)
    ds_sort = ds.sortby('time')
    ds_sliced = ds_sort.sel(time=slice(sdate, edate))
    # make dict and start with stdvarname for varalias
    stdvarname = variable_info[varalias]['standard_name']
    var_sliced = ds_sliced[ncvar]
    vardict = {}
    vardict[stdvarname] = list(var_sliced.values)
    # add coords to vardict
    # 1. retrieve list of coordinates
    coords_lst = list(var_sliced.coords.keys())
    # 2. iterate over coords_lst
    for varname in coords_lst:
        stdcoordname = ds_sliced[varname].attrs['standard_name']
        if stdcoordname == 'longitude':
            vardict[stdcoordname] = \
                list(((ds_sliced[varname].values - 180) % 360) - 180)
        elif stdcoordname == 'time':
            # convert to unixtime
            df_time = ds_sliced[varname].to_dataframe()
            unxt = (pd.to_datetime(df_time[varname]).view(int) / 10**9)
            vardict[stdcoordname] = unxt.values
            vardict['time_unit'] = variable_info[stdcoordname]['units']
        else:
            vardict[stdcoordname] = list(ds_sliced[varname].values)
    return vardict

def read_local_ncfiles_swim(**kwargs):
    """
    Wrapping function to read swim netcdf files.

    param:
        pathlst - list of paths to be parsed
        product - product as specified in satellite_specs.yaml
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
                'cmems_L3_NRT':read_local_ncfiles,
                'cmems_L3_s6a':read_local_ncfiles,
                'cmems_L3_MY':read_local_ncfiles,
                'cci_L2P':read_local_ncfiles,
                'cci_L3':read_local_ncfiles,
                'eumetsat_L2':read_local_files_eumetsat,
                'cfo_swim_L2P':read_local_ncfiles_swim
                }
    product = kwargs.get('product')
    vardict = dispatch_reader[product](**kwargs)
    return vardict
