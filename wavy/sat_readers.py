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
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import read_netcdfs, get_filevarname
from wavy.wconfig import load_or_default
from wavy.utils import parse_date
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_specs.yaml')
variable_info = load_or_default('variable_info.yaml')

def unzip_eumetsat(pathlst,tmpdir):
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

    # --- find variable cf names --- #
    print ("Processing " + str(int(len(pathlst))) + " files")
    print (pathlst[0])
    print (pathlst[-1])
    # --- find ncvar cf names --- #
    tmpdir = tempfile.TemporaryDirectory()
    zipped = zipfile.ZipFile(pathlst[0])
    enhanced_measurement = zipped.namelist()[-1]
    extracted = zipped.extract(enhanced_measurement, path=tmpdir.name)
    stdname = variable_info[varalias]['standard_name']
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
    tmpdir = tempfile.TemporaryDirectory()
    print('tmp directory is established:',tmpdir.name)
    for f in tqdm(pathlst):
        #path = f[0:-len(f.split('/')[-1])]
        zipped = zipfile.ZipFile(f)
        enhanced_measurement = zipped.namelist()[-1]
        extracted = zipped.extract(enhanced_measurement,
                                   path=tmpdir.name)
        nc = netCDF4.Dataset(extracted)
        if stdname in vardict.keys():
            vardict[stdname] += list(nc.variables[ncvar][:])
        else:
            vardict[stdname] = list(nc.variables[ncvar][:])
        if 'longitude' in vardict.keys():
            vardict['longitude'] += list(nc.variables[lonname][:])
        else:
            vardict['longitude'] = list(nc.variables[lonname][:])
        if 'latitude' in vardict.keys():
            vardict['latitude'] += list(nc.variables[latname][:])
        else:
            vardict['latitude'] = list(nc.variables[latname][:])
        time_obj = nc.variables[timename]
        dtime_tmp = netCDF4.num2date(time_obj[:],time_obj.units)
        dtime = [parse_date(str(d)) for d in dtime_tmp]
        time = netCDF4.date2num(dtime,variable_info['time']['units'])
        if timename in vardict.keys():
            vardict['time'] += list(time)
            vardict['datetime'] += list(dtime)
        else:
            vardict['time'] = list(time)
            vardict['datetime'] = list(dtime)
        # transform to datetime
    # transform to -180 to 180 degrees
    tmp = np.array(vardict['longitude'])
    vardict['longitude'] = list(((tmp - 180) % 360) - 180)
    vardict['time_unit'] = variable_info['time']['units']
    tmpdir.cleanup()
    return vardict

def read_local_files_eumetsat_old(**kwargs):
    '''
    Read and concatenate all data to one timeseries for each variable.
    Fct is tailored to EUMETSAT files.
    '''
    pathlst = kwargs.get('pathlst')
    product = kwargs.get('product')
    varalias = kwargs.get('varalias')

    # --- find variable cf names --- #
    print ("Processing " + str(int(len(pathlst))) + " files")
    print (pathlst[0])
    print (pathlst[-1])
    # --- find ncvar cf names --- #
    tmpdir = tempfile.TemporaryDirectory()
    zipped = zipfile.ZipFile(pathlst[0])
    enhanced_measurement = zipped.namelist()[-1]
    extracted = zipped.extract(enhanced_measurement, path=tmpdir.name)
    stdname = variable_info[varalias]['standard_name']
    ncmeta = ncdumpMeta(extracted)
    ncvar = get_filevarname(varalias,variable_info,
                            satellite_dict[product],ncmeta)
    tmpdir.cleanup()
    # --- create vardict --- #
    vardict = {}
    tmpdir = tempfile.TemporaryDirectory()
    print('tmp directory is established:',tmpdir.name)
    for f in tqdm(pathlst):
        #path = f[0:-len(f.split('/')[-1])]
        zipped = zipfile.ZipFile(f)
        enhanced_measurement = zipped.namelist()[-1]
        extracted = zipped.extract(enhanced_measurement,
                                   path=tmpdir.name)
        ds = xr.open_dataset(extracted)
        ds_var = ds[ncvar]
        if stdname in vardict.keys():
            vardict[stdname] += list(ds[ncvar])
        else:
            vardict[stdname] = list(ds[ncvar])
        coords_lst = list(ds_var.coords.keys())
        for varname in coords_lst:
            stdcoordname = ds[varname].attrs['standard_name']
            if stdcoordname == 'longitude':
                if stdcoordname in vardict.keys():
                    vardict[stdcoordname] += list(ds[varname].values)
                else:
                    vardict[stdcoordname] = list(ds[varname].values)
            elif stdcoordname == 'time':
                # convert to unixtime
                df_time = ds[varname].to_dataframe()
                unxt = (pd.to_datetime(df_time[varname])\
                        .view(int) / 10**9)
                if stdcoordname in vardict.keys():
                    vardict[stdcoordname] += list(unxt.values)
                else:
                    vardict[stdcoordname] = list(unxt.values)
            else:
                if stdcoordname in vardict.keys():
                    vardict[stdcoordname] += list(ds[varname].values)
                else:
                    vardict[stdcoordname] = list(ds[varname].values)
    # transform to -180 to 180 degrees
    tmp = np.array(vardict['longitude'])
    vardict['longitude'] = list(((tmp - 180) % 360) - 180)
    vardict['time_unit'] = variable_info['time']['units']
    tmpdir.cleanup()
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

def read_local_files(**kwargs):
    '''
    wrapping function to read altimetry files

    return:
        vardict - dictionary of variables for altimeter data
    '''
    dispatch_reader = {
                'cmems_L3_NRT':read_local_ncfiles,
                'cmems_L3_MY':read_local_ncfiles,
                'cci_L2P':read_local_ncfiles,
                'cci_L3':read_local_ncfiles,
                'eumetsat_L2':read_local_files_eumetsat,
                }
    product = kwargs.get('product')
    vardict = dispatch_reader[product](**kwargs)
    return vardict
