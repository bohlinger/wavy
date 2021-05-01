#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and write to netcdf 
files from model, station, or satellite output. I try to mostly follow 
the PEP convention for python code style. Constructive comments on style 
and effecient programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# read files
from netCDF4 import Dataset
import netCDF4

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os
import yaml

import sys

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/buoy_specs.yaml'))
with open(moddir,'r') as stream:
    buoy_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/station_specs.yaml'))
with open(moddir,'r') as stream:
    station_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/d22_var_dicts.yaml'))
with open(moddir,'r') as stream:
    d22_dict=yaml.safe_load(stream)

# --- global functions ------------------------------------------------#
"""
definition of some global functions
"""
# currently None
# ---------------------------------------------------------------------#

def get_nc_time(pathtofile):
    """
    timestep: "first" or "last" time step in nc-file
    pathtofile: complete path to file
    """
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
    else:
        nc = netCDF4.Dataset(
                    pathtofile,mode='r',
                    )
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:],time_var.units)
        nc.close()
    return dtime

def get_sat_alt_coll_var(pathtofile,varname):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        sys.exit('File does not exist')
    else:
        nc = netCDF4.Dataset(
            pathtofile,mode='r',
            )
        if varname == 'dtime':
            time_var = nc.variables['time']
            var = netCDF4.num2date(time_var[:],time_var.units)
        else:
            var = nc.variables[varname][:]
        nc.close()
    return var

def get_arcmfc_ts(pathtofile):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        sys.exit('File does not exist')
    else:
        nc = netCDF4.Dataset(
            pathtofile,mode='r',
            )
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:],time_var.units)
        sHs = nc.variables['obs_values'][:]
        mHs = nc.variables['model_values'][:]
        nc.close()
    return dtime,sHs,mHs

def get_arcmfc_stats(pathtofile):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        print('File does not exist')
        return
    else:
        nc = netCDF4.Dataset(
            pathtofile,mode='r',
            )
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:],time_var.units)
        mop = nc.variables['mop'][:]
        mor = nc.variables['mor'][:]
        rmsd = nc.variables['rmsd'][:]
        msd = nc.variables['msd'][:]
        corr = nc.variables['corr'][:]
        mad = nc.variables['mad'][:]
        bias = nc.variables['bias'][:]
        SI = nc.variables['SI'][:]
        nov = nc.variables['nov'][:]
        nc.close()
        valid_dict = {
            'mop':mop,
            'mor':mor,
            'msd':msd,
            'nov':nov,
            'rmsd':rmsd,
            'msd':msd,
            'corr':corr,
            'mad':mad,
            'bias':bias,
            'SI':SI}
        return valid_dict, dtime

def get_nc_ts(pathtofile,varlst):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        sys.exit('File does not exist')
    else:
        vardict = {}
        for name in varlst:
            nc = netCDF4.Dataset(
                pathtofile,mode='r',
                )
            var = nc.variables[name][:]
            vardict[name] = var
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:],time_var.units)
        vardict['dtime'] = dtime
        vardict['time'] = time_var[:]
        vardict['time_unit'] = time_var.units
        nc.close()
    return vardict

def get_nc_1D(pathtofile,varlst):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        sys.exit('File does not exist')
    else:
        vardict = {}
        for name in varlst:
            nc = netCDF4.Dataset(
                pathtofile,mode='r',
                )
            var = nc.variables[name][:]
            vardict[name]=var
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:],time_var.units)
        vardict['dtime']=dtime
        nc.close()
    return vardict

def dumptonc_ts(outpath,filename,title,model_time_unit,results_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time_dt = results_dict['date_matches']
    mHs = results_dict['model_matches']
    mlons = results_dict['model_lons_matches']
    mlats = results_dict['model_lats_matches']
    sHs = results_dict['sat_matches']
    slons = results_dict['sat_lons_matches']
    slats = results_dict['sat_lats_matches']
    dists = results_dict['dist_matches']
    
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False
                        )
        # create time vector in seconds since first date
        model_time_unit = nc.variables['time'].units
        time = netCDF4.date2num(time_dt, model_time_unit)
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['mHs'][startidx:endidx] = mHs[:]
        nc.variables['mlons'][startidx:endidx] = mlons[:]
        nc.variables['mlats'][startidx:endidx] = mlats[:]
        nc.variables['sHs'][startidx:endidx] = sHs[:]
        nc.variables['slons'][startidx:endidx] = slons[:]
        nc.variables['slats'][startidx:endidx] = slats[:]
        nc.variables['dists'][startidx:endidx] = dists[:]
    else:
        # create time vector in seconds since first date
        time = netCDF4.date2num(time_dt, model_time_unit)
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        )
        nc.title = title
        dimsize = None
        # dimensions
        dimtime = nc.createDimension(
                                'time',
                                size=dimsize
                                )
        # variables
        nctime = nc.createVariable(
                               'time',
                               np.float64,
                               dimensions=('time')
                               )
        ncmlats = nc.createVariable(
                               'mlats',
                               np.float64,
                               dimensions=('time')
                               )
        ncmlons = nc.createVariable(
                               'mlons',
                               np.float64,
                               dimensions=('time')
                               )
        ncmHs = nc.createVariable(
                               'mHs',
                               np.float64,
                               dimensions=('time')
                               )
        ncslats = nc.createVariable(
                               'slats',
                               np.float64,
                               dimensions=('time')
                               )
        ncslons = nc.createVariable(
                               'slons',
                               np.float64,
                               dimensions=('time')
                               )
        ncsHs = nc.createVariable(
                               'sHs',
                               np.float64,
                               dimensions=('time')
                               )
        ncdists = nc.createVariable(
                               'dists',
                               np.float64,
                               dimensions=('time')
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time matches'
        nctime.long_name = 'associated time steps between model and observation'
        nctime.units = model_time_unit
        nctime[:] = time
        # mHs
        ncmHs.standard_name = 'model Hs'
        ncmHs.long_name = 'significant wave height from wave model'
        ncmHs.units = 'm'
        ncmHs[:] = mHs
        # mlons
        ncmlons.standard_name = 'model lons'
        ncmlons.long_name = 'longitudes of associated model grid points'
        ncmlons.units = 'degrees east'
        ncmlons[:] = mlons
        # mlats
        ncmlats.standard_name = 'model lats'
        ncmlats.long_name = 'latitudes of associated model grid points'
        ncmlats.units = 'degrees north'
        ncmlats[:] = mlats
        # sHs
        ncsHs.standard_name = 'observed Hs'
        ncsHs.long_name = 'significant wave height from wave observation'
        ncsHs.units = 'm'
        ncsHs[:] = sHs
        # slons
        ncslons.standard_name = 'obs lons'
        ncslons.long_name = 'longitudes of observations'
        ncslons.units = 'degrees east'
        ncslons[:] = slons
        # slats
        ncslats.standard_name = 'obs lats'
        ncslats.long_name = 'latitudes of observations'
        ncslats.units = 'degrees north'
        ncslats[:] = slats
        # dists
        ncdists.standard_name = 'dists'
        ncdists.long_name = 'distances between observations and model grids'
        ncdists.units = 'km'
        ncdists[:] = dists
    nc.close()

def dumptonc_ts_station(st_obj,pathtofile,title):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file and folder structure
    """
    print('Dump data to netCDF4 file')
    stdvarname = st_obj.stdvarname
    time = st_obj.vars['time']
    lon = st_obj.vars['longitude']
    lat = st_obj.vars['latitude']
    var = np.array(st_obj.vars[stdvarname])
    var[var<variable_info[st_obj.varalias]['valid_range'][0]] = -999.
    var[var>variable_info[st_obj.varalias]['valid_range'][1]] = -999.
    var = list(var)
    print ('Dump data to file: ' + pathtofile)
    if os.path.isfile(pathtofile):
        nc = netCDF4.Dataset(
                        pathtofile,mode='a',
                        clobber=False
                        )
        # compare existing times in input time and existing time
        startin = time[0]
        timeex = list(nc.variables['time'][:])
        if startin in timeex:
            print('Time already detected in ncfile')
            print('Find correct index to start from there')
            print('Overwrite double time stamps')
            startidx = timeex.index(startin)
        else:
            startidx = len(nc['time'])
        endidx = startidx+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['longitude'][startidx:endidx] = lon[:]
        nc.variables['latitude'][startidx:endidx] = lat[:]
        nc.variables[st_obj.varname][startidx:endidx] = var[:]
        nc.close()
    else:
        outpath = os.path.dirname(pathtofile)
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        pathtofile,mode='w',
                        )
        # dimensions
        dimsize = None
        dimtime = nc.createDimension(
                                'time',
                                size=dimsize
                                )
        # variables
        nclon = nc.createVariable(
                               'longitude',
                               np.float64,
                               dimensions=('time')
                               )
        nclat = nc.createVariable(
                               'latitude',
                               np.float64,
                               dimensions=('time')
                               )
        nctime = nc.createVariable(
                               'time',
                               np.float64,
                               dimensions=('time')
                               )
        ncvar = nc.createVariable(
                               st_obj.varname,
                               np.float64,
                               dimensions=('time'),
                               fill_value=-999.
                               )
        # generate time for netcdf file
        # time
        nctime[:] = time
        nctime.units = str(st_obj.vars['time_unit'])
        nctime.setncatts(variable_info['time'])
        # longitude
        nclon[:] = lon
        nclon.setncatts(variable_info['lons'])
        # latitude
        nclat[:] = lat
        nclat.setncatts(variable_info['lats'])
        # var
        ncvar[:] = var
        ncvar.setncatts(variable_info[st_obj.varalias])
        # coordinate system info
        nc_crs = nc.createVariable('latlon',np.int32)
        nc_crs.proj4_string = "+proj=latlong +R=6370997.0 +ellps=WGS84"
        nc_crs.grid_mapping_name = 'latitude_longitude'
        # close file
        nc.close()
        #add global attributes
        nc = netCDF4.Dataset(pathtofile,mode='r+')
        nowstr = datetime.utcnow().isoformat()
        globalAttribs = {}
        globalAttribs['title'] = title
        globalAttribs['Conventions'] = "CF-1.6"
        globalAttribs['institution'] = \
                                "Norwegian Meteorological Institute"
        globalAttribs['history'] = nowstr + ". Created."
        globalAttribs['netcdf_version'] = "NETCDF4"
        if 'superob' in vars(st_obj).keys():
            globalAttribs['processing_level'] = 'post-processing performed' 
            globalAttribs['method_superobbing'] = str(st_obj.superob)
            globalAttribs['method_outlier_detection'] = (\
                                                st_obj.outlier_detection)
            globalAttribs['missing_value_handling'] = (st_obj.missing_data)
        else:
            globalAttribs['processing_level'] = \
                                "No post-processing performed"
        globalAttribs['static_position_station'] =  ("Latitude: "
                            + str(station_dict['platform']\
                                            [st_obj.platform]\
                                            ['coords']['lat'])
                            + ", Longitude: "
                            + str(station_dict['platform']\
                                            [st_obj.platform]\
                                            ['coords']['lon']))
        globalAttribs['platform_operator'] = station_dict['platform']\
                                           [st_obj.platform]\
                                           ['operator']
        nc.setncatts(globalAttribs)
        nc.sync()
        nc.close()

def dumptonc_ts_collocation(col_obj,pathtofile,title):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file and folder structure
    """
    print('Dump data to netCDF4 file')
    stdvarname = col_obj.stdvarname
    time = col_obj.vars['time']
    modlon = col_obj.vars['model_lons']
    modlat = col_obj.vars['model_lats']
    obslon = col_obj.vars['obs_lons']
    obslat = col_obj.vars['obs_lats']
    colidx_x = col_obj.vars['collocation_idx_x']
    colidx_y = col_obj.vars['collocation_idx_y']
    leadtime = col_obj.leadtime
    varobs = np.array(col_obj.vars['obs_values'])
    varobs[varobs<variable_info[col_obj.varalias]['valid_range'][0]] = -999.
    varobs[varobs>variable_info[col_obj.varalias]['valid_range'][1]] = -999.
    varobs = list(varobs)
    varmod = np.array(col_obj.vars['model_values'])
    varmod[varmod<variable_info[col_obj.varalias]['valid_range'][0]] = -999.
    varmod[varmod>variable_info[col_obj.varalias]['valid_range'][1]] = -999.
    varmod = list(varmod)
    dists = np.array(col_obj.vars['distance'])*1000 # from km to m
    print ('Dump data to file: ' + pathtofile)
    if os.path.isfile(pathtofile):
        nc = netCDF4.Dataset(
                        pathtofile,mode='a',
                        clobber=False
                        )
        # compare existing times in input time and existing time
        startin = time[0]
        timeex = list(nc.variables['time'][:])
        if startin in timeex:
            print('Time already detected in ncfile')
            print('Find correct index to start from there')
            print('Overwrite double time stamps')
            startidx = timeex.index(startin)
        else:
            startidx = len(nc['time'])
        endidx = startidx+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['model_lons'][startidx:endidx] = modlon[:]
        nc.variables['model_lats'][startidx:endidx] = modlat[:]
        nc.variables['model_values'][startidx:endidx] = varmod[:]
        nc.variables['obs_lons'][startidx:endidx] = obslon[:]
        nc.variables['obs_lats'][startidx:endidx] = obslat[:]
        nc.variables['obs_values'][startidx:endidx] = varobs[:]
        nc.variables['dist'][startidx:endidx] = dists[:]
        nc.variables['colidx'][startidx:endidx] = colidx[:]
        nc.close()
    else:
        outpath = os.path.dirname(pathtofile)
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        pathtofile,mode='w',
                        )
        # dimensions
        dimsize = None
        dimtime = nc.createDimension(
                                'time',
                                size=dimsize
                                )
        # variables
        ncmodlon = nc.createVariable(
                               'model_lons',
                               np.float64,
                               dimensions=('time')
                               )
        ncmodlat = nc.createVariable(
                               'model_lats',
                               np.float64,
                               dimensions=('time')
                               )
        ncobslon = nc.createVariable(
                               'obs_lons',
                               np.float64,
                               dimensions=('time')
                               )
        ncobslat = nc.createVariable(
                               'obs_lats',
                               np.float64,
                               dimensions=('time')
                               )
        nccolidx_x = nc.createVariable(
                               'colidx_x',
                               np.float64,
                               dimensions=('time')
                               )
        nccolidx_y = nc.createVariable(
                               'colidx_y',
                               np.float64,
                               dimensions=('time')
                               )
        nctime = nc.createVariable(
                               'time',
                               np.float64,
                               dimensions=('time')
                               )
        ncvarobs = nc.createVariable(
                               'obs_values',
                               np.float64,
                               dimensions=('time'),
                               fill_value=-999.
                               )
        ncvarmod = nc.createVariable(
                               'model_values',
                               np.float64,
                               dimensions=('time'),
                               fill_value=-999.
                               )
        ncdist = nc.createVariable(
                               'dist',
                               np.float64,
                               dimensions=('time'),
                               fill_value=-999.
                               )
        # generate time for netcdf file
        # time
        nctime[:] = time
        nctime.units = str(col_obj.vars['time_unit'])
        nctime.setncatts(variable_info['time'])
        # observation longitude
        ncobslon[:] = obslon
        ncobslon.setncatts(variable_info['lons'])
        # observation latitude
        ncobslat[:] = obslat
        ncobslat.setncatts(variable_info['lats'])
        # model longitude
        ncmodlon[:] = modlon
        ncmodlon.setncatts(variable_info['lons'])
        # model latitude
        ncmodlat[:] = modlat
        ncmodlat.setncatts(variable_info['lats'])
        # varobs
        ncvarobs[:] = varobs
        ncvarobs.setncatts(variable_info[col_obj.varalias])
        ncvarobs.observation_name = col_obj.obsname
        # varmod
        ncvarmod[:] = varmod
        ncvarmod.setncatts(variable_info[col_obj.varalias])
        ncvarmod.model_name = col_obj.model
        # dists
        ncdist[:] = dists
        ncdist.setncatts(variable_info['dist'])
        # colidx
        nccolidx_x[:] = colidx_x
        nccolidx_x.setncatts(variable_info['colidx_x'])
        nccolidx_y[:] = colidx_y
        nccolidx_y.setncatts(variable_info['colidx_y'])
        # coordinate system info
        nc_crs = nc.createVariable('latlon',np.int32)
        nc_crs.proj4_string = "+proj=latlong +R=6370997.0 +ellps=WGS84"
        nc_crs.grid_mapping_name = 'latitude_longitude'
        # close file
        nc.close()
        #add global attributes
        nc = netCDF4.Dataset(pathtofile,mode='r+')
        nowstr = datetime.utcnow().isoformat()
        globalAttribs = {}
        globalAttribs['title'] = title
        globalAttribs['Conventions'] = "CF-1.6"
        globalAttribs['institution'] = \
                                "Norwegian Meteorological Institute"
        globalAttribs['history'] = nowstr + ". Created."
        globalAttribs['netcdf_version'] = "NETCDF4"
        if 'superob' in vars(col_obj).keys():
            globalAttribs['processing_level'] = ('superob: '
                                                + str(col_obj.superob)
                                                + '; '
                                                + 'outlier_detection: '
                                                + col_obj.outlier_detection
                                                + '; '
                                                + 'missing_data: '
                                                + col_obj.missing_data)
        else:
            globalAttribs['processing_level'] = \
                                "No post-processing performed"
        globalAttribs['leadtime'] = str(leadtime) + 'h'
        nc.setncatts(globalAttribs)
        nc.sync()
        nc.close()

def dumptonc_stats(pathtofile,title,time_dt,time_unit,valid_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    # create time vector in seconds since first date
    time = netCDF4.date2num(time_dt,time_unit)
    mop = np.array(valid_dict['mop'])
    mor = np.array(valid_dict['mor'])
    rmsd = np.array(valid_dict['rmsd'])
    msd = np.array(valid_dict['msd'])
    corr = np.array(valid_dict['corr'])
    mad = np.array(valid_dict['mad'])
    bias = np.array(valid_dict['bias'])
    SI = np.array(valid_dict['SI'][1])
    nov = np.array(valid_dict['nov'])
    print ('Dump data to file: ' + pathtofile)
    if os.path.isfile(pathtofile):
        nc = netCDF4.Dataset(
                        pathtofile,mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+1
        nc.variables['time'][startidx:endidx] = time
        nc.variables['mop'][startidx:endidx] = mop
        nc.variables['mor'][startidx:endidx] = mor
        nc.variables['rmsd'][startidx:endidx] = rmsd
        nc.variables['msd'][startidx:endidx] = msd
        nc.variables['corr'][startidx:endidx] = corr
        nc.variables['mad'][startidx:endidx] = mad
        nc.variables['bias'][startidx:endidx] = bias
        nc.variables['SI'][startidx:endidx] = SI
        nc.variables['nov'][startidx:endidx] = nov
    else:
        outpath = pathtofile[0:-len(pathtofile.split('/')[-1])]
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        pathtofile,mode='w',
#                        format='NETCDF4'
                        )
        nc.title = title
        dimsize = None
        # dimensions
        dimtime = nc.createDimension(
                                'time',
                                size=dimsize
                                )
        # variables
        nctime = nc.createVariable(
                               'time',
                               np.float64,
                               dimensions=('time')
                               )
        ncmop = nc.createVariable(
                               'mop',
                               np.float64,
                               dimensions=('time')
                               )
        ncmor = nc.createVariable(
                               'mor',
                               np.float64,
                               dimensions=('time')
                               )
        ncrmsd = nc.createVariable(
                               'rmsd',
                               np.float64,
                               dimensions=('time')
                               )
        ncmsd = nc.createVariable(
                               'msd',
                               np.float64,
                               dimensions=('time')
                               )
        nccorr = nc.createVariable(
                               'corr',
                               np.float64,
                               dimensions=('time')
                               )
        ncmad = nc.createVariable(
                               'mad',
                               np.float64,
                               dimensions=('time')
                               )
        ncbias = nc.createVariable(
                               'bias',
                               np.float64,
                               dimensions=('time')
                               )
        ncSI = nc.createVariable(
                               'SI',
                               np.float64,
                               dimensions=('time')
                               )
        ncnov = nc.createVariable(
                               'nov',
                               np.float64,
                               dimensions=('time')
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time matches'
        nctime.long_name = 'associated time steps between model and observation'
        nctime.units = time_unit
        nctime[:] = time
        # mop
        ncmop.standard_name = 'mop'
        ncmop.long_name = 'mean of product (wave model)'
        ncmop.units = 'm'
        ncmop[:] = mop
        # mor
        ncmor.standard_name = 'mor'
        ncmor.long_name = 'mean of reference (observations)'
        ncmor.units = 'm'
        ncmor[:] = mor
        # rmsd
        ncrmsd.standard_name = 'rmsd'
        ncrmsd.long_name = 'root mean square deviation'
        ncrmsd.units = 'm'
        ncrmsd[:] = rmsd
        # msd
        ncmsd.standard_name = 'msd'
        ncmsd.long_name = 'mean square deviation'
        ncmsd.units = 'm^2'
        ncmsd[:] = msd
        # corr
        nccorr.standard_name = 'corr'
        nccorr.long_name = 'correlation coefficient'
        nccorr.units = 'none'
        nccorr[:] = corr
        # mad
        ncmad.standard_name = 'mad'
        ncmad.long_name = 'mean absolute deviation'
        ncmad.units = 'm'
        ncmad[:] = mad
        # bias
        ncbias.standard_name = 'bias'
        ncbias.long_name = 'Bias (mean error)'
        ncbias.units = 'm'
        ncbias[:] = bias
        # SI
        ncSI.standard_name = 'SI'
        ncSI.long_name = 'scatter index'
        ncSI.units = 'none'
        ncSI[:] = SI
        # nov
        ncnov.standard_name = 'nov'
        ncnov.long_name = 'number of values'
        ncnov.units = 'none'
        ncnov[:] = nov
    nc.close()

def dumptonc_sat(sa_obj,outpath,mode=None):
    """
    dump satellite altimetry data to netcdf-file
    """
    sdate=sa_obj.sdate
    edate=sa_obj.edate
    filename = (sa_obj.sat
                + "_"
                + sa_obj.region
                + "_"
                + sdate.strftime("%Y%m%d%H%M%S")
                + "_"
                + edate.strftime("%Y%m%d%H%M%S")
                + ".nc")
    fullpath = outpath + filename
    os.system('mkdir -p ' + outpath)
    print ('Dump altimeter wave data from '
            + sa_obj.sat
            + ' to file: ' + fullpath)
    nc = netCDF4.Dataset(
                    fullpath,mode='w',
                    )
    nc.title = sa_obj.sat + ' altimeter significant wave height'
    timerange=len(sa_obj.vars['time'])
    dimsize = None
    # dimensions
    dimtime = nc.createDimension(
                            'time',
                            size=dimsize
                            )
    # variables
    nctime = nc.createVariable(
                           'time',
                           np.float64,
                           dimensions=('time')
                           )
    nclatitude = nc.createVariable(
                           'latitude',
                           np.float64,
                           dimensions=('time')
                           )
    nclongitude = nc.createVariable(
                           'longitude',
                           np.float64,
                           dimensions=('time')
                           )
    ncHs = nc.createVariable(
                           'Hs',
                           np.float64,
                           dimensions=('time')
                           )

    # generate time for netcdf file
    nctime.units = sa_obj.vars['time_unit']
    nctime[:] = sa_obj.vars['time']
    nctime.standard_name = 'time'
    ncHs.units = 'm'
    ncHs[:] = sa_obj.vars['sea_surface_wave_significant_height']
    ncHs.standard_name = 'sea_surface_wave_significant_height'
    ncHs.long_name = \
        'Significant wave height estimate from altimeter wave form'
    ncHs.valid_range = 0., 25.
    nclongitude.units = 'degree_east'
    nclongitude[:] = sa_obj.vars['longitude']
    nclongitude.standard_name = 'longitude'
    nclongitude.valid_min = -180.
    nclongitude.valid_max = 180.
    nclatitude[:] = sa_obj.vars['latitude']
    nclatitude.standard_name = 'latitude'
    nclatitude.units = 'degree_north'
    nclatitude.valid_min = -90.
    nclatitude.valid_max = 90.
    nc.close()

def dumptonc_pointsat(sa_obj,outpath,mode=None):
    """
    dump satellite altimetry data to netcdf-file
    """
    sdate=sa_obj.sdate
    edate=sa_obj.edate
    filename = (sa_obj.sat
                + "_"
                + sa_obj.region
                + "_"
                + sdate.strftime("%Y%m%d%H%M%S")
                + "_"
                + edate.strftime("%Y%m%d%H%M%S")
                + ".nc")
    fullpath = outpath + filename
    os.system('mkdir -p ' + outpath)
    print ('Dump altimeter wave data from '
            + sa_obj.sat
            + ' to file: ' + fullpath)
    nc = netCDF4.Dataset(
                    fullpath,mode='w',
                    )
    nc.title = (sa_obj.sat + 
                ' altimeter significant wave height close to ' 
                + sa_obj.region)
    timerange=len(sa_obj.Hs)
    dimsize = None
    # dimensions
    dimtime = nc.createDimension(
                            'time',
                            size=dimsize
                            )
    # variables
    nctime = nc.createVariable(
                           'time',
                           np.float64,
                           dimensions=('time')
                           )
    nclatitude = nc.createVariable(
                           'latitude',
                           np.float64,
                           dimensions=('time')
                           )
    nclongitude = nc.createVariable(
                           'longitude',
                           np.float64,
                           dimensions=('time')
                           )
    ncHs = nc.createVariable(
                           'Hs',
                           np.float64,
                           dimensions=('time')
                           )
    ncdist = nc.createVariable(
                           'dist',
                           np.float64,
                           dimensions=('time')
                           )
    # generate time for netcdf file
    basetime=sa_obj.basetime
    nctime.units = 'seconds since 2000-01-01 00:00:00'
    nctime[:] = sa_obj.time
    ncHs[:] = sa_obj.Hs
    ncHs.units = 'm'
    ncHs.standard_name = 'sea_surface_wave_significant_height'
    ncHs.long_name = \
        'Significant wave height estimate from altimeter wave form'
    ncHs.valid_range = 0., 25.
    nclongitude.units = 'degree_east'
    nclongitude[:] = sa_obj.loc[1]
    nclongitude.standard_name = 'longitude'
    nclongitude.valid_min = -180.
    nclongitude.valid_max = 180.
    nclatitude[:] = sa_obj.loc[0]
    nclatitude.standard_name = 'latitude'
    nclatitude.units = 'degree_north'
    nclatitude.valid_min = -90.
    nclatitude.valid_max = 90.
    ncdist.units = 'km'
    ncdist[:] = sa_obj.dist
    ncdist.long_name = ('distance from footprint ' 
                    + 'to location according '
                    + 'to haversine')
    nc.close()

def dumptonc_ts_pos(outpath,filename,title,coll_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    # extract dict
    model = coll_dict['model']
    varname = coll_dict['varname']
    basetime = coll_dict['basetime']
    time = coll_dict['time']
    var_model = coll_dict[varname]
    lons_model = coll_dict['lons_model']
    lats_model = coll_dict['lats_model']
    lons_pos = coll_dict['lons_pos']
    lats_pos = coll_dict['lats_pos']
    dist = coll_dict['hdist']
    idx = coll_dict['idx']
    idy = coll_dict['idy']
    # writing/appending
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(fullpath,mode='a',clobber=False)
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables[varname][startidx:endidx] = var_model[:]
    else:
        os.system('mkdir -p ' + outpath)
        # create nc-file
        nc = netCDF4.Dataset(fullpath,mode='w')
        # create dimension time
        dimtime = nc.createDimension('time',size=None)
        # add time
        nctime = nc.createVariable('time',np.float64,dimensions=('time'))
        nctime.standard_name = 'time'
        nctime.units = 'seconds since ' + str(basetime)
        nctime[:] = time
        # coordinate system info
        nc_crs = nc.createVariable('latlon',np.int32)
        nc_crs.proj4_string = "+proj=latlong +R=6370997.0 +ellps=WGS84"
        nc_crs.grid_mapping_name = 'latitude_longitude'
        # close file
        nc.close()
        #add global attributes
        nc = netCDF4.Dataset(fullpath,mode='r+')
        nowstr = datetime.utcnow().isoformat()
        globalAttribs = {}
        globalAttribs['title'] = title
        globalAttribs['Conventions'] = "CF-1.6"
        globalAttribs['institution'] = \
                                "Norwegian Meteorological Institute"
        globalAttribs['history'] = nowstr + ". Created."
        globalAttribs['netcdf_version'] = "NETCDF4"
        globalAttribs['processing_level'] = \
                                "No post-processing performed"
        globalAttribs['static_position_station'] =  ("Latitude: "
                                + "{:.4f}".format(lats_pos[0])
                                + ", Longitude: "
                                + "{:.4f}".format(lons_pos[0]))
        globalAttribs['static_position_model'] =  ("Latitude: "
                                + "{:.4f}".format(lats_model[0])
                                + ", Longitude: "
                                + "{:.4f}".format(lons_model[0]))
        globalAttribs['static_collocation_idx'] =  ("idx: "
                                + str(idx[0])
                                + ", idy: "
                                + str(idy[0]))
        globalAttribs['static_collocation_distance'] =  \
                                ("{:.4f}".format(dist[0]) + " km")
        nc.setncatts(globalAttribs)
        nc.sync()
        nc.close()
        # append all other variables
        for varstr in coll_dict:
            if varstr in [varname]:
                nc = netCDF4.Dataset(fullpath,mode='r+')
                ncvar = nc.createVariable(varstr,
                            np.float64,dimensions=('time'))
                # add variable attributes
                varAttribs = {}
                varAttribs['standard_name'] = variable_info[varname]\
                                                ['standard_name']
                varAttribs['units'] = variable_info[varname]['units']
                varAttribs['valid_range'] = variable_info[varname]\
                                                ['valid_range'][0], \
                                            variable_info[varname]\
                                                ['valid_range'][1]
                varAttribs['convention'] = variable_info[varname]\
                                                    ['convention']
                ncvar.setncatts(varAttribs)
                ncvar[:] = coll_dict[varstr][:]
                nc.close()

def check_vals_in_nc(filestr,varname,pytime_in):
    print('check for time: ', pytime_in)
    if os.path.exists(filestr):
        nc = netCDF4.Dataset(filestr,mode='r')
        var = nc.variables[varname][:]
        time = nc.variables['time'][:]
        unit = nc.variables['time'].units
        pytime_file = netCDF4.num2date(time,units = unit)
        try:
            idx = list(pytime_file).index(pytime_in)
        except ValueError:
            idx = None
        nc.close()
    else:
        idx = None
    return idx

def ncdump(nc_fid, verb=True):
    '''
    Function from:
    http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
    #
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.
    #
    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed
    #
    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print ("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print ('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print ("WARNING: %s does not contain variable attributes" % key)
    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print ("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print ('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print ("NetCDF dimension information:")
        for dim in nc_dims:
            print ("\tName:", dim)
            print ("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print ("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print ('\tName:', var)
                print ("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print ("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

def ncdumpMeta(pathtofile):
    '''
    Returns dict of netcdf-file content
    Input: str pointing to netcdf-file
    Output: dict of attributes
    '''
    # init netCDF4 instance
    nc = netCDF4.Dataset(pathtofile,mode='r')
    # init empty dict
    ncdict = {}
    # retrieve variable attributes
    for var in nc.variables:
        ncattrs = nc.variables[var].ncattrs()
        ncdict[var] = {}
        for ncattr in ncattrs:
            ncdict[var][ncattr] = nc.variables[var].getncattr(ncattr)
    # retrieve global attributes
    ncdict['global'] = {}
    nc_attrs = nc.ncattrs()
    for nc_attr in nc_attrs:
        ncdict['global'][nc_attr] = nc.getncattr(nc_attr)
    return ncdict

def find_attr_in_nc(attrstr,pathtofile=None,ncdict=None,subattrstr=None):
    """
    fct to find a specific attribute with its value in a netcdf-file
    when only a fraction of attribute name is given, can also search 
    in a deeper layer of attribute hierarchy.
    input:  path to the nc-file or ncdict
            string of desired attribute
            string of desired sub-attribute (optional)
    output: dict
    """
    if pathtofile is not None and ncdict is None:
        # get content of netcdf-file
        ncdict = ncdumpMeta(pathtofile)
    # browse for str using list comprehension
    res1 = [i for i in ncdict.keys() if attrstr in i]
    if subattrstr is None:
        return ncdict[res1[0]]
    else:
        res2 = [i for i in ncdict[res1[0]] if subattrstr in i]
        return ncdict[res1[0]][res2[0]]

def get_varname_for_cf_stdname_in_ncfile(ncdict,stdname):
    lst = [i for i in ncdict.keys() \
            if ('standard_name' in ncdict[i].keys() \
            and stdname in ncdict[i]['standard_name']) \
            ]
    if len(lst) >= 1:
        return lst
    else: return None
