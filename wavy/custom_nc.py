#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and write to netcdf 
files from model, station, or sentinel output. I try to mostly follow 
the PEP convention for python code style. Constructive comments on style 
and effecient programming are most welcome!

I wish it to work something like this in future:
 1: get_model for given time period
 2: dumptonc based on model (e.g. MWAM4, ARCMFC, ARCMFCnew)
 3: choose create or append based on the existence of the file
    Must have one unlimited dimension (time), and two spatial dimensions
    (latitude, longitude, which depend on rlat,rlon)
'''
__version__ = "0.5.0"
__author__="Patrik Bohlinger, Norwegian Meteorological Institute"
__maintainer__ = "Patrik Bohlinger"
__email__ = "patrikb@met.no"
__status__ = "under development with operation ARCMFC branch"

# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# ignore irrelevant warnings from matplotlib for stdout
import warnings
warnings.filterwarnings("ignore")

# plotting
import matplotlib.pyplot as plt

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

# progress bar
import sys

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# get necessary paths for module
import pathfinder

# --- global functions ------------------------------------------------#
"""
definition of some global functions
"""
# Currently None
# ---------------------------------------------------------------------#


class custom_nc():
    '''
    class to write to netcdf files from sentinel, station, or model data
    sentinel: level 3 data i.e. Hs[time], lat[time], lon[time] 
    station: e.g. Hs[time], lat, lon
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the sentinel, model, and 
    station classes.
    '''
    satpath_lustre = pathfinder.satpath_lustre
    satpath_copernicus = pathfinder.satpath_copernicus
    satpath_ftp_008_052 = pathfinder.satpath_ftp_008_052
    satpath_ftp_014_001 = pathfinder.satpath_ftp_014_001
    
    from region_specs import region_dict
    from model_specs import model_dict

    def __init__(self,sdate,edate=None,model=None,timewin=None,region=None):
        print ('# ----- ')
        print (" ### Initializing custom_nc instance ###")
        print ('# ----- ')
        if region is None:
            model='ARCMFC'
        if timewin is None:
            timewin = int(30)
        if edate is None:
            edate=sdate
            if timewin is None:
                timewin = int(30)
            print ("Requested time: ", str(sdate))
            print (" with timewin: ", str(timewin))
        else:
            print ("Requested time frame: " +
                str(sdate) + " - " + str(edate))
        self.sdate = sdate
        self.edate = edate
        self.model = model
        self.basetime = model_dict[model]['basetime']

def get_nc_time(timestep,pathtofile):
    """
    timestep: "first" or "last" time step in nc-file
    pathtofile: complete path to file
    """
    nc = netCDF4.Dataset(
                    pathtofile,mode='r',
                    clobber=False
                    )
    time_var = nc.variables['time']
    dtime = netCDF4.num2date(time_var[:],time_var.units)
    return dtime[timestep]

def dumptonc_ts(outpath,filename,title,basetime,results_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time_dt = results_dict['date_matches']
    # create time vector in seconds since first date
    time = []
    for dt in time_dt:
        time.append((dt-basetime).total_seconds())
    time = np.array(time)
    mHs = results_dict['model_Hs_matches']
    mlons = results_dict['model_lons_matches']
    mlats = results_dict['model_lats_matches']
    sHs = results_dict['sat_Hs_matches']
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
        # variables
        #ncrtime=nc.variables['rtime'][:]
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
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        format='NETCDF4'
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
        nctime.units = 'seconds since ' + str(basetime)
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

def dumptonc_stats(outpath,filename,title,basetime,time_dt,valid_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    # create time vector in seconds since first date
    time = np.array((time_dt-basetime).total_seconds())
    mop = np.array(valid_dict['mop'])
    mor = np.array(valid_dict['mor'])
    rmsd = np.array(valid_dict['rmsd'])
    corr = np.array(valid_dict['corr'])
    mad = np.array(valid_dict['mad'])
    bias = np.array(valid_dict['bias'])
    SI = np.array(valid_dict['SI'][1])
    nov = np.array(valid_dict['nov'])
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+1
        #print(time)
        #startidx = len(nc['time'])-1
        print(startidx)
        #endidx = len(nc['time'])
        print(endidx)
        nc.variables['time'][startidx:endidx] = time
        nc.variables['mop'][startidx:endidx] = mop
        nc.variables['mor'][startidx:endidx] = mor
        nc.variables['rmsd'][startidx:endidx] = rmsd
        nc.variables['corr'][startidx:endidx] = corr
        nc.variables['mad'][startidx:endidx] = mad
        nc.variables['bias'][startidx:endidx] = bias
        nc.variables['SI'][startidx:endidx] = SI
        nc.variables['nov'][startidx:endidx] = nov
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        format='NETCDF4'
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
        nctime.units = 'seconds since ' + str(basetime)
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
