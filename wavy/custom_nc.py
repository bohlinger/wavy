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

# specs
from buoy_specs import buoy_dict

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
        sHs = nc.variables['sHs'][:]
        mHs = nc.variables['mHs'][:]
        nc.close()
    return dtime,sHs,mHs

def get_arcmfc_stats(pathtofile):
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
            vardict[name]=var
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:],time_var.units)
        vardict['dtime']=dtime
        nc.close()
    return vardict

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

def dumptonc_ts_Tennholmen(outpath,filename,title,basetime,obs_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time_dt = obs_dict['time_dt']
    time = obs_dict['time_s']
    Hm0 = obs_dict['Hm0']
    Tm02 = obs_dict['Tm02']
    lons = obs_dict['lons']
    lats = obs_dict['lats']
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['Hm0'][startidx:endidx] = Hm0[:]
        nc.variables['Tm02'][startidx:endidx] = Tm02[:]
        nc.variables['longitude'][startidx:endidx] = lons[:]
        nc.variables['latitude'][startidx:endidx] = lats[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = "Tennholmen"
        nc.buoy_type = "Directional Waverider DWR MkIII"
        nc.buoy_specs = "http://www.datawell.nl/products/buoys.aspx"
        nc.buoy_manufacturer = "Datawell"
        nc.netcdf_version = "4"
        nc.data_owner = ("Norwegian Coastal Administration, " 
                        + "Institute of Marine Research, "
                        + "and Norwegian Meteorological Institute")
        nc.licence = ("Data and products are licensed under Norwegian" 
                    + "license for public data (NLOD) and " 
                    + "Creative Commons Attribution 3.0 Norway. "
                    + "See https://www.met.no/en/"
                    + "free-meteorological-data/Licensing-and-crediting")
        # dimensions
        dimsize = None
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
        ncHm0 = nc.createVariable(
                               'Hm0',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        ncTm02 = nc.createVariable(
                               'Tm02',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        nclons = nc.createVariable(
                               'longitude',
                               np.float64,
                               dimensions=('time'),
                               )
        nclats = nc.createVariable(
                               'latitude',
                               np.float64,
                               dimensions=('time'),
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time'
        nctime.long_name = 'Time of measurement'
        nctime.units = 'seconds since ' + str(basetime)
#        time.comment = "hourly values" ;
        nctime[:] = time
        # Hm0
        ncHm0.standard_name = 'sea_surface_wave_significant_height'
        ncHm0.long_name = 'Significant wave height estimate from spectrum'
        ncHm0.units = 'm'
        ncHm0.valid_range = 0., 25.
        ncHm0[:] = Hm0
        # Tm02
        ncTm02.standard_name = ('sea_surface_wave_mean_period'
                                + '_from_variance_spectral_density'
                                + '_second_frequency_moment')
        ncTm02.long_name = ('Mean wave period estimated from 0th' 
                            + 'and 2nd moment of spectrum')
        ncTm02.units = 's'
        ncTm02.valid_range = 0., 30.
        ncTm02[:] = Tm02
        # lons
        nclons.standard_name = ('longitude')
        nclons.units = 'degree_east'
        nclons.valid_min = -180.
        nclons.valid_max = 180.
        nclons[:] = lons
        # lats
        nclats.standard_name = ('latitude')
        nclats.units = 'degree_north'
        nclats.valid_min = -90.
        nclats.valid_max = 90.
        nclats[:] = lats
    nc.close()

def dumptonc_ts_Tennholmen_ext(outpath,filename,title,basetime,obs_dict):
    """
    Extended version of dumptonc_ts_Tennholmen
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time_dt = obs_dict['time_dt']
    time = obs_dict['time_s']
    Hm0 = obs_dict['Hm0']
    Tm02 = obs_dict['Tm02']
    lons = obs_dict['lons']
    lats = obs_dict['lats']
    Hs = obs_dict['Hs']
    TI = obs_dict['TI']
    TE = obs_dict['TE']
    T1 = obs_dict['T1']
    TZ = obs_dict['TZ']
    T3 = obs_dict['T3']
    Tc = obs_dict['Tc']
    Tdw = obs_dict['Tdw']
    Tp = obs_dict['Tp']
    Qp = obs_dict['Qp']
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['Hm0'][startidx:endidx] = Hm0[:]
        nc.variables['Tm02'][startidx:endidx] = Tm02[:]
        nc.variables['longitude'][startidx:endidx] = lons[:]
        nc.variables['latitude'][startidx:endidx] = lats[:]
        nc.variables['Hs'][startidx:endidx] = Hs[:]
        nc.variables['TI'][startidx:endidx] = TI[:]
        nc.variables['TE'][startidx:endidx] = TE[:]
        nc.variables['T1'][startidx:endidx] = T1[:]
        nc.variables['TZ'][startidx:endidx] = TZ[:]
        nc.variables['T3'][startidx:endidx] = T3[:]
        nc.variables['Tc'][startidx:endidx] = Tc[:]
        nc.variables['Tdw'][startidx:endidx] = Tdw[:]
        nc.variables['Tp'][startidx:endidx] = Tp[:]
        nc.variables['Qp'][startidx:endidx] = Qp[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = "Tennholmen"
        nc.buoy_type = "Directional Waverider DWR MkIII"
        nc.buoy_specs = "http://www.datawell.nl/products/buoys.aspx"
        nc.buoy_manufacturer = "Datawell"
        nc.netcdf_version = "4"
        nc.data_owner = ("Norwegian Coastal Administration, "
                        + "Institute of Marine Research, "
                        + "and Norwegian Meteorological Institute")
        nc.licence = ("Data and products are licensed under Norwegian"
                    + "license for public data (NLOD) and "
                    + "Creative Commons Attribution 3.0 Norway. "
                    + "See https://www.met.no/en/"
                    + "free-meteorological-data/Licensing-and-crediting")
        # dimensions
        dimsize = None
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
        ncHm0 = nc.createVariable(
                               'Hm0',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        ncTm02 = nc.createVariable(
                               'Tm02',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        nclons = nc.createVariable(
                               'longitude',
                               np.float64,
                               dimensions=('time'),
                               )
        nclats = nc.createVariable(
                               'latitude',
                               np.float64,
                               dimensions=('time'),
                               )
        ncHs = nc.createVariable(
                               'Hs',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        ncTI = nc.createVariable(
                               'TI',
                               np.float64,
                               dimensions=('time'),
                               )
        ncTE = nc.createVariable(
                               'TE',
                               np.float64,
                               dimensions=('time'),
                               )
        ncT1 = nc.createVariable(
                               'T1',
                               np.float64,
                               dimensions=('time'),
                               )
        ncTZ = nc.createVariable(
                               'TZ',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        ncT3 = nc.createVariable(
                               'T3',
                               np.float64,
                               dimensions=('time'),
                               )
        ncTc = nc.createVariable(
                               'Tc',
                               np.float64,
                               dimensions=('time'),
                               )
        ncTdw = nc.createVariable(
                               'Tdw',
                               np.float64,
                               dimensions=('time'),
                               )
        ncTp = nc.createVariable(
                               'Tp',
                               np.float64,
                               dimensions=('time'),
                               )
        ncQp = nc.createVariable(
                               'Qp',
                               np.float64,
                               dimensions=('time'),
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time'
        nctime.long_name = 'Time of measurement'
        nctime.units = 'seconds since ' + str(basetime)
#        time.comment = "hourly values" ;
        nctime[:] = time
        # Hm0
        ncHm0.standard_name = 'sea_surface_wave_significant_height'
        ncHm0.long_name = 'Significant wave height estimate from spectrum'
        ncHm0.units = 'm'
        ncHm0.valid_range = 0., 25.
        ncHm0[:] = Hm0
        # Tm02
        ncTm02.standard_name = ('sea_surface_wave_mean_period'
                                + '_from_variance_spectral_density'
                                + '_second_frequency_moment')
        ncTm02.long_name = ('Mean wave period estimated from 0th'
                            + 'and 2nd moment of spectrum')
        ncTm02.units = 's'
        ncTm02.valid_range = 0., 30.
        ncTm02[:] = Tm02
        # lons
        nclons.standard_name = ('longitude')
        nclons.units = 'degree_east'
        nclons.valid_min = -180.
        nclons.valid_max = 180.
        nclons[:] = lons
        # lats
        nclats.standard_name = ('latitude')
        nclats.units = 'degree_north'
        nclats.valid_min = -90.
        nclats.valid_max = 90.
        nclats[:] = lats
        # Hs
        ncHs.standard_name = ('Hs')
        ncHs.units = 'm'
        ncHs[:] = Hs
        # TI
        ncTI.standard_name = ('TI')
        ncTI.units = 's'
        ncTI[:] = TI
        # TE
        ncTE.standard_name = ('TE')
        ncTE.units = 's'
        ncTE[:] = TE
        # T1
        ncT1.standard_name = ('T1')
        ncT1.units = 's'
        ncT1[:] = T1
        # TZ
        ncTZ.standard_name = ('TZ')
        ncTZ.units = 's'
        ncTZ[:] = TZ
        # T3
        ncT3.standard_name = ('T3')
        ncT3.units = 's'
        ncT3[:] = T3
        # Tc
        ncTc.standard_name = ('Tc')
        ncTc.units = 's'
        ncTc[:] = Tc
        # Tdw
        ncTdw.standard_name = ('Tdw')
        ncTdw.units = 's'
        ncTdw[:] = Tdw
        # Tp
        ncTp.standard_name = ('Tp')
        ncTp.units = 's'
        ncTp[:] = Tp
        # Qp
        ncQp.standard_name = ('Qp')
        ncQp.units = 'None'
        ncQp[:] = Qp
    nc.close()

def dumptonc_coll_ts_Tennholmen(outpath,filename,title,basetime,obs_dict,model):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = obs_dict['time']
    Hm0_model = obs_dict['Hm0_model']
#    Tm02 = obs_dict['Tm02']
    lons_model = obs_dict['lons_model']
    lats_model = obs_dict['lats_model']
    Hm0_buoy = obs_dict['Hm0_buoy']
    lons_buoy = obs_dict['lons_buoy']
    lats_buoy = obs_dict['lats_buoy']
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['Hm0_model'][startidx:endidx] = Hm0_model[:]
#        nc.variables['Tm02'][startidx:endidx] = Tm02[:]
        nc.variables['longitude_model'][startidx:endidx] = lons_model[:]
        nc.variables['latitude_model'][startidx:endidx] = lats_model[:]
        nc.variables['Hm0_buoy'][startidx:endidx] = Hm0_buoy[:]
        nc.variables['longitude_buoy'][startidx:endidx] = lons_buoy[:]
        nc.variables['latitude_buoy'][startidx:endidx] = lats_buoy[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = "Tennholmen"
        nc.buoy_type = "Directional Waverider DWR MkIII"
        nc.buoy_specs = "http://www.datawell.nl/products/buoys.aspx"
        nc.buoy_manufacturer = "Datawell"
        nc.netcdf_version = "4"
        nc.data_owner = ("Norwegian Coastal Administration, "
                        + "Institute of Marine Research, "
                        + "and Norwegian Meteorological Institute")
        nc.licence = ("Data and products are licensed under Norwegian"
                    + "license for public data (NLOD) and "
                    + "Creative Commons Attribution 3.0 Norway. "
                    + "See https://www.met.no/en/"
                    + "free-meteorological-data/Licensing-and-crediting")
#        nc.processing_level = "Missing data has been filled with fillValue."
#        nc.processing_level = "No imputation for missing values."
        nc.static_position =  ("Latitude: "
                            + str(buoy_dict['Tennholmen']['lat']) 
                            + ", Longitude: " 
                            + str(buoy_dict['Tennholmen']['lon']))
        # dimensions
        dimsize = None
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
        ncHm0_model = nc.createVariable(
                               'Hm0_model',
                               np.float64,
                               dimensions=('time'),
                               fill_value=-999.,
                               )
#        ncTm02 = nc.createVariable(
#                               'Tm02',
#                               np.float64,
#                               dimensions=('time'),
#                               fill_value=9999.,
#                               )
        nclons_model = nc.createVariable(
                               'longitude_model',
                               np.float64,
                               dimensions=('time'),
                               )
        nclats_model = nc.createVariable(
                               'latitude_model',
                               np.float64,
                               dimensions=('time'),
                               )
        ncHm0_buoy = nc.createVariable(
                               'Hm0_buoy',
                               np.float64,
                               dimensions=('time'),
                               fill_value=9999.,
                               )
        nclons_buoy = nc.createVariable(
                               'longitude_buoy',
                               np.float64,
                               dimensions=('time'),
                               )
        nclats_buoy = nc.createVariable(
                               'latitude_buoy',
                               np.float64,
                               dimensions=('time'),
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time'
        nctime.long_name = 'Time of measurement'
        nctime.units = 'seconds since ' + str(basetime)
        nctime[:] = time
        # Hm0_model
        ncHm0_model.standard_name = (
                          'sea_surface_wave_significant_height '\
                        + 'from wave model'
                                    )
        ncHm0_model.long_name = ( 
                          'Significant wave height estimate '\
                        + 'from spectrum from wave model'
                                )
        ncHm0_model.units = 'm'
        ncHm0_model.valid_range = 0., 25.
        ncHm0_model[:] = Hm0_model
#        # Tm02
#        ncTm02.standard_name = ('sea_surface_wave_mean_period'
#                                + '_from_variance_spectral_density'
#                                + '_second_frequency_moment')
#        ncTm02.long_name = ('Mean wave period estimated from 0th'
#                            + 'and 2nd moment of spectrum')
#        ncTm02.units = 's'
#        ncTm02.valid_range = 0., 30.
#        ncTm02[:] = Tm02
        # lons_model
        nclons_model.standard_name = ('longitude_model')
        nclons_model.units = 'degree_east'
        nclons_model.valid_min = -180.
        nclons_model.valid_max = 180.
        nclons_model[:] = lons_model
        # lats_model
        nclats_model.standard_name = ('latitude_model')
        nclats_model.units = 'degree_north'
        nclats_model.valid_min = -90.
        nclats_model.valid_max = 90.
        nclats_model[:] = lats_model
        # Hm0_buoy
        ncHm0_buoy.standard_name = (
                          'sea_surface_wave_significant_height '\
                        + 'from buoy'
                                    )
        ncHm0_buoy.long_name = ( 'Significant wave height estimate'\
                                + 'from spectrum from buoy')
        ncHm0_buoy.units = 'm'
        ncHm0_buoy.valid_range = 0., 25.
        ncHm0_buoy[:] = Hm0_buoy
        # lons_buoy
        nclons_buoy.standard_name = ('longitude_buoy')
        nclons_buoy.units = 'degree_east'
        nclons_buoy.valid_min = -180.
        nclons_buoy.valid_max = 180.
        nclons_buoy[:] = lons_buoy
        # lats_buoy
        nclats_buoy.standard_name = ('latitude_buoy')
        nclats_buoy.units = 'degree_north'
        nclats_buoy.valid_min = -90.
        nclats_buoy.valid_max = 90.
        nclats_buoy[:] = lats_buoy
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
    msd = np.array(valid_dict['msd'])
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
        print(startidx)
        print(endidx)
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
