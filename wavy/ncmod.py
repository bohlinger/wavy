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

# progress bar
import sys

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# read yaml config files:
with open("../wavy/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)
with open("../wavy/buoy_specs.yaml", 'r') as stream:
    buoy_dict=yaml.safe_load(stream)
with open("../wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)
with open("../wavy/variable_info.yaml", 'r') as stream:
    var_dict=yaml.safe_load(stream)
with open("../wavy/pathfinder.yaml", 'r') as stream:
    pathfinder=yaml.safe_load(stream)

# --- global functions ------------------------------------------------#
"""
definition of some global functions
"""
# currently None
# ---------------------------------------------------------------------#


class ncmod():
    '''
    class to write to netcdf files from satellite, station, or model data
    satellite: level 3 data i.e. Hs[time], lat[time], lon[time] 
    station: e.g. Hs[time], lat, lon
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and 
    station classes.
    '''
    satpath_lustre = pathfinder['satpath_lustre']
    satpath_copernicus = pathfinder['satpath_copernicus']
    satpath_ftp_014_001 = pathfinder['satpath_ftp_014_001']
    
    def __init__(self,sdate,edate=None,model=None,timewin=None,region=None):
        print ('# ----- ')
        print (" ### Initializing ncmod instance ###")
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
#                        format='NETCDF4'
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
#                        format='NETCDF4'
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
#                        format='NETCDF4'
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

def dumptonc_coll_ts_buoy(outpath,filename,title,basetime,obs_dict,model):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = obs_dict['time']
    Hm0_model = obs_dict['Hm0_model']
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
        nc.variables['longitude_model'][startidx:endidx] = lons_model[:]
        nc.variables['latitude_model'][startidx:endidx] = lats_model[:]
        nc.variables['Hm0_buoy'][startidx:endidx] = Hm0_buoy[:]
        nc.variables['longitude_buoy'][startidx:endidx] = lons_buoy[:]
        nc.variables['latitude_buoy'][startidx:endidx] = lats_buoy[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
#                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = obs_dict['buoy']
        nc.buoy_type = buoy_dict[obs_dict['buoy']]['type']
#        nc.buoy_specs = "http://www.datawell.nl/products/buoys.aspx"
        nc.buoy_manufacturer = buoy_dict[obs_dict['buoy']]['manufacturer']
        nc.netcdf_version = "4"
        nc.data_owner = buoy_dict[obs_dict['buoy']]['data_owner']
        nc.licence = buoy_dict[obs_dict['buoy']]['licence']
        nc.processing_level = "No imputation for missing values."
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
    del nc

def dumptonc_coll_ts_Tp_station(outpath,filename,title,basetime,\
                        obs_dict,model,statname,sensorname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = np.array(obs_dict['time'])
    Tp_model = np.array(obs_dict['Tp_model'])
    lons_model = np.array(obs_dict['lons_model'])
    lats_model = np.array(obs_dict['lats_model'])
    Tp_stat_10min = np.array(obs_dict['Tp_stat_10min'])
    idx = np.array(obs_dict['idx'])
    idy = np.array(obs_dict['idy'])
    hdist = np.array(obs_dict['hdist'])
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    if os.path.isfile(fullpath):
        print(fullpath)
        print('--- file exists, appending to nc-file ---')
        nc = netCDF4.Dataset(
                        fullpath,mode='a',
                        clobber=False,
#                        format='NETCDF4'
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+len(time)
        nc.variables['time'][startidx:endidx] = time[:]
        nc.variables['Tp_model'][startidx:endidx] = Tp_model[:]
        nc.variables['Tp_station_10min'][startidx:endidx] = Tp_stat_10min[:]
    else:
        print(fullpath)
        print('--- file does not exist, creating new nc-file ---')
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
#                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = statname
        nc.instrument_type = sensorname
        nc.instrument_specs = "NA"
        nc.instrument_manufacturer = \
                        station_dict[statname]['manufacturer'][sensorname]
        nc.netcdf_version = "4"
        nc.data_owner = ("NA")
        nc.licence = ("NA")
        nc.processing_level = "No imputation for missing or erroneous values."
        nc.static_position_station =  ("Latitude: "
                            + str(station_dict[statname]['coords']['lat'])
                            + ", Longitude: "
                            + str(station_dict[statname]['coords']['lon']))
        nc.static_position_model =  ("Latitude: "
                            + "{:.2f}".format(lats_model[0])
                            + ", Longitude: "
                            + "{:.2f}".format(lons_model[0]))
        nc.collocation_distance = "{:.2f}".format(hdist[0][0]) + " km"
        nc.static_collocation_idx =  ("idx: "
                            + str(idx[0])
                            + ", idy: "
                            + str(idy[0]))
        # dimensions
        dimsize = None
        #dimsize = len(time)
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
        ncTp_model = nc.createVariable(
                               'Tp_model',
                               np.float64,
                               dimensions=('time')
                               )
        ncTp_stat_10min = nc.createVariable(
                               'Tp_station_10min',
                               np.float64,
                               dimensions=('time'),
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time'
        nctime.long_name = 'Time of measurement'
        nctime.units = ('seconds since ' + str(basetime))
        nctime[:] = time
        # Tp_model
        ncTp_model.standard_name = (
                          'Wave_Peak_Period'
                                    )
        ncTp_model.long_name = (
                          'Wave Peak Period '\
                        + 'from wave model'
                                )
        ncTp_model.units = 's'
        ncTp_model.valid_range = 0., 30.
        ncTp_model[:] = Tp_model
        # Tp_stat_10min
        ncTp_stat_10min.standard_name = (
                          'Wave_Peak_Period'
                                    )
        ncTp_stat_10min.long_name = ( 'Wave Peak Period estimate '\
                                + 'from ' + statname)
        ncTp_stat_10min.units = 's'
        ncTp_stat_10min.valid_range = 0., 30.
        ncTp_stat_10min[:] = Tp_stat_10min
    nc.close()

def dumptonc_coll_ts_station(outpath,filename,title,basetime,\
                        obs_dict,model,statname,sensorname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = obs_dict['time']
    Hm0_model = obs_dict['Hm0_model']
    lons_model = obs_dict['lons_model']
    lats_model = obs_dict['lats_model']
    Hs_stat_1h = obs_dict['Hs_stat_1h']
    Hs_stat_10min = obs_dict['Hs_stat_10min']
    idx = obs_dict['idx']
    idy = obs_dict['idy']
    hdist = obs_dict['hdist']
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
        nc.variables['Hs_station_1h'][startidx:endidx] = Hs_stat_1h[:]
        nc.variables['Hs_station_10min'][startidx:endidx] = Hs_stat_10min[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
#                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = statname
        nc.instrument_type = sensorname
        nc.instrument_specs = "NA"
        nc.instrument_manufacturer = \
                        station_dict[statname]['manufacturer'][sensorname]
        nc.netcdf_version = "4"
        nc.data_owner = ("NA")
        nc.licence = ("NA")
        nc.processing_level = "No imputation for missing or erroneous values."
        nc.static_position_station =  ("Latitude: "
                            + str(station_dict[statname]['coords']['lat'])
                            + ", Longitude: "
                            + str(station_dict[statname]['coords']['lon']))
        nc.static_position_model =  ("Latitude: "
                            + "{:.2f}".format(lats_model[0])
                            + ", Longitude: "
                            + "{:.2f}".format(lons_model[0]))
        nc.collocation_distance = "{:.2f}".format(hdist[0][0]) + " km"
        nc.static_collocation_idx =  ("idx: "
                            + str(idx[0])
                            + ", idy: "
                            + str(idy[0]))
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
                               dimensions=('time')
                               )
        ncHs_stat_1h = nc.createVariable(
                               'Hs_station_1h',
                               np.float64,
                               dimensions=('time'),
                               )
        ncHs_stat_10min = nc.createVariable(
                               'Hs_station_10min',
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
        # Hs_stat_1h
        ncHs_stat_1h.standard_name = (
                          'sea_surface_wave_significant_height_hourly '\
                        + 'from ' + statname
                                    )
        ncHs_stat_1h.long_name = ( 'Significant wave height hourly estimate'
                                + 'from spectrum (not sure) from ' 
                                + statname 
                                + '. Estimate is computed using a ' 
                                + 'convolution with window of 3 ' 
                                + 'effectively smoothing the 10min ' 
                                + 'values as retrieved with the '  
                                + 'plotd22 program')
        ncHs_stat_1h.units = 'm'
        ncHs_stat_1h.valid_range = 0., 25.
        ncHs_stat_1h[:] = Hs_stat_1h
        # Hs_stat_10min
        ncHs_stat_10min.standard_name = (
                          'sea_surface_wave_significant_height_10min '\
                        + 'from ' + statname
                                    )
        ncHs_stat_10min.long_name = ( 'Significant wave height 10 min estimate'\
                                + 'from spectrum (not sure) from ' + statname)
        ncHs_stat_10min.units = 'm'
        ncHs_stat_10min.valid_range = 0., 25.
        ncHs_stat_10min[:] = Hs_stat_10min
    nc.close()

def dumptonc_ts_station(outpath,filename,title,basetime,\
                        obs_dict,statname,sensorname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = obs_dict['time']
    Hs_stat_10min = obs_dict['Hs_stat_10min']
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
        nc.variables['Hs'][startidx:endidx] = Hs_stat_10min[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
#                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = statname
        nc.instrument_type = sensorname
        nc.instrument_specs = "?"
        nc.instrument_manufacturer = station_dict[statname]\
                                    ['manufacturer'][sensorname]
        nc.netcdf_version = "4"
        nc.data_owner = "?"
        nc.licence = "?"
        nc.processing_level = "No imputation for missing or erroneous values."
        nc.static_position_station =  ("Latitude: "
                            + str(station_dict[statname]['coords']['lat'])
                            + ", Longitude: "
                            + str(station_dict[statname]['coords']['lon']))
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
        ncHs_stat_10min = nc.createVariable(
                               'Hs',
                               np.float64,
                               dimensions=('time'),
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time'
        nctime.long_name = 'Time of measurement'
        nctime.units = 'seconds since ' + str(basetime)
        nctime.delta_t = '10 min'
        nctime[:] = time
        # Hs_stat_10min
        ncHs_stat_10min.standard_name = (
                          'sea_surface_wave_significant_height_10min'                                       )
        ncHs_stat_10min.long_name = ( 'Significant wave height retrieved '
                                    + 'at imposed 10 min interval, '
                                    + 'estimation method currently unknown '
                                    + '(spectrum or zero crossings?)')
        ncHs_stat_10min.units = 'm'
        ncHs_stat_10min.valid_range = 0., 25.
        ncHs_stat_10min[:] = Hs_stat_10min
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

    # generate time for netcdf file
    basetime=sa_obj.basetime
    nctime.units = 'seconds since 2000-01-01 00:00:00'
    nctime[:] = sa_obj.time
    ncHs.units = 'm'
    ncHs[:] = sa_obj.Hs
    ncHs.standard_name = 'sea_surface_wave_significant_height'
    ncHs.long_name = \
        'Significant wave height estimate from altimeter wave form'
    ncHs.valid_range = 0., 25.
    nclongitude.units = 'degree_east'
    nclongitude[:] = sa_obj.rloc[1]
    nclongitude.standard_name = 'longitude'
    nclongitude.valid_min = -180.
    nclongitude.valid_max = 180.
    nclatitude[:] = sa_obj.rloc[0]
    nclatitude.standard_name = 'latitude'
    nclatitude.units = 'degree_north'
    nclatitude.valid_min = -90.
    nclatitude.valid_max = 90.
    nc.close()

def dumptonc_ts_pos(outpath,filename,title,basetime,\
                    coll_dict,model,varname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = coll_dict['time']
    var_model = coll_dict[varname]
    lons_model = coll_dict['lons_model']
    lats_model = coll_dict['lats_model']
    lons_pos = coll_dict['lons_pos']
    lats_pos = coll_dict['lats_pos']
    dist = coll_dict['hdist']
    idx = coll_dict['idx']
    idy = coll_dict['idy']
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
#        nc.variables['dist'][startidx:endidx] = dist[:]
        nc.variables[varname][startidx:endidx] = var_model[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
#                        format='NETCDF4_CLASSIC'
                        )
        # global attributes
        nc.title = title
#        nc.netcdf_version = "NETCDF4_CLASSIC"
        nc.netcdf_version = "NETCDF4"
        nc.processing_level = "No post-processing performed"
        nc.static_position_station =  ("Latitude: "
                            + "{:.4f}".format(lats_pos[0])
                            + ", Longitude: "
                            + "{:.4f}".format(lons_pos[0]))
        nc.static_position_model =  ("Latitude: "
                            + "{:.4f}".format(lats_model[0])
                            + ", Longitude: "
                            + "{:.4f}".format(lons_model[0]))
        nc.static_collocation_idx =  ("idx: "
                            + str(idx[0])
                            + ", idy: "
                            + str(idy[0]))
        nc.static_collocation_distance =  ("{:.4f}".format(dist[0]) + " km")
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
#        ncdist = nc.createVariable(
#                               'dist',
#                               np.float64,
#                               dimensions=('time')
#                               )
        ncvar_model = nc.createVariable(
                               varname,
                               np.float64,
                               dimensions=('time')
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = var_dict['time']['standard_name']
        nctime.units = var_dict['time']['units'] + ' ' + str(basetime)
        nctime[:] = time
#        # dist
#        ncdist.standard_name = 'collocation_distance'
#        ncdist.units = 'km'
#        ncdist[:] = dist
        # var_model
        ncvar_model.standard_name = var_dict[varname]['standard_name']
        ncvar_model.units = var_dict[varname]['units']
        ncvar_model.valid_range = var_dict[varname]['valid_range'][0], \
                                  var_dict[varname]['valid_range'][1]
        ncvar_model[:] = var_model
    nc.close()

def dumptonc_ts_pos_wind(outpath,filename,title,basetime,\
                    obs_dict,model,statname,sensorname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = obs_dict['time']
    u10_model = obs_dict['u10_model']
    v10_model = obs_dict['v10_model']
    lons_model = obs_dict['lons_model']
    lats_model = obs_dict['lats_model']
    lons_stat = obs_dict['lons_stat']
    lats_stat = obs_dict['lats_stat']
    idx = obs_dict['idx']
    idy = obs_dict['idy']
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
        nc.variables['u10_model'][startidx:endidx] = u10_model[:]
        nc.variables['v10_model'][startidx:endidx] = v10_model[:]
        nc.variables['longitude_model'][startidx:endidx] = lons_model[:]
        nc.variables['latitude_model'][startidx:endidx] = lats_model[:]
    else:
        os.system('mkdir -p ' + outpath)
        nc = netCDF4.Dataset(
                        fullpath,mode='w',
#                        format='NETCDF4'
                        )
        # global attributes
        nc.title = title
        nc.station_name = statname
        nc.instrument_type = "? " + sensorname + " ?"
        nc.instrument_specs = "?"
        nc.instrument_manufacturer = "?"
        nc.netcdf_version = "4"
        nc.data_owner = ("?")
        nc.licence = ("?")
        nc.processing_level = "No imputation for missing or erroneous values."
        nc.static_position_station =  ("Latitude: "
                            + "{:.4f}".format(lats_stat[0])
                            + ", Longitude: "
                            + "{:.4f}".format(lons_stat[0]))
        nc.static_position_model =  ("Latitude: "
                            + "{:.4f}".format(lats_model[0])
                            + ", Longitude: "
                            + "{:.4f}".format(lons_model[0]))
        nc.static_collocation_idx =  ("idx: "
                            + str(idx[0])
                            + ", idy: "
                            + str(idy[0]))
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
        ncu10_model = nc.createVariable(
                               'u10_model',
                               np.float64,
                               dimensions=('time')
                               )
        ncv10_model = nc.createVariable(
                               'v10_model',
                               np.float64,
                               dimensions=('time')
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time'
        nctime.long_name = 'Time of measurement'
        nctime.units = 'seconds since ' + str(basetime)
        nctime[:] = time
        # u10_model
        ncu10_model.standard_name = 'u10'
        ncu10_model.long_name = (
                            '10m wind speed in x-direction from ' 
                            + model + ' forcing file' )
        ncu10_model.units = 'm/s'
        ncu10_model[:] = u10_model
        # v10_model
        ncv10_model.standard_name = 'v10'
        ncv10_model.long_name = (
                            '10m wind speed in y-direction from '
                            + model + ' forcing file' )
        ncv10_model.units = 'm/s'
        ncv10_model[:] = v10_model
    nc.close()

def dumptonc_LCWVF(outpath,filename,title,basetime,coll_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = coll_dict['time']
    lt = coll_dict['lt']
    stationid = coll_dict['stationid']
    Hs = coll_dict['Hs']
    Tp = coll_dict['Tp']
    Tz = coll_dict['Tz']
    thq = coll_dict['thq']
    u10 = coll_dict['u10']
    v10 = coll_dict['v10']
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    os.system('mkdir -p ' + outpath)
    nc = netCDF4.Dataset(
                    fullpath,mode='w',
                    )
    # global attributes
    nc.title = title
    # dimensions
    dimtime = nc.createDimension(
                            'time',
                            size=len(time)
                            )
    dimlt = nc.createDimension(
                            'leadtime',
                            size=len(lt)
                            ) 
    dimstationid = nc.createDimension(
                            'stationid',
                            size=len(stationid)
                            )
    dimstr = nc.createDimension(
                            'strdim',
                            size=1
                            )
    # base variables
    nctime = nc.createVariable(
                           'time',
                           np.float64,
                           dimensions=('time')
                           )
    nclt = nc.createVariable(
                           'leadtime',
                           np.int,
                           dimensions=('leadtime')
                           )
    ncstationid = nc.createVariable(
                           'stationid',
                           np.str,
                           dimensions=('stationid','strdim')
                           )
    # other variables
    ncHs = nc.createVariable(
                           'Hs',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncTp = nc.createVariable(
                           'Tp',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncTz = nc.createVariable(
                           'Tz',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncthq = nc.createVariable(
                           'thq',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncu10 = nc.createVariable(
                           'u10',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncv10 = nc.createVariable(
                           'v10',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    # time
    nctime.standard_name = 'time'
    nctime.units = 'seconds since ' + str(basetime)
    nctime[:] = time
    # lt
    nclt.standard_name = 'leadtime'
    nclt.units = 'hours'
    nclt[:] = lt
    # stationid
    ncstationid.standard_name = 'stationid'
    ncstationid.units = 'None'
    ncstationid[:] = stationid
    # Hs
    ncHs.standard_name = 'sea_surface_wave_significant_height'
    ncHs.units = 'm'
    ncHs[:] = Hs
    # Tp
    ncTp.standard_name = 'sea_surface_wave_peak_period_from_variance_spectral_density'
    ncTp.units = 's'
    ncTp[:] = Tp
    # Tz
    ncTz.standard_name = 'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
    ncTz.units = 's'
    ncTz[:] = Tz
    # thq
    ncthq.standard_name = 'sea_surface_wave_to_direction'
    ncthq.units = 'degree'
    ncthq[:] = thq
    # u10
    ncu10.standard_name = 'x_wind'
    ncu10.units = 'm/s'
    ncu10[:] = u10
    # v10
    ncv10.standard_name = 'y_wind'
    ncv10.units = 'm/s'
    ncv10[:] = v10
    nc.close()

def dumptonc_LCWVF(outpath,filename,title,basetime,coll_dict,varname):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    time = coll_dict['time']
    lt = coll_dict['lt']
    stationid = coll_dict['stationid']
    var = coll_dict[varname]
    fullpath = outpath + filename
    print ('Dump data to file: ' + fullpath)
    os.system('mkdir -p ' + outpath)
    nc = netCDF4.Dataset(
                    fullpath,mode='w',
                    )
    # global attributes
    nc.title = title
    # dimensions
    dimtime = nc.createDimension(
                            'time',
                            size=len(time)
                            )
    dimlt = nc.createDimension(
                            'leadtime',
                            size=len(lt)
                            )
    dimstationid = nc.createDimension(
                            'stationid',
                            size=len(stationid)
                            )
    dimstr = nc.createDimension(
                            'strdim',
                            size=1
                            )
    # base variables
    nctime = nc.createVariable(
                           'time',
                           np.float64,
                           dimensions=('time')
                           )
    nclt = nc.createVariable(
                           'leadtime',
                           np.int,
                           dimensions=('leadtime')
                           )
    ncstationid = nc.createVariable(
                           'stationid',
                           np.str,
                           dimensions=('stationid','strdim')
                           )
    # other variables
    ncHs = nc.createVariable(
                           'Hs',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncTp = nc.createVariable(
                           'Tp',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncTz = nc.createVariable(
                           'Tz',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncthq = nc.createVariable(
                           'thq',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncu10 = nc.createVariable(
                           'u10',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    ncv10 = nc.createVariable(
                           'v10',
                           np.float64,
                           dimensions=('time','leadtime','stationid')
                           )
    # time
    nctime.standard_name = 'time'
    nctime.units = 'seconds since ' + str(basetime)
    nctime[:] = time
    # lt
    nclt.standard_name = 'leadtime'
    nclt.units = 'hours'
    nclt[:] = lt
    # stationid
    ncstationid.standard_name = 'stationid'
    ncstationid.units = 'None'
    ncstationid[:] = stationid
    # Hs
    ncHs.standard_name = 'sea_surface_wave_significant_height'
    ncHs.units = 'm'
    ncHs[:] = Hs
    # Tp
    ncTp.standard_name = 'sea_surface_wave_peak_period_from_variance_spectral_density'
    ncTp.units = 's'
    ncTp[:] = Tp
    # Tz
    ncTz.standard_name = 'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
    ncTz.units = 's'
    ncTz[:] = Tz
    # thq
    ncthq.standard_name = 'sea_surface_wave_to_direction'
    ncthq.units = 'degree'
    ncthq[:] = thq
    # u10
    ncu10.standard_name = 'x_wind'
    ncu10.units = 'm/s'
    ncu10[:] = u10
    # v10
    ncv10.standard_name = 'y_wind'
    ncv10.units = 'm/s'
    ncv10[:] = v10
    nc.close()
