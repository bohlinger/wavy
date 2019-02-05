#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from stations. I try to mostly follow the PEP 
convention for python code style. Constructive comments on style and 
effecient programming are most welcome!
'''
__version__ = "0.5.0"
__author__="Patrik Bohlinger, Norwegian Meteorological Institute"
__maintainer__ = "Patrik Bohlinger"
__email__ = "patrikb@met.no"
__status__ = "operational ARCMFC branch"

# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# ignore irrelevant warnings from matplotlib for stdout
import warnings
warnings.filterwarnings("ignore")

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os

# get_altim
import urllib
import gzip
import ftplib
from ftplib import FTP

# read_altim
import netCDF4 as netCDF4

# libraries for quip and quim imported in fct for speedup

# create_file
import calendar

# libraries for parallel computing
from joblib import Parallel, delayed
import multiprocessing as mp

# bintime
import math

# progress bar
import sys

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# get d22 files
from pathfinder import station_d22_starc, station_d22_opdate
from pathfinder import stationpath_lustre_om, stationpath_lustre_hi
import pylab as pl
from d22_var_dicts import dat,  WM1, WM2, WM3, \
                                WL1, WL2, WL3, \
                                WIA, WIB, WIC
from datetime import datetime
import scipy as sp

# --- global functions ------------------------------------------------#

# define flatten function for lists
''' fct does the following:
flat_list = [item for sublist in TIME for item in sublist]
or:
for sublist in TIME:
for item in sublist:
flat_list.append(item)
'''
flatten = lambda l: [item for sublist in l for item in sublist]

# ---------------------------------------------------------------------#


class station_class():
    '''
    class to handle station files on  Hs[time], lat[time], lon[time] 
    This class offers the following added functionality:
     - get the closest time stamp(s)
     - get Hs value for this time
    '''
    #stationpath_lustre_om = '/lustre/storeA/project/fou/om/' + \
    #                        'waveverification/data/'
    #stationpath_lustre_hi = '/lustre/storeA/project/fou/hi/' + \
    #                        'waveverification/data/'
    basedate = datetime(1970,1,1)
    from stationlist import locations

    def __init__(self,statname,sdate,edate,mode=None):
        print ('# ----- ')
        print (" ### Initializing station_class instance ###")
        print ('# ----- ')
        if mode is None:
            mode = 'nc' # or if 'd22' in future
        hs, hs_obs, time, timedt = self.get_station(statname,sdate,edate,mode)
        self.hs = hs
        self.hs_obs = hs_obs
        self.time = time
        self.timedt = timedt
        self.basedate = self.basedate
        self.lat = self.locations[statname][0]
        self.lon = self.locations[statname][1]

    def get_station(self,statname,sdate,edate,mode):
        if mode == 'nc':
            basedate = self.basedate
            tmpdate = sdate
            hs = []
            hs_obs = []
            time = []
            timedt = []
            while (tmpdate <= edate):
                filepath = stationpath_lustre_hi + statname + \
                    tmpdate.strftime('_%Y%m') + '.nc'
                nc = netCDF4.Dataset(filepath,'r')
                #ncg = nc.groups['OBS_d22']
                #var = np.mean(ncg.variables['Hs'][:],axis=0)
                Hs_OBS = nc.variables['Hs_OBS'][:]
                var = np.nanmean(Hs_OBS,axis=0)
                hs.append(var)
                hs_obs.append(Hs_OBS)
                var = nc.variables['time'][:]
                time.append(var)
                for s in var:
                    timedt.append(basedate + timedelta(seconds=s))
                nc.close()
                tmpdate = tmpdate + relativedelta(months = +1)
            hs = flatten(hs)
            time = flatten(time)
        elif mode == 'd22':
            print ("no d22 fileformat included yet")
        return hs, hs_obs, time, timedt

def read_Tennholmen(date):
    from buoy_specs import buoy_dict
    url=date.strftime(buoy_dict['Tennholmen']['url_template'])
    filename=date.strftime(buoy_dict['Tennholmen']['file_template'])
    filename_heave_spec=date.strftime(buoy_dict['Tennholmen']
                                    ['file_template_heave_spec'])
    filename_prim_dir_spec=date.strftime(buoy_dict['Tennholmen']
                                    ['file_template_prim_dir_spec'])
    filename_sec_dir_spec=date.strftime(buoy_dict['Tennholmen']
                                    ['file_template_sec_dir_spec'])
    filename_pos=date.strftime(buoy_dict['Tennholmen']
                                    ['file_template_pos'])
    basetime=buoy_dict['Tennholmen']['basetime']
    tmpdir = 'tmp_Tennholmen/'
    t=os.system('mkdir -p ' + tmpdir)
    t=os.system('wget ' + url + filename + ' -P ' + tmpdir)
    time_s, Hm0, Tm02 = np.loadtxt(tmpdir + filename, skiprows=1, \
                        usecols=(buoy_dict['Tennholmen']['time'],
                                    buoy_dict['Tennholmen']['Hm0'],
                                    buoy_dict['Tennholmen']['Tm02']), 
                        unpack=True)
    t=os.system('wget ' + url + filename_pos + ' -P ' + tmpdir)
    lons, lats        = np.loadtxt(tmpdir + filename_pos, skiprows=1, \
                        usecols=(buoy_dict['Tennholmen']['lons'],
                                buoy_dict['Tennholmen']['lats']),
                        unpack=True)
    time_dt = [basetime + timedelta(seconds=time_s[i]) \
                for i in range(len(time_s))]
    print('cleaning up ...')
    t=os.system('rm -r tmp_Tennholmen')
    # convert lons, lats from radians to degree
    lons = lons * 180. / np.pi
    lats = lats * 180. / np.pi
    return time_s, time_dt, Hm0, Tm02, lons, lats
    
def get_buoy(sdate,edate,buoyname=None,mode=None):
    if (sdate.month == edate.month and sdate.year == edate.year):
        time_s, time_dt, Hm0, Tm02, lons, lats \
                             = read_Tennholmen(sdate)
        sidx=time_dt.index(sdate)
        eidx=time_dt.index(edate) + 1
        time_s_lst = time_s[sidx:eidx]
        time_dt_lst = time_dt[sidx:eidx]
        Hm0_lst = Hm0[sidx:eidx]
        Tm02_lst = Tm02[sidx:eidx]
        lons_lst = lons[sidx:eidx]
        lats_lst = lats[sidx:eidx]
    else:
        tmpdate = deepcopy(sdate)
        time_s_lst = []
        time_dt_lst = []
        Hm0_lst = []
        Tm02_lst = []
        lons_lst = []
        lats_lst = []
        while tmpdate <= edate:
            time_s, time_dt, Hm0, Tm02, lons, lats \
                            = read_Tennholmen(tmpdate)
            time_s_lst.append(time_s)
            time_dt_lst.append(time_dt)
            Hm0_lst.append(Hm0)
            Tm02_lst.append(Tm02)
            lons_lst.append(lons)
            lats_lst.append(lats)
            tmpdate = tmpdate + relativedelta(months = +1)
        del tmpdate
        time_s_lst = flatten(time_s_lst)
        time_dt_lst = flatten(time_dt_lst)
        Hm0_lst = flatten(Hm0_lst)
        Tm02_lst = flatten(Tm02_lst)
        lons_lst = flatten(lons_lst)
        lats_lst = flatten(lats_lst)
        sidx=time_dt_lst.index(sdate)
        eidx=time_dt_lst.index(edate) + 1
        time_s_lst = time_s_lst[sidx:eidx]
        time_dt_lst = time_dt_lst[sidx:eidx]
        Hm0_lst = Hm0_lst[sidx:eidx]
        Tm02_lst = Tm02_lst[sidx:eidx]
        lons_lst = lons_lst[sidx:eidx]
        lats_lst = lats_lst[sidx:eidx]
    return time_s_lst, time_dt_lst, Hm0_lst, Tm02_lst, lons_lst, lats_lst

def parse_d22(statname,sdate,edate):
    # Read all lines in file and append to searchlines
    searchlines=[]
    for d in range(int(pl.date2num(sdate)),int(pl.date2num(edate))+1): 
        dy = pl.num2date(d).strftime("%Y%m%d")
        YY = pl.num2date(d).strftime("%Y")
        ifile_starc = (station_d22_starc + statname 
                        + "/d22/" + YY + "/" + dy + ".d22")
        ifile_opdata = (station_d22_opdate + statname 
                        + "/d22/" + dy + ".d22")
        if os.path.isfile(ifile_opdata):
            f = open(ifile_opdata, "r")
            searchlines = searchlines + f.readlines()
            f.close()
        elif os.path.isfile(ifile_starc):
            f = open(ifile_starc, "r")
            searchlines = searchlines + f.readlines()
            f.close()
    return searchlines

# Function that converts 's' to float32 or nan if floater throws exception 
def floater(s):
    try:
        x = np.float32(s)
    except:
        x = np.nan
    return x

def extract_d22(searchlines):
    #Extract data of choice - reading searchlines
    tseries=[]
    for i, line in enumerate(searchlines):
        if "!!!!" in line: 
            tseriesl = []
            for l in searchlines[i+3:i+5]:
                tseriesl.append(l.strip())
            tseries.append(' '.join(tseriesl))
            date_object = datetime.strptime(' '.join(tseriesl),'%d-%m-%Y %H:%M')
            dat['10min'].append(date_object)
            #legger inn nan for alle variable
            for WM in [WM1,WM2,WM3,WIA,WIB,WIC,WL1,WL2,WL3]:
                for var in WM.keys():
                    if var != 'name':
                        WM[var].append(sp.nan)
        for WM in [WM1,WM2,WM3]:  
            if str(WM['name']) in line:
                for l in searchlines[i+2:i+3]:
                    WM['Hs_10min'][-1]=floater(l.strip())
                for l in searchlines[i+4:i+5]:
                    WM['Hmax_10min'][-1]=floater(l.strip())
                for l in searchlines[i+5:i+6]:
                    WM['Tp_10min'][-1]=floater(l.strip())  
                for l in searchlines[i+11:i+12]:
                    WM['Tm_10min'][-1]=floater(l.strip())
        for WI in [WIA,WIB,WIC]:
            if str(WI['name']) in line:
                for l in searchlines[i+10:i+11]:
                    WI['FF_10min'][-1]=floater(l.strip())
                for l in searchlines[i+13:i+14]:
                    WI['DD_10min'][-1]=floater(l.strip()) 
        for WL in [WL1,WL2,WL3]:
            if str(WL['name']) in line:
                for l in searchlines[i+2:i+3]:
                    WL['Hlat'][-1]=floater(l.strip())
    # From here ....
    N = len(dat['10min'])
    #Convert data to arrays
    dat['10min']=sp.array(dat['10min'])
    for WM in [WM1,WM2,WM3,WIA,WIB,WIC,WL1,WL2,WL3]:
        for var in WM.keys():
            if var != 'name':
                WM[var] = sp.array(WM[var])
                #Set all negative values to nan
                #WM[var][WM[var]<0] = sp.nan
                #If variable does not exist
                #if len(sp.array(WM[var]))==0:
                #print len(dat['10min'])
                if len(sp.array(WM[var]))!=N:
                    WM[var] = np.empty(len(dat['10min']))
                    WM[var][:] = sp.nan
                    os.system("pause")  
                #Create 1 hour averages
                if var == 'Hs_10min':
                    WINDOW = 3
                    wmt = sp.power(WM[var],2)
                    weights = sp.ones((WINDOW))/WINDOW
                    wmt[2:-2:2] = np.convolve(wmt[0::2], weights,mode='valid')
                    wmt[3:-2:2] = np.convolve(wmt[1::2], weights,mode='valid')
                    wmt[[0,1,-2,-1]]= sp.nan
                    WM['Hs_1hr'] = sp.sqrt(wmt)
                if var == 'Tm_10min':
                    WINDOW = 3
                    wmt = sp.power(WM[var],2)
                    weights = sp.ones((WINDOW))/WINDOW
                    wmt[2:-2:2] = np.convolve(wmt[0::2], weights,mode='valid')
                    wmt[3:-2:2] = np.convolve(wmt[1::2], weights,mode='valid')
                    wmt[[0,1,-2,-1]]= sp.nan
                    WM['Tm_1hr'] = sp.sqrt(wmt)
    # CAUTION: 10min data is extracted for entire days only 00:00h - 23:50h
    return WM1, WM2, WM3, dat

def matchtime(sdate,edate,time,basetime,timewin=None):
    '''
    fct to obtain the index of the time step closest to the 
    requested time including the respective time stamp(s). 
    Similarily, indices are chosen for the time and defined region.
    '''
    if timewin is None:
        timewin = 0
    # create list of datetime instances
    timelst=[]
    ctime=[]
    cidx=[]
    idx=0
    print ('Time window is: ', timewin)
    if (edate is None or sdate==edate):
        for element in time:
            tmp=basetime + timedelta(seconds=element)
            timelst.append(tmp)
            # choose closest match within window of win[minutes]
            if (tmp >= sdate-timedelta(minutes=timewin)
            and tmp <= sdate+timedelta(minutes=timewin)):
                ctime.append(tmp)
                cidx.append(idx)
            del tmp
            idx=idx+1
    if (edate is not None and edate!=sdate):
        for element in time:
            tmp=basetime + timedelta(seconds=element)
            timelst.append(tmp)
            if (tmp >= sdate-timedelta(minutes=timewin)
            and tmp < edate+timedelta(minutes=timewin)):
                ctime.append(tmp)
                cidx.append(idx)
            del tmp
            idx=idx+1
    return ctime, cidx

def get_loc_idx(init_lats,init_lons,target_lat,target_lon,mask=None):
    from utils import haversine
    distM = np.zeros(init_lats.shape)*np.nan
    for i in range(init_lats.shape[0]):
        for j in range(init_lats.shape[1]):
            if mask is None:
                distM[i,j] = haversine(init_lons[i,j],init_lats[i,j],
                                    target_lon,target_lat)
            else:
                if isinstance(mask[i,j],(np.float32)):
                    distM[i,j] = haversine(init_lons[i,j],init_lats[i,j],
                                    target_lon,target_lat)
    idx,idy = np.where(distM==np.nanmin(distM))
    return idx, idy, distM, init_lats[idx,idy], init_lons[idx,idy]

def get_model(model,fc_date,init_date=None,leadtime=None):
    """ 
    fct to get model data
    if model ARCMFC you need fc_date, init_date
    if model mwam4 you need fc_date, leadtime
    """
    from utils import haversine
    from model_specs import model_dict
    print ("Get model data according to selected date ....")
    if init_date is None:
        print ("leadtime:",leadtime,"h")
    else:
        print ("init_date:",init_date)
    print ("fc_date:",fc_date)
    if model == 'ARCMFC':
        filestr = (model_dict[model]['path']
          + fc_date.strftime('%Y%m%d')
          + init_date.strftime(model_dict[model]['file_template']))
    elif (model == 'mwam4' or model=='mwam8'):
        if fc_date == init_date:
            filestr = (init_date.strftime(model_dict[model]['path_template'])
                + init_date.strftime(model_dict[model]['file_template']))
        else:
            if leadtime%6!=0:
                print ("leadtime needs to be multiple of 6h")
                print ("exit loop ...")
                #sys.exit()
            else:
                tmpdate = fc_date - timedelta(hours=leadtime)
                filedate = tmpdate
                filestr = (filedate.strftime(model_dict[model]['path_template'])
                    + filedate.strftime(model_dict[model]['file_template']))
            del tmpdate
    print (filestr)
    f = netCDF4.Dataset(filestr,'r')
    model_lons = f.variables[model_dict[model]['lons']][:]
    model_lats = f.variables[model_dict[model]['lats']][:]
    model_time = f.variables[model_dict[model]['time']][:]
    # Hs [time,lat,lon]
    model_Hs = f.variables[model_dict[model]['Hs']][:].squeeze()
    f.close()
    model_basetime = model_dict[model]['basetime']
    model_time_dt=[]
    for element in model_time:
        model_time_dt.append(model_basetime
                    + timedelta(seconds=element))
    model_time_dt_valid = [model_time_dt[model_time_dt.index(fc_date)]]
    model_hs_valid = model_Hs[model_time_dt.index(fc_date),:,:]
    return model_time_dt, model_hs_valid, model_lons, model_lats

def dumptonc(time,model,obs,outpath,filename):
    """
    create a simple netcdf
    """
    # 1. check if nc file already exists
    fullpath = outpath + filename
    os.system('mkdir -p ' + outpath)
    print ('Dump data to file: ' + fullpath)
    nc = netCDF4.Dataset(
                    fullpath,mode='w',
                    format='NETCDF4'
                    )
    nc.title = 'quick data dump'
    rtimerange=len(time)
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
    ncmodelHs = nc.createVariable(
                           'modelHs',
                           np.float64,
                           dimensions=('time')
                           )
    ncobsHs = nc.createVariable(
                           'obsHs',
                           np.float64,
                           dimensions=('time')
                           )
    # generate time for netcdf file
    basetime=datetime(1970,1,1)
    nctime.units = 'seconds since 1970-01-01 00:00:00'
    nctime[:] = time
    ncmodelHs.units = 'm'
    ncmodelHs[:] = model
    ncobsHs.units = 'm'
    ncobsHs[:] = obs
    nc.close()

# --- help ------------------------------------------------------------#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
This module encompasses classes and methods to read and process wave
field related data from stations.\n
Usage:
from stationmod import station_class as sc
from datetime import datetime, timedelta
sc_obj = sc('ekofiskL',sdate,edate)
        """,
        formatter_class = RawTextHelpFormatter
        )
    args = parser.parse_args()
