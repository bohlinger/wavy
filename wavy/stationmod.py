#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from stations. I try to mostly follow the PEP 
convention for python code style. Constructive comments on style and 
effecient programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# read_altim
import os
import sys
import netCDF4

# ignore irrelevant warnings from matplotlib for stdout
#import warnings
#warnings.filterwarnings("ignore")

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import yaml
import os

# get_altim
import urllib
import gzip
import ftplib
from ftplib import FTP

# create_file
import calendar

import sys

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# get d22 files
import pylab as pl
from datetime import datetime
import scipy as sp

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/buoy_specs.yaml'))
with open(moddir,'r') as stream:
    buoy_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/stationlist.yaml'))
with open(moddir,'r') as stream:
    locations=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/station_specs.yaml'))
with open(moddir,'r') as stream:
    station_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/pathfinder.yaml'))
with open(moddir,'r') as stream:
    pathfinder=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                        '..', 'config/d22_var_dicts.yaml'))
with open(moddir, 'r') as stream:
    d22_var_dicts=yaml.safe_load(stream)

station_d22_starc = pathfinder['station_d22_starc']
station_d22_opdate = pathfinder['station_d22_opdate']
stationpath_lustre_om = pathfinder['stationpath_lustre_om']
stationpath_lustre_hi = pathfinder['stationpath_lustre_hi']

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
    Class to handle platform based time series.
    '''
    basedate = datetime(1970,1,1)
    time_unit = 'seconds since 1970-01-01 00:00:00.0'
    def __init__(self,statname,sensorname,sdate,edate,
                mode='d22',deltat=10,varalias='Hs_10min'):
        print ('# ----- ')
        print (" ### Initializing station_class object ###")
        print ('Chosen period: ' + str(sdate) + ' - ' + str(edate))
        print (" Please wait ...")
        stdvarname = d22_var_dicts['standard_name'][varalias]
        var, time, timedt = self.get_station(
                                    statname, # change to stdvarname in future
                                    sdate,edate,
                                    mode,deltat,
                                    sensorname,
                                    varalias)
        vardict = {
                    stdvarname:var,
                    'time':time,
                    'datetime':timedt,
                    'time_unit':self.time_unit,
                    'longitude':[locations[statname][1]]*len(var),
                    'latitude':[locations[statname][0]]*len(var)
                    }
        # in future coordinates need to be properly 
        # defined with names longitude and latitude 
        # in yaml file not as it is now in stationlist.yaml
        self.vars = vardict
        #self.varname = varname # if d22: Hs_10m; if nc: variable name
        self.stdvarname = d22_var_dicts['standard_name'][varalias]
        self.varname = varalias
        self.varalias = varalias
        self.sdate = sdate
        self.edate = edate
        self.lat = locations[statname][0]
        self.lon = locations[statname][1]
        self.sensorname = sensorname
        self.statname = statname
        print (" ### station_class object initialized ###")
        print ('# ----- ')

    def get_station(self,statname,sdate,edate,mode,deltat,sensorname,varalias):
        if mode == 'nc':
            basedate = self.basedate
            tmpdate = sdate
            var = []
            var_obs = []
            time = []
            timedt = []
            while (tmpdate <= edate):
                filepath = stationpath_lustre_hi + statname + \
                    tmpdate.strftime('_%Y%m') + '.nc'
                nc = netCDF4.Dataset(filepath,'r')
                Hs_OBS = nc.variables[varname][:]
                tmp = np.nanmean(Hs_OBS,axis=0)
                var.append(tmp)
                var_obs.append(Hs_OBS)
                tmp = nc.variables['time'][:]
                time.append(tmp)
                for s in tmp:
                    timedt.append(basedate + timedelta(seconds=s))
                nc.close()
                tmpdate = tmpdate + relativedelta(months = +1)
            var = flatten(var)
            time = flatten(time)
        elif mode == 'd22':
            sdatetmp = sdate - timedelta(days=1)
            edatetmp = edate + timedelta(days=1)
            sl = parse_d22(statname,sdatetmp,edatetmp)
            sensor_lst, dates = extract_d22(sl)
            tmpdate = sdate
            var = []
            var_obs = []
            time = []
            timedt = []
            while (tmpdate <= edate):
                ctime, idxtmp = matchtime(tmpdate,tmpdate,dates['10min'],
                                          timewin=2)
                try:
                    tmp = sensor_lst[station_dict['platform'][statname]\
                                        ['sensor'][sensorname]]\
                                        [varalias][idxtmp][0]
                    time.append((tmpdate-self.basedate).total_seconds())
                    var.append(np.real(tmp))
                    timedt.append(tmpdate)
                except Exception as e:
                    print(e)
                tmpdate = tmpdate + timedelta(minutes=deltat)
        return var, time, timedt

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
    with open("../config/d22_var_dicts.yaml", 'r') as stream:
        d22_var_dicts=yaml.safe_load(stream)
    dat=d22_var_dicts['dat']
    WM1=d22_var_dicts['WM1']
    WM2=d22_var_dicts['WM2']
    WM3=d22_var_dicts['WM3']
    WL1=d22_var_dicts['WL1']
    WL2=d22_var_dicts['WL2']
    WL3=d22_var_dicts['WL3']
    WIA=d22_var_dicts['WIA']
    WIB=d22_var_dicts['WIB']
    WIC=d22_var_dicts['WIC']
    tseries=[]
    for i, line in enumerate(searchlines):
        if "!!!!" in line: 
            tseriesl = [l.strip() for l in searchlines[i+3:i+5]]
            tseries.append(' '.join(tseriesl))
            date_object = datetime.strptime(' '.join(tseriesl),'%d-%m-%Y %H:%M')
            dat['10min'].append(date_object)
            # nan for all variables
            for WM in [WM1,WM2,WM3,WIA,WIB,WIC,WL1,WL2,WL3]:
                for var in WM.keys():
                    if var != 'name':
                        WM[var].append(sp.nan)
        for WM in [WM1,WM2,WM3]:  
            if str(WM['name']) in line:
                for l in searchlines[i+2:i+3]:
                    WM['Hs_10min'][-1]=floater(l.strip())
                for l in searchlines[i+5:i+6]:
                    WM['Tp_10min'][-1]=floater(l.strip())  
                for l in searchlines[i+11:i+12]:
                    WM['Tm02_10min'][-1]=floater(l.strip())
                for l in searchlines[i+12:i+13]:
                    WM['Tm01_10min'][-1]=floater(l.strip())
                for l in searchlines[i+9:i+10]:
                    WM['Tm10_10min'][-1]=floater(l.strip())
                for l in searchlines[i+17:i+18]:
                    WM['Mdir_10min'][-1]=floater(l.strip())
                for l in searchlines[i+16:i+17]:
                    WM['Pdir_10min'][-1]=floater(l.strip())
        for WI in [WIA,WIB,WIC]:
            if str(WI['name']) in line:
                for l in searchlines[i+10:i+11]:
                    WI['FF_10min_10m'][-1]=floater(l.strip())
                for l in searchlines[i+16:i+17]:
                    WI['FF_10min_sensor'][-1]=floater(l.strip())
                for l in searchlines[i+13:i+14]:
                    WI['DD_10min_sensor'][-1]=floater(l.strip()) 
        for WL in [WL1,WL2,WL3]:
            if str(WL['name']) in line:
                for l in searchlines[i+2:i+3]:
                    WL['Hlat'][-1]=floater(l.strip())
    N = len(dat['10min'])
    #Convert data to arrays
    dat['10min']=sp.array(dat['10min'])
    for WM in [WM1,WM2,WM3,WIA,WIB,WIC,WL1,WL2,WL3]:
        for var in WM.keys():
            if var != 'name':
                WM[var] = sp.array(WM[var])
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
                if var == 'Tm02_10min':
                    WINDOW = 3
                    wmt = sp.power(WM[var],2)
                    weights = sp.ones((WINDOW))/WINDOW
                    wmt[2:-2:2] = np.convolve(wmt[0::2], weights,mode='valid')
                    wmt[3:-2:2] = np.convolve(wmt[1::2], weights,mode='valid')
                    wmt[[0,1,-2,-1]]= sp.nan
                    WM['Tm02_1hr'] = sp.sqrt(wmt)
                if var == 'Tm01_10min':
                    WINDOW = 3
                    wmt = sp.power(WM[var],2)
                    weights = sp.ones((WINDOW))/WINDOW
                    wmt[2:-2:2] = np.convolve(wmt[0::2], weights,mode='valid')
                    wmt[3:-2:2] = np.convolve(wmt[1::2], weights,mode='valid')
                    wmt[[0,1,-2,-1]]= sp.nan
                    WM['Tm01_1hr'] = sp.sqrt(wmt)
                if var == 'Tm10_10min':
                    WINDOW = 3
                    wmt = sp.power(WM[var],2)
                    weights = sp.ones((WINDOW))/WINDOW
                    wmt[2:-2:2] = np.convolve(wmt[0::2], weights,mode='valid')
                    wmt[3:-2:2] = np.convolve(wmt[1::2], weights,mode='valid')
                    wmt[[0,1,-2,-1]]= sp.nan
                    WM['Tm10_1hr'] = sp.sqrt(wmt)
    # CAUTION: 10min data is extracted for entire days only 00:00h - 23:50h
    sensor_lst = [WM1,WM2,WM3]
    dates = deepcopy(dat)
    del dat
    return sensor_lst, dates

def matchtime(sdate,edate,time,time_unit=None,timewin=None):
    '''
    fct to obtain the index of the time step closest to the 
    requested time including the respective time stamp(s). 
    Similarily, indices are chosen for the time and defined region.
    time_win is in [minutes]
    '''
    if timewin is None:
        timewin = 0
    # create list of datetime instances
    timelst=[]
    ctime=[]
    cidx=[]
    idx=0
    if (edate is None or sdate==edate):
        for element in time:
            if time_unit is None:
                tmp = element
            else:
                tmp = netCDF4.num2date(element,time_unit)
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
            if basetime is None:
                tmp = element
            else:
                tmp = netCDF4.num2date(element,time_unit)
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
