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
import yaml

# get_altim
import urllib
import gzip
import ftplib
from ftplib import FTP

# read_altim
import netCDF4 as netCDF4

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
import pylab as pl
from datetime import datetime
import scipy as sp

# read yaml config files:
with open("/home/patrikb/wavy/wavy/buoy_specs.yaml", 'r') as stream:
    buoy_dict=yaml.safe_load(stream)
with open("/home/patrikb/wavy/wavy/pathfinder.yaml", 'r') as stream:
    pathfinder=yaml.safe_load(stream)
with open("/home/patrikb/wavy/wavy/stationlist.yaml", 'r') as stream:
    locations=yaml.safe_load(stream)

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
    class to handle station files on  Hs[time], lat[time], lon[time] 
    This class offers the following added functionality:
     - get the closest time stamp(s)
     - get Hs value for this time
    '''
    basedate = datetime(1970,1,1)
    def __init__(self,statname,sdate,edate,
                mode=None,deltat=None,sensorname=None,varname=None):
        print ('# ----- ')
        print (" ### Initializing station_class object ###")
        print ('Chosen period: ' + str(sdate) + ' - ' + str(edate))
        print (" Please wait ...")
        print ('# ----- ')
        if mode is None:
            mode = 'nc' # mode: 'nc', 'd22'
        if varname is None:
            varname = 'Hs_10min'
        hs, time, timedt = self.get_station(
                                    statname,
                                    sdate,edate,
                                    mode,deltat,
                                    sensorname,
                                    varname
                                )
        self.hs = hs
        self.time = time
        self.timedt = timedt
        self.basedate = self.basedate
        self.lat = locations[statname][0]
        self.lon = locations[statname][1]
        self.sensorname = sensorname
        self.statname = statname
        print (" ### station_class object initialized ###")

    def get_station(self,statname,sdate,edate,mode,deltat,sensorname,varname):
        if mode == 'nc':
            if varname is None:
                varname = 'Hs_OBS'
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
#            from station_specs import station_dict
            with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
                station_dict=yaml.safe_load(stream)
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
                    tmp = sensor_lst[station_dict[statname]['sensor']
                                            [sensorname]][varname][idxtmp][0]
                    time.append((tmpdate-self.basedate).total_seconds())
                    var.append(np.real(tmp))
                    timedt.append(tmpdate)
                except:
#                    print('no entry --> pass')
                    pass
                tmpdate = tmpdate + timedelta(minutes=deltat)
        return var, time, timedt

def read_Tennholmen(date):
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

def read_Tennholmen_ext(date):
    # extended version of read_Tennholmen
    url=date.strftime(buoy_dict['Tennholmen']['url_template'])
    filename=date.strftime(buoy_dict['Tennholmen']['file_template'])
    filename_spec_ext=date.strftime(buoy_dict['Tennholmen']
                                    ['file_template_spec_ext'])
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
    t=os.system('wget ' + url + filename_spec_ext + ' -P ' + tmpdir)
    Hs,TI,TE,T1,TZ,T3,Tc,Tdw,Tp,Qp = \
                        np.loadtxt(tmpdir + filename_spec_ext, skiprows=1, \
                        usecols=(buoy_dict['Tennholmen']['Hs'],
                                buoy_dict['Tennholmen']['TI'],
                                buoy_dict['Tennholmen']['TE'],
                                buoy_dict['Tennholmen']['T1'],
                                buoy_dict['Tennholmen']['TZ'],
                                buoy_dict['Tennholmen']['T3'],
                                buoy_dict['Tennholmen']['Tc'],
                                buoy_dict['Tennholmen']['Tdw'],
                                buoy_dict['Tennholmen']['Tp'],
                                buoy_dict['Tennholmen']['Qp']),
                        unpack=True)
    # time conversion
    time_dt = [basetime + timedelta(seconds=time_s[i]) \
                for i in range(len(time_s))]
    print('cleaning up ...')
    t=os.system('rm -r tmp_Tennholmen')
    # convert lons, lats from radians to degree
    lons = lons * 180. / np.pi
    lats = lats * 180. / np.pi
    return time_s,time_dt,Hm0,Tm02,lons,lats,\
            Hs,TI,TE,T1,TZ,T3,Tc,Tdw,Tp,Qp

def read_buoy(buoyname,date):
    path=date.strftime(buoy_dict[buoyname]['path_template'])
    filename=date.strftime(buoy_dict[buoyname]['file_template'])
    basetime=buoy_dict[buoyname]['basetime']
    nc = netCDF4.Dataset(path+filename,mode='r')
    time_s = nc.variables[buoy_dict[buoyname]['time']][:]
    # time conversion
    time_dt = [basetime + timedelta(seconds=time_s[i]) \
                for i in range(len(time_s))]
    Hm0 = nc.variables[buoy_dict[buoyname]['Hm0']][:]
    Tm02 = nc.variables[buoy_dict[buoyname]['Tm02']][:]
    lats = nc.variables[buoy_dict[buoyname]['lats']][:]
    lons = nc.variables[buoy_dict[buoyname]['lons']][:]
    nc.close()
    ctime, cidx = matchtime(date,date,time_s,basetime=basetime)
    idx = cidx[0]
    Hm0 = Hm0[idx]
    Tm02 = Tm02[idx]
    lats = lats[idx]
    lons = lons[idx]
    time_s = time_s[idx]
    time_dt = time_dt[idx]
    return time_s, time_dt, Hm0, Tm02, lons, lats
    
def get_buoy(sdate,edate,buoyname=None,mode=None):
    if (buoyname == 'Tennholmen' or buoyname is None):
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
            try:
                sidx=time_dt_lst.index(sdate)
                eidx=time_dt_lst.index(edate) + 1
                time_s_lst = time_s_lst[sidx:eidx]
                time_dt_lst = time_dt_lst[sidx:eidx]
                Hm0_lst = Hm0_lst[sidx:eidx]
                Tm02_lst = Tm02_lst[sidx:eidx]
                lons_lst = lons_lst[sidx:eidx]
                lats_lst = lats_lst[sidx:eidx]
            except ValueError:
                print('Date not available')
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
                        = read_buoy(buoyname,tmpdate)
            time_s_lst.append(time_s)
            time_dt_lst.append(time_dt)
            Hm0_lst.append(Hm0)
            Tm02_lst.append(Tm02)
            lons_lst.append(lons)
            lats_lst.append(lats)
            tmpdate = tmpdate + timedelta(hours = 1)
        del tmpdate
        time_s_lst = np.array(time_s_lst)
        time_dt_lst = np.array(time_dt_lst)
        Hm0_lst = np.array(Hm0_lst)
        Tm02_lst = np.array(Tm02_lst)
        lons_lst = np.array(lons_lst)
        lats_lst = np.array(lats_lst)
    return time_s_lst, time_dt_lst, Hm0_lst, Tm02_lst, lons_lst, lats_lst

def get_buoy_ext(sdate,edate,buoyname=None,mode=None):
    # extended version of get_buoy
    if (sdate.month == edate.month and sdate.year == edate.year):
        time_s, time_dt, Hm0, Tm02, lons, lats, \
        Hs, TI, TE, T1, TZ, T3, Tc, Tdw, Tp, Qp \
                             = read_Tennholmen_ext(sdate)
        sidx=time_dt.index(sdate)
        eidx=time_dt.index(edate) + 1
        time_s_lst = time_s[sidx:eidx]
        time_dt_lst = time_dt[sidx:eidx]
        Hm0_lst = Hm0[sidx:eidx]
        Tm02_lst = Tm02[sidx:eidx]
        lons_lst = lons[sidx:eidx]
        lats_lst = lats[sidx:eidx]
        Hs_lst = Hs[sidx:eidx]
        TI_lst = TI[sidx:eidx]
        TE_lst = TE[sidx:eidx]
        T1_lst = T1[sidx:eidx]
        TZ_lst = TZ[sidx:eidx]
        T3_lst = T3[sidx:eidx]
        Tc_lst = Tc[sidx:eidx]
        Tdw_lst = Tdw[sidx:eidx]
        Tp_lst = Tp[sidx:eidx]
        Qp_lst = Qp[sidx:eidx]
    else:
        tmpdate = deepcopy(sdate)
        time_s_lst = []
        time_dt_lst = []
        Hm0_lst = []
        Tm02_lst = []
        lons_lst = []
        lats_lst = []
        Hs_lst = []
        TI_lst = []
        TE_lst = []
        T1_lst = []
        TZ_lst = []
        T3_lst = []
        Tc_lst = []
        Tdw_lst = []
        Tp_lst = []
        Qp_lst = []
        while tmpdate <= edate:
            time_s, time_dt, Hm0, Tm02, lons, lats, \
            Hs, TI, TE, T1, TZ, T3, Tc, Tdw, Tp, Qp \
                            = read_Tennholmen_ext(tmpdate)
            time_s_lst.append(time_s)
            time_dt_lst.append(time_dt)
            Hm0_lst.append(Hm0)
            Tm02_lst.append(Tm02)
            lons_lst.append(lons)
            lats_lst.append(lats)
            Hs_lst.append(Hs)
            TI_lst.append(TI)
            TE_lst.append(TE)
            T1_lst.append(T1)
            TZ_lst.append(TZ)
            T3_lst.append(T3)
            Tc_lst.append(Tc)
            Tdw_lst.append(Tdw)
            Tp_lst.append(Tp)
            Qp_lst.append(Qp)
            tmpdate = tmpdate + relativedelta(months = +1)
        del tmpdate
        time_s_lst = flatten(time_s_lst)
        time_dt_lst = flatten(time_dt_lst)
        Hm0_lst = flatten(Hm0_lst)
        Tm02_lst = flatten(Tm02_lst)
        lons_lst = flatten(lons_lst)
        lats_lst = flatten(lats_lst)
        Hs_lst = flatten(Hs_lst)
        TI_lst = flatten(TI_lst)
        TE_lst = flatten(TE_lst)
        T1_lst = flatten(T1_lst)
        TZ_lst = flatten(TZ_lst)
        T3_lst = flatten(T3_lst)
        Tc_lst = flatten(Tc_lst)
        Tdw_lst = flatten(Tdw_lst)
        Tp_lst = flatten(Tp_lst)
        Qp_lst = flatten(Qp_lst)
        sidx = time_dt_lst.index(sdate)
        eidx = time_dt_lst.index(edate) + 1
        time_s_lst = time_s_lst[sidx:eidx]
        time_dt_lst = time_dt_lst[sidx:eidx]
        Hm0_lst = Hm0_lst[sidx:eidx]
        Tm02_lst = Tm02_lst[sidx:eidx]
        lons_lst = lons_lst[sidx:eidx]
        lats_lst = lats_lst[sidx:eidx]
        Hs_lst = Hs_lst[sidx:eidx]
        TI_lst = TI_lst[sidx:eidx]
        TE_lst = TE_lst[sidx:eidx]
        T1_lst = T1_lst[sidx:eidx]
        TZ_lst = TZ_lst[sidx:eidx]
        T3_lst = T3_lst[sidx:eidx]
        Tc_lst = Tc_lst[sidx:eidx]
        Tdw_lst = Tdw_lst[sidx:eidx]
        Tp_lst = Tp_lst[sidx:eidx]
        Qp_lst = Qp_lst[sidx:eidx]
    return time_s_lst, time_dt_lst, Hm0_lst, Tm02_lst,\
            lons_lst, lats_lst, Hs_lst,TI_lst,TE_lst,T1_lst,\
            TZ_lst,T3_lst,Tc_lst,Tdw_lst,Tp_lst,Qp_lst

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
#    import d22_var_dicts
    with open("/home/patrikb/wavy/wavy/d22_var_dicts.yaml", 'r') as stream:
        d22_var_dicts=yaml.safe_load(stream)
#    try:
#        if sys.version_info <= (3, 0):
#            reload(d22_var_dicts)
#        else:
#            import importlib
#            importlib.reload(module)
#    except:
#        pass
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
                if var == 'Tm_10min':
                    WINDOW = 3
                    wmt = sp.power(WM[var],2)
                    weights = sp.ones((WINDOW))/WINDOW
                    wmt[2:-2:2] = np.convolve(wmt[0::2], weights,mode='valid')
                    wmt[3:-2:2] = np.convolve(wmt[1::2], weights,mode='valid')
                    wmt[[0,1,-2,-1]]= sp.nan
                    WM['Tm_1hr'] = sp.sqrt(wmt)
    # CAUTION: 10min data is extracted for entire days only 00:00h - 23:50h
    sensor_lst = [WM1,WM2,WM3]
    dates = deepcopy(dat)
    del dat
    return sensor_lst, dates

def matchtime(sdate,edate,time,basetime=None,timewin=None):
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
    if (edate is None or sdate==edate):
        for element in time:
            if basetime is None:
                tmp = element
            else:
                tmp = basetime + timedelta(seconds=element)
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
                tmp = basetime + timedelta(seconds=element)
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
