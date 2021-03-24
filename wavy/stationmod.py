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
# standard library imports
import os
import sys
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import datetime
import argparse
from argparse import RawTextHelpFormatter
import yaml
import os
import urllib
import gzip
import ftplib
from ftplib import FTP
import calendar
import sys
from dateutil.relativedelta import relativedelta
from copy import deepcopy
import time
import pylab as pl
from datetime import datetime
import scipy as sp
# own imports
from ncmod import ncdumpMeta, get_varname_for_cf_stdname_in_ncfile
from ncmod import dumptonc_ts_station
from utils import collocate_times
# ---------------------------------------------------------------------#

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

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)

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
    def __init__(self,platform,sensor,sdate,edate,
                mode='d22',varalias='Hs'):
        print ('# ----- ')
        print (" ### Initializing station_class object ###")
        print ('Chosen period: ' + str(sdate) + ' - ' + str(edate))
        print (" Please wait ...")
        stdvarname = variable_info[varalias]['standard_name']
        try:
            var, time, timedt, \
            pathtofile = self.get_station( platform,
                                           sdate,edate,
                                           mode,
                                           sensor,
                                           varalias )
            vardict = {
                    stdvarname:var,
                    'time':time,
                    'datetime':timedt,
                    'time_unit':self.time_unit,
                    'longitude':[ station_dict['platform'][platform]\
                                  ['coords']['lon'] ]*len(var),
                    'latitude':[ station_dict['platform'][platform]\
                                 ['coords']['lat'] ]*len(var)
                    }
            self.vars = vardict
            self.stdvarname = stdvarname
            if mode == 'd22':
                self.varname = varalias
            elif mode == 'nc':
                model_meta = ncdumpMeta(pathtofile)
                self.vars['model_meta'] = model_meta
                self.varname = get_varname_for_cf_stdname_in_ncfile( 
                                model_meta,stdvarname)
            self.varalias = varalias
            self.sdate = sdate
            self.edate = edate
            self.lat = station_dict['platform'][platform]\
                                   ['coords']['lat']
            self.lon = station_dict['platform'][platform]\
                                   ['coords']['lon']
            self.platform = platform
            self.sensor = sensor
            print (" ### station_class object initialized ###")
        except Exception as e:
            print(e)
            self.error = True
            print ("! no station_class object initialized !")
        print ('# ----- ')

    def get_station(self,platform,sdate,edate,mode,sensor,varalias):
        stdvarname = variable_info[varalias]['standard_name']
        path_template = station_dict['path']['platform']['local']\
                                    [mode]['path_template']
        file_template = station_dict['path']['platform']['local']\
                                    [mode]['file_template']
        pathlst = [p + ('/' + file_template) for p in path_template]
        strsublst = station_dict['path']['platform']['local'][mode]['strsub']
        pathtofile = None
        if mode == 'nc':
            tmpdate = sdate
            var = []
            time = []
            timedt = []
            while (tmpdate <= edate):
                pathtofile = get_pathtofile( platform,sensor,varalias,
                                             pathlst,strsublst,tmpdate )
                ncdict = ncdumpMeta(pathtofile)
                varname = get_varname_for_cf_stdname_in_ncfile(ncdict,
                                                               stdvarname)
                if len(varname) > 1:
                    # overwrite like for satellite altimetry files
                    varname = station_dict['platform']['misc']\
                                          ['vardef'][stdvarname]
                else:
                    varname = varname[0]
                nc = netCDF4.Dataset(pathtofile,'r')
                var.append(nc.variables[varname][:])
                timeobj = nc.variables['time']
                time.append(nc.variables['time'][:])
                timedt.append(netCDF4.num2date(timeobj[:],timeobj.units))
                nc.close()
                tmpdate = tmpdate + relativedelta(months = +1)
            var = np.array(var).flatten()
            time = np.array(time).flatten()
            timedt = np.array(timedt).flatten()
        elif mode == 'd22':
            sdatetmp = sdate - timedelta(days=1)
            edatetmp = edate + timedelta(days=1)
            sl = parse_d22(platform,sensor,varalias,sdatetmp,edatetmp,
                          pathlst,strsublst,mode)
            var, timedt = extract_d22(sl,varalias,platform,sensor)
            time = np.array(
                    [(t-self.basedate).total_seconds() for t in timedt]
                    )
        idxtmp = collocate_times(unfiltered_t=timedt,
                                sdate=sdate,edate=edate,twin=1)
        # convert to list for consistency with other classes
        var = list(np.real(var[idxtmp]))
        time = list(time[idxtmp])
        timedt = list(timedt[idxtmp])
        return var, time, timedt, pathtofile
    
    def write_to_monthly_nc(self,path=None,filename=None):
        # divide time into months by loop over months from sdate to edate
        if 'error' in vars(self):
            print('Erroneous station_class file detected')
            print('--> dump to netCDF not possible !')
        else:
            tmpdate = self.sdate
            edate = self.edate
            while tmpdate <= edate:
                idxtmp = collocate_times(unfiltered_t=self.vars['datetime'],
                                     sdate = datetime(tmpdate.year,
                                                      tmpdate.month,1),
                                     edate = datetime(tmpdate.year,
                                                      tmpdate.month,
                                                      calendar.monthrange(
                                                        tmpdate.year,
                                                        tmpdate.month)[1],
                                                        23,59) )
                if (path is not None and filename is not None):
                    pathtofile = path + '/' + filename
                else:
                    if path is None:
                        path_template = station_dict['path']['platform']\
                                                    ['local']['nc']\
                                                    ['path_template'][0]
                    if filename is None:
                        file_template = station_dict['path']['platform']\
                                                    ['local']['nc']\
                                                    ['file_template']
                    strsublst = station_dict['path']['platform']\
                                                    ['local']['nc']\
                                                    ['strsub']
                    tmppath = path_template + '/' + file_template
                    pathtofile = make_pathtofile(self.platform,self.sensor,
                                          self.varalias,tmppath,strsublst,
                                          tmpdate)
                title = ( self.varname + ' observations from ' 
                        + self.platform + ' ' + self.sensor )
                dumptonc_ts_station(self,pathtofile,title)
                tmpdate = tmpdate + relativedelta(months = +1)
        return

def get_pathtofile(platform,sensor,varalias,pathlst,strsublst,date):
    i = 0
    pathtofile = date.strftime(pathlst[i])
    for strsub in strsublst:
        pathtofile = pathtofile.replace(strsub,locals()[strsub])
    while os.path.isfile(pathtofile) is False:
        i += 1
        pathtofile = date.strftime(pathlst[i])
        for strsub in strsublst:
            pathtofile = pathtofile.replace(strsub,locals()[strsub])
    return pathtofile

def make_pathtofile(platform,sensor,varalias,tmppath,strsublst,date):
    pathtofile = date.strftime(tmppath)
    for strsub in strsublst:
        pathtofile = pathtofile.replace(strsub,locals()[strsub])
    return pathtofile

def compute_superobs(st_obj,smoother='running_mean',**kwargs):
    """
    Applies a smoothing filter to create a super-observed ts
    **kwargs includes method specific input for chosen smoother
    Smoother on wish list are:
            block-average
            running mean using convolution
            GP
            GAM
            ...
    Caution:    for some smoothers much more of time series has 
                to be included.
    """
    print('under construction')
    return

def parse_d22(platform,sensor,varalias,sdate,edate,pathlst,strsublst,mode):
    """
    Read all lines in file and append to sl
    """
    sl=[]
    for d in range(int(pl.date2num(sdate)),int(pl.date2num(edate))): 
        pathtofile = get_pathtofile(platform,sensor,varalias,\
                                    pathlst,strsublst,pl.num2date(d))
        print('Parsing:', pathtofile)
        f = open(pathtofile, "r")
        sl = sl + f.readlines()
        f.close()
    return sl

def floater(s):
    """
    Function that converts 's' to float32 or nan if floater throws exception 
    """
    try:
        x = np.float32(s)
    except:
        x = np.nan
    return x

def find_category_for_variable(varalias):
    lst = [ i for i in d22_var_dicts.keys() \
            if (varalias in d22_var_dicts[i]) ]
    if len(lst) == 1:
        return lst[0]
    else: return None

def get_revised_categories(sl,category):
    """
    finds number of occurences of string (category) to determine
    revised_categories (type: list)
    """
    revised_categories = []
    count = 1
    newsl = sl[1:(sl[1:-1].index('!!!!\n'))]
    for element in newsl: 
        if category in element:
            revised_categories.append(category+str(count))
            count+=1
    return revised_categories

def extract_d22(sl,varalias,platform,sensor):
    """
    Extracting data of choice - reading sl from parse_d22
    CAUTION: 10min data is extracted for entire days only 00:00h - 23:50h
    Returns values of chosen variable (ts) and corresponding datetimes (dt)
    as type: np.array
    """
    print('Extracting data from parsed .d22-files')
    category = find_category_for_variable(varalias)
    revised_categories = get_revised_categories(sl,category)
    ts = []
    dt = []
    for i, line in enumerate(sl):
        # get ts for date and time
        if "!!!!" in line:
            datestr = sl[  i
                         + d22_var_dicts['datetime']['date']['idx']
                        ].strip()
            timestr = sl[  i
                         + d22_var_dicts['datetime']['time']['idx']
                        ].strip()
            date_object = datetime.strptime(datestr 
                                            + ' ' 
                                            + timestr,
                                            '%d-%m-%Y %H:%M')
            dt.append(date_object)
        # get ts for variable of interest
        revised_category_for_sensor = revised_categories[
                                        station_dict['platform']\
                                        [platform]['sensor']\
                                        [sensor] ]
        if revised_category_for_sensor in line:
            value = sl[  i
                       + d22_var_dicts[category][varalias]['idx']
                       ].strip()
            ts.append(floater(value))
    #Convert data to arrays
    dt = np.array(dt)
    ts = np.array(ts)
    return ts, dt

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
