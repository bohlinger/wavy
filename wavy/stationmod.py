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
from wavy.ncmod import ncdumpMeta, get_varname_for_cf_stdname_in_ncfile
from wavy.ncmod import dumptonc_ts_station, get_varlst_from_nc_1D
from wavy.ncmod import get_filevarname_from_nc
from wavy.utils import collocate_times
from wavy.utils import make_pathtofile, get_pathtofile
from wavy.utils import convert_meteorologic_oceanographic
from wavy.superobmod import superobbing
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#

# read yaml config files:
buoy_dict = load_or_default('buoy_specs.yaml')
insitu_dict = load_or_default('insitu_specs.yaml')
variable_info = load_or_default('variable_info.yaml')
d22_dict = load_or_default('d22_var_dicts.yaml')
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
    Class to handle insitu based time series.
    '''
    basedate = datetime(1970,1,1)
    time_unit = 'seconds since 1970-01-01 00:00:00.0'
    def __init__(self,nID,sensor,sdate,edate,fifo='d22',
    varalias='Hs',superobserve = False, **kwargs):
        print ('# ----- ')
        print (" ### Initializing station_class object ###")
        print ('Chosen period: ' + str(sdate) + ' - ' + str(edate))
        print (" Please wait ...")
        stdvarname = variable_info[varalias]['standard_name']
        try:
#        for i in range(1):
            self.stdvarname = stdvarname
            self.varalias = varalias
            self.sdate = sdate
            self.edate = edate
            self.lat = insitu_dict[fifo][nID]['coords']\
                                   [sensor]['lat']
            self.lon = insitu_dict[fifo][nID]['coords']\
                                   [sensor]['lon']
            self.nID = nID
            self.sensor = sensor
            if ('tags' in insitu_dict[fifo].keys and
            len(insitu_dict[fifo]['tags']>0)):
                self.tags = insitu_dict[fifo]['tags']
            if superobserve == False:
                var, time, timedt, lon, lat, pathtofile = \
                    self.get_insitu_ts(\
                                nID, sensor,sdate,edate,fifo,
                                varalias,**kwargs)
                vardict = {
                    stdvarname:var,
                    'time':time,
                    'datetime':timedt,
                    'time_unit':self.time_unit,
                    'longitude':lon,
                    'latitude':lat
                    }
            elif superobserve == True:
                print(kwargs)
                # determine start and end date
                if 'stwin' not in kwargs.keys():
                    kwargs['stwin'] = 3
                if 'etwin' not in kwargs.keys():
                    kwargs['etwin'] = 0
                sdate_new = sdate - timedelta(hours=kwargs['stwin'])
                edate_new = edate + timedelta(hours=kwargs['etwin'])
                var, time, timedt, lon, lat, pathtofile = \
                    self.get_insitu_ts(nID,sensor,
                                sdate_new,edate_new,fifo,
                                varalias,**kwargs)
                tmp_vardict = {
                    stdvarname:var,
                    'time':time,
                    'datetime':timedt,
                    'time_unit':self.time_unit,
                    'longitude':lon,
                    'latitude':lat
                    }
                vardict = superobbing(varalias,tmp_vardict,**kwargs)
                # cut to original sdate and edate
                time_cut = np.array(vardict['time'])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ]
                var_cut = np.array(vardict[stdvarname])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ]
                lon_cut = np.array(vardict['longitude'])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ]
                lat_cut = np.array(vardict['latitude'])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ]
                dtime_cut = np.array(vardict['datetime'])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)) ]
                vardict['time'] = list(time_cut)
                vardict['datetime'] = list(dtime_cut)
                vardict[stdvarname] = list(var_cut)
                vardict['longitude'] = list(lon_cut)
                vardict['latitude'] = list(lat_cut)
                self.superob = kwargs['superob']
                self.outlier_detection = kwargs['outlier_detection']
                self.missing_data = kwargs['missing_data']
            self.vars = vardict
            if fifo == 'nc':
                meta = ncdumpMeta(pathtofile)
                self.vars['meta'] = meta
                self.varname = get_varname_for_cf_stdname_in_ncfile(
                                meta,stdvarname)
            else:
                self.varname = varalias
            print (" ### station_class object initialized ###")
        except Exception as e:
            print(e)
            self.error = e
            print ("! No station_class object initialized !")
        print ('# ----- ')

    def get_insitu_ts(self,nID,sensor,sdate,edate,fifo,
    varalias,basedate,**kwargs):
        stdvarname = variable_info[varalias]['standard_name']
        path_template = insitu_dict[fifo][nID]['src']['path_template']
        file_template = insitu_dict[fifo][nID]['src']['file_template']
        pathlst = [p + ('/' + file_template) for p in path_template]
        strsublst = insitu_dict[fifo][nID]['src']['strsub']
        pathtofile = None
        if mode == 'nc':
            tmpdate = sdate
            var = []
            time = []
            timedt = []
            while (datetime(tmpdate.year,tmpdate.month,1) \
            <= datetime(edate.year,edate.month,1)):
                pathtofile = get_pathtofile(pathlst,strsublst,tmpdate,
                                            platform=platform,
                                            sensor=sensor,
                                            varalias=varalias)
                print('Parsing:',pathtofile)
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
                print(tmpdate)
            var = flatten(var)
            time = flatten(time)
            timedt = flatten(timedt)
            # convert to datetime object
            timedt = [datetime(t.year,t.month,t.day,t.hour,t.minute,t.second)
                    for t in timedt]
            # remove duplicates
            timedt, \
            unique_time_idx = list(np.unique(timedt,return_index=True)[0]),\
                              list(np.unique(timedt,return_index=True)[1])
            var = list(np.array(var)[unique_time_idx])
            time = list(np.array(time)[unique_time_idx])
        elif mode == 'd22':
            sdatetmp = sdate
            edatetmp = edate
            sl = parse_d22(ID,sensor,varalias,sdatetmp,edatetmp,
                          pathlst,strsublst)
            var, timedt = extract_d22(sl,varalias,project,ID,sensor)
            time = np.array(
                    [(t-self.basedate).total_seconds() for t in timedt]
                    )
        if 'twin' in insitu_dict[project][ID]:
            idxtmp = collocate_times(unfiltered_t=timedt,\
                                sdate=sdate,edate=edate,
                                twin=insitu_dict[project][ID]['twin'])
        else:
            # default to allow for a 1 min variation
            idxtmp = collocate_times(unfiltered_t=timedt,\
                                sdate=sdate,edate=edate,
                                twin=1)
        # convert to list for consistency with other classes
        # and make sure that only steps with existing obs are included
        time = [time[i] for i in idxtmp if i < len(var)]
        timedt = [timedt[i] for i in idxtmp if i < len(var)]
        var = [np.real(var[i]) for i in idxtmp if i < len(var)]
        # rm double entries due to 10min spacing
        if ('unique' in kwargs.keys() and kwargs['unique'] is True):
            # delete 10,30,50 min times, keep 00,20,40
            # 1. create artificial time vector for collocation
            tmpdate = deepcopy(sdate)
            tmpdatelst = []
            while tmpdate<edate:
                tmpdatelst.append(tmpdate)
                tmpdate += timedelta(minutes=20)
            # 2. collocate times
            if 'twin' in insitu_dict[project][ID]:
                idxtmp = collocate_times(\
                            unfiltered_t=timedt,
                            target_t=tmpdatelst,
                            twin=insitu_dict[project][ID]['twin'])
            else:
                idxtmp = collocate_times(unfiltered_t=timedt,\
                                target_t=tmpdatelst,
                                twin=1)
            time = list(np.array(time)[idxtmp])
            timedt = list(np.array(timedt)[idxtmp])
            var = list(np.array(var)[idxtmp])
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
                    if 'superob' in vars(self).keys():
                        file_template = 'superobbed_' + file_template
                    tmppath = path_template + '/' + file_template
                    pathtofile = make_pathtofile(tmppath,strsublst,
                                                tmpdate,
                                                platform=self.platform,
                                                sensor=self.sensor,
                                                varalias=self.varalias)
                title = ( self.varname + ' observations from '
                        + self.platform + ' ' + self.sensor )
                dumptonc_ts_station(self,pathtofile,title)
                tmpdate = tmpdate + relativedelta(months = +1)
        return


def get_insitu_ts(nID,sensor,sdate,edate,fifo,varalias,
basedate,**kwargs):
    stdvarname = variable_info[varalias]['standard_name']
    path_template = insitu_dict[fifo][nID]['src']['path_template']
    file_template = insitu_dict[fifo][nID]['src']['file_template']
    pathlst = [p + ('/' + file_template) for p in path_template]
    strsublst = insitu_dict[fifo][nID]['src']['strsub']
    pathtofile = None
    if fifo == 'nc':
        var, time, timedt, lons, lats = \
            get_nc_ts(nID,sensor,varalias,sdate,edate,pathlst,strsublst)
    elif fifo == 'd22':
        var, time, timedt = \
            get_d22_ts(sdate,edate,basedate,fifo,nID,sensor,varalias,\
                        pathlst,strsublst)
        if 'twin' in insitu_dict[fifo][nID]:
            idxtmp = collocate_times(unfiltered_t=timedt,\
                                sdate=sdate,edate=edate,
                                twin=insitu_dict[fifo][nID]['twin'])
        else:
            # default to allow for a 1 min variation
            idxtmp = collocate_times(unfiltered_t=timedt,\
                                sdate=sdate,edate=edate,
                                twin=1)
        # convert to list for consistency with other classes
        # and make sure that only steps with existing obs are included
        time = [time[i] for i in idxtmp if i < len(var)]
        timedt = [timedt[i] for i in idxtmp if i < len(var)]
        var = [np.real(var[i]) for i in idxtmp if i < len(var)]
        # rm double entries due to 10min spacing
        if ('unique' in kwargs.keys() and kwargs['unique'] is True):
            # delete 10,30,50 min times, keep 00,20,40
            # 1. create artificial time vector for collocation
            tmpdate = deepcopy(sdate)
            tmpdatelst = []
            while tmpdate<edate:
                tmpdatelst.append(tmpdate)
                tmpdate += timedelta(minutes=20)
            # 2. collocate times
            if 'twin' in insitu_dict[fifo][nID]:
                idxtmp = collocate_times(\
                            unfiltered_t=timedt,
                            target_t=tmpdatelst,
                            twin=insitu_dict[fifo][nID]['twin'])
            else:
                idxtmp = collocate_times(unfiltered_t=timedt,\
                                target_t=tmpdatelst,
                                twin=1)
            time = list(np.array(time)[idxtmp])
            timedt = list(np.array(timedt)[idxtmp])
            var = list(np.array(var)[idxtmp])
        lons = [insitu_dict[fifo][nID]['coords'][sensor]['lon']]\
               *len(var)
        lats = [insitu_dict[fifo][nID]['coords'][sensor]['lat']]\
               *len(var)
    return var, time, timedt, lons, lats

def get_d22_ts(sdate,edate,basedate,fifo,nID,sensor,varalias,
pathlst,strsublst):
    sdatetmp = sdate
    edatetmp = edate
    sl = parse_d22(nID,sensor,varalias,sdatetmp,edatetmp,
                   pathlst,strsublst)
    var, timedt = extract_d22(sl,varalias,fifo,nID,sensor)
    time = np.array(
           [(t-basedate).total_seconds() for t in timedt]
           )
    return var, time, timedt

def get_nc_ts(nID,sensor,varalias,sdate,edate,pathlst,strsublst):
    # loop from sdate to edate with dateincr
    tmpdate = deepcopy(sdate)
    varlst = []
    lonlst = []
    latlst = []
    timelst = []
    dtimelst = []
    count = 0
    while tmpdate <= edate:
        flatten_switch += 1
        # make pathtofile
        pathtofile = get_pathtofile(pathlst,strsublst,tmpdate,nID=nID)
        # get ncdump
        ncdict = ncdumpMeta(pathtofile)
        # retrieve filevarname for varalias
        filevarname = get_filevarname_from_nc(varalias,
                                          variable_info,
                                          insitu_dict['nc'][nID],
                                          ncdict)
        varstrlst = [filevarname,'longitude','latitude','time']
        # query
        vardict = get_varlst_from_nc_1D(pathtofile,
                                        varstrlst,
                                        sdate,edate)
        varlst.append(list(vardict[filevarname]))
        lonlst.append(list(vardict['longitude']))
        latlst.append(list(vardict['latitude']))
        timelst.append(list(vardict['time']))
        dtimelst.append(list(vardict['dtime']))
        # determine date increment
        date_incr = insitu_dict['nc'][nID]['src']['date_incr']
        if date_incr == 'm':
            tmpdate += relativedelta(months = +1)
        elif date_incr == 'Y':
            tmpdate += relativedelta(years = +1)
        elif date_incr == 'd':
            tmpdate += timedelta(days = +1)
    if count > 1:
        # flatten lists
    #turn timedt into datetime objects
    return varlst, timelst, dtimelst, lonlst, latlst

def parse_d22(nID,sensor,varalias,sdate,edate,pathlst,strsublst):
    """
    Read all lines in file and append to sl
    """
    sl=[]
    for d in range(int(pl.date2num(sdate))-1,int(pl.date2num(edate))+2):
        try:
            pathtofile = get_pathtofile(pathlst,strsublst,
                                        pl.num2date(d),
                                        varalias=varalias,
                                        nID=nID,sensor=sensor)
            print('Parsing:', pathtofile)
            f = open(pathtofile, "r")
            sl = sl + f.readlines()
            f.close()
        except Exception as e:
            print('Error in parse_d22:')
            print(e)
    return sl

# flatten all lists before returning them
# define flatten function for lists
''' fct does the following:
flat_list = [item for sublist in TIME for item in sublist]
or:
for sublist in TIME:
for item in sublist:
flat_list.append(item)
'''
flatten = lambda l: [item for sublist in l for item in sublist]

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
    lst = [ i for i in d22_dict.keys() \
            if (varalias in d22_dict[i]) ]
    if len(lst) == 1:
        return lst[0]
    else: return None

def get_revised_categories(sl,category):
    """
    finds number of occurences of string (category) to determine
    revised_categories (type: list)
    """
    revised_categories = []
    idxlst = []
    count = 1
    searching = True
    while (searching is True or count<10):
        revised_category = category+str(count)
        if find_category(sl,revised_category) is True:
            revised_categories.append(revised_category)
            idxlst.append(count-1)
            count += 1
        else:
            searching = False
            count += 1
    return revised_categories,idxlst

def find_category(sl,category):
    for element in sl:
        if category in element:
            return True

def check_sensor_availability(revised_categories,idxlst,fifo,nID,sensor):
    idxyaml = insitu_dict[fifo][nID]['sensor'][sensor]
    if idxyaml in idxlst:
        return idxlst.index(idxyaml)
    else:
        return None

def extract_d22(sl,varalias,fifo,nID,sensor):
    """
    Extracting data of choice - reading sl from parse_d22
    CAUTION: 10min data is extracted for entire days only 00:00h - 23:50h
    Returns values of chosen variable (ts) and corresponding datetimes (dt)
    as type: np.array
    """
    print('Extracting data from parsed .d22-files')
    category = find_category_for_variable(varalias)
    revised_categories,idxlst = get_revised_categories(sl,category)
    print( 'Consistency check: \n'
           ' --> compare found #sensors against defined in insitu_specs.yaml')
    sensornr = len(insitu_dict[fifo][nID]['sensor'].keys())
    if len(revised_categories) == sensornr:
        print('Consistency check: OK!')
    else:
        print('Consistency check: Failed!')
        print(    '!!! Caution:\n'
                + 'found #sensor ('
                + str(len(revised_categories))
                + ') is not equal to defined ' +
                '#sensors ('
                + str(sensornr)
                + ') in insitu_specs.yaml')
    # check that the defined sensors are actually the ones being found
    check = check_sensor_availability(revised_categories,\
                                      idxlst,fifo,nID,sensor)
    ts = []
    dt = []
    if check is not None:
        print('Sensor is available and defined in insitu_specs.yaml')
        for i, line in enumerate(sl):
            # get ts for date and time
            if "!!!!" in line:
                datestr = sl[  i
                             + d22_dict['datetime']['date']['idx']
                            ].strip()
                timestr = sl[  i
                             + d22_dict['datetime']['time']['idx']
                            ].strip()
                date_object = datetime.strptime(datestr
                                                + ' '
                                                + timestr,
                                                '%d-%m-%Y %H:%M')
                dt.append(date_object)
            # get ts for variable of interest
            revised_category_for_sensor = revised_categories[check]
            #print(revised_category_for_sensor)
            if revised_category_for_sensor in line:
                value = sl[  i
                            + d22_dict[category][varalias]['idx']
                           ].strip()
                ts.append(floater(value))
    else:
        print('Caution: Sensor is not defined or available')
    #Convert data to arrays
    dt = np.array(dt)
    ts = np.array(ts)
    # adjust conventions
    if ('convention' in d22_dict[category][varalias].keys() and
    d22_dict[category][varalias]['convention'] == 'meteorologic'):
        print('Convert from meteorologic to oceanographic convention')
        ts = convert_meteorologic_oceanographic(ts)
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
