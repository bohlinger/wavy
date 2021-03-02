#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field from model output. I try to mostly follow the PEP convention for 
python code style. Constructive comments on style and effecient 
programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
import sys
import yaml
import os

# read files
import netCDF4

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
#import os
import math

# progress bar
from utils import progress, hour_rounder

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# matchtime fct
from stationmod import matchtime

# netcdf specific
from ncmod import ncdumpMeta, get_varname_for_cf_stdname_in_ncfile

# 1: get_model for given time period
# 2: dumptonc based on model (e.g. MWAM4, ARCMFC, ARCMFCnew)
# 3: choose create or append based on the existence of the file
# Must have one unlimited dimension (time), and two spatial dimensions
#   (latitude, longitude, which depend on rlat,rlon)

# --- global functions ------------------------------------------------#
"""
definition of some global functions
"""
def get_model_filedate(model,fc_date,leadtime):
    '''
    get init_date for latest model output file and checks if available
    '''
    if ('init_times' in model_dict[model].keys() and 
    model_dict[model]['init_times'] is not None):
        init_times = \
            np.array(model_dict[model]['init_times']).astype('float')
    else:
        print('init_times for chosen model not specified in config file')
        print('Assuming continuous simulation with hourly values')
        init_times = np.array(range(25)).astype('float')
    date = fc_date - timedelta(hours=leadtime)
    if date.hour in init_times:
        init_diffs = date.hour - init_times
        init_diffs[init_diffs<0] = np.nan
        h_idx = np.where(init_diffs==\
                        np.min(init_diffs[~np.isnan(init_diffs)]))
        h = int(init_times[h_idx[0][0]])
        return datetime(date.year,date.month,date.day,h)
    else:
        raise ValueError('!!! leadtime not available !!!')

def make_model_filename(model=None,fc_date=None,leadtime=None):
    """
    creates/returns filename based on fc_date,init_date,leadtime
    """
    if model in model_dict:
        if 'xtra_h' in model_dict[model]:
            filedate = get_model_filedate(model,fc_date,leadtime)
            tmpstr = model_dict[model]['file_template']
            for i in range(model_dict[model]['nr_filedates']):
                filedatestr = model_dict[model]['filedate_formats'][i]
                replacestr = (filedate \
                            + timedelta(hours = \
                                        model_dict[model]['xtra_h'][i])).\
                            strftime(filedatestr)
                tmpstr = tmpstr.replace('filedate',replacestr,1)
            filename = (filedate.strftime(model_dict[model]['path_template'])
                        + tmpstr)
        else:
            filedate = get_model_filedate(model,fc_date,leadtime)
            filename = (
                filedate.strftime(model_dict[model]['path_template'])
                + filedate.strftime(model_dict[model]['file_template'])
                )
    else:
        raise ValueError("chosen model is not specified in model_specs.yaml")
    return filename

def get_model_fc_mode(filestr,model,fc_date,leadtime=None,varalias=None):
    """ 
    fct to retrieve model data for correct time
    """
    vardict = {}
    print ("Get model data according to selected date ....")
    print(filestr)
    t0=time.time()
    model_meta = ncdumpMeta(filestr)
    t1=time.time()
    print('------- TIMER t1-t0: ',t1-t0,'----------')
    f = netCDF4.Dataset(filestr,'r')
    stdvarname = variable_info[varalias]['standard_name']
    # get coordinates and time
    lonsname = get_filevarname(model,'lons',variable_info,
                                model_dict,model_meta)
    latsname = get_filevarname(model,'lats',variable_info,
                                model_dict,model_meta)
    timename = get_filevarname(model,'time',variable_info,
                                model_dict,model_meta)
    t2=time.time()
    print('------- TIMER t2-t1: ',t2-t1,'----------')
    model_lons = f.variables[lonsname][:]
    t3=time.time()
    print('------- TIMER t3-t2: ',t3-t2,'----------')
    vardict[variable_info['lons']['standard_name']]=model_lons
    t4=time.time()
    print('------- TIMER t4-t3: ',t4-t3,'----------')
    model_lats = f.variables[latsname][:]
    vardict[variable_info['lats']['standard_name']]=model_lats
    model_time = f.variables[timename]
    t5=time.time()
    print('------- TIMER t5-t4: ',t5-t4,'----------')
    # get other variables e.g. Hs [time,lat,lon]
    filevarname = get_filevarname(model,varalias,variable_info,
                                    model_dict,model_meta)
    t6=time.time()
    print('------- TIMER t6-t5: ',t6-t5,'----------')
    if (type(filevarname) is dict):
        print('Target variable can be computed from vector \n' 
              'components with the following aliases: ', filevarname)
        model_time_dt = list( netCDF4.num2date(model_time[:],
                                           units = model_time.units) )
        model_time_dt_valid = [model_time_dt[model_time_dt.index(fc_date)]]
        model_time_valid = [model_time[model_time_dt.index(fc_date)]]
        model_time_unit = model_time.units
        vardict[variable_info['time']['standard_name']]=model_time_valid
        vardict['datetime']=model_time_dt_valid
        vardict['time_unit']=model_time_unit
        for key in filevarname.keys():
            filevarname_dummy = get_filevarname(model,
                                    filevarname[key][0],
                                    variable_info,
                                    model_dict,model_meta)
            if filevarname_dummy is not None:
                print(filevarname[key][0], 'exists')
                break
        print('Use aliases:',filevarname[key])
        model_var_dummy = f.variables[filevarname_dummy]
        if len(model_var_dummy.shape)>2: # for multiple time steps
            model_var_valid_tmp = \
                model_var_dummy[model_time_dt.index(fc_date),:,:].squeeze()**2
            for i in range(1,len(filevarname[key])):
                filevarname_dummy = get_filevarname(model,
                                                    filevarname[key][i],
                                                    variable_info,
                                                    model_dict,
                                                    model_meta)
                model_var_valid_tmp += \
                    f.variables[filevarname_dummy][
                        model_time_dt.index(fc_date),:,:
                        ].squeeze()**2
            model_var_valid = np.sqrt(model_var_valid_tmp)
        else:# if only one time step
            model_var_valid_tmp = model_var_dummy[:,:].squeeze()**2
            for i in range(1,len(filevarname[key])):
                filevarname_dummy = get_filevarname(model,
                                                    filevarname[key][i],
                                                    variable_info,
                                                    model_dict,
                                                    model_meta)
                model_var_valid_tmp += \
                    f.variables[filevarname_dummy][:,:].squeeze()**2
            model_var_valid = np.sqrt(model_var_valid_tmp)
        vardict[stdvarname] = model_var_valid
    else:
        model_time_dt = list( netCDF4.num2date(model_time[:],
                                           units = model_time.units) )
        model_time_dt_valid = [model_time_dt[model_time_dt.index(fc_date)]]
        model_time_valid = [model_time[model_time_dt.index(fc_date)]]
        model_time_unit = model_time.units
        vardict[variable_info['time']['standard_name']]=model_time_valid
        vardict['datetime']=model_time_dt_valid
        vardict['time_unit']=model_time_unit
        model_var_link = f.variables[filevarname]
        if len(model_var_link.shape)>2: # for multiple time steps
            model_var_valid = \
                model_var_link[model_time_dt.index(fc_date),:,:].squeeze()
        else:# if only one time step
            model_var_valid = model_var_link[:,:].squeeze()
        vardict[variable_info[varalias]['standard_name']] = \
                                                    model_var_valid
    f.close()
    vardict['model_meta'] = model_meta
    t7=time.time()
    print('------- TIMER t7-t6: ',t7-t6,'----------')
    print('------- TIMER t7-t0: ',t7-t0,'----------')
    return vardict, filevarname

def make_dates_and_lt(fc_date,init_date=None,leadtime=None):
    if (init_date is None and leadtime is None):
        leadtime = 0
        init_date = fc_date
    elif (init_date is None):
        init_date = fc_date - timedelta(hours = leadtime)
    elif (leadtime is None):
        leadtime = int(np.abs(((fc_date - init_date).total_seconds()))/60/60)
    return fc_date, init_date, leadtime

def get_filevarname(model,varalias,variable_info,model_dict,ncdict):
    stdname = variable_info[varalias]['standard_name']
    print('Get filevarname for \n' +
          'stdvarname:',stdname,'\n' + 
          'varalias:',varalias)
    filevarname = get_varname_for_cf_stdname_in_ncfile(ncdict,stdname)
    print(filevarname)
    if (filevarname is None and 'alias' in variable_info[varalias]):
        filevarname = get_varname_for_cf_stdname_in_ncfile(
                        ncdict,variable_info[varalias]['alias'])
    if (filevarname is not None and len(filevarname) > 1):
        print('!!! standard_name: ',stdname,' is not unique !!!',
            '\nThe following variables have the same standard_name:\n',
            filevarname)
        print('Searching model_specs.yaml config file for definition')
        filevarname = None
    if filevarname is not None:
        return filevarname[0]
    if (filevarname is None 
    and varalias in model_dict[model]['vars'].keys()):
        filevarname = model_dict[model]['vars'][varalias]
        print('Variable defined in model_specs.yaml is:')
        print(varalias, '=', filevarname)
        return filevarname
    elif (filevarname is None
    and varalias not in model_dict[model]['vars'].keys()
    and 'aliases_of_vector_components' in variable_info[varalias]):
        print('Checking variable_info if variable can be ' + 
        'computed from vector components')
        filevarname = variable_info[varalias]['aliases_of_vector_components']
        return filevarname
    else:
        print('!!! variable not defined or '
            + 'available in model output file !!!')

def get_model(model=None,sdate=None,edate=None,
    fc_date=None,init_date=None,leadtime=None,varalias=None):
    """ 
    toplevel function to get model data
    """
    fc_date, init_date, leadtime = make_dates_and_lt(
                                            fc_date=fc_date,
                                            init_date=init_date,
                                            leadtime=leadtime)
    filestr = make_model_filename(model=model,fc_date=fc_date,
                                    leadtime=leadtime)
    vardict, \
    filevarname = get_model_fc_mode(filestr=filestr,model=model,
                    fc_date=fc_date,leadtime=leadtime,varalias=varalias)
    return vardict, fc_date, init_date, leadtime, filestr, filevarname
# ---------------------------------------------------------------------#

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                        '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)


class model_class():

    '''
    class to read and process model data 
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and 
    station classes.
    '''

    def __init__(self,model='mwam4',sdate=None,edate=None,fc_date=None,
                init_date=None,leadtime=None,varalias='Hs'):
        print ('# ----- ')
        print (" ### Initializing model_class object ###")
        if fc_date is not None:
            print ("Requested time: ", str(fc_date))
        elif (edate is None and fc_date is None and sdate is not None):
            fc_date = sdate
            print ("Requested time: ", str(fc_date))
        elif (sdate is None and edate is None and fc_date is None):
            now = datetime.now()
            fc_date = datetime(now.year,now.month,now.day,now.hour)
            print ("Requested time: ", str(fc_date))
        else:
            # time frame function to access a temporal subset
            # --> not yet in use
            print ("Requested time frame: " +
                str(sdate) + " - " + str(edate))

        print('Please wait ...')

        vardict, \
        fc_date, init_date, \
        leadtime, filestr, \
        filevarname = get_model(model=model,sdate=sdate,edate=edate,
                            fc_date=fc_date,init_date=init_date,
                            leadtime=leadtime,varalias=varalias)

        stdname = variable_info[varalias]['standard_name']
        varname = filevarname
        # define class variables
        self.fc_date = fc_date
        self.init_date = init_date
        self.sdate = sdate # not yet in use
        self.edate = edate # not yet in use
        self.model = model
        self.varalias = varalias
        self.varname = varname
        self.stdvarname = stdname
        self.vars = vardict
        self.filestr = filestr
        print (" ### model_class object initialized ###")
        print ('# ----- ')
