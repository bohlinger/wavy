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
    init_times = np.array(model_dict[model]['init_times']).astype('float')
    date = fc_date - timedelta(hours=leadtime)
    if date.hour in init_times:
        init_diffs = date.hour - init_times
        init_diffs[init_diffs<0] = np.nan
        h_idx = np.where(init_diffs==np.min(init_diffs[~np.isnan(init_diffs)]))
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

def get_model_fc_mode(filestr,model,fc_date,leadtime=None,varname=None):
    """ 
    fct to retrieve model data for correct time
    """
    print ("Get model data according to selected date ....")
    print(filestr)
    f = netCDF4.Dataset(filestr,'r')
    model_lons = f.variables[model_dict[model]['coords']['lons']][:]
    model_lats = f.variables[model_dict[model]['coords']['lats']][:]
    model_time = f.variables[model_dict[model]['vars']['time']]
    # Hs [time,lat,lon]
    model_var_link = f.variables[model_dict[model]['vars'][varname]]
    model_time_dt = list(netCDF4.num2date(model_time[:],
                                          units = model_time.units) )
    model_time_dt_valid = [model_time_dt[model_time_dt.index(fc_date)]]
    model_time_valid = [model_time[model_time_dt.index(fc_date)]]
    model_time_unit = model_time.units
    if len(model_var_link.shape)>2: # for mulitple time steps
        model_var_valid = \
            model_var_link[model_time_dt.index(fc_date),:,:].squeeze()
    else:# if only one time step
        model_var_valid = model_var_link[:,:].squeeze()
    f.close()
    return model_var_valid, model_lats, model_lons, model_time_valid,\
         model_time_dt_valid, model_time_unit

def make_dates_and_lt(fc_date,init_date=None,leadtime=None):
    if (init_date is None and leadtime is None):
        leadtime = 0
        init_date = fc_date
    elif (init_date is None):
        init_date = fc_date - timedelta(hours=leadtime)
    elif (leadtime is None):
        leadtime = int(np.abs(((fc_date - init_date).total_seconds()))/60/60)
    return fc_date, init_date, leadtime

def get_model(model=None,sdate=None,edate=None,
    fc_date=None,init_date=None,leadtime=0,varname='Hs'):
    """ 
    toplevel function to get model data
    """
    fc_date, init_date, leadtime = make_dates_and_lt(
                                            fc_date=fc_date,
                                            init_date=init_date,
                                            leadtime=leadtime)
    filename = make_model_filename(model=model,fc_date=fc_date,
                                    leadtime=leadtime)
    model_var, \
    model_lats, \
    model_lons, \
    model_time, \
    model_time_dt, \
    model_time_unit = get_model_fc_mode(filestr=filename,model=model,
                    fc_date=fc_date,leadtime=leadtime,varname=varname)
    model_var_dict = {
                    'model_var':model_var,
                    'model_lats':model_lats,
                    'model_lons':model_lons,
                    'model_time':model_time,
                    'model_time_dt':model_time_dt,
                    'model_time_unit':model_time_unit,
                    }
    return model_var_dict, fc_date, init_date, leadtime
# ---------------------------------------------------------------------#

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

class model_class():
    '''
    class to read and process model data 
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and 
    station classes.
    '''

    def __init__(self,model=None,sdate=None,edate=None,fc_date=None,
                init_date=None,leadtime=0,varname='Hs'):
        print ('# ----- ')
        print (" ### Initializing modelmod instance ###")
        print ('# ----- ')
        if fc_date is not None:
            print ("Requested time: ", str(fc_date))
        elif (edate is None and fc_date is None):
            fc_date=sdate
            print ("Requested time: ", str(fc_date))
        else:
            # time frame function to access a temporal subset
            # --> not yet in use
            print ("Requested time frame: " +
                str(sdate) + " - " + str(edate))

        model_var_dict, \
        fc_date, init_date, \
        leadtime = get_model(model=model,sdate=sdate,edate=edate,
                            fc_date=fc_date,init_date=init_date,
                            leadtime=leadtime,varname=varname)

        self.fc_date = fc_date
        self.init_date = init_date
        self.sdate = sdate
        self.edate = edate
        self.model = model
        self.model_var_dict = model_var_dict
