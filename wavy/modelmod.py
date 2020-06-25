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
# currently None
# ---------------------------------------------------------------------#

# read yaml config files:
with open("/home/patrikb/wavy/wavy/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)
with open("/home/patrikb/wavy/wavy/pathfinder.yaml", 'r') as stream:
    pathfinder=yaml.safe_load(stream)

class model_class():
    '''
    class to read and process model data 
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and 
    station classes.
    '''
    satpath_lustre = pathfinder['satpath_lustre']
    satpath_copernicus = pathfinder['satpath_copernicus']
    satpath_ftp_014_001 = pathfinder['satpath_ftp_014_001']
    

    def __init__(self,sdate,edate=None,model=None,timewin=None,region=None):
        print ('# ----- ')
        print (" ### Initializing modelmod instance ###")
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

def check_date(model,fc_date=None,init_date=None,leadtime=None):
    """
    mwam4    6 hourly (00h, 06h, 12h, 18h)
    mwam8   12 hourly (06h, 18h)
    ecwam   12 hourly (00h, 12h)
    ARCMFC  24 hourly (00h)
    ARCMFC3 12 hourly (00h, 12h); naming convention +6h for bulletin
    Erin1W   3 hourly (00h, 03h, 06h, 09h, 12h, 15h, 18h)
    ww3        hourly
    """
    if (model == 'mwam4' or model == 'mwam4force' or model == 'ww3'):
        multsix = int(leadtime/6)
        restsix = leadtime%6
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and 
            (fc_date - timedelta(hours=leadtime)).hour != 6 and 
            (fc_date - timedelta(hours=leadtime)).hour !=12 and 
            (fc_date - timedelta(hours=leadtime)).hour !=18):
            print('error: --> leadtime is not available') 
        if leadtime>60:
            print('error: --> Leadtime must be less than 60')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date 
                       - timedelta(hours=multsix*6) 
                       - timedelta(hours=restsix))
    elif model == 'mwam8':
        multsix = int(leadtime/12)
        restsix = leadtime%12
        dummy_date = fc_date + timedelta(hours=6)
        if ((dummy_date - timedelta(hours=leadtime)).hour != 0 and 
            (dummy_date - timedelta(hours=leadtime)).hour !=12):
            print('error: --> leadtime is not available')
        if leadtime>228:
            print('error: --> Leadtime must be less than 228')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*12)
                       - timedelta(hours=restsix))
    elif model == 'ARCMFC':
        multsix = int(leadtime/24)
        restsix = leadtime%24
        if ((fc_date - timedelta(hours=leadtime)).hour != 0):
            print('error: --> leadtime is not available')
        if leadtime>228:
            print('error: --> Leadtime must be less than 228')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*24)
                       - timedelta(hours=restsix))
    elif (model == 'ARCMFC3'):
        multsix = int(leadtime/12)
        restsix = leadtime%12
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and
            (fc_date - timedelta(hours=leadtime)).hour !=12):
            print('error: --> leadtime is not available')
        if leadtime>228:
            print('error: --> Leadtime must be less than 228')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*12)
                       - timedelta(hours=restsix))
    elif (model == 'mwam3force' or model == 'mwam3'):
        multsix = int(leadtime/12)
        restsix = leadtime%12
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and
            (fc_date - timedelta(hours=leadtime)).hour !=12):
            print('error: --> leadtime is not available')
        if leadtime>228:
            print('error: --> Leadtime must be less than 228')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*12)
                       - timedelta(hours=restsix)
                       + timedelta(hours=6))
    elif (model == 'ecwam' or model == 'mwam800c3'):
        multsix = int(leadtime/12)
        restsix = leadtime%12
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and
            (fc_date - timedelta(hours=leadtime)).hour !=12):
            print('error: --> leadtime is not available')
        if leadtime>228:
            print('error: --> Leadtime must be less than 228')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*12)
                       - timedelta(hours=restsix))
    elif (model == 'Erin1W' or model == 'Erin2W'):
        multsix = int(leadtime/3)
        restsix = leadtime%3
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and
            (fc_date - timedelta(hours=leadtime)).hour != 3 and
            (fc_date - timedelta(hours=leadtime)).hour != 6 and
            (fc_date - timedelta(hours=leadtime)).hour != 9 and
            (fc_date - timedelta(hours=leadtime)).hour !=12 and
            (fc_date - timedelta(hours=leadtime)).hour !=15 and
            (fc_date - timedelta(hours=leadtime)).hour !=18 and
            (fc_date - timedelta(hours=leadtime)).hour !=21):
            print('error: --> leadtime is not available')
        if leadtime>60:
            print('error: --> Leadtime must be less than 60')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date 
                       - timedelta(hours=multsix*3)
                       - timedelta(hours=restsix))
    elif (model=='ww3' or model == 'swan_karmoy250'):
        tmp_date = fc_date
    return tmp_date

def get_latest_output_init_date(model):
    '''
    get init_date for latest model output file
    '''
    now = datetime.now()
    init_times = np.array(model_dict[model]['init_times']).astype('float')
    init_diffs = now.hour - init_times
    init_diffs[init_diffs<0] = np.nan
    h_idx = np.where(init_diffs==np.min(init_diffs[~np.isnan(init_diffs)]))
    h = int(init_times[h_idx[0][0]])
    return datetime(now.year,now.month,now.day,h)

def make_filename(simmode=None,model=None,datein=None,
    expname=None,fc_date=None,init_date=None,leadtime=None):
    filetemplate = 'file_template'
    if simmode == 'fc':
        if model == 'ARCMFC':
            filename = (model_dict[model]['path']
              + fc_date.strftime('%Y%m%d')
              + init_date.strftime(model_dict[model][filetemplate]))
        elif model == 'ARCMFC3':
            init_date = init_date + timedelta(hours=6)
            filename = (fc_date.strftime(model_dict[model]['path'])
              + fc_date.strftime('%Y%m%d%H')
              + init_date.strftime(model_dict[model][filetemplate]))
        elif (model == 'Erin1W' or model == 'Erin2W'):
            filename = (model_dict[model]['path']
              + fc_date.strftime(model_dict[model][filetemplate]))
        elif (model == 'ErinFix'):
            filename = (model_dict[model]['path'] 
                        + model_dict[model][filetemplate])
        elif (model == 'mwam4' or model=='mwam8' or model=='ecwam' or\
            model=='mwam800c3' or model == 'mwam4force' or \
            model=='mwam8force' or model=='ecwamforce' or \
            model == 'ww3' or model=='mwam3' or model=='mwam3force'):
            if (fc_date == init_date or leadtime == 0):
                filename = (fc_date.strftime(
                            model_dict[model]['path_template'])
                            + fc_date.strftime(
                            model_dict[model][filetemplate])
                            )
            else:
                filedate = check_date(model,
                                fc_date=fc_date,leadtime=leadtime)
                filename = (filedate.strftime(
                            model_dict[model]['path_template'])
                            + filedate.strftime(
                            model_dict[model][filetemplate])
                            )
        elif (model == 'MoskNC' or model == 'MoskWC' or \
            model=='swanKC' or model=='ww3' or model=='swan_karmoy250'):
            filename = (model_dict[model]['path'] 
                    + fc_date.strftime(model_dict[model][filetemplate]))
    elif simmode == 'cont':
        """
        explst was in model_specs.py for continuous simulations 
        is now removed, this part of the code will be removed as well
        """
        explst = []
        if expname in explst:
            days = [1,10,20]
            tmp = np.abs(np.array(days)-datein.day)
            idx = list(tmp).index(np.min(tmp))
            date = datetime(datein.year,datein.month,days[idx])
            filepath = (model_dict[model]['path']
                + expname
                + date.strftime(model_dict[model]['file_template']))
            filename = (expname + 
                date.strftime(model_dict[model]['file_template']))
    return filename

def get_model_filepathlst(simmode=None,model=None,sdate=None,edate=None,
    expname=None,fc_date=None,init_date=None,leadtime=None):
    if model in model_dict:
        filestr = make_filename(simmode=simmode,model=model,
                        fc_date=fc_date,init_date=init_date,
                        leadtime=leadtime)
        filepathlst = [filestr]
    else:
        print("chosen model is not specified in model_specs.yaml")
    return filepathlst

def get_model_cont_mode(model,sdate,edate,filestr,expname,
    timewin):
    print ("Read model output file: ")
    print (filestr)
    # read the file
    f = netCDF4.Dataset(filestr,'r')
    model_lons = f.variables[model_dict[model]['vars']['lons']][:]
    model_lats = f.variables[model_dict[model]['vars']['lats']][:]
    model_time = f.variables[model_dict[model]['vars']['time']][:]
    # Hs [time,lat,lon]
    model_Hs = f.variables[model_dict[model]['vars']['Hs']][:].squeeze()
    f.close()
    # create datetime objects
    model_basetime = model_dict[model]['basetime']
    model_time_dt=[]
    for element in model_time:
        model_time_dt.append(model_basetime
                        + timedelta(seconds=element))
    # adjust to sdate and edate
    ctime,cidx = matchtime(sdate,edate,model_time,model_basetime,timewin)
    model_Hs = model_Hs[cidx,:,:]
    model_time = model_time[cidx]
    model_time_dt = np.array(model_time_dt)[cidx]
    return model_Hs, model_lats, model_lons, model_time, model_time_dt

def get_model_fc_mode(filestr=None,model=None,fc_date=None,
    init_date=None,leadtime=None,varname=None):
    """ 
    fct to get model data
    if model ARCMFC you need fc_date, init_date
    if model mwam4 you need fc_date, leadtime
    """
    from utils import haversine
    import glob
    print ("Get model data according to selected date ....")
    print(filestr)
    f = netCDF4.Dataset(filestr,'r')
    model_lons = f.variables[model_dict[model]['vars']['lons']][:]
    model_lats = f.variables[model_dict[model]['vars']['lats']][:]
    model_time = f.variables[model_dict[model]['vars']['time']][:]
    # Hs [time,lat,lon]
    if (varname == 'Hs' or varname is None):
        model_var_link = f.variables[model_dict[model]['vars']['Hs']]
    else: 
        model_var_link = f.variables[model_dict[model]['vars'][varname]]
    model_basetime = model_dict[model]['basetime']
    if (model == 'ww3' or model == 'ErinFix'):
        model_time_dt=[]
        for element in model_time:
            # hour_rounder used because ww3 deviates slightly from hours
            model_time_dt.append(hour_rounder(model_basetime
                    + timedelta(days=element)))
    elif model == 'swanKC':
        model_time_dt=[]
        for element in model_time:
            minutes = element/60.
            model_time_dt.append(model_basetime
                    + timedelta(minutes=minutes))
    else:
        model_time_dt=[]
        for element in model_time:
            model_time_dt.append(model_basetime
                    + timedelta(seconds=element))
    model_time_dt_valid = [model_time_dt[model_time_dt.index(fc_date)]]
    model_time_valid = [model_time[model_time_dt.index(fc_date)]]
    print(model_var_link.shape)
    if len(model_var_link.shape)>2:
        model_var_valid = \
            model_var_link[model_time_dt.index(fc_date),:,:].squeeze()
    else:
        model_var_valid = model_var_link[:,:].squeeze()
    f.close()
    return model_var_valid, model_lats, model_lons, model_time_valid,\
         model_time_dt_valid

def get_model(simmode=None,model=None,sdate=None,edate=None,
    fc_date=None,init_date=None,leadtime=None,expname=None,
    sa_obj=None,timewin=None,varname=None):
    """ 
    Get model data.
    """
    if sa_obj is not None:
        sdate = sa_obj.sdate
        edate = sa_obj.edate
    if timewin is None:
        timewin = int(30)
    if (simmode == 'cont'):
        filepathlst = get_model_filepathlst(simmode=simmode,
                        model=model,sdate=sdate,edate=edate,
                        expname=expname)
        model_Hs_lst, \
        model_time_lst, \
        model_time_dt_lst = [],[],[]
        for element in filepathlst:
            model_Hs, \
            model_lats, \
            model_lons, \
            model_time, \
            model_time_dt = \
            get_model_cont_mode(model,\
                                sdate,edate,\
                                element,expname,\
                                timewin)
            for i in range(len(model_time)):
                model_Hs_lst.append(model_Hs[i,:,:])
                model_time_lst.append(model_time[i])
                model_time_dt_lst.append(model_time_dt[i])
        model_Hs = np.array(model_Hs_lst)
    elif (simmode == 'fc'):
        filepathlst = get_model_filepathlst(simmode=simmode,model=model,
                            fc_date=fc_date,init_date=init_date,
                            leadtime=leadtime)
        for element in filepathlst:
            model_Hs, \
            model_lats, \
            model_lons, \
            model_time, \
            model_time_dt = \
            get_model_fc_mode(filestr=element,model=model,
                    fc_date=fc_date,init_date=init_date,
                    leadtime=leadtime,varname=varname)
            model_Hs_lst, \
            model_time_lst, \
            model_time_dt_lst = [],[],[]
            for i in range(len(model_time)):
                model_Hs_lst.append(model_Hs)
                model_time_lst.append(model_time[i])
                model_time_dt_lst.append(model_time_dt[i])
        model_Hs = np.array(model_Hs_lst)
    return model_Hs, model_lats, model_lons, model_time_lst, \
           model_time_dt_lst
