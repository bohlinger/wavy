#!/usr/bin/env python
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

# read files
import netCDF4
from netCDF4 import Dataset

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os
import math

# progress bar
from utils import progress, hour_rounder

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# get necessary paths for module
import pathfinder

# import outsorced specs
from model_specs import model_dict

# matchtime fct
from stationmod import matchtime

# module to dump satellite_class object into nc-file
# should also treat other similar type data

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


class model_class():
    '''
    class to read and process model data 
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and 
    station classes.
    '''
    satpath_lustre = pathfinder.satpath_lustre
    satpath_copernicus = pathfinder.satpath_copernicus
    satpath_ftp_014_001 = pathfinder.satpath_ftp_014_001
    
    from region_specs import region_dict
    from model_specs import model_dict

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
    Erin1W   3 hourly (00h, 03h, 06h, 09h, 12h, 15h, 18h)
    ww3        hourly
    """
    if (model == 'mwam4' or model == 'mwam4force'):
        multsix = int(leadtime/6)
        restsix = leadtime%6
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and 
            (fc_date - timedelta(hours=leadtime)).hour != 6 and 
            (fc_date - timedelta(hours=leadtime)).hour !=12 and 
            (fc_date - timedelta(hours=leadtime)).hour !=18):
            sys.exit('error: --> leadtime is not available') 
        if leadtime>60:
            sys.exit('error: --> Leadtime must be less than 60')
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
            sys.exit('error: --> leadtime is not available')
        if leadtime>60:
            sys.exit('error: --> Leadtime must be less than 60')
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
            sys.exit('error: --> leadtime is not available')
        if leadtime>228:
            sys.exit('error: --> Leadtime must be less than 228')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*24)
                       - timedelta(hours=restsix))
    elif (model == 'ecwam' or model == 'mwam800c3'):
        multsix = int(leadtime/12)
        restsix = leadtime%12
        if ((fc_date - timedelta(hours=leadtime)).hour != 0 and
            (fc_date - timedelta(hours=leadtime)).hour !=12):
            sys.exit('error: --> leadtime is not available')
        if leadtime>60:
            sys.exit('error: --> Leadtime must be less than 60')
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
            sys.exit('error: --> leadtime is not available')
        if leadtime>60:
            sys.exit('error: --> Leadtime must be less than 60')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date 
                       - timedelta(hours=multsix*3)
                       - timedelta(hours=restsix))
    elif (model=='ww3' or model == 'swan_karmoy250'):
        tmp_date = fc_date
    return tmp_date

def make_filename(simmode=None,model=None,datein=None,
    expname=None,fc_date=None,init_date=None,leadtime=None):
    from model_specs import explst
    filetemplate = 'file_template'
    if simmode == 'fc':
        if model == 'ARCMFC':
            filename = (model_dict[model]['path']
              + fc_date.strftime('%Y%m%d')
              + init_date.strftime(model_dict[model][filetemplate]))
        elif (model == 'Erin1W' or model == 'Erin2W'):
            filename = (model_dict[model]['path']
              + fc_date.strftime(model_dict[model][filetemplate]))
        elif (model == 'mwam4' or model=='mwam8' or model=='ecwam' or\
            model=='mwam800c3' or model == 'mwam4force' or \
            model=='mwam8force' or model=='ecwamforce'):
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
    if (model in model_dict.keys() and model is not 'ARCNFCnew'):
        filestr = make_filename(simmode=simmode,model=model,
                        fc_date=fc_date,init_date=init_date,
                        leadtime=leadtime)
        filepathlst = [filestr]
    elif model == 'ARCMFCnew':
        filestr_start = make_filename(simmode=simmode,model=model,
                            datein=sdate,expname=expname)
        filestr_end = make_filename(simmode=simmode,model=model,
                            datein=edate,expname=expname)
        filelst = list(np.sort(os.listdir(model_dict[model]['path'])))
        filepathlst = []
        for element in filelst:
            filepathlst.append(model_dict[model]['path'] + element)
        sidx = filelst.index(filestr_start)
        eidx = filelst.index(filestr_end)
        filepathlst = filepathlst[sidx:eidx+1]
    return filepathlst

def get_model_cont_mode(model,sdate,edate,filestr,expname,
    timewin):
    print ("Read model output file: ")
    print (filestr)
    # read the file
    f = netCDF4.Dataset(filestr,'r')
    model_lons = f.variables[model_dict[model]['lons']][:]
    model_lats = f.variables[model_dict[model]['lats']][:]
    model_time = f.variables[model_dict[model]['time']][:]
    # Hs [time,lat,lon]
    model_Hs = f.variables[model_dict[model]['Hs']][:].squeeze()
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
    if (model == 'MoskNC' or model == 'MoskWC'):
        f = netCDF4.Dataset(model_dict[model]['file_coords'],'r')
        model_lons = f.variables['longitude'][:]
        model_lats = f.variables['latitude'][:]
        f.close()
        filestr = glob.glob(filestr)[0]
        print (filestr)
        f = netCDF4.Dataset(filestr,'r')
    else:
        print(filestr)
        f = netCDF4.Dataset(filestr,'r')
        model_lons = f.variables[model_dict[model]['lons']][:]
        model_lats = f.variables[model_dict[model]['lats']][:]
    model_time = f.variables[model_dict[model]['time']][:]
    # Hs [time,lat,lon]
    if (varname == 'Hs' or varname is None):
        model_Hs = f.variables[model_dict[model]['Hs']][:].squeeze()
    else: 
        model_Hs = f.variables[model_dict[model][varname]][:].squeeze()
    f.close()
    model_basetime = model_dict[model]['basetime']
    if (model == 'Erin1W' or model == 'Erin2W' or model == 'ww3'):
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
    print(model_Hs.shape)
    if len(model_Hs.shape)>2:
        model_Hs_valid = model_Hs[model_time_dt.index(fc_date),:,:].squeeze()
    else:
        model_Hs_valid = model_Hs[:,:].squeeze()
    return model_Hs_valid, model_lats, model_lons, model_time_valid,\
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
