#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field from model output. I try to mostly follow the PEP convention for 
python code style. Constructive comments on style and effecient 
programming are most welcome!

Currently only mode for continuous model runs, forecast mode is missing
'''
__version__ = "0.5.0"
__author__="Patrik Bohlinger, Norwegian Meteorological Institute"
__maintainer__ = "Patrik Bohlinger"
__email__ = "patrikb@met.no"
__status__ = "under development with operation ARCMFC branch"

# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# ignore irrelevant warnings from matplotlib for stdout
#import warnings
#warnings.filterwarnings("ignore")

# plotting
import matplotlib.pyplot as plt

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
import sys

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

# collocate
from utils import haversine

# module to dump sentinel_class object into nc-file
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
def progress(count, total, status=''):
    "from: https://gist.github.com/vladignatyev/06860ec2040cb497f0f3"
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()
# ---------------------------------------------------------------------#


class model_class():
    '''
    class to read and process model data 
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the sentinel, model, and 
    station classes.
    '''
    satpath_lustre = pathfinder.satpath_lustre
    satpath_copernicus = pathfinder.satpath_copernicus
    satpath_ftp_008_052 = pathfinder.satpath_ftp_008_052
    satpath_ftp_014_001 = pathfinder.satpath_ftp_014_001
    
    from region_specs import region_dict
    from model_specs import model_dict

    def __init__(self,sdate,edate=None,model=None,timewin=None,region=None):
        print ('# ----- ')
        print (" ### Initializing custom_nc instance ###")
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
    """
    if model == 'mwam4':
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
        if leadtime>60:
            sys.exit('error: --> Leadtime must be less than 60')
        if leadtime is None:
            pass
        else:
            tmp_date = (fc_date
                       - timedelta(hours=multsix*24)
                       - timedelta(hours=restsix))
    elif model == 'ecwam':
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
    return tmp_date

def make_filename(simmode=None,model=None,datein=None,
    expname=None,fc_date=None,init_date=None,leadtime=None):
    from model_specs import explst
    if simmode == 'fc':
        if model == 'ARCMFC':
            filename = (model_dict[model]['path']
              + fc_date.strftime('%Y%m%d')
              + init_date.strftime(model_dict[model]['file_template']))
        elif (model == 'mwam4' or model=='mwam8' or model=='ecwam'):
            if (fc_date == init_date or leadtime == 0):
                if ((fc_date.hour != 6 and fc_date.hour != 18) 
                    and model == 'mwam8'):
                    sys.exit('error: --> Check your date. mwam8 is only' 
                            + ' initiated at 06 and 18')
                elif (fc_date.hour%6 != 0 and model == 'mwam4'):
                    sys.exit('error: --> Check your date. mwam4 is only' 
                            + ' initiated at 06, 12, 18, and 00')
                elif (fc_date.hour%12 != 0 and model == 'ecwam'):
                    sys.exit('error: --> Check your date. ecwam is only' 
                            + ' initiated at 00, 12')
                else:
                    filename = (fc_date.strftime(
                            model_dict[model]['path_template'])
                            + fc_date.strftime(
                            model_dict[model]['file_template'])
                            )
            else:
                filedate = check_date(model,
                                fc_date=fc_date,leadtime=leadtime)
                filename = (filedate.strftime(
                            model_dict[model]['path_template'])
                            + filedate.strftime(
                            model_dict[model]['file_template'])
                            )
        elif (model == 'MoskNC' or model == 'MoskWC'):
            filename = (model_dict[model]['path'] 
                    + fc_date.strftime(model_dict[model]['file_template']))
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
    if (model == 'ARCMFC' or model == 'MoskNC' or model == 'MoskWC' or \
        model == 'mwam4' or model=='mwam8' or model=='ecwam'):
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
    init_date=None,leadtime=None):
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
    model_Hs = f.variables[model_dict[model]['Hs']][:].squeeze()
    f.close()
    model_basetime = model_dict[model]['basetime']
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
    #return model_time_dt, model_hs_valid, model_lons, model_lats
    return model_Hs_valid, model_lats, model_lons, model_time_valid,\
         model_time_dt_valid


def tmploop_get_model(
    j,sat_time_dt,model_time_dt_valid,distlim,model,
    sat_rlats,sat_rlons,sat_rHs,
    model_rlats,model_rlons,model_rHs,
    moving_win
    ):
    lat_win = 0.1
    if model in model_dict:
        sat_rlat=sat_rlats[j]
        sat_rlon=sat_rlons[j]
        # constraints to reduce workload
        if model == 'ecwam':
            print"ECWAM not yet implemented, caution ecwam has other dimensions!"
            pass
        else:
            model_rlats_new=model_rlats[
                    (model_rlats>=sat_rlat-lat_win)
                    &
                    (model_rlats<=sat_rlat+lat_win)
                    &
                    (model_rlons>=sat_rlon-moving_win)
                    &
                    (model_rlons<=sat_rlon+moving_win)
                    ]
            model_rlons_new=model_rlons[
                    (model_rlats>=sat_rlat-lat_win)
                    &
                    (model_rlats<=sat_rlat+lat_win)
                    &
                    (model_rlons>=sat_rlon-moving_win)
                    &
                    (model_rlons<=sat_rlon+moving_win)
                    ]
            tmp=range(len(model_rlats))
            tmp_idx=np.array(tmp)[
                    (model_rlats>=sat_rlat-lat_win)
                    &
                    (model_rlats<=sat_rlat+lat_win)
                    &
                    (model_rlons>=sat_rlon-moving_win)
                    &
                    (model_rlons<=sat_rlon+moving_win)
                    ]
            # compute distances
            distlst=map(
                haversine,
                [sat_rlon]*len(model_rlons_new),
                [sat_rlat]*len(model_rlons_new),
                model_rlons_new,model_rlats_new
                )
            tmp_idx2 = distlst.index(np.min(distlst))
            idx_valid = tmp_idx[tmp_idx2]
        if (distlst[tmp_idx2]<=distlim and model_rHs[idx_valid]>=0):
            nearest_all_dist_matches=distlst[tmp_idx2]
            nearest_all_date_matches=sat_time_dt[j]
            nearest_all_model_Hs_matches=\
                           model_rHs[idx_valid]
            nearest_all_sat_Hs_matches=sat_rHs[j]
            nearest_all_sat_lons_matches=sat_rlon
            nearest_all_sat_lats_matches=sat_rlat
            nearest_all_model_lons_matches=\
                            model_rlons[idx_valid]
            nearest_all_model_lats_matches=\
                            model_rlats[idx_valid]
            return nearest_all_date_matches,nearest_all_dist_matches,\
                nearest_all_model_Hs_matches,nearest_all_sat_Hs_matches,\
                nearest_all_sat_lons_matches, nearest_all_sat_lats_matches,\
                nearest_all_model_lons_matches, nearest_all_model_lats_matches
        else:
           return

def collocate(model,model_Hs,model_lats,model_lons,model_time_dt,\
    sa_obj,datein,distlim=None):
    """
    get stellite time steps close to model time step. 
    """
    if distlim is None:
        distlim = int(6)
    timewin = sa_obj.timewin
    ctime, cidx = matchtime(datein,datein,sa_obj.rTIME,
                    sa_obj.basetime,timewin=sa_obj.timewin)
    sat_time_dt=np.array(sa_obj.rtime)[cidx]
    model_time_idx = model_time_dt.index(datein)
    model_time_dt_valid=[model_time_dt[model_time_idx]]
    print ("date matches found:")
    print (model_time_dt_valid)
    # Constrain to region
    if ((model in model_dict)
    and (sa_obj.region=='ecwam' 
    and (sa_obj.region=='Global' or sa_obj.region=='ecwam'))):
        model_rlats = model_lats
        model_rlons = model_lons
        model_rHs = model_Hs.squeeze()
    elif ((model in model_dict)
    and (sa_obj.region!='Arctic' and sa_obj.region!='ARCMFC')):
        model_rlats = model_lats[
                    (model_lats
                    >= sa_obj.region_dict[sa_obj.region]["llcrnrlat"])
                  & (model_lats
                    <= sa_obj.region_dict[sa_obj.region]["urcrnrlat"])
                  & (model_lons
                    >= sa_obj.region_dict[sa_obj.region]["llcrnrlon"])
                  & (model_lons
                    <= sa_obj.region_dict[sa_obj.region]["urcrnrlon"])
                            ]
        model_rlons = model_lons[
                    (model_lats
                    >= sa_obj.region_dict[sa_obj.region]["llcrnrlat"])
                  & (model_lats
                    <= sa_obj.region_dict[sa_obj.region]["urcrnrlat"])
                  & (model_lons
                    >= sa_obj.region_dict[sa_obj.region]["llcrnrlon"])
                  & (model_lons
                    <= sa_obj.region_dict[sa_obj.region]["urcrnrlon"])
                            ]
        tmpA=model_Hs.squeeze()
        tmpB = tmpA[
                (model_lats
                >= sa_obj.region_dict[sa_obj.region]["llcrnrlat"])
                & (model_lats
                <= sa_obj.region_dict[sa_obj.region]["urcrnrlat"])
              & (model_lons
                >= sa_obj.region_dict[sa_obj.region]["llcrnrlon"])
              & (model_lons
                <= sa_obj.region_dict[sa_obj.region]["urcrnrlon"])                            ]
        model_rHs=tmpB
        del tmpA, tmpB
    elif ((model in model_dict)
    and (sa_obj.region=='Arctic' or sa_obj.region=='ARCMFC')):
        model_rlats = model_lats[
                    (model_lats
                    >= sa_obj.region_dict[sa_obj.region]["boundinglat"])
                            ]
        model_rlons = model_lons[
                    (model_lats
                    >= sa_obj.region_dict[sa_obj.region]["boundinglat"])
                            ]
        model_rHs=[]
        tmpA = model_Hs[model_time_idx,:]
        tmpB = tmpA[
                (model_lats
                >= sa_obj.region_dict[sa_obj.region]["boundinglat"])
                        ]
        model_rHs = tmpB
        del tmpA, tmpB
    # Compare wave heights of satellite with model with 
    # constraint on distance and time frame
    nearest_all_date_matches=[]
    nearest_all_dist_matches=[]
    nearest_all_model_Hs_matches=[]
    nearest_all_sat_Hs_matches=[]
    nearest_all_sat_lons_matches=[]
    nearest_all_sat_lats_matches=[]
    nearest_all_model_lons_matches=[]
    nearest_all_model_lats_matches=[]
    # create local variables before loop
    sat_rlats=sa_obj.rloc[0][cidx]
    sat_rlons=sa_obj.rloc[1][cidx]
    sat_rHs=np.array(sa_obj.rHs)[cidx]
    # moving window compensating for increasing latitudes
    try:
        moving_win = round(
                (distlim /
                 haversine(0,
                    np.max(np.abs(sat_rlats)),
                    1,
                    np.max(np.abs(sat_rlats)))
                ),
                2)
    except (ValueError):
        moving_win = .6
    print ("Searching for matches with moving window of degree:",\
            moving_win)
    for j in range(len(sat_time_dt)):
        progress(j,str(int(len(sat_time_dt))),'')
#        for k in range(1):
        try:
            resultlst = tmploop_get_model(\
                j,sat_time_dt,model_time_dt_valid,distlim,model,\
                sat_rlats,sat_rlons,sat_rHs,\
                model_rlats,model_rlons,model_rHs,moving_win)
            nearest_all_date_matches.append(resultlst[0])
            nearest_all_dist_matches.append(resultlst[1])
            nearest_all_model_Hs_matches.append(resultlst[2])
            nearest_all_sat_Hs_matches.append(resultlst[3])
            nearest_all_sat_lons_matches.append(resultlst[4])
            nearest_all_sat_lats_matches.append(resultlst[5])
            nearest_all_model_lons_matches.append(resultlst[6])
            nearest_all_model_lats_matches.append(resultlst[7])
        except:
#            print "Unexpected error:", sys.exc_info()[0]
            pass
    results_dict = {
        'valid_date':np.array(model_time_dt_valid),
        'date_matches':np.array(nearest_all_date_matches),
        'dist_matches':np.array(nearest_all_dist_matches),
        'model_Hs_matches':np.array(nearest_all_model_Hs_matches),
        'sat_Hs_matches':np.array(nearest_all_sat_Hs_matches),
        'sat_lons_matches':np.array(nearest_all_sat_lons_matches),
        'sat_lats_matches':np.array(nearest_all_sat_lats_matches),
        'model_lons_matches':np.array(nearest_all_model_lons_matches),
        'model_lats_matches':np.array(nearest_all_model_lats_matches)
        }
    return results_dict

def get_model(simmode=None,model=None,sdate=None,edate=None,
    fc_date=None,init_date=None,leadtime=None,expname=None,
    sa_obj=None,timewin=None):
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
                    leadtime=leadtime)
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
