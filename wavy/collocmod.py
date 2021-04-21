#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
"""
- Module that should take care of collocation of points or swaths
- Needs input from modules that retrieve from observational platforms
  and models
"""
# --- import libraries ------------------------------------------------#
# standard library imports
import sys
import numpy as np
import yaml
import netCDF4
from datetime import datetime, timedelta
import os
import time
import calendar
from dateutil.relativedelta import relativedelta

# own imports
from utils import haversine, haversine_new, collocate_times
from utils import progress, make_fc_dates
from utils import make_pathtofile
from utils import hour_rounder
from modelmod import model_class, make_model_filename_wrapper
from modelmod import get_model_filedate
from ncmod import dumptonc_ts_collocation
# ---------------------------------------------------------------------#

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/station_specs.yaml'))
with open(moddir,'r') as stream:
    station_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/collocation_specs.yaml'))
with open(moddir,'r') as stream:
    collocation_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)

flatten = lambda l: [item for sublist in l for item in sublist]

def get_collocation_idx(distlst,tmp_idx):
    tmp_idx2 = distlst.index(np.nanmin(distlst))
    collocation_idx = tmp_idx[tmp_idx2]
    return collocation_idx,tmp_idx2

def collocation_loop(j,distlim,obs_lats,obs_lons,model_lats,
model_lons,model_vals,lon_win,lat_win):
    obs_lat = obs_lats[j]
    obs_lon = obs_lons[j]
    # constraints to reduce workload
    model_lats_new = model_lats[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-lon_win)
                    &
                    (model_lons<=obs_lon+lon_win)
                    ]
    model_lons_new = model_lons[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-lon_win)
                    &
                    (model_lons<=obs_lon+lon_win)
                    ]
    tmp = range(len(model_lats))
    tmp_idx = np.array(tmp)[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-lon_win)
                    &
                    (model_lons<=obs_lon+lon_win)
                    ]
    # compute distances
    #distlst = list(map(
    #                haversine,
    #                [obs_lon]*len(model_lons_new),
    #                [obs_lat]*len(model_lons_new),
    #                model_lons_new,model_lats_new
    #                ))
    distlst = haversine_new([obs_lon]*len(model_lons_new),
                            [obs_lat]*len(model_lons_new),
                            model_lons_new,model_lats_new)
    collocation_idx,tmp_idx2 = get_collocation_idx(distlst,tmp_idx)
    dist = distlst[tmp_idx2]
    if dist <= distlim:
        while(model_vals[collocation_idx]<=0):
            distlst[tmp_idx2] = np.nan
            collocation_idx,tmp_idx2 = get_collocation_idx(distlst,tmp_idx)
            dist = distlst[tmp_idx2]
        return collocation_idx,dist
    else:
        return

def find_valid_fc_dates_for_model_and_leadtime(fc_dates,model,leadtime):
    '''
    Finds valid dates that are close to desired dates at a precision
    of complete hours
    '''
    if (leadtime is None or leadtime == 'best'):
        fc_dates_new = [hour_rounder(d) for d in fc_dates]
    else:
        fc_dates_new = [hour_rounder(d) for d in fc_dates \
                    if get_model_filedate(model,d,leadtime) != False]
    return fc_dates_new

def check_if_file_is_valid(fc_date,model,leadtime):
    fname = make_model_filename_wrapper(model,fc_date,leadtime)
    print('Check if requested file:\n',fname,'\nis available and valid')
    try:
        nc = netCDF4.Dataset(fname,mode='r')
        time = nc.variables['time']
        dt = netCDF4.num2date(time[:],time.units)
        if fc_date in list(dt):
            print('File is available')
            return True
        else:
            print('Desired date ' + str(fc_date) +  ' is not in', fname)
            return False
    except FileNotFoundError as e:
        print('File is not available')
        print(e)
        return False

def get_closest_date(overdetermined_lst,target_lst):
    idx = []
    for i in range(len(target_lst)):
        diffs=np.abs( [ ( target_lst[i]
                        - overdetermined_lst[j] ).total_seconds()
                        for j in range(len(overdetermined_lst)) ] )
        mindiff= np.min(diffs)
        idx.append(list(diffs).index(mindiff))
    return idx

def collocate_station_ts(obs_obj=None,model=None,distlim=None,\
    leadtime=None,date_incr=None):
    """
    Some info
    """
    fc_date = make_fc_dates(obs_obj.sdate,obs_obj.edate,date_incr)
    # get coinciding date between fc_date and dates in obs_obj
    idx1 = collocate_times( unfiltered_t = obs_obj.vars['datetime'],
                                target_t = fc_date, twin = obs_obj.twin )
    # find valid/coinciding fc_dates
    if len(idx1) > len(fc_date):
        print('Muliple assignments within given time window')
        print('--> only closest to time stamp is chosen')
        idx_closest = get_closest_date(\
                    list(np.array(obs_obj.vars['datetime'])[idx1]),\
                    fc_date)
        idx1 = list(np.array(idx1)[idx_closest])
    # adjust obs_obj according to valid dates
    for key in obs_obj.vars.keys():
        if (key != 'time_unit' and key !='model_meta'):
            obs_obj.vars[key] = list(np.array(obs_obj.vars[key])[idx1])
    # adjust again assumed fc_dates by filtered obs dates
    fc_date = obs_obj.vars['datetime']
    # find valid dates for given leadtime and model
    fc_date = find_valid_fc_dates_for_model_and_leadtime(\
                                    fc_date,model,leadtime)
    # check if file exists and if it includes desired time
    # if not check next possible file
    check = False
    for d in range(len(fc_date)):
        check = check_if_file_is_valid(fc_date[d],model,leadtime)
        if check == True:
            break
    if check == True:
        mc_obj = model_class( model=model,
                              fc_date=fc_date[d],
                              leadtime=leadtime,
                              varalias=obs_obj.varalias)
        col_obj = collocation_class( mc_obj=mc_obj,
                                     st_obj=obs_obj,
                                     distlim=distlim )
        model_vals = [col_obj.vars['model_values'][0]]
        tmpdate = hour_rounder(col_obj.vars['datetime'][0])
        model_datetime = [ datetime(tmpdate.year,
                                    tmpdate.month,
                                    tmpdate.day,
                                    tmpdate.hour) ]
        model_time = [netCDF4.date2num(model_datetime[0],
                        units=col_obj.vars['time_unit'])]
    if check == False:
        print('No valid model file available!')
    else:
        for i in range(d+1,len(fc_date)):
#            for t in range(1):
            try:
                check = check_if_file_is_valid(fc_date[i],model,leadtime)
                if check == False:
                    raise FileNotFoundError
                mc_obj = model_class( model=model,
                                      fc_date=fc_date[i],
                                      leadtime=leadtime,
                                      varalias=obs_obj.varalias )
                model_vals.append(list(mc_obj.vars[\
                                        mc_obj.stdvarname ].flatten()\
                                    [ col_obj.vars['collocation_idx']\
                                    ] )[0])
                model_time.append(mc_obj.vars['time'][0])
                model_datetime.append( datetime(\
                                    mc_obj.vars['datetime'][0].year,
                                    mc_obj.vars['datetime'][0].month,
                                    mc_obj.vars['datetime'][0].day,
                                    mc_obj.vars['datetime'][0].hour ) )
            except FileNotFoundError as e:
                print(e)
        # potentially there are different number of values 
        # for obs and model
        # double check and use only coherent datetimes
        idx2 = collocate_times( model_datetime,
                                target_t=obs_obj.vars['datetime'],
                                twin = obs_obj.twin)
        col_obj.vars['model_values'] = list(np.array(\
                                                    model_vals)[idx2])
        col_obj.vars['time'] = list(np.array(model_time)\
                                            [idx2])
        col_obj.vars['datetime'] = list(np.array(\
                                            model_datetime)[idx2])
        idx3 = collocate_times(  \
                            unfiltered_t = obs_obj.vars['datetime'],
                            target_t = col_obj.vars['datetime'],
                            twin = obs_obj.twin)
        col_obj.vars['obs_values'] = list(np.array(
                                        obs_obj.vars[
                                            obs_obj.stdvarname
                                                    ])[idx3])
    # valid_date is meaningless for ts application and set to None
    col_obj.vars['valid_date'] = None
    # inflate length of constant sized variables
    col_obj.vars['distance'] = col_obj.vars['distance']*\
                                       len(col_obj.vars['datetime'])
    col_obj.vars['obs_lats'] = col_obj.vars['obs_lats']*\
                                        len(col_obj.vars['datetime'])
    col_obj.vars['obs_lons'] = col_obj.vars['obs_lons']*\
                                        len(col_obj.vars['datetime'])
    col_obj.vars['collocation_idx'] = col_obj.vars['collocation_idx']*\
                                        len(col_obj.vars['datetime'])
    col_obj.vars['model_lats'] = col_obj.vars['model_lats']*\
                                       len(col_obj.vars['datetime'])
    col_obj.vars['model_lons'] = col_obj.vars['model_lons']*\
                                        len(col_obj.vars['datetime'])
    results_dict = col_obj.vars
    return results_dict

def collocate_field(mc_obj=None,obs_obj=None,col_obj=None,distlim=None):
    """
    Some info
    """
    dtime = netCDF4.num2date(obs_obj.vars['time'],obs_obj.vars['time_unit'])
    datein = netCDF4.num2date(mc_obj.vars['time'],mc_obj.vars['time_unit'])
    if isinstance(dtime,np.ndarray):
        dtime = list(dtime)
    if isinstance(datein,np.ndarray):
        datein = list(datein)
    cidx = collocate_times(dtime,target_t=datein,twin=obs_obj.twin)
    obs_time_dt = np.array(dtime)[cidx]
    obs_time = np.array(obs_obj.vars['time'])[cidx]
    obs_time_unit = obs_obj.vars['time_unit']
    # Compare wave heights of satellite with model with 
    # constraint on distance and time frame
    collocation_idx_lst = []
    dist_lst = []
    time_idx_lst = []
    # create local variables before loop
    obs_lats = np.array(obs_obj.vars['latitude'])[cidx]
    obs_lons = np.array(obs_obj.vars['longitude'])[cidx]
    obs_vals = np.array(obs_obj.vars[obs_obj.stdvarname])[cidx]
    # flatten numpy arrays
    model_lats = mc_obj.vars['latitude'].flatten()
    model_lons = mc_obj.vars['longitude'].flatten()
    model_vals = mc_obj.vars[mc_obj.stdvarname].flatten()
    # moving window compensating for increasing latitudes
    if distlim == None:
        distlim = 6
    lon_win = round(
                (distlim /
                 haversine(0,
                    np.max(np.abs(obs_lats)),
                    1,
                    np.max(np.abs(obs_lats))) ) 
                + 0.01, 2)
    lat_win = round(distlim/111.+0.01,2)
    if (col_obj is None):
        print ("No collocation idx available")
        print (len(obs_time_dt),"footprints to be collocated")
        print ("Perform collocation with moving window of degree\n",\
            "lon:",lon_win,"lat:",lat_win)
        for j in range(len(obs_time_dt)):
            progress(j,str(int(len(obs_time_dt))),'')
            try:
#            for i in range(1):
                collocation_idx,dist = collocation_loop(\
                    j,distlim,obs_lats,obs_lons,\
                    model_lats,model_lons,model_vals,\
                    lon_win,lat_win)
                collocation_idx_lst.append(collocation_idx)
                dist_lst.append(dist)
                time_idx_lst.append(j)
            except:
                print ("Collocation error -> no collocation:",
                        sys.exc_info()[0])
        results_dict = {
            'valid_date':np.array(datein),
            'time':list(np.array(obs_time[time_idx_lst])),
            'time_unit':obs_time_unit,
            'datetime':list(np.array(obs_time_dt[time_idx_lst])),
            'distance':dist_lst,
            'model_values':list(model_vals[collocation_idx_lst]),
            'model_lons':list(model_lons[collocation_idx_lst]),
            'model_lats':list(model_lats[collocation_idx_lst]),
            'obs_values':list(obs_vals[time_idx_lst]),
            'obs_lons':list(obs_lons[time_idx_lst]),
            'obs_lats':list(obs_lats[time_idx_lst]),
            'collocation_idx':collocation_idx_lst
            }
    elif (col_obj is not None and len(col_obj.vars['collocation_idx']) > 0):
        print("Collocation idx given through collocation_class object")
        results_dict = col_obj.vars
        results_dict['model_values'] = model_vals[\
                                        col_obj.vars['collocation_idx']]
    return results_dict


def collocate(mc_obj=None,obs_obj=None,col_obj=None,
    model=None,distlim=None,leadtime=None,date_incr=None):
    """
    get obs value for model value for given 
        temporal and spatial constraints
    """
    if (len(obs_obj.vars[obs_obj.stdvarname]) < 1):
        raise Exception ( '\n###\n'
                        + 'Collocation not possible, '
                        + 'no observation values for collocation!'
                        + '\n###'
                        )
    if ((mc_obj is None or len(mc_obj.vars[mc_obj.stdvarname]) < 1) 
    and model is None):
        raise Exception ( '\n###\n'
                        + 'Collocation not possible, '
                        + 'no model values available for collocation!'
                        + '\n###'
                        )
    if (mc_obj is None and model is not None and obs_obj is not None):
        results_dict = collocate_station_ts(obs_obj=obs_obj,
                                            model=model,\
                                            distlim=distlim,\
                                            leadtime=leadtime,\
                                            date_incr=date_incr)
        return results_dict
    else:
        results_dict = collocate_field(mc_obj=mc_obj,\
                                              obs_obj=obs_obj,\
                                              col_obj=col_obj,\
                                              distlim=distlim)
    return results_dict


class collocation_class():
    '''
    draft of envisioned collocation class object
    '''

    def __init__(self,mc_obj=None,sa_obj=None,st_obj=None,col_obj=None,
        model=None,distlim=None,leadtime=None,date_incr=None):
        print ('# ----- ')
        print (" ### Initializing collocation_class object ###")
        print (" Please wait ...")
        if sa_obj is not None:
            obs_obj = sa_obj
            obs_obj.twin = sa_obj.twin
            self.obsname = sa_obj.sat
            self.sat = sa_obj.sat
            self.obstype = "satellite_altimeter"
            self.region = sa_obj.region
        if st_obj is not None:
            obs_obj = st_obj
            if 'twin' in station_dict['platform'][st_obj.platform].keys():
                obs_obj.twin =  station_dict['platform']\
                                [st_obj.platform]['twin']
            else:
                obs_obj.twin = None
            self.obsname = st_obj.platform + '_' +  st_obj.sensor
            self.obstype = 'platform'
            self.platform = st_obj.platform
            self.sensor = st_obj.sensor
        if mc_obj is not None:
            model = mc_obj.model
        # define class variables
        self.sdate = obs_obj.sdate
        self.edate = obs_obj.edate
        self.model = model
        self.varalias = obs_obj.varalias
        self.stdvarname = obs_obj.stdvarname
        self.leadtime = leadtime
        if leadtime is None:
            self.leadtime = 'best'
        # get vars dictionary
        try:
#        for i in range(1):
            t0=time.time()
            results_dict = collocate(mc_obj=mc_obj,
                                    obs_obj=obs_obj,
                                    col_obj=col_obj,
                                    model=model,
                                    distlim=distlim,
                                    leadtime=leadtime,
                                    date_incr=date_incr)
            t1=time.time()
            print("time used for collocation:",round(t1-t0,2),"seconds")
            self.vars = results_dict
            self.fc_date = results_dict['datetime']
            print(len(self.vars['time'])," values collocated")
            if 'superob' in vars(obs_obj).keys():
                self.superob = obs_obj.superob
                self.outlier_detection = obs_obj.outlier_detection
                self.missing_data = obs_obj.missing_data
            print (" ### Collocation_class object initialized ###")
        except Exception as e:
            print(e)
            self.error = e
            print ("! No collocation_class object initialized !")
        print ('# ----- ')

    def write_to_monthly_nc(self,path=None,filename=None):
        # divide time into months by loop over months from sdate to edate
        if 'error' in vars(self):
            print('Erroneous collocation_class file detected')
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
                        path_template = collocation_dict['path']\
                                                    [self.obstype]\
                                                    ['local']['nc']\
                                                    ['path_template'][0]
                    if filename is None:
                        file_template = collocation_dict['path'][self.obstype]\
                                                    ['local']['nc']\
                                                    ['file_template']
                    strsublst = collocation_dict['path'][self.obstype]\
                                                    ['local']['nc']\
                                                    ['strsub']
                    if 'superob' in vars(self).keys():
                        file_template = 'superobbed_' + file_template
                    tmppath = path_template + '/' + file_template
                    if isinstance(self.leadtime,str):
                        leadtimestr=self.leadtime
                    else:
                        leadtimestr="{:0>3d}h".format(self.leadtime)
                    if self.obstype=='platform':
                        pathtofile = make_pathtofile(tmppath,strsublst,
                                            tmpdate,
                                            varalias=self.varalias,
                                            model=self.model,
                                            platform=self.platform,
                                            sensor=self.sensor,
                                            leadtime=leadtimestr)
                        title = ( 'Collocation of ' + self.stdvarname
                                + ' observations from '
                                + self.platform + ' ' + self.sensor
                                + ' vs ' + self.model)
                    elif self.obstype=='satellite_altimeter':
                        pathtofile = make_pathtofile(tmppath,strsublst,
                                                tmpdate,
                                                varalias=self.varalias,
                                                model=self.model,
                                                satellite=self.sat,
                                                region=self.region,
                                                leadtime=leadtimestr)
                        title = ( 'Collocation of ' + self.stdvarname 
                                + ' observations from ' + self.sat
                                + ' vs ' + self.model)
                dumptonc_ts_collocation(self,pathtofile,title)
                tmpdate = tmpdate + relativedelta(months = +1)
        return

    def validate_collocated_values(self,date=None,twin=None,path=None,
    filename=None):
        validation_dict = validate_collocated_valuesa(col_obj=self,date=date,\
                                    twin=twin,path=path,filename=filename)
        return validation_dict

def validate_collocated_values(col_obj=None,mods=None,obs=None,\
dtime=None,date=None,twin=None,path=None,filename=None):
    if col_obj is not None:
        mods = col_obj.vars['model_values']
        obs = col_obj.vars['obs_values']
        dtime = col_obj.vars['datetime']
    # get idx for date and twin
    from utils import find_included_times
    idx = find_included_times(dtime,target_t=date,twin=twin)
    mods = np.array(mods)[idx]
    obs = np.array(obs)[idx]
    results_dict = {'model_values':mods,'obs_values':obs}
    # validate
    from validationmod import validate
    validation_dict = validate(results_dict)
    return validation_dict

