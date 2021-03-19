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

# own imports
from utils import haversine, haversine_new, collocate_times
from utils import progress
# ---------------------------------------------------------------------#

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

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

def collocate(mc_obj,obs_obj=None,col_obj=None,collocation_idx=None,
            distlim=None):
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
    if len(mc_obj.vars[mc_obj.stdvarname]) < 1:
        raise Exception ( '\n###\n'
                        + 'Collocation not possible, '
                        + 'no model values available for collocation!'
                        + '\n###'
                        )
    if (len(mc_obj.vars['time'])>1 and len(obs_obj.vars['time'])>1): 
        # time collocation
        idx = collocate_times(  unfiltered_t = obs_obj.vars['datetime'],
                                target_t = mc_obj.vars['datetime'] )
        dist = haversine_new( mc_obj.vars['longitude'][0],
                              mc_obj.vars['latitude'][0],
                              obs_obj.lon,obs_obj.lat )
        results_dict = {
                'time':mc_obj.vars['time'],
                'time_unit':mc_obj.vars['time_unit'],
                'datetime':list(np.array(obs_obj.vars['datetime'])[idx]),
                'distance':dist*len(mc_obj.vars['time']),
                'model_values':mc_obj.vars[mc_obj.stdvarname],
                'model_lons':mc_obj.vars['longitude'],
                'model_lats':mc_obj.vars['latitude'],
                'obs_values':list(np.array(obs_obj.vars[
                                            obs_obj.stdvarname])[idx]),
                'obs_lons':obs_obj.vars['longitude'],
                'obs_lats':obs_obj.vars['latitude'],
                'collocation_idx':None
                }
        
    else: # space collocation
        dtime = netCDF4.num2date(obs_obj.vars['time'],obs_obj.vars['time_unit'])
        datein = netCDF4.num2date(mc_obj.vars['time'],mc_obj.vars['time_unit'])
        if isinstance(dtime,np.ndarray):
            dtime = list(dtime)
        if isinstance(datein,np.ndarray):
            datein = list(datein)
        cidx = collocate_times(dtime,target_t=datein,twin=obs_obj.timewin)
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
        if (collocation_idx is None and col_obj is None):
            print ("No collocation idx available")
            print ("Perform collocation with moving window of degree\n",\
                "lon:",lon_win,"lat:",lat_win)
            for j in range(len(obs_time_dt)):
                progress(j,str(int(len(obs_time_dt))),'')
                try:
#                for i in range(1):
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
                #'valid_date':np.array(datein),
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
            results_dict['model_values'] = model_vals[collocation_idx]
    return results_dict


class collocation_class():
    '''
    draft of envisioned collocation class object
    '''

    def __init__(self,mc_obj,sa_obj=None,st_obj=None,col_obj=None,
    distlim=None):
        print ('# ----- ')
        print (" ### Initializing collocation_class object ###")
        print (" Please wait ...")
        if sa_obj is not None:
            obs_obj = sa_obj
            obs_obj.timewin = sa_obj.timewin
            obsname = sa_obj.sat
        if st_obj is not None:
            obs_obj = st_obj
            obs_obj.timewin = None
            obsname = st_obj.platform + '_' +  st_obj.sensor
        t0=time.time()
        results_dict = collocate(mc_obj,
                                obs_obj=obs_obj,
                                col_obj=col_obj,
                                distlim=distlim)
        t1=time.time()
        print("time used for collocation:",round(t1-t0,2),"seconds")
        # define class variables
        self.fc_date = mc_obj.fc_date
        self.sdate = obs_obj.sdate
        self.edate = obs_obj.edate
        self.model = mc_obj.model
        self.obsname = obsname
        self.varalias = mc_obj.varalias
        self.stdvarname = mc_obj.stdvarname
        #self.vars = vardict # divided into model and obs
        self.vars = results_dict
        print(len(self.vars['time'])," values collocated")
        print (" ### Collocation_class object initialized ###")
        print ('# ----- ')
