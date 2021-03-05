"""
- Module that should take care of collocation of points or swaths
- Needs input from modules that retrieve from observational platforms
  and models
"""
# standard library import
import sys
import numpy as np
from utils import progress
import yaml
import netCDF4
from datetime import datetime, timedelta
import os
import time
# own imports
from utils import haversine

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)


def matchtime(sdate,edate,dtime,timewin=None):
    '''
    fct to obtain the index of the time step closest to the 
    requested time including the respective time stamp(s). 
    Similarily, indices are chosen for the time and defined region.
    '''
    if timewin is None:
        timewin = 0
    # create list of datetime instances
    ctime=[]
    cidx=[]
    idx=0
    if (edate is None or sdate==edate):
        for element in dtime:
            # choose closest match within window of win[minutes]
            if (element >= sdate-timedelta(minutes=timewin)
            and element <= sdate+timedelta(minutes=timewin)):
                ctime.append(element)
                cidx.append(idx)
            idx=idx+1
    if (edate is not None and edate!=sdate):
        for element in dtime:
            if (element >= sdate-timedelta(minutes=timewin)
            and element < edate+timedelta(minutes=timewin)):
                ctime.append(element)
                cidx.append(idx)
            idx=idx+1
    return ctime, cidx

def collocation_loop(j,distlim,obs_lats,obs_lons,model_lats,
model_lons,model_vals,lon_win,lat_win):
    obs_lat=obs_lats[j]
    obs_lon=obs_lons[j]
    # constraints to reduce workload
    model_lats_new=model_lats[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-lon_win)
                    &
                    (model_lons<=obs_lon+lon_win)
                    ]
    model_lons_new=model_lons[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-lon_win)
                    &
                    (model_lons<=obs_lon+lon_win)
                    ]
    tmp=range(len(model_lats))
    tmp_idx=np.array(tmp)[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-lon_win)
                    &
                    (model_lons<=obs_lon+lon_win)
                    ]
    # compute distances
    distlst=list(map(
                    haversine,
                    [obs_lon]*len(model_lons_new),
                    [obs_lat]*len(model_lons_new),
                    model_lons_new,model_lats_new
                    ))
    tmp_idx2 = distlst.index(np.min(distlst))
    collocation_idx = tmp_idx[tmp_idx2]
    dist = distlst[tmp_idx2]
    if (distlst[tmp_idx2]<=distlim and model_vals[collocation_idx]>=0):
        return collocation_idx,distlst[tmp_idx2]
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
    dtime = netCDF4.num2date(obs_obj.vars['time'],obs_obj.vars['time_unit'])
    datein = netCDF4.num2date(mc_obj.vars['time'],mc_obj.vars['time_unit'])
    ctime, cidx = matchtime(datein,datein,dtime,timewin=obs_obj.timewin)
    obs_time_dt = dtime[cidx]
    obs_time = np.array(obs_obj.vars['time'])[cidx]
    obs_time_unit = obs_obj.vars['time_unit']

    # Compare wave heights of satellite with model with 
    # constraint on distance and time frame
    collocation_idx_lst=[]
    dist_lst=[]
    time_idx_lst=[]
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
                    np.max(np.abs(obs_lats)))
                ) + 0.01,
                2)
    lat_win = round(distlim/111.+0.01,2)
    if (collocation_idx is None and col_obj is None):
        print ("No collocation idx available")
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
                'time':np.array(obs_time[time_idx_lst]),
                'time_unit':obs_time_unit,
                'datetime':np.array(obs_time_dt[time_idx_lst]),
                'distance':np.array(dist_lst),
                'model_values':model_vals[collocation_idx_lst],
                'model_lons':model_lons[collocation_idx_lst],
                'model_lats':model_lats[collocation_idx_lst],
                'obs_values':obs_vals[time_idx_lst],
                'obs_lons':obs_lons[time_idx_lst],
                'obs_lats':obs_lats[time_idx_lst],
                'collocation_idx':np.array(collocation_idx_lst)
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
            obsname = st_obj.statname + '_' +  st_obj.sensorname
        t0=time.time()
        results_dict = collocate(mc_obj,
                                obs_obj=obs_obj,
                                col_obj=col_obj,
                                distlim=distlim)
        t1=time.time()
        print("time used for collocation:",round(t1-t0,2),"seconds")
        # define class variables
        self.fc_date = mc_obj.fc_date
        self.init_date = mc_obj.init_date
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
