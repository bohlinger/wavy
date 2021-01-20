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

def collocation_loop(
    j,obs_time_dt,model_time_dt_valid,distlim,model,
    obs_lats,obs_lons,obs_val,model_lats,model_lons,model_val,
    moving_win):
    from utils import haversine
    lat_win = 0.1
    if model in model_dict:
        obs_lat=obs_lats[j]
        obs_lon=obs_lons[j]
        # constraints to reduce workload
        model_lats_new=model_lats[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-moving_win)
                    &
                    (model_lons<=obs_lon+moving_win)
                    ]
        model_lons_new=model_lons[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-moving_win)
                    &
                    (model_lons<=obs_lon+moving_win)
                    ]
        tmp=range(len(model_lats))
        tmp_idx=np.array(tmp)[
                    (model_lats>=obs_lat-lat_win)
                    &
                    (model_lats<=obs_lat+lat_win)
                    &
                    (model_lons>=obs_lon-moving_win)
                    &
                    (model_lons<=obs_lon+moving_win)
                    ]
        # compute distances
        if sys.version_info <= (3, 0):
            distlst=map(
                        haversine,
                        [obs_lon]*len(model_lons_new),
                        [obs_lat]*len(model_lons_new),
                        model_lons_new,model_lats_new
                        )
        else:
            distlst=list(map(
                        haversine,
                        [obs_lon]*len(model_lons_new),
                        [obs_lat]*len(model_lons_new),
                        model_lons_new,model_lats_new
                        ))
        tmp_idx2 = distlst.index(np.min(distlst))
        collocation_idx = tmp_idx[tmp_idx2]
        if (distlst[tmp_idx2]<=distlim and model_val[collocation_idx]>=0):
            nearest_all_dist_matches=distlst[tmp_idx2]
            nearest_all_date_matches=obs_time_dt[j]
            nearest_all_model_matches=\
                           model_val[collocation_idx]
            nearest_all_obs_matches=obs_val[j]
            nearest_all_obs_lons_matches=obs_lon
            nearest_all_obs_lats_matches=obs_lat
            nearest_all_model_lons_matches=\
                            model_lons[collocation_idx]
            nearest_all_model_lats_matches=\
                            model_lats[collocation_idx]
            return nearest_all_date_matches,nearest_all_dist_matches,\
                nearest_all_model_matches,nearest_all_obs_matches,\
                nearest_all_obs_lons_matches, nearest_all_obs_lats_matches,\
                nearest_all_model_lons_matches, \
                nearest_all_model_lats_matches, collocation_idx
        else:
           return

def collocate(mc_obj,sa_obj=None,st_obj=None,collocation_idx=None,
            distlim=None):
    """
    get obs value for model value for given 
        temporal and spatial constraints
    """
    if len(sa_obj.vars[sa_obj.stdvarname]) < 1:
        raise Exception ( '\n###\n'
                        + 'Collocation not possible, '
                        + 'no satellite values for collocation!'
                        + '\n###'
                        )
    if len(mc_obj.vars[mc_obj.stdvarname]) < 1:
        raise Exception ( '\n###\n'
                        + 'Collocation not possible, '
                        + 'no model values available for collocation!'
                        + '\n###'
                        )
    dtime = netCDF4.num2date(sa_obj.vars['time'],sa_obj.vars['time_unit'])
    datein = netCDF4.num2date(mc_obj.vars['time'],mc_obj.vars['time_unit'])
    ctime, cidx = matchtime(datein,datein,dtime,timewin=sa_obj.timewin)
    obs_time_dt = dtime[cidx]

    #model_time_idx = model_time_dt.index(datein)
    #model_time_dt_valid = [model_time_dt[model_time_idx]]
    #print ("date matches found:")
    #print (model_time_dt_valid)

    # Compare wave heights of satellite with model with 
    # constraint on distance and time frame
    nearest_all_date_matches=[]
    nearest_all_dist_matches=[]
    nearest_all_model_matches=[]
    nearest_all_obs_matches=[]
    nearest_all_obs_lons_matches=[]
    nearest_all_obs_lats_matches=[]
    nearest_all_model_lons_matches=[]
    nearest_all_model_lats_matches=[]
    collocation_idx_lst=[]
    # create local variables before loop
    obs_lats = np.array(sa_obj.vars['latitude'])[cidx]
    obs_lons = np.array(sa_obj.vars['longitude'])[cidx]
    obs_val = np.array(sa_obj.vars[sa_obj.stdvarname])[cidx]
    # flatten numpy arrays
    model_lats = mc_obj.vars['latitude'].flatten()
    model_lons = mc_obj.vars['longitude'].flatten()
    model_val = mc_obj.vars[mc_obj.stdvarname].flatten()
    # moving window compensating for increasing latitudes
    try:
        moving_win = round(
                (distlim /
                 haversine(0,
                    np.max(np.abs(obs_lats)),
                    1,
                    np.max(np.abs(obs_lats)))
                ),
                2)
        if moving_win == 0.0:
            raise ValueError
    except (ValueError):
        moving_win = .6
    print ("Searching for matches with moving window of degree:",\
            moving_win)
    if collocation_idx is None:
        for j in range(len(obs_time_dt)):
            progress(j,str(int(len(obs_time_dt))),'')
#            try:
            for i in range(1):
                resultlst = collocation_loop(\
                    j,obs_time_dt,datein,distlim,
                    mc_obj.model,\
                    obs_lats,obs_lons,obs_val,\
                    model_lats,model_lons,model_val,\
                    moving_win)
                nearest_all_date_matches.append(resultlst[0])
                nearest_all_dist_matches.append(resultlst[1])
                nearest_all_model_matches.append(resultlst[2])
                nearest_all_obs_matches.append(resultlst[3])
                nearest_all_obs_lons_matches.append(resultlst[4])
                nearest_all_obs_lats_matches.append(resultlst[5])
                nearest_all_model_lons_matches.append(resultlst[6])
                nearest_all_model_lats_matches.append(resultlst[7])
                collocation_idx_lst.append(resultlst[8])
#            except:
#                print ("Collocation error -> no collocation:", 
#                        sys.exc_info()[0])
        results_dict = {
                'valid_date':np.array(datein),
                'time':np.array(nearest_all_date_matches),
                'distance':np.array(nearest_all_dist_matches),
                'model_values':np.array(nearest_all_model_matches),
                'model_lons':np.array(nearest_all_model_lons_matches),
                'model_lats':np.array(nearest_all_model_lats_matches),
                'obs_values':np.array(nearest_all_obs_matches),
                'obs_lons':np.array(nearest_all_obs_lons_matches),
                'obs_lats':np.array(nearest_all_obs_lats_matches),
                'collocation_idx':np.array(collocation_idx_lst)
                }
    else:
        results_dict = {'model_values':model_val[collocation_idx]}
    return results_dict


class collocation_class():
    '''
    draft of envisioned collocation class object
    '''

    def __init__(self,mc_obj,sa_obj=None,st_obj=None,collocation_idx=None,
    distlim=None):
        print ('# ----- ')
        print (" ### Initializing collocation_class instance ###")
        print ('# ----- ')
        results_dict = collocate(mc_obj,
                                sa_obj=sa_obj,
                                st_obj=st_obj,
                                collocation_idx=collocation_idx,
                                distlim=distlim)
        if sa_obj is not None:
            obs = sa_obj.sat
        if st_obj is not None:
            obs = sa_obj.stat
        # define class variables
        self.fc_date = mc_obj.fc_date
        self.init_date = mc_obj.init_date
        self.sdate = sa_obj.sdate
        self.edate = sa_obj.edate
        self.model = mc_obj.model
        self.obs = obs
        self.varalias = mc_obj.varalias
        self.stdvarname = mc_obj.stdvarname
        #self.vars = vardict # divided into model and obs
        self.vars = results_dict
