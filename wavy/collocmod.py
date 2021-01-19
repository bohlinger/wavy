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
    sat_rlats,sat_rlons,sat_rval,
    model_rlats,model_rlons,model_rval,
    moving_win):
    from utils import haversine
    lat_win = 0.1
    if model in model_dict:
        sat_rlat=sat_rlats[j]
        sat_rlon=sat_rlons[j]
        # constraints to reduce workload
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
        if sys.version_info <= (3, 0):
            distlst=map(
                        haversine,
                        [sat_rlon]*len(model_rlons_new),
                        [sat_rlat]*len(model_rlons_new),
                        model_rlons_new,model_rlats_new
                        )
        else:
            distlst=list(map(
                        haversine,
                        [sat_rlon]*len(model_rlons_new),
                        [sat_rlat]*len(model_rlons_new),
                        model_rlons_new,model_rlats_new
                        ))
        tmp_idx2 = distlst.index(np.min(distlst))
        idx_valid = tmp_idx[tmp_idx2]
        if (distlst[tmp_idx2]<=distlim and model_rval[idx_valid]>=0):
            nearest_all_dist_matches=distlst[tmp_idx2]
            nearest_all_date_matches=obs_time_dt[j]
            nearest_all_model_matches=\
                           model_rval[idx_valid]
            nearest_all_sat_matches=sat_rval[j]
            nearest_all_sat_lons_matches=sat_rlon
            nearest_all_sat_lats_matches=sat_rlat
            nearest_all_model_lons_matches=\
                            model_rlons[idx_valid]
            nearest_all_model_lats_matches=\
                            model_rlats[idx_valid]
            return nearest_all_date_matches,nearest_all_dist_matches,\
                nearest_all_model_matches,nearest_all_sat_matches,\
                nearest_all_sat_lons_matches, nearest_all_sat_lats_matches,\
                nearest_all_model_lons_matches, \
                nearest_all_model_lats_matches, idx_valid
        else:
           return

#def collocate(model,model_val,model_lats,model_lons,model_time_dt,\
#    sa_obj,var,datein,distlim=6,idx_valid=None):
def collocate(mc_obj,sa_obj=None,st_obj=None,idx_valid=None,
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
    idx_valid_lst=[]
    # create local variables before loop
    obs_lats = np.array(sa_obj.vars['latitude'])[cidx]
    obs_lons = np.array(sa_obj.vars['longitude'])[cidx]
    obs_val = np.array(sa_obj.vars[sa_obj.stdvarname])[cidx]
    # flatten numpy arrays
    model_lats = mc_obj.vars['latitude'].flatten()
    model_lons = mc_obj.vars['longitude'].flatten()
    model_val = mc_obj.vars[mc_obj.stdvarname].flatten()
    #model_rval = model_val.squeeze().flatten()
    #model_rlons = model_lons.flatten()
    #model_rlats = model_lats.flatten()
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
    if idx_valid is None:
        for j in range(len(obs_time_dt)):
            progress(j,str(int(len(obs_time_dt))),'')
            try:
#            for i in range(1):
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
                idx_valid_lst.append(resultlst[8])
            except:
                print ("Collocation error -> no collocation:", 
                        sys.exc_info()[0])
#                pass
        results_dict = {
                'valid_date':np.array(datein),
                'date_matches':np.array(nearest_all_date_matches),
                'dist_matches':np.array(nearest_all_dist_matches),
                'model_matches':np.array(nearest_all_model_matches),
                'obs_matches':np.array(nearest_all_obs_matches),
                'obs_lons_matches':np.array(nearest_all_obs_lons_matches),
                'obs_lats_matches':np.array(nearest_all_obs_lats_matches),
                'model_lons_matches':np.array(nearest_all_model_lons_matches),
                'model_lats_matches':np.array(nearest_all_model_lats_matches),
                'idx_valid':np.array(idx_valid_lst)
                }
    else:
        results_dict = {'model_matches':model_rval[idx_valid]}
    return results_dict


class collocation_class():
    '''
    draft of envisioned collocation class object
    '''

    def __init__(self,mc_obj,sa_obj=None,st_obj=None,idx_valid=None,
    distlim=None):
        print ('# ----- ')
        print (" ### Initializing collocation_class instance ###")
        print ('# ----- ')
        results_dict = collocate(mc_obj,
                                sa_obj=sa_obj,
                                st_obj=st_obj,
                                idx_valid=idx_valid,
                                distlim=distlim)
#                    model,model_val,model_lats,model_lons,model_time_dt
#                    sa_obj,var,datein,distlim=6,idx_valid=None)
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
