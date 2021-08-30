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
import pyresample
import xarray as xa
import pyproj
from tqdm import tqdm
from copy import deepcopy

# own imports
from wavy.utils import haversine, haversine_new, collocate_times
from wavy.utils import progress, make_fc_dates
from wavy.utils import make_pathtofile
from wavy.utils import hour_rounder
from wavy.utils import NoStdStreams
from wavy.modelmod import model_class, make_model_filename_wrapper
from wavy.modelmod import get_model_filedate, get_filevarname
from wavy.modelmod import model_class,get_model
from wavy.ncmod import dumptonc_ts_collocation
from wavy.ncmod import find_attr_in_nc, ncdumpMeta
from wavy.satmod import satellite_class
from wavy.stationmod import station_class
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

def collocation_fct(obs_lons,obs_lats,model_lons,model_lats):
    grid = pyresample.geometry.GridDefinition(\
                                lats=model_lats, \
                                lons=model_lons)
    # Define some sample points
    swath = pyresample.geometry.SwathDefinition(lons=obs_lons,
                                                lats=obs_lats)
    # Determine nearest (great circle distance) neighbour in the grid.
    valid_input_index, valid_output_index, index_array, distance_array = \
                            pyresample.kd_tree.get_neighbour_info(
                            source_geo_def=grid,
                            target_geo_def=swath,
                            radius_of_influence=1000000000,
                            neighbours=1)
    # get_neighbour_info() returns indices in the
    # flattened lat/lon grid. Compute the 2D grid indices:
    index_array_2d = np.unravel_index(index_array, grid.shape)
    return  index_array_2d, distance_array, valid_output_index,

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
    except (FileNotFoundError, OSError) as e:
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
        if (key != 'time_unit' and key !='meta'):
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
        col_obj = collocation_class( mc_obj_in=mc_obj,
                                     obs_obj_in=obs_obj,
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
                model_vals.append( mc_obj.vars[\
                                   mc_obj.stdvarname][ \
                                   col_obj.vars['collocation_idx_x'],\
                                   col_obj.vars['collocation_idx_y']\
                                   ][0] )
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
                                target_t = obs_obj.vars['datetime'],
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
    col_obj.vars['collocation_idx_x'] = col_obj.vars['collocation_idx_x']*\
                                        len(col_obj.vars['datetime'])
    col_obj.vars['collocation_idx_y'] = col_obj.vars['collocation_idx_y']*\
                                        len(col_obj.vars['datetime'])
    col_obj.vars['model_lats'] = col_obj.vars['model_lats']*\
                                       len(col_obj.vars['datetime'])
    col_obj.vars['model_lons'] = col_obj.vars['model_lons']*\
                                        len(col_obj.vars['datetime'])
    results_dict = col_obj.vars
    return results_dict

def collocate_satellite_ts(obs_obj=None,model=None,distlim=None,\
    leadtime=None,date_incr=None):
    """
    Some info
    """
    fc_date = make_fc_dates(obs_obj.sdate,obs_obj.edate,date_incr)
    fc_date = find_valid_fc_dates_for_model_and_leadtime(\
                            fc_date,model,leadtime)
    results_dict = {
            'valid_date':[],
            'time':[],
            'time_unit':obs_obj.vars['time_unit'],
            'datetime':[],
            'distance':[],
            'model_values':[],
            'model_lons':[],
            'model_lats':[],
            'obs_values':[],
            'obs_lons':[],
            'obs_lats':[],
            'collocation_idx_x':[],
            'collocation_idx_y':[],
            }
    # get model coordinates
    #filestr = make_model_filename_wrapper(model=model,
    #                                      fc_date=fc_date[0],
    #                                      leadtime=leadtime)
    #meta = ncdumpMeta(filestr)
    #flon = get_filevarname(model,'lons',variable_info,\
    #                                model_dict,meta)
    #flat = get_filevarname(model,'lats',variable_info,\
    #                                model_dict,meta)
    #M = xa.open_dataset(filestr, decode_cf=True)
    #model_lons = M[flon].data
    #model_lats = M[flat].data
    for i in tqdm(range(len(fc_date))):
#        for j in range(1):
        with NoStdStreams():
            try:
                # filter needed obs within time period
                idx = collocate_times( obs_obj.vars['datetime'],
                                       target_t = [fc_date[i]],
                                       twin = obs_obj.twin )
                # make tmp obs_obj with filtered data
                obs_obj_tmp = deepcopy(obs_obj)
                obs_obj_tmp.vars['time'] = list(\
                        np.array(obs_obj.vars['time'])[idx] )
                obs_obj_tmp.vars['latitude'] = list(\
                        np.array(obs_obj.vars['latitude'])[idx] )
                obs_obj_tmp.vars['longitude'] = list(\
                        np.array(obs_obj.vars['longitude'])[idx] )
                obs_obj_tmp.vars[obs_obj.stdvarname] = \
                        list(np.array(\
                        obs_obj.vars[obs_obj.stdvarname])[idx] )
                # collocate
                #filestr = make_model_filename_wrapper(model=model,
                #                              fc_date=fc_date[i],
                #                              leadtime=leadtime)
                #M = xa.open_dataset(filestr, decode_cf=True)
                #filevarname = get_filevarname(model,obs_obj.varalias,
                #                    variable_info,model_dict,meta)
                ## need function if variable consists of vector components
                #model_vals = M.sel(time=fc_date[i])[filevarname].data
                vardict,_,_,_,_ = get_model(model=model,
                                    fc_date=fc_date[i],
                                    varalias=obs_obj.varalias,
                                    leadtime=leadtime)
                results_dict_tmp = collocate_field(\
                                datein=fc_date[i],\
#                                model_lats=model_lats,\
#                               model_lons=model_lons,\
#                               model_vals=model_vals,\
                                model_lats=vardict['latitude'],\
                                model_lons=vardict['longitude'],\
                                model_vals=vardict[obs_obj.stdvarname],\
                                obs_obj=obs_obj_tmp,\
                                distlim=distlim )
                # append to dict
                results_dict['valid_date'].append(fc_date[i])
                results_dict['time'].append(results_dict_tmp['time'])
                results_dict['datetime'].append(results_dict_tmp['datetime'])
                results_dict['distance'].append(results_dict_tmp['distance'])
                results_dict['model_values'].append(results_dict_tmp['model_values'])
                results_dict['model_lons'].append(results_dict_tmp['model_lons'])
                results_dict['model_lats'].append(results_dict_tmp['model_lats'])
                results_dict['obs_values'].append(results_dict_tmp['obs_values'])
                results_dict['obs_lats'].append(results_dict_tmp['obs_lats'])
                results_dict['obs_lons'].append(results_dict_tmp['obs_lons'])
                results_dict['collocation_idx_x'].append(\
                                results_dict_tmp['collocation_idx_x'])
                results_dict['collocation_idx_y'].append(\
                                results_dict_tmp['collocation_idx_y'])
                if 'results_dict_tmp' in locals():
                    del results_dict_tmp
            except ValueError as e:
                print(e)
    # flatten all aggregated entries
    #results_dict['valid_date'] = flatten(results_dict['valid_date'])
    results_dict['time'] = flatten(results_dict['time'])
    results_dict['datetime'] = flatten(results_dict['datetime'])
    results_dict['distance'] = flatten(results_dict['distance'])
    results_dict['model_values'] = flatten(results_dict['model_values'])
    results_dict['model_lons'] = flatten(results_dict['model_lons'])
    results_dict['model_lats'] = flatten(results_dict['model_lats'])
    results_dict['obs_values'] = flatten(results_dict['obs_values'])
    results_dict['obs_lats'] = flatten(results_dict['obs_lats'])
    results_dict['obs_lons'] = flatten(results_dict['obs_lons'])
    results_dict['collocation_idx_x'] = flatten(\
                                results_dict['collocation_idx_x'])
    results_dict['collocation_idx_y'] = flatten(\
                                results_dict['collocation_idx_y'])
    return results_dict

def collocate_field(mc_obj=None,obs_obj=None,col_obj=None,distlim=None,
                    datein=None,model_lats=None,model_lons=None,
                    model_vals=None):
    """
    Some info
    """
    if mc_obj is not None:
        datein = netCDF4.num2date(mc_obj.vars['time'],mc_obj.vars['time_unit'])
        model_lats = mc_obj.vars['latitude']
        model_lons = mc_obj.vars['longitude']
        model_vals = mc_obj.vars[mc_obj.stdvarname]
        if len(model_lats.shape)==1:
            model_lons,model_lats = np.meshgrid(model_lons,model_lats)
    dtime = netCDF4.num2date(obs_obj.vars['time'],
                             obs_obj.vars['time_unit'])
    if isinstance(dtime,np.ndarray):
        dtime = list(dtime)
    if isinstance(datein,np.ndarray):
        datein = list(datein)
    if isinstance(datein,datetime):
        datein = [datein]
    cidx = collocate_times(dtime,target_t=datein,twin=obs_obj.twin)
    obs_time_dt = np.array(dtime)[cidx]
    obs_time_dt = [datetime(t.year,t.month,t.day,t.hour,t.minute,t.second)\
                   for t in obs_time_dt]
    datein = [datetime(t.year,t.month,t.day,t.hour,t.minute,t.second)\
                   for t in datein]
    obs_time = np.array(obs_obj.vars['time'])[cidx]
    obs_time_unit = obs_obj.vars['time_unit']
    # Compare wave heights of satellite with model with
    #   constraint on distance and time frame
    # 1. time constraint
    obs_lats = np.array(obs_obj.vars['latitude'])[cidx]
    obs_lons = np.array(obs_obj.vars['longitude'])[cidx]
    obs_vals = np.array(obs_obj.vars[obs_obj.stdvarname])[cidx]
    if distlim == None:
        distlim = 6
    if (col_obj is None):
        print ("No collocation idx available")
        print (len(obs_time_dt),"footprints to be collocated")
        print ("Perform collocation with distance limit\n",\
                "distlim:",distlim)
        index_array_2d, distance_array, valid_output_index =\
                                collocation_fct(
                                obs_lons, obs_lats,
                                model_lons, model_lats)
        # caution: index_array_2d is tuple
        # impose distlim
        dist_idx = np.where( (distance_array<distlim*1000)&\
                             (~np.isnan(\
                                model_vals[index_array_2d[0],\
                                index_array_2d[1]])) )[0]
        idx_x = index_array_2d[0][dist_idx]
        idx_y = index_array_2d[1][dist_idx]
        results_dict = {
            'valid_date':datein,
            'time':list(obs_time[dist_idx]),
            'time_unit':obs_time_unit,
            'datetime':list(np.array(obs_time_dt)[dist_idx]),
            'distance':list(distance_array[dist_idx]),
            'model_values':list(model_vals[idx_x,\
                                           idx_y]),
            'model_lons':list(model_lons[idx_x,\
                                         idx_y]),
            'model_lats':list(model_lats[idx_x,\
                                         idx_y]),
            'obs_values':list(obs_vals[dist_idx]),
            'obs_lons':list(obs_lons[dist_idx]),
            'obs_lats':list(obs_lats[dist_idx]),
            'collocation_idx_x':list(idx_x),
            'collocation_idx_y':list(idx_y),
            }
    elif (col_obj is not None and \
    len(col_obj.vars['collocation_idx'][0]) > 0):
        print("Collocation idx given through collocation_class object")
        results_dict = col_obj.vars
        results_dict['model_values'] = list(model_vals[\
                                        col_obj.vars['collocation_idx_x'],
                                        col_obj.vars['collocation_idx_y']])
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
    if (mc_obj is None and model is not None and obs_obj is not None\
    and isinstance(obs_obj,station_class)):
        results_dict = collocate_station_ts(obs_obj=obs_obj,
                                            model=model,\
                                            distlim=distlim,\
                                            leadtime=leadtime,\
                                            date_incr=date_incr)
    elif (mc_obj is None and model is not None and obs_obj is not None\
    and isinstance(obs_obj,satellite_class)):
        results_dict = collocate_satellite_ts(obs_obj=obs_obj,
                                            model=model,\
                                            distlim=distlim,\
                                            leadtime=leadtime,\
                                            date_incr=date_incr)
    else:
        datein = mc_obj.fc_date
        results_dict = collocate_field( mc_obj=mc_obj,\
                                        obs_obj=obs_obj,\
                                        col_obj=col_obj,\
                                        distlim=distlim )
    return results_dict


class collocation_class():
    '''
    draft of envisioned collocation class object
    '''

    def __init__(self,mc_obj_in=None,obs_obj_in=None,
        col_obj_in=None,model=None,distlim=None,leadtime=None,
        date_incr=None):
        print ('# ----- ')
        print (" ### Initializing collocation_class object ###")
        print (" Please wait ...")
        # make clones to prevent overwriting
        mc_obj = deepcopy(mc_obj_in)
        obs_obj = deepcopy(obs_obj_in)
        col_obj = deepcopy(col_obj_in)
        if isinstance(obs_obj,satellite_class):
            self.obsname = obs_obj.sat
            self.sat = obs_obj.sat
            self.obstype = "satellite_altimeter"
            self.region = obs_obj.region
        if isinstance(obs_obj,station_class):
            if 'twin' in station_dict['platform'][obs_obj.platform].keys():
                obs_obj.twin =  station_dict['platform']\
                                [obs_obj.platform]['twin']
            else:
                obs_obj.twin = None
            self.obsname = obs_obj.platform + '_' +  obs_obj.sensor
            self.obstype = 'platform'
            self.platform = obs_obj.platform
            self.sensor = obs_obj.sensor
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
#        try:
        for i in range(1):
            t0=time.time()
            results_dict = collocate(mc_obj=mc_obj,
                                    obs_obj=obs_obj,
                                    col_obj=col_obj,
                                    model=model,
                                    distlim=distlim,
                                    leadtime=leadtime,
                                    date_incr=date_incr)
            t1=time.time()
            print("Time used for collocation:",round(t1-t0,2),"seconds")
            self.vars = results_dict
            self.fc_date = results_dict['datetime']
            print(len(self.vars['time'])," values collocated")
            if 'superob' in vars(obs_obj).keys():
                self.superob = obs_obj.superob
                self.outlier_detection = obs_obj.outlier_detection
                self.missing_data = obs_obj.missing_data
            print (" ### Collocation_class object initialized ###")
#        except Exception as e:
#            print(e)
#            self.error = e
#            print ("! No collocation_class object initialized !")
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
                                                mission=self.sat,
                                                region=self.region,
                                                leadtime=leadtimestr)
                        title = ( 'Collocation of ' + self.stdvarname
                                + ' observations from ' + self.sat
                                + ' vs ' + self.model)
                dumptonc_ts_collocation(self,pathtofile,title)
                tmpdate = tmpdate + relativedelta(months = +1)
        return

    def validate_collocated_values(self,**kwargs):
        dtime = self.vars['datetime']
        mods = self.vars['model_values']
        obs = self.vars['obs_values']
        sdate = self.vars['datetime'][0]
        edate = self.vars['datetime'][-1]
        validation_dict = validate_collocated_values(
                                dtime,obs,mods,\
                                sdate=sdate,edate=edate,\
                                **kwargs)
        return validation_dict

def validate_collocated_values(dtime,obs,mods,**kwargs):
    target_t, sdate, edate, twin = None, None, None, None
    if ('col_obj' in kwargs.keys() and kwargs['col_obj'] is not None):
        mods = col_obj.vars['model_values']
        obs = col_obj.vars['obs_values']
        dtime = col_obj.vars['datetime']
    # get idx for date and twin
    from utils import collocate_times
    if 'target_t' in kwargs.keys():
        target_t = kwargs['target_t']
    if 'sdate' in kwargs.keys():
        sdate = kwargs['sdate']
    if 'edate' in kwargs.keys():
        edate = kwargs['edate']
    if 'twin' in kwargs.keys():
        twin = kwargs['twin']
    idx = collocate_times(dtime,
                          target_t=target_t,
                          sdate=sdate,
                          edate=edate,
                          twin=twin)
    mods = np.array(mods)[idx]
    obs = np.array(obs)[idx]
    results_dict = {'model_values':mods,'obs_values':obs}
    # validate
    from validationmod import validate,disp_validation
    validation_dict = validate(results_dict)
    disp_validation(validation_dict)
    return validation_dict

