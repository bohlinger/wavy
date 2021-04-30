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
# standard library imports
import sys
import yaml
import os
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import math
from dateutil.relativedelta import relativedelta
from copy import deepcopy
import time

# own imports
from utils import progress, hour_rounder, collocate_times
from utils import make_fc_dates
#from collocmod import collocation_class
from ncmod import ncdumpMeta, get_varname_for_cf_stdname_in_ncfile

# --- global functions ------------------------------------------------#
"""
definition of some global functions
"""
def get_model_filedate(model,fc_date,leadtime):
    '''
    get init_date for latest model output file and checks if available
    '''
    if ('init_times' in model_dict[model].keys() and 
    model_dict[model]['init_times'] is not None):
        init_times = \
            np.array(model_dict[model]['init_times']).astype('float')
    else:
        print('init_times for chosen model not specified in config file')
        print('Assuming continuous simulation with hourly values')
        init_times = np.array(range(25)).astype('float')
    date = fc_date - timedelta(hours=leadtime)
    date_hour = hour_rounder(date).hour
    if date_hour in init_times:
        init_diffs = date_hour - init_times
        init_diffs[init_diffs<0] = np.nan
        h_idx = np.where(init_diffs==\
                        np.min(init_diffs[~np.isnan(init_diffs)]))
        h = int(init_times[h_idx[0][0]])
        return datetime(date.year,date.month,date.day,h)
    else:
        print('Leadtime', leadtime , \
              'not available for date' ,fc_date, '!')
        return False

def make_model_filename(model,fc_date,leadtime):
    """
    creates/returns filename based on fc_date,leadtime
    """
    if model in model_dict:
        if 'xtra_h' in model_dict[model]:
            filedate = get_model_filedate(model,fc_date,leadtime)
            pathdate = filedate + timedelta(hours=leadtime) \
                                * model_dict[model]['lt_switch_p']
            tmpstr = model_dict[model]['file_template']
            for i in range(model_dict[model]['nr_filedates']):
                filedatestr = model_dict[model]['filedate_formats'][i]
                replacestr = (filedate \
                            + timedelta(hours = leadtime \
                                        - (leadtime % \
                                           model_dict[model]['init_step']))
                                    * model_dict[model]['lt_switch_f'][i]
                            + timedelta(hours = \
                                    model_dict[model]['xtra_h'][i])).\
                            strftime(filedatestr)
                tmpstr = tmpstr.replace('filedate',replacestr,1)
            filename = (pathdate.strftime(model_dict[model]['path_template'])
                        + tmpstr)
        else:
            filedate = get_model_filedate(model,fc_date,leadtime)
            filename = (
                filedate.strftime(model_dict[model]['path_template'])
                + filedate.strftime(model_dict[model]['file_template'])
                )
    else:
        raise ValueError("Chosen model is not specified in model_specs.yaml")
    return filename

def make_model_filename_wrapper(model,fc_date,leadtime):
    if leadtime is None:
        leadtime = 'best'
    if (isinstance(fc_date,datetime) and leadtime != 'best'):
        filename = make_model_filename(model,fc_date,leadtime)
    elif (isinstance(fc_date,datetime) and leadtime == 'best'):
        leadtime = generate_bestguess_leadtime(model,fc_date)
        filename = make_model_filename(model,fc_date,leadtime)
    elif (isinstance(fc_date,list) and isinstance(leadtime,int)):
        filename = [make_model_filename(model,date,leadtime) \
                    for date in fc_date]
    elif (isinstance(fc_date,list) and leadtime == 'best'):
        leadtime = generate_bestguess_leadtime(model,fc_date)
        filename = [make_model_filename(model,fc_date[i],leadtime[i]) \
                    for i in range(len(fc_date))]
    return filename

def get_model_fc_mode(filestr,model,fc_date,varalias=None):
    """ 
    fct to retrieve model data for correct time
    """
    vardict = {}
    print ("Get model data according to selected date ....")
    print(filestr)
    model_meta = ncdumpMeta(filestr)
    f = netCDF4.Dataset(filestr,'r')
    stdvarname = variable_info[varalias]['standard_name']
    # get coordinates and time
    lonsname = get_filevarname(model,'lons',variable_info,
                                model_dict,model_meta)
    latsname = get_filevarname(model,'lats',variable_info,
                                model_dict,model_meta)
    timename = get_filevarname(model,'time',variable_info,
                                model_dict,model_meta)
    model_lons = f.variables[lonsname][:]
    vardict[variable_info['lons']['standard_name']]=model_lons
    model_lats = f.variables[latsname][:]
    vardict[variable_info['lats']['standard_name']]=model_lats
    model_time = f.variables[timename]
    model_time_dt = list( netCDF4.num2date(model_time[:],
                        units = model_time.units) )
    # get other variables e.g. Hs [time,lat,lon]
    filevarname = get_filevarname(model,varalias,variable_info,
                                    model_dict,model_meta)
    if (type(filevarname) is dict):
        print('Target variable can be computed from vector \n' 
              'components with the following aliases: ', filevarname)
        model_time_dt_valid = model_time_dt[model_time_dt.index(fc_date)]
        model_time_valid = float(model_time[model_time_dt.index(fc_date)])
        model_time_unit = model_time.units
        vardict[variable_info['time']['standard_name']] = model_time_valid
        vardict['datetime'] = model_time_dt_valid
        vardict['time_unit'] = model_time_unit
        for key in filevarname.keys():
            filevarname_dummy = get_filevarname(model,
                                    filevarname[key][0],
                                    variable_info,
                                    model_dict,model_meta)
            if filevarname_dummy is not None:
                print(filevarname[key][0], 'exists')
                break
        print('Use aliases:',filevarname[key])
        model_var_dummy = f.variables[filevarname_dummy]
        if len(model_var_dummy.shape)>2: # for multiple time steps
            model_var_valid_tmp = \
                model_var_dummy[model_time_dt.index(fc_date),:,:].squeeze()**2
            for i in range(1,len(filevarname[key])):
                filevarname_dummy = get_filevarname(model,
                                                    filevarname[key][i],
                                                    variable_info,
                                                    model_dict,
                                                    model_meta)
                model_var_valid_tmp += \
                    f.variables[filevarname_dummy][
                        model_time_dt.index(fc_date),:,:
                        ].squeeze()**2
            model_var_valid = np.sqrt(model_var_valid_tmp)
        else:# if only one time step
            model_var_valid_tmp = model_var_dummy[:,:].squeeze()**2
            for i in range(1,len(filevarname[key])):
                filevarname_dummy = get_filevarname(model,
                                                    filevarname[key][i],
                                                    variable_info,
                                                    model_dict,
                                                    model_meta)
                model_var_valid_tmp += \
                    f.variables[filevarname_dummy][:,:].squeeze()**2
            model_var_valid = np.sqrt(model_var_valid_tmp)
        vardict[stdvarname] = model_var_valid
    else:
        model_time_dt_valid = model_time_dt[model_time_dt.index(fc_date)]
        model_time_valid = float(model_time[model_time_dt.index(fc_date)])
        model_time_unit = model_time.units
        vardict[variable_info['time']['standard_name']] = model_time_valid
        vardict['datetime'] = model_time_dt_valid
        vardict['time_unit'] = model_time_unit
        model_var_link = f.variables[filevarname]
        if len(model_var_link.shape)>2: # for multiple time steps
            model_var_valid = \
                model_var_link[model_time_dt.index(fc_date),:,:].squeeze()
        else:# if only one time step
            model_var_valid = model_var_link[:,:].squeeze()
        vardict[variable_info[varalias]['standard_name']] = \
                                                    model_var_valid
    # transform masked array to numpy array with NaNs
    f.close()
    vardict[variable_info[varalias]['standard_name']] = \
        vardict[variable_info[varalias]['standard_name']].filled(np.nan)
    vardict['model_meta'] = model_meta
    return vardict, filevarname

def generate_bestguess_leadtime(model,fc_date):
    """
    fct to return leadtimes for bestguess
    """
    if isinstance(fc_date,list):
        leadtime = \
            [generate_bestguess_leadtime(model,date) for date in fc_date]
    else:
        init_times = \
            np.array(model_dict[model]['init_times']).astype('float')
        init_step = \
            np.array(model_dict[model]['init_step']).astype('int')
        diffs = fc_date.hour - np.array(init_times)
        leadtime = int(np.min(diffs[diffs>=0]))
    return leadtime

def get_filevarname(model,varalias,variable_info,model_dict,ncdict):
    stdname = variable_info[varalias]['standard_name']
    print('Get filevarname for \n' +
          'stdvarname:',stdname,'\n' + 
          'varalias:',varalias)
    filevarname = get_varname_for_cf_stdname_in_ncfile(ncdict,stdname)
    if (filevarname is None and 'alias' in variable_info[varalias]):
        filevarname = get_varname_for_cf_stdname_in_ncfile(
                        ncdict,variable_info[varalias]['alias'])
    if (filevarname is not None and len(filevarname) > 1):
        print('!!! standard_name: ',stdname,' is not unique !!!',
            '\nThe following variables have the same standard_name:\n',
            filevarname)
        print('Searching model_specs.yaml config file for definition')
        filevarname = None
    if filevarname is not None:
        return filevarname[0]
    if (filevarname is None 
    and varalias in model_dict[model]['vars'].keys()):
        filevarname = model_dict[model]['vars'][varalias]
        print('Variable defined in model_specs.yaml is:')
        print(varalias, '=', filevarname)
        return filevarname
    elif (filevarname is None
    and varalias not in model_dict[model]['vars'].keys()
    and 'aliases_of_vector_components' in variable_info[varalias]):
        print('Checking variable_info if variable can be ' + 
        'computed from vector components')
        filevarname = variable_info[varalias]['aliases_of_vector_components']
        return filevarname
    else:
        print('!!! variable not defined or '
            + 'available in model output file !!!')

def get_model(model=None,sdate=None,edate=None,date_incr=None,
    fc_date=None,leadtime=None,varalias=None,
    st_obj=None,distlim=None):
    """ 
    toplevel function to get model data
    """
    if st_obj is not None:
        sdate = st_obj.sdate
        edate = st_obj.edate
        varalias = st_obj.varalias
    if (sdate is not None and edate is not None 
    and date_incr is not None):
        fc_date = make_fc_dates(sdate,edate,date_incr)
    
    filestr = make_model_filename_wrapper(model=model,
                                          fc_date=fc_date,
                                          leadtime=leadtime)
    
    if (isinstance(filestr,list) and st_obj is None):
        vardict, \
        filevarname = get_model_fc_mode(filestr=filestr[0],model=model,
                                    fc_date=fc_date[0],varalias=varalias)
        vardict[variable_info[varalias]['standard_name']]=\
                    [vardict[variable_info[varalias]['standard_name']]]
        vardict['time'] = [vardict['time']]
        vardict['datetime'] = [vardict['datetime']]
        for i in range(1,len(filestr)):
            tmpdict, \
            filevarname = get_model_fc_mode(filestr=filestr[i],model=model,
                                    fc_date=fc_date[i],varalias=varalias)
            vardict[variable_info[varalias]['standard_name']].append(
                    tmpdict[variable_info[varalias]['standard_name']])
            vardict['time'].append(tmpdict['time'])
            vardict['datetime'].append(tmpdict['datetime'])
        vardict[variable_info[varalias]['standard_name']]=\
            np.array(vardict[variable_info[varalias]['standard_name']])
    else:
        vardict, \
        filevarname = get_model_fc_mode(filestr=filestr,model=model,
                                    fc_date=fc_date,varalias=varalias)
        vardict['time'] = [vardict['time']]
        vardict['datetime'] = [vardict['datetime']]
    return vardict, fc_date, leadtime, filestr, filevarname
# ---------------------------------------------------------------------#

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                        '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)


class model_class():

    '''
    class to read and process model data 
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and 
    station classes.
    '''

    def __init__(self,model='mwam4',sdate=None,edate=None,date_incr=1,
    fc_date=None,leadtime=None,varalias='Hs',st_obj=None,distlim=6):
        print ('# ----- ')
        print (" ### Initializing model_class object ###")
        if st_obj is not None:
            sdate = st_obj.sdate
            edate = st_obj.edate
            varalias = st_obj.varalias
        if fc_date is not None:
            print ("Requested time: ", str(fc_date))
        elif (edate is None and fc_date is None and sdate is not None):
            fc_date = sdate
            print ("Requested time: ", str(fc_date))
        elif (sdate is None and edate is None and fc_date is None):
            now = datetime.now()
            fc_date = datetime(now.year,now.month,now.day,now.hour)
            print ("Requested time: ", str(fc_date))
        elif (sdate is not None and edate is not None 
        and date_incr is not None):
            # time frame function to access a temporal subset
            # --> not yet in use
            print ("Requested time frame: " +
                str(sdate) + " - " + str(edate))

        print('Please wait ...')
        t0=time.time()
        vardict, \
        fc_date, leadtime, \
        filestr, \
        filevarname = get_model(model=model,sdate=sdate,edate=edate,
                            date_incr=date_incr,fc_date=fc_date,
                            leadtime=leadtime,varalias=varalias,
                            st_obj=st_obj)

        stdname = variable_info[varalias]['standard_name']
        varname = filevarname
        # define class variables
        self.fc_date = fc_date
        self.sdate = sdate
        self.edate = edate
        self.model = model
        self.varalias = varalias
        self.varname = varname
        self.stdvarname = stdname
        self.vars = vardict
        self.filestr = filestr
        t1=time.time()
        print("Time used for retrieving model data:",round(t1-t0,2),"seconds")
        print (" ### model_class object initialized ###")
        print ('# ----- ')
