#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to read insitu obs for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import os
import numpy as np
from datetime import datetime, timedelta
import datetime
import os
from dateutil.relativedelta import relativedelta
from copy import deepcopy
import pylab as pl
from datetime import datetime
import pandas as pd
import netCDF4

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import get_varlst_from_nc_1D
from wavy.ncmod import get_filevarname
from wavy.utils import collocate_times
from wavy.utils import get_pathtofile
from wavy.utils import convert_meteorologic_oceanographic
from wavy.utils import make_subdict
from wavy.utils import parse_date
from wavy.utils import flatten
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#
# read yaml config files:
insitu_dict = load_or_default('insitu_specs.yaml')
variable_info = load_or_default('variable_info.yaml')
d22_dict = load_or_default('d22_var_dicts.yaml')
# ---------------------------------------------------------------------#

def get_d22_ts(sdate,edate,basedate,nID,sensor,varalias,
pathlst,strsublst,dict_for_sub):
    sdatetmp = sdate
    edatetmp = edate
    sl = parse_d22(sdatetmp,edatetmp,pathlst,strsublst,dict_for_sub)
    var, timedt = extract_d22(sl,varalias,nID,sensor)
    time = np.array(
           [(t-basedate).total_seconds() for t in timedt]
           )
    return var, time, timedt

def get_d22_dict(**kwargs):
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    basedate = kwargs.get('basedate')
    nID = kwargs.get('nID')
    sensor = kwargs.get('sensor')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    strsublst = kwargs.get('strsublst')
    dict_for_sub = kwargs.get('dict_for_sub')
    stdvarname = variable_info[varalias]['standard_name']
    var, time, timedt = \
        get_d22_ts(sdate,edate,basedate,nID,sensor,varalias,\
                    pathlst,strsublst,dict_for_sub)
    if 'twin' in insitu_dict[nID]:
        idxtmp = collocate_times(unfiltered_t=timedt,\
                            sdate=sdate,edate=edate,
                            twin=insitu_dict[nID]['twin'])
    else:
        # default to allow for a 1 min variation
        idxtmp = collocate_times(unfiltered_t=timedt,\
                            sdate=sdate,edate=edate,
                            twin=1)
    # convert to list for consistency with other classes
    # and make sure that only steps with existing obs are included
    time = [time[i] for i in idxtmp if i < len(var)]
    timedt = [timedt[i] for i in idxtmp if i < len(var)]
    var = [np.real(var[i]) for i in idxtmp if i < len(var)]
    # rm double entries due to 10min spacing
    if ('unique' in kwargs.keys() and kwargs['unique'] is True):
        # delete 10,30,50 min times, keep 00,20,40
        # 1. create artificial time vector for collocation
        tmpdate = deepcopy(sdate)
        tmpdatelst = []
        while tmpdate<edate:
            tmpdatelst.append(tmpdate)
            tmpdate += timedelta(minutes=20)
        # 2. collocate times
        if 'twin' in insitu_dict[nID]:
            idxtmp = collocate_times(\
                        unfiltered_t=timedt,
                        target_t=tmpdatelst,
                        twin=insitu_dict[nID]['twin'])
        else:
            idxtmp = collocate_times(unfiltered_t=timedt,\
                            target_t=tmpdatelst,
                            twin=1)
        time = list(np.array(time)[idxtmp])
        timedt = list(np.array(timedt)[idxtmp])
        var = list(np.array(var)[idxtmp])
    lons = [insitu_dict[nID]['coords'][sensor]['lon']]\
            *len(var)
    lats = [insitu_dict[nID]['coords'][sensor]['lat']]\
            *len(var)
    vardict = {
                stdvarname:var,
                'time':time,
                'datetime':timedt,
                'time_unit':variable_info['time']['units'],
                'longitude':lons,
                'latitude':lats
                }
    return vardict

def make_frost_reference_time_period(sdate,edate):
    sdate = parse_date(sdate)
    edate = parse_date(edate)
    formatstr = '%Y-%m-%dT%H:%M:00.000Z'
    refstr = '{}/{}'.format(sdate.strftime(formatstr),
                            edate.strftime(formatstr))
    return refstr

def call_frost_api_v1(nID, varalias, frost_reference_time):
    import requests
    import dotenv
    dotenv.load_dotenv()
    client_id = os.getenv('CLIENT_ID', None)
    if client_id is None:
        print("No Frost CLIENT_ID given!")
    stdvarname = variable_info[varalias]['standard_name']
    ID = insitu_dict[nID]['ID']
    endpoint = 'https://frost-prod.met.no/api/v1/obs/met.no/filter/get?'
    #endpoint = 'https://frost-prod.met.no/api/v1/obs/met.no/kvkafka/get?'
    parameters = {
                'stationids': ID,
                'elementids': stdvarname,
                'time': frost_reference_time,
                'timeoffsets': 'default', # handled by filter
                'levels': 0,
                'incobs': 'true'
                }
    return requests.get(endpoint, parameters, auth=(client_id, client_id))

def get_frost_df_v1(r,nID,sensor,varalias):
    # base df
    df = pd.json_normalize(r.json()['data']['tseries'])
    # df to be concatenated initialized with time
    dfc = pd.json_normalize(r.json()
      ['data']['tseries'][0]['observations'])['time'].to_frame()
    stdvarname = variable_info[varalias]['standard_name']
    sensoridx = insitu_dict[nID]['sensor'][sensor]
    idx = df['header.extra.element.id']\
            [df['header.extra.element.id']==stdvarname]\
            .index.to_list()
    idxs = df['header.id.sensor']\
            [df['header.id.sensor']==sensoridx]\
            .index.to_list()
    dftmp = pd.json_normalize(r.json()\
                    ['data']['tseries'][idxs[sensoridx]]['observations'])\
                    ['body.data'].to_frame()
    dftmp = dftmp.rename(columns={ dftmp.columns[0]: stdvarname })
    dftmp[stdvarname] = dftmp[stdvarname].mask(dftmp[stdvarname] < 0, np.nan)
    dfc = pd.concat([dfc, dftmp.reindex(dfc.index)], axis=1)
    var = list(dfc[stdvarname].values)
    timedt = [parse_date(d) for d in dfc['time'].values]
    timedt = [dt.replace(tzinfo=None) for dt in timedt]
    time = netCDF4.date2num(timedt,variable_info['time']['units'])
    sensor = list(insitu_dict[nID]['sensor'].keys())[0]
    lons = [insitu_dict[nID]['coords'][sensor]['lon']]*len(var)
    lats = [insitu_dict[nID]['coords'][sensor]['lat']]*len(var)
    vardict = {
                stdvarname:var,
                'time':time,
                'datetime':timedt,
                'time_unit':variable_info['time']['units'],
                'longitude':lons,
                'latitude':lats
                }
    #return dfc
    return vardict

def get_frost_dict(**kwargs):
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    nID = kwargs.get('nID')
    sensor = kwargs.get('sensor')
    varalias = kwargs.get('varalias')
    reftp = make_frost_reference_time_period(sdate,edate)
    r = call_frost_api_v1(nID,varalias,reftp)
    vardict = get_frost_df_v1(r,nID,sensor,varalias)
    return vardict

def get_nc_dict(**kwargs):
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    nID = kwargs.get('nID')
    #sensor = kwargs.get('sensor')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    strsublst = kwargs.get('strsublst')
    dict_for_sub = kwargs.get('dict_for_sub')
    # loop from sdate to edate with dateincr
    stdvarname = variable_info[varalias]['standard_name']
    tmpdate = deepcopy(sdate)
    varlst = []
    lonlst = []
    latlst = []
    timelst = []
    dtimelst = []
    # make subdict
    subdict = make_subdict(strsublst,class_object_dict=dict_for_sub)
    while datetime(tmpdate.year,tmpdate.month,1)\
    <= datetime(edate.year,edate.month,1):
        # get pathtofile
        pathtofile = get_pathtofile(pathlst,strsublst,\
                                        subdict,tmpdate)
        # get ncdump
        ncdict = ncdumpMeta(pathtofile)
        # retrieve filevarname for varalias
        filevarname = get_filevarname(varalias,
                                      variable_info,
                                      insitu_dict[nID],
                                      ncdict)
        varstrlst = [filevarname,'longitude','latitude','time']
        # query
        vardict = get_varlst_from_nc_1D(pathtofile,
                                        varstrlst,
                                        sdate,edate)
        varlst.append(list(vardict[filevarname]))
        lonlst.append(list(vardict['longitude']))
        latlst.append(list(vardict['latitude']))
        timelst.append(list(vardict['time']))
        dtimelst.append(list(vardict['dtime']))
        # determine date increment
        file_date_incr = insitu_dict[nID]['src'].get('file_date_incr','m')
        if file_date_incr == 'm':
            tmpdate += relativedelta(months = +1)
        elif file_date_incr == 'Y':
            tmpdate += relativedelta(years = +1)
        elif file_date_incr == 'd':
            tmpdate += timedelta(days = +1)
    varlst = flatten(varlst)
    lonlst = flatten(lonlst)
    latlst = flatten(latlst)
    timelst = flatten(timelst)
    dtimelst = flatten(dtimelst)
    vardict = {
                stdvarname:varlst,
                'time':timelst,
                'datetime':dtimelst,
                'time_unit':variable_info['time']['units'],
                'longitude':lonlst,
                'latitude':latlst
                }
    # adjust conventions
    # check if variable is one with conventions
    if 'convention' in variable_info[varalias].keys():
        convention_set = False
        print('Chosen variable is defined with conventions')
        print('... checking if correct convention is used ...')
        # 1. check if clear from standard_name
        file_stdvarname = find_direction_convention(filevarname,ncdict)
        if "to_direction" in file_stdvarname:
            print('Convert from oceanographic to meteorologic convention')
            vardict[variable_info[varalias]['standard_name']] = \
                    convert_meteorologic_oceanographic(\
                        vardict[variable_info[varalias]['standard_name']])
            convention_set = True
        elif "from_direction" in file_stdvarname:
            print('standard_name indicates meteorologic convention')
            convention_set = True
            pass
        # 2. overwrite convention from config file
        if ('convention' in insitu_dict[nID].keys() and
        insitu_dict[model]['convention'] == 'oceanographic' and
        convention_set is False):
            print('Convention is set in config file')
            print('This will overwrite conventions from standard_name in file!')
            print('\n')
            print('Convert from oceanographic to meteorologic convention')
            vardict[variable_info[varalias]['standard_name']] = \
                    convert_meteorologic_oceanographic(\
                        vardict[variable_info[varalias]['standard_name']])
            convention_set = True
    return vardict

def parse_d22(sdate,edate,pathlst,strsublst,dict_for_sub):
    """
    Read all lines in file and append to sl
    """
    subdict = make_subdict(strsublst,class_object_dict=dict_for_sub)
    sl=[]
    for d in range(int(pl.date2num(sdate))-1,int(pl.date2num(edate))+2):
        try:
            pathtofile = get_pathtofile(pathlst,strsublst,
                                        subdict,pl.num2date(d))
            print('Parsing:', pathtofile)
            f = open(pathtofile, "r")
            sl = sl + f.readlines()
            f.close()
        except Exception as e:
            print('Error in parse_d22:')
            print(e)
    return sl

def floater(s):
    """
    Function that converts 's' to float32 or nan if floater throws exception
    """
    try:
        x = np.float32(s)
    except:
        x = np.nan
    return x

def find_category_for_variable(varalias):
    lst = [ i for i in d22_dict.keys() \
            if (varalias in d22_dict[i]) ]
    if len(lst) == 1:
        return lst[0]
    else: return None

def get_revised_categories(sl,category):
    """
    finds number of occurences of string (category) to determine
    revised_categories (type: list)
    """
    revised_categories = []
    idxlst = []
    count = 1
    searching = True
    while (searching is True or count<10):
        revised_category = category+str(count)
        if find_category(sl,revised_category) is True:
            revised_categories.append(revised_category)
            idxlst.append(count-1)
            count += 1
        else:
            searching = False
            count += 1
    return revised_categories,idxlst

def find_category(sl,category):
    for element in sl:
        if category in element:
            return True

def check_sensor_availability(idxlst,nID,sensor):
    idxyaml = insitu_dict[nID]['sensor'][sensor]
    if idxyaml in idxlst:
        return idxlst.index(idxyaml)
    else:
        return None

def extract_d22(sl,varalias,nID,sensor):
    """
    Extracting data of choice - reading sl from parse_d22
    CAUTION: 10min data is extracted for entire days only 00:00h - 23:50h
    Returns values of chosen variable (ts) and corresponding datetimes (dt)
    as type: np.array
    """
    print('Extracting data from parsed .d22-files')
    category = find_category_for_variable(varalias)
    revised_categories,idxlst = get_revised_categories(sl,category)
    print( 'Consistency check: \n'
           ' --> compare found #sensors against defined in insitu_specs.yaml')
    sensornr = len(insitu_dict[nID]['sensor'].keys())
    if len(revised_categories) == sensornr:
        print('Consistency check: OK!')
    else:
        print('Consistency check: Failed!')
        print(    '!!! Caution:\n'
                + 'found #sensor ('
                + str(len(revised_categories))
                + ') is not equal to defined ' +
                '#sensors ('
                + str(sensornr)
                + ') in insitu_specs.yaml')
    # check that the defined sensors are actually the ones being found
    check = check_sensor_availability(idxlst,nID,sensor)
    ts = []
    dt = []
    if check is not None:
        print('Sensor is available and defined in insitu_specs.yaml')
        for i, line in enumerate(sl):
            # get ts for date and time
            if "!!!!" in line:
                datestr = sl[  i
                             + d22_dict['datetime']['date']['idx']
                            ].strip()
                timestr = sl[  i
                             + d22_dict['datetime']['time']['idx']
                            ].strip()
                date_object = datetime.strptime(datestr
                                                + ' '
                                                + timestr,
                                                '%d-%m-%Y %H:%M')
                dt.append(date_object)
            # get ts for variable of interest
            revised_category_for_sensor = revised_categories[check]
            #print(revised_category_for_sensor)
            if revised_category_for_sensor in line:
                value = sl[  i
                            + d22_dict[category][varalias]['idx']
                           ].strip()
                ts.append(floater(value))
    else:
        print('Caution: Sensor is not defined or available')
    #Convert data to arrays
    dt = np.array(dt)
    ts = np.array(ts)
    # adjust conventions
    if ('convention' in d22_dict[category][varalias].keys() and
    d22_dict[category][varalias]['convention'] == 'oceanographic'):
        print('Convert from oceanographic to meteorologic convention')
        ts = convert_meteorologic_oceanographic(ts)
    return ts, dt

def insitu_reader(**kwargs):
    '''
    wrapping function to read insitu files

    return:
        vardict - dictionary of variables for insitu data
    '''
    dispatch_reader = {
                'd22':get_d22_dict,
                'frost':get_frost_dict,
                'nc':get_nc_dict,
                'thredds':get_nc_dict,
                }
    product = kwargs.get('fifo') # change to reader
    vardict = dispatch_reader[product](**kwargs)
    return vardict
