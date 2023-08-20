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
from dateutil.relativedelta import relativedelta
from copy import deepcopy
import pylab as pl
import pandas as pd
import netCDF4
import requests
import dotenv
import xarray as xr

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
from wavy.utils import find_direction_convention
from wavy.utils import build_xr_ds
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#
# read yaml config files:
insitu_dict = load_or_default('insitu_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')
variables_frost = load_or_default('variables_frost.yaml')
# ---------------------------------------------------------------------#

def get_typeid(insitu_dict: dict, s: str) -> str:
    typeid = insitu_dict[s].get('typeids', 22)
    return typeid

def make_frost_reference_time_period(sdate, edate):
    sdate = parse_date(sdate)
    edate = parse_date(edate)
    formatstr = '%Y-%m-%dT%H:%M:00.000Z'
    refstr = '{}/{}'.format(sdate.strftime(formatstr),
                            edate.strftime(formatstr))
    return refstr

def call_frost_api(
    sdate: datetime, edate: datetime,
    nID: str, varstr: str, sensor: str) -> 'requests.models.Response':
    """
    make frost api call
    """
    print('Make frost api call ...')
    dotenv.load_dotenv()
    client_id = os.getenv('CLIENT_ID', None)
    frost_reference_time = make_frost_reference_time_period(sdate, edate)
    if client_id is None:
        print("No Frost CLIENT_ID given!")
    sensor = insitu_dict[nID]['name'][sensor]
    r = call_frost_api_v1(nID, varstr,
                          frost_reference_time,
                          client_id, sensor)
    print(r.url)
    print('\nr.status_code:', r.status_code, '\n')
    if r.status_code != 200:
        print('Error! Returned status code %s' % r.status_code)
        error = r.json()['error']
        for part in ['message', 'reason', 'help']:
            if part in error:
                print(part.upper(), ': ', error[part])
    else:
        return r

def call_frost_api_v1(
        nID: str, varstr: str, frost_reference_time: str,
        client_id: str, sensor: str)\
    -> 'requests.models.Response':
    """
    frost call, retrieve data from frost v1
    """
    ID = insitu_dict[nID]['misc']['ID']
    endpoint = 'https://frost-beta.met.no/api/v1/obs/met.no/kvkafka/get?'
    parameters = {
                'stationids': ID,
                'elementids': varstr,
                'time': frost_reference_time,
                'levels': 'all',
                'incobs': 'true',
                # 'sensors': '0,1,2,3,4,5',
                'sensors': sensor,  # limit to one sensor
                'typeids': str(get_typeid(insitu_dict, nID))
                }
    print('parameters forst api call: ', parameters)
    return requests.get(endpoint, parameters, auth=(client_id, client_id))

def find_preferred(idx, sensors, refs, pref):
    sensorsU = np.unique(sensors)
    preferred_idx = []
    for s in sensorsU:
        no = len(refs[sensors == s])
        idx_1 = idx[sensors == s]
        if no > 1:
            idx_2 = np.where(refs[sensors == s] == pref)
            idx_3 = idx_1[idx_2]
            preferred_idx.append(list(idx_3)[0])
        else:
            preferred_idx.append(list(idx_1)[0])
    return preferred_idx

def get_frost_df_v1(r: 'requests.models.Response')\
    -> 'pandas.core.frame.DataFrame':
    """
    create pandas dataframe from frost call for v1
    """
    # empy sensor id lst
    # base df
    df = pd.json_normalize(r.json()['data']['tseries'])
    # coordinates for static station (sensor 0)
    lon = float(df['header.extra.station.location'][0][0]['value']['longitude'])
    lat = float(df['header.extra.station.location'][0][0]['value']['latitude'])
    # df to be concatenated initialized with time
    # select time index, some ts have less than others
    # choose the one with most values
    no_of_ts = len(pd.json_normalize(r.json()['data']['tseries'][:]))
    no_of_ts = min(4, no_of_ts)
    lenlst = []
    for t in range(no_of_ts):
        lenlst.append(len(pd.json_normalize(r.json()
                      ['data']['tseries'][t]['observations'])['time']
                      .to_frame()))
    time_idx = lenlst.index(max(lenlst))
    dfc = pd.json_normalize(r.json()
            ['data']['tseries'][time_idx]['observations'])['time'].to_frame()
    dinfo = {'sensor': {}, 'level': {}, 'parameterid': {},
             'geometric height': {}, 'masl': {}}
    for vn in variables_frost:
        frostvar = variables_frost[vn]['frost_name']
        idx = np.array(df['header.extra.element.id']
                [df['header.extra.element.id'] == frostvar].index.to_list())
        sensors = df['header.id.sensor'][idx].values
        parameterids = df['header.id.parameterid'][idx].values
        levels = df['header.id.level'][idx].values
        if len(sensors) != len(np.unique(sensors)):
            print("-> id.sensor was not unique " 
                    + "selecting according to variable_def.yaml")
            print("   affected variable: ", frostvar)
            # 1. prioritize according to parameterid
            if len(np.unique(parameterids)) > 1:
                print('multiple parameterids (',
                        len(np.unique(parameterids)), ')')
                print('parameterids:', np.unique(parameterids))
                idx = find_preferred(
                        idx, sensors, parameterids,
                        variables_frost[vn]['prime_parameterid'])
                sensors = df['header.id.sensor'][idx].values
                parameterids = df['header.id.parameterid'][idx].values
                levels = df['header.id.level'][idx].values
            # 2. prioritize according to level
            if len(np.unique(levels)) > 1:
                print('multiple levels (', len(np.unique(levels)), ')')
                print('unique(levels):', np.unique(levels))
                idx = find_preferred(
                        idx, sensors, levels,
                        variables_frost[vn]['prime_level'])
                sensors = df['header.id.sensor'][idx].values
                parameterids = df['header.id.parameterid'][idx].values
                levels = df['header.id.level'][idx].values
        for n, i in enumerate(idx):
            dftmp = pd.json_normalize(r.json()\
                        ['data']['tseries'][i]['observations'])\
                        ['body.value'].to_frame()
            vns = vn
            #vns = vn + '_'a \
            #            + str(df['header.id.sensor'][i])
            dftmp = dftmp.rename(columns={ dftmp.columns[0]: vns }).\
                            astype(float)
            dftmp[vns] = dftmp[vns].mask(dftmp[vns] < 0, np.nan)
            dfc = pd.concat([dfc, dftmp.reindex(dfc.index)], axis=1)
            # sensor
            dinfo['sensor'][vns] = sensors[n]
            # level
            if levels[n] == 0:
                dinfo['level'][vns] = variables_frost[vn]['default_level']
            else:
                dinfo['level'][vns] = levels[n]
            # parameterid
            dinfo['parameterid'][vns] = parameterids[n]
    return dfc, dinfo, lon, lat


def get_frost_dict(**kwargs):
    sdate = kwargs.get('sd')
    edate = kwargs.get('ed')
    nID = kwargs.get('nID')
    varalias = kwargs.get('varalias')
    varstr = [variables_frost[varalias]['frost_name']]
    sensor = kwargs.get('sensor', 0)
    r = call_frost_api(sdate, edate, nID, varstr, sensor)
    df, dinfo, lon, lat = get_frost_df_v1(r)
    var = df[varalias].values
    timevec = df['time'].values
    timedt = [parse_date(str(d)) for d in timevec]

    # rm datetime timezone info
    timedt = [d.replace(tzinfo=None) for d in timedt]
    lons = len(var)*[lon]
    lats = len(var)*[lat]
    varnames = (varalias, 'lons', 'lats', 'time')
    var_tuple = (var, lons, lats, timedt)

    # build xarray ds
    ds = build_xr_ds(var_tuple, varnames)
    return ds

def get_thredds_dict(**kwargs):
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    nID = kwargs.get('nID')
    #sensor = kwargs.get('sensor')
    varalias = kwargs.get('varalias')
    src_tmplt = kwargs.get('cfg').wavy_input['src_tmplt']
    print('HERE')
    ds = None
    return ds

def get_nc_dict(**kwargs):
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    nID = kwargs.get('nID')
    #sensor = kwargs.get('sensor')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    print(pathlst)
    strsublst = kwargs.get('strsublst')
    print(strsublst)
    dict_for_sub = kwargs.get('dict_for_sub')
    print(dict_for_sub)
    # loop from sdate to edate with dateincr
    tmpdate = deepcopy(sd)
    varlst = []
    lonlst = []
    latlst = []
    timelst = []
    dtimelst = []
    # make subdict
    subdict = make_subdict(strsublst, class_object_dict=dict_for_sub)
    while datetime(tmpdate.year, tmpdate.month, 1)\
    <= datetime(ed.year, ed.month, 1):
        # get pathtofile
        pathtofile = get_pathtofile(pathlst, strsublst,
                                    subdict, tmpdate)
        # get ncdump
        ncdict = ncdumpMeta(pathtofile)
        # retrieve filevarname for varalias
        filevarname = get_filevarname(varalias,
                                      variable_info,
                                      insitu_dict[nID],
                                      ncdict)
        varstrlst = [filevarname, 'longitude', 'latitude', 'time']
        # query
        vardict = get_varlst_from_nc_1D(pathtofile,
                                        varstrlst,
                                        sd, ed)
        varlst.append(list(vardict[filevarname]))
        lonlst.append(list(vardict['longitude']))
        latlst.append(list(vardict['latitude']))
        timelst.append(list(vardict['time']))
        dtimelst.append(list(vardict['dtime']))
        # determine date increment
        file_date_incr = insitu_dict[nID]['src'].get('file_date_incr', 'm')
        if file_date_incr == 'm':
            tmpdate += relativedelta(months=+1)
        elif file_date_incr == 'Y':
            tmpdate += relativedelta(years=+1)
        elif file_date_incr == 'd':
            tmpdate += timedelta(days=+1)
    varlst = flatten(varlst)
    lonlst = flatten(lonlst)
    latlst = flatten(latlst)
    timelst = flatten(timelst)
    dtimelst = flatten(dtimelst)

    varnames = (varalias, 'lons', 'lats', 'time')
    var_tuple = (varlst, lonlst, latlst, dtimelst)

    # build xarray ds
    ds = build_xr_ds(var_tuple, varnames)
    return ds

def insitu_reader(**kwargs):
    '''
    wrapping function to read insitu files

    return:
        vardict - dictionary of variables for insitu data
    '''
    dispatch_reader = {
                'frost': get_frost_dict,
                'nc': get_nc_dict,
                'thredds': get_nc_dict,
                }
    product = kwargs.get('fifo')  # change to reader
    vardict = dispatch_reader[product](**kwargs)
    return vardict
