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
import logging
# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import get_varlst_from_nc_1D
from wavy.ncmod import get_filevarname
from wavy.ncmod import read_netcdfs
from wavy.utils import collocate_times
from wavy.utils import get_pathtofile
from wavy.utils import convert_meteorologic_oceanographic
from wavy.utils import make_subdict
from wavy.utils import parse_date
from wavy.utils import flatten
from wavy.utils import find_direction_convention
from wavy.utils import build_xr_ds
from wavy.utils import date_dispatcher
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#
# read yaml config files:
insitu_dict = load_or_default('insitu_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')
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
    r = requests.get(endpoint, parameters, auth=(client_id, client_id))
    # print(r.status_code, r.text)
    return r


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
    # empty sensor id lst
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


def get_frost(**kwargs):
    sdate = kwargs.get('sd')
    edate = kwargs.get('ed')
    nID = kwargs.get('nID')
    varalias = kwargs.get('varalias')
    varstr = [variables_frost[varalias]['frost_name']]
    sensor = insitu_dict[nID]['name'][kwargs.get('name', 0)]
    r = call_frost_api(sdate, edate, nID, varstr, sensor)
    df, _, lon, lat = get_frost_df_v1(r)
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


def get_nc_thredds(**kwargs):
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    cfg = vars(kwargs['cfg'])

    # determine ncvarname
    meta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias, variable_def,
                            cfg, meta)
    lonstr = get_filevarname('lons', variable_def,
                             cfg, meta)
    latstr = get_filevarname('lats', variable_def,
                             cfg, meta)
    timestr = get_filevarname('time', variable_def,
                              cfg, meta)

    # read all paths
    ds = read_netcdfs(pathlst, dim=timestr)
    ds_sort = ds.sortby(timestr)
    ds_sliced = ds_sort.sel({timestr: slice(sd, ed)})
    var_sliced = ds_sliced[[ncvar, lonstr, latstr]]
    return var_sliced


def get_nc_thredds_static_coords(**kwargs):
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    cfg = vars(kwargs['cfg'])

    # determine ncvarname
    meta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias, variable_def,
                            cfg, meta)
    lonstr = get_filevarname('lons', variable_def,
                             cfg, meta)
    latstr = get_filevarname('lats', variable_def,
                             cfg, meta)
    timestr = get_filevarname('time', variable_def,
                              cfg, meta)

    # read all paths
    ds = read_netcdfs(pathlst, dim=timestr)
    ds_sort = ds.sortby(timestr)

    ds_sliced = ds_sort.sel({timestr: slice(sd, ed)})

    # if lonstr is None try static from cfg
    if (lonstr is None and latstr is None):
        lons = np.ones(ds_sliced[timestr].shape)\
                *cfg['misc']['coords'][kwargs.get('name')]['lon']
        lats = np.ones(ds_sliced[timestr].shape)\
                *cfg['misc']['coords'][kwargs.get('name')]['lat']

    var_sliced = ds_sliced[[ncvar]]

    # combine and create dataset
    ds_combined = xr.Dataset({
            ncvar: xr.DataArray(data=var_sliced[ncvar].data,
                                dims=[timestr],
                                coords={timestr: var_sliced.time.data}),
            "lons": xr.DataArray(data=lons, dims=[timestr],
                                 coords={timestr: var_sliced.time.data}),
            "lats": xr.DataArray(data=lats, dims=[timestr],
                                 coords={timestr: var_sliced.time.data})
            })

    return ds_combined


def get_nc_thredds_static_coords_single_file(**kwargs):
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    cfg = vars(kwargs['cfg'])

    # determine ncvarname
    meta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias, variable_def,
                            cfg, meta)
    lonstr = get_filevarname('lons', variable_def,
                             cfg, meta)
    latstr = get_filevarname('lats', variable_def,
                             cfg, meta)
    timestr = get_filevarname('time', variable_def,
                              cfg, meta)

    # read all paths
    ds = xr.open_dataset(pathlst[0])
    ds_sort = ds.sortby(timestr)

    ds_sliced = ds_sort.sel({timestr: slice(sd, ed)})

    # if lonstr is None try static from cfg
    if (lonstr is None and latstr is None):
        lons = np.ones(ds_sliced[timestr].shape)\
                *cfg['misc']['coords'][kwargs.get('name')]['lon']
        lats = np.ones(ds_sliced[timestr].shape)\
                *cfg['misc']['coords'][kwargs.get('name')]['lat']

    var_sliced = ds_sliced[[ncvar]]

    # combine and create dataset
    ds_combined = xr.Dataset({
            ncvar: xr.DataArray(data=var_sliced[ncvar].data,
                                dims=[timestr],
                                coords={timestr: var_sliced.time.data}),
            "lons": xr.DataArray(data=lons, dims=[timestr],
                                 coords={timestr: var_sliced.time.data}),
            "lats": xr.DataArray(data=lats, dims=[timestr],
                                 coords={timestr: var_sliced.time.data})
            })

    return ds_combined


def get_cmems(**kwargs):
    sd = kwargs.get('sd')
    ed = kwargs.get('ed')
    varalias = kwargs.get('varalias')
    pathlst = kwargs.get('pathlst')
    cfg = vars(kwargs['cfg'])
    # check if dimensions are fixed
    fixed_dim_str = list(cfg['misc']['fixed_dim'].keys())[0]
    fixed_dim_idx = cfg['misc']['fixed_dim'][fixed_dim_str]

    # determine ncvarname
    meta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias, variable_def,
                            cfg, meta)
    lonstr = get_filevarname('lons', variable_def,
                             cfg, meta)
    latstr = get_filevarname('lats', variable_def,
                             cfg, meta)
    timestr = get_filevarname('time', variable_def,
                              cfg, meta)

    var_list = [ncvar, lonstr, latstr]

    ds_list = []

    # build a list of datasets using files that matches given dates
    for p in pathlst:
        try:
            ds = xr.open_dataset(p)    
            ds = ds[var_list]

            # builds the dictionary given as an argument to
            dict_var = {coord: ds.coords[coord].values
                        for coord in list(ds.coords) if coord 
                        in [lonstr, latstr, timestr]}

            dict_var.update({var: rebuild_split_variable(ds,
                                          fixed_dim_str, var) 
                             for var in list(ds.data_vars)})

            # build an xr.dataset with timestr as the only coordinate
            # using build_xr_ds function
            ds_list.append(build_xr_ds_cmems(dict_var, timestr))

        except Exception as e:
            logger.exception(e)
    
    ds_combined = xr.concat(ds_list, timestr,
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')

    ds_sort = ds_combined.sortby(timestr)
    ds_sliced = ds_sort.sel({timestr: slice(sd, ed)})

    return ds_sliced


def rebuild_split_variable(ds, fixed_dim_str, var):
    '''
    Gather values of a given variable, for which 
    values are split between several levels of
    a given dimension of a dataset.
    
    Args:
        ds (xarray dataset): dataset
        fixed_dim_str (string): name of the dimension
        var (string): name of the variable
    
    Returns:
        1D numpy array, returns the complete variable
        serie of values on a single dimension  
    '''
    lvl_nb = len(ds[fixed_dim_str].data)

    if lvl_nb==1:
        ts = list(ds.isel({fixed_dim_str: 0})[var].data)

    elif lvl_nb > 1:

        lvl_not_nan = []
        for i in range(lvl_nb):

            if not np.isnan(ds.isel({fixed_dim_str: i})[var].data).all():
                lvl_not_nan.append(i) 

        if len(lvl_not_nan)==1:
            ts= list(ds.isel({fixed_dim_str: lvl_not_nan[0]})[var].data)

        else:

            ts = ds.isel({fixed_dim_str: 0})[var].data
            dict_not_nan = {}
            for i in range(1, lvl_nb):

                nan_val_tmp = np.isnan(ds.isel({fixed_dim_str: i})[var].data)
                not_nan_idx = [j for j in range(len(nan_val_tmp)) 
                               if not nan_val_tmp[j]]
                ts[not_nan_idx] = ds.isel({fixed_dim_str: i})[var].data[not_nan_idx]

    return np.array(ts, dtype='f')


def build_xr_ds_cmems(dict_var, var_name_ref):

    ds = xr.Dataset({
        var_name: xr.DataArray(
            data=dict_var[var_name],
            dims=[var_name_ref],
            coords={var_name_ref: dict_var[var_name_ref]}
            ) for var_name in dict_var.keys()},
                attrs={'title': 'wavy dataset'})

    return ds
