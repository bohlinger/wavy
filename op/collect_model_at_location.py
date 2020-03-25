#!/usr/bin/env python3
'''
    - retrieve data from station
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf
!!!

'''
import os
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
#sys.path.append(os.getenv("HOME") + "/met-ecflow-support/lib/python")
#import python
#python.module('load', 'compiler/intelPE2018')
#python.module('load', 'hdf5/1.10.5-intel2018')
#python.module('load', 'netcdf/4.7.0-intel2018')

from modelmod import get_model, check_date
from datetime import datetime, timedelta
#from station_specs import station_dict
from ncmod import dumptonc_ts_pos,check_vals_in_nc
import numpy as np
from stationmod import matchtime, get_loc_idx
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter
import yaml

with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieves data from a station and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./collect_model_at_location.py -sd 2019010100 -ed 2019020200 -station ekofiskL -mod mwam4 -var Hs
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-station", metavar='statname',
    help="stationname")
parser.add_argument("-mod", metavar='modelname',
    help="modelname")
parser.add_argument("-var", metavar='varname',
    help="varname")

args = parser.parse_args()

if args.mod is None:
    args.mod = 'mwam4'

now = datetime.now()

if args.sd is None:
    sdate = datetime(now.year,now.month,now.day)-timedelta(days=1)
else:
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
if args.ed is None:
    edate = datetime(now.year,now.month,now.day)-timedelta(hours=1)
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))

# retrieve PID
grab_PID()

# settings
station = args.station
model = args.mod
varname = args.var
basetime = datetime(1970,1,1)

if model == 'mwam4':
    init_step = 6
    leadtimes = [0, 6, 12, 18, 24, 36, 48, 60]
if (model == 'mwam8' or model == 'ecwam' or model == 'mwam3'):
    init_step = 12
    leadtimes = [0, 12, 24, 36, 48, 60, 72, 96, 120, 144]

tmpdate = deepcopy(sdate)
idx, idy = np.nan, np.nan

mod_lst = []
dist_lst = []
mod_lon_lst = []
mod_lat_lst = []
stat_lat_lst = []
stat_lon_lst = []
time_dt_lst = []
time_s_lst = []
idx_lst = []
idy_lst = []

print('---')
print(sdate)
print(edate)
print('---')

while tmpdate <= edate:
    print(tmpdate)
    if np.isnan(idx):
        for i in range(len(leadtimes)):
            element = leadtimes[i]
            print('leadtime: ', element)
            fc_date = tmpdate
            print('fc_date: ', fc_date)
            init_date = tmpdate - timedelta(hours=element)
            print('init_date: ', init_date)
            outpath = fc_date.strftime('/lustre/storeB/project/fou/om/'
                            + 'waveverification/' + model + '/stations/'
                            + 'CollocationFiles/'
                            + '%Y/%m/')
            os.system('mkdir -p ' + outpath)
            filename_ts=fc_date.strftime(model
                                    + "_"
                                    + varname
                                    + "_at_"
                                    + station
                                    + "_ts_lt" 
                                    + "{:0>3d}".format(element)
                                    + "h_%Y%m.nc")
            try:
                # check if values for this date already exist
                vidx = check_vals_in_nc(outpath+filename_ts,varname,fc_date)
                if vidx is None:
                    print('time does not yet exist, filling slot...')
                    # get wave model
                    check_date(model,fc_date=fc_date,leadtime=element)
                    print('Read model file')
                    model_var,model_lats,model_lons,model_time,model_time_dt = \
                        get_model(simmode="fc",model=model,fc_date=fc_date,
                        init_date=init_date,leadtime=element,varname=varname)
                    # collocate with wave model
                    print('Collocate')
                    if model=='ecwam': # inflate coords
                        model_lons, model_lats = \
                                        np.meshgrid(model_lons, model_lats)
                    idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
                            model_lats,model_lons,\
                            station_dict[station]['coords']['lat'],\
                            station_dict[station]['coords']['lon'],\
                            mask=None)
                    # append values
                    idx_lst.append(idx)
                    idy_lst.append(idy)
                    dist_lst.append(distM[idx,idy])
                    mod_lst.append(model_var.squeeze()[idx,idy][0])
                    mod_lat_lst.append(picked_lat[0])
                    mod_lon_lst.append(picked_lon[0])
                    stat_lat_lst.append(station_dict[station]['coords']['lat'])
                    stat_lon_lst.append(station_dict[station]['coords']['lon'])
                    time_dt_lst.append(fc_date)
                    time_s = (fc_date - basetime).total_seconds()
                    time_s_lst.append(time_s)
                    # dump tp nc-file
                    coll_dict = {'basetime':basetime,
                         'time':[time_s_lst[-1]],
                         varname:[mod_lst[-1]],
                         'lats_model':[mod_lat_lst[-1]],
                         'lons_model':[mod_lon_lst[-1]],
                         'lats_pos':[stat_lat_lst[-1]],
                         'lons_pos':[stat_lon_lst[-1]],
                         'hdist':[dist_lst[-1][0]],
                         'idx':[idx_lst[-1][0]],
                         'idy':[idy_lst[-1][0]],
                         'model':model,
                         'station':station,
                         'varname':varname
                        }
                    print(coll_dict)
                    title_ts=(
                        model + ' ' + varname + ' at location ' + station 
                        + ' with leadtime '
                        + "{:0>3d}".format(element)
                        + ' h')
                    dumptonc_ts_pos(outpath,
                                filename_ts,
                                title_ts,
                                coll_dict,
                                )
            except Exception as e: print(e)
    else:
        for i in range(len(leadtimes)):
            element = leadtimes[i]
            print('leadtime: ', element)
            fc_date = tmpdate
            print('fc_date: ', fc_date)
            init_date = tmpdate - timedelta(hours=element)
            print('init_date: ', init_date)
            outpath = fc_date.strftime('/lustre/storeB/project/fou/om/'
                            + 'waveverification/' + model + '/stations/'
                            + 'CollocationFiles/'
                            + '%Y/%m/')
            os.system('mkdir -p ' + outpath)
            filename_ts=fc_date.strftime(model
                                    + "_"
                                    + varname
                                    + "_at_"
                                    + station
                                    + "_ts_lt"
                                    + "{:0>3d}".format(element)
                                    + "h_%Y%m.nc")
            try:
                # check if values for this date already exist
                vidx = check_vals_in_nc(outpath+filename_ts,varname,fc_date)
                if vidx is None:
                    print('time does not yet exist, filling slot...')
                    # get wave model
                    check_date(model,fc_date=fc_date,leadtime=element)
                    print('Read model file')
                    model_var,model_lats,model_lons,model_time,model_time_dt = \
                        get_model(simmode="fc",model=model,fc_date=fc_date,
                        init_date=init_date,leadtime=element,varname=varname)
                    # append values
                    idx_lst.append(idx)
                    idy_lst.append(idy)
                    dist_lst.append(distM[idx,idy])
                    mod_lst.append(model_var.squeeze()[idx,idy][0])
                    mod_lat_lst.append(picked_lat[0])
                    mod_lon_lst.append(picked_lon[0])
                    stat_lat_lst.append(station_dict[station]['coords']['lat'])
                    stat_lon_lst.append(station_dict[station]['coords']['lon'])
                    time_dt_lst.append(fc_date)
                    time_s = (fc_date - basetime).total_seconds()
                    time_s_lst.append(time_s)
                    # dump tp nc-file
                    coll_dict = {'basetime':basetime,
                         'time':[time_s_lst[-1]],
                         varname:[mod_lst[-1]],
                         'lats_model':[mod_lat_lst[-1]],
                         'lons_model':[mod_lon_lst[-1]],
                         'lats_pos':[stat_lat_lst[-1]],
                         'lons_pos':[stat_lon_lst[-1]],
                         'hdist':[dist_lst[-1][0]],
                         'idx':[idx_lst[-1][0]],
                         'idy':[idy_lst[-1][0]],
                         'model': model,
                         'station':station,
                         'varname':varname
                        }
                    print(coll_dict)
                    title_ts=(
                        model + ' ' + varname + ' at location ' + station
                        + ' with leadtime '
                        + "{:0>3d}".format(element)
                        + ' h')
                    dumptonc_ts_pos(outpath,
                                filename_ts,
                                title_ts,
                                coll_dict,
                                )
            except Exception as e: print(e)
    tmpdate = tmpdate + timedelta(hours=init_step)
