#!/usr/bin/env python
'''
    - retrieve data from station
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf
!!!

'''
import sys
import os
#from stationmod import get_buoy, get_loc_idx
from stationmod import parse_d22, extract_d22, matchtime, get_loc_idx
from modelmod import get_model, check_date
from collocmod import collocate
#from buoy_specs import buoy_dict
from station_specs import station_dict
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np

#from ncmod import dumptonc_coll_ts_station
from ncmod import dumptonc_ts_pos
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieves data from a station and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./get_bestEstimate.py -sd 2018110112 -ed 2018110118 -m mwam4 -lat latitude -lon longitude -o outpath/
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-m", metavar='model',
    help="model")
parser.add_argument("-lat", metavar='latitude', type=float,
    help="latitude")
parser.add_argument("-lon", metavar='longitude', type=float,
    help="longitude")
parser.add_argument("-o", metavar='outpath',
    help="outpath")

args = parser.parse_args()

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

edate = edate + timedelta(hours=1)

# retrieve PID
grab_PID()
model = args.m
lat = args.lat
lon = args.lon
location = str('lat: ' + str(lat) + ', lon: ' + str(lon))
sensorname = '?'
basetime = datetime(1970,1,1)

stat_hourly_lst = []
stat_10min_lst = []
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
if model == 'mwam4':
    init_step = 6
    init_start = np.min([0,6,12,18])
if model == 'mwam8':
    init_step = 12
    init_start = 12
tdeltas = range(1,init_step+1)[::-1]
leadtimes = range(init_step)

tmpdate = deepcopy(sdate) + timedelta(hours=init_step)
idx, idy = np.nan, np.nan
while tmpdate <= edate:
    print('tmpdate: ', tmpdate)
    if np.isnan(idx):
        for i in range(len(leadtimes)):
            try:
                element = leadtimes[i]
                print('leadtimes: ', element)
                fc_date = ( tmpdate 
                            - timedelta(hours=tdeltas[i]) 
#                            + timedelta(hours=init_step)
                            )
                print('fc_date: ', fc_date)
                print('init_date: ', tmpdate)
                # get wave model
                check_date(model,fc_date=fc_date,leadtime=element)
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc",model=model,fc_date=fc_date,
                    init_date=tmpdate,leadtime=element)
                # collocate with wave model
                idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
                            model_lats,model_lons,\
                            lat,\
                            lon,\
                            mask=None)
                # append values
                idx_lst.append(idx)
                idy_lst.append(idy)
                dist_lst.append(distM[idx,idy])
                mod_lst.append(model_Hs.squeeze()[idx,idy][0])
                mod_lat_lst.append(picked_lat[0])
                mod_lon_lst.append(picked_lon[0])
                time_dt_lst.append(fc_date)
                time_s = (fc_date - basetime).total_seconds()
                time_s_lst.append(time_s)
            except SystemExit:
                print('error: --> leadtime is not available')
            except (ValueError,IOError,KeyError) as e:
                print(e)
    else:
        for i in range(len(leadtimes)):
            try:
                element = leadtimes[i]
                print('leadtimes: ', element)
                fc_date = tmpdate - timedelta(hours=tdeltas[i])
                print('fc_date: ', fc_date)
                # get wave model
                check_date(model,fc_date=fc_date,leadtime=element)
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc",model=model,fc_date=fc_date,
                    init_date=tmpdate,leadtime=element)
                # append values
                idx_lst.append(idx)
                idy_lst.append(idy)
                dist_lst.append(distM[idx,idy])
                mod_lst.append(model_Hs.squeeze()[idx,idy][0])
                mod_lat_lst.append(picked_lat[0])
                mod_lon_lst.append(picked_lon[0])
                time_dt_lst.append(fc_date)
                time_s = (fc_date - basetime).total_seconds()
                time_s_lst.append(time_s)
            except SystemExit:
                print('error: --> leadtime is not available')
            except (ValueError,IOError,KeyError) as e:
                print(e)
    tmpdate = tmpdate + timedelta(hours=init_step)


coll_dict = {'basetime':basetime,
             'time':time_s_lst,
             'Hm0_model':mod_lst,
             'lats_model':mod_lat_lst,
             'lons_model':mod_lon_lst,
             'Hs_stat_1h':stat_hourly_lst,
             'Hs_stat_10min':stat_10min_lst,
             'lats_stat':[lat],
             'lons_stat':[lon],
             'hdist':dist_lst,
             'idx':idx_lst,
             'idy':idy_lst,
            }

# dump tp nc-file
outpath = args.o
locationstr = ('lat_' + str(lat) + '_' + 'lon_' + str(lon))
os.system('mkdir -p ' + outpath)
filename_ts=tmpdate.strftime(model 
                            + "_" 
                            + locationstr 
                            + "_" + args.sd + "_" + args.ed 
                            + ".nc")
title_ts=('Hs from ' + model + ' for ' + location)
dumptonc_ts_pos(outpath,filename_ts,title_ts,basetime,coll_dict,model,location,sensorname)
