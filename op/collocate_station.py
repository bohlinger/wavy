#!/usr/bin/env python
'''
    - retrieve data from station
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf
!!!

'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from stationmod import station_class as sc
from stationmod import parse_d22, extract_d22, matchtime, get_loc_idx
from modelmod import get_model, check_date
from collocmod import collocate
from station_specs import station_dict
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np

from ncmod import dumptonc_coll_ts_station
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
./collocate_station.py -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-station", metavar='statname',
    help="stationname")
parser.add_argument("-sensor", metavar='sensorname',
    help="sensorname")

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

# retrieve PID
grab_PID()
#sdate=datetime(2019,1,1)
#edate=datetime(2019,1,2) # only entire days
#edate=datetime.now()
model='mwam4'
station = args.station
sensorname = args.sensor
title_ts = (model + ' vs ' + station + ' (' + sensorname + ')')
basetime = datetime(1970,1,1)

tmpdate = deepcopy(sdate)
element = 0
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
    print(tmpdate)
    if np.isnan(idx):
        try:
            sc_obj10 = sc(station,tmpdate,tmpdate+timedelta(minutes=10),
                            mode='d22',sensorname=sensorname,deltat=10)
            hs_obs_10min = sc_obj10.hs_obs[0:-1]
            sc_obj60 = sc(station,tmpdate,tmpdate+timedelta(minutes=60),
                            mode='d22',sensorname=sensorname,deltat=60)
            hs_obs_hourly = sc_obj60.hs[0:-1]
            # get wave model
            check_date(model,fc_date=tmpdate,leadtime=element)
            model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=tmpdate,
                init_date=tmpdate,leadtime=element)
            # collocate with wave model
            idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
                            model_lats,model_lons,\
                            station_dict[station]['coords']['lat'],\
                            station_dict[station]['coords']['lon'],\
                            mask=None)
            # append values
            idx_lst.append(idx)
            idy_lst.append(idy)
            dist_lst.append(distM[idx,idy])
            mod_lst.append(model_Hs.squeeze()[idx,idy][0])
            mod_lat_lst.append(picked_lat[0])
            mod_lon_lst.append(picked_lon[0])
            stat_10min_lst.append(hs_obs_10min[0])
            stat_hourly_lst.append(hs_obs_hourly[0])
            stat_lat_lst.append(station_dict[station]['coords']['lat'])
            stat_lon_lst.append(station_dict[station]['coords']['lon'])
            time_dt_lst.append(tmpdate)
            time_s = (tmpdate - basetime).total_seconds()
            time_s_lst.append(time_s)
            # dump tp nc-file
            coll_dict = {'basetime':basetime,
                 'time':[time_s_lst[-1]],
                 'Hm0_model':[mod_lst[-1]],
                 'lats_model':[mod_lat_lst[-1]],
                 'lons_model':[mod_lon_lst[-1]],
                 'Hs_stat_1h':[stat_hourly_lst[-1]],
                 'Hs_stat_10min':[stat_10min_lst[-1]],
                 'lats_stat':[stat_lat_lst[-1]],
                 'lons_stat':[stat_lon_lst[-1]],
                 'hdist':[dist_lst[-1]],
                 'idx':[idx_lst[-1]],
                 'idy':[idy_lst[-1]],
                }
            print(coll_dict)
            outpath = tmpdate.strftime('/lustre/storeB/project/fou/om/'
                            + 'waveverification/mwam4/stations/'
                            + 'CollocationFiles/'
                            + '%Y/%m/')
            os.system('mkdir -p ' + outpath)
            filename_ts=tmpdate.strftime(model 
                                    + "_" 
                                    + station 
                                    + "_" 
                                    + sensorname 
                                    + "_coll_ts_lt000h_%Y%m.nc")
            title_ts=(model + ' vs ' + station)
            dumptonc_coll_ts_station(outpath,
                                filename_ts,
                                title_ts,
                                basetime,
                                coll_dict,
                                model,
                                station,sensorname)
        except SystemExit:
            print('error: --> leadtime is not available')
        except (ValueError,IOError,KeyError) as e:
            print(e)
        except (IndexError) as e:
            print('no measurement available')
            print(e)
    else:
        try:
            sc_obj10 = sc(station,tmpdate,tmpdate+timedelta(minutes=10),
                            mode='d22',sensorname=sensorname,deltat=10)
            hs_obs_10min = sc_obj10.hs_obs[0:-1]
            sc_obj60 = sc(station,tmpdate,tmpdate+timedelta(minutes=60),
                            mode='d22',sensorname=sensorname,deltat=60)
            hs_obs_hourly = sc_obj60.hs[0:-1]
            # get wave model
            check_date(model,fc_date=tmpdate,leadtime=element)
            model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=tmpdate,
                init_date=tmpdate,leadtime=element)
            # append values
            idx_lst.append(idx)
            idy_lst.append(idy)
            dist_lst.append(distM[idx,idy])
            mod_lst.append(model_Hs.squeeze()[idx,idy][0])
            mod_lat_lst.append(picked_lat[0])
            mod_lon_lst.append(picked_lon[0])
            stat_10min_lst.append(hs_obs_10min[0])
            stat_hourly_lst.append(hs_obs_hourly[0])
            stat_lat_lst.append(station_dict[station]['coords']['lat'])
            stat_lon_lst.append(station_dict[station]['coords']['lon'])
            time_dt_lst.append(tmpdate)
            time_s = (tmpdate - basetime).total_seconds()
            time_s_lst.append(time_s)
            # dump tp nc-file
            coll_dict = {'basetime':basetime,
                 'time':[time_s_lst[-1]],
                 'Hm0_model':[mod_lst[-1]],
                 'lats_model':[mod_lat_lst[-1]],
                 'lons_model':[mod_lon_lst[-1]],
                 'Hs_stat_1h':[stat_hourly_lst[-1]],
                 'Hs_stat_10min':[stat_10min_lst[-1]],
                 'lats_stat':[stat_lat_lst[-1]],
                 'lons_stat':[stat_lon_lst[-1]],
                 'hdist':[dist_lst[-1]],
                 'idx':[idx_lst[-1]],
                 'idy':[idy_lst[-1]],
                }
            print(coll_dict)
            outpath = tmpdate.strftime('/lustre/storeB/project/fou/om/'
                            + 'waveverification/mwam4/stations/'
                            + 'CollocationFiles/'
                            + '%Y/%m/')
            os.system('mkdir -p ' + outpath)
            filename_ts=tmpdate.strftime(model
                                    + "_"
                                    + station
                                    + "_"
                                    + sensorname
                                    + "_coll_ts_lt000h_%Y%m.nc")
            title_ts=(model + ' vs ' + station)
            dumptonc_coll_ts_station(outpath,
                                filename_ts,
                                title_ts,
                                basetime,
                                coll_dict,
                                model,
                                station,sensorname)
        except SystemExit:
            print('error: --> leadtime is not available')
        except (ValueError,IOError,KeyError) as e:
            print(e)
        except (IndexError) as e:
            print('no measurement available')
            print(e)
    tmpdate = tmpdate + timedelta(hours=6)
