#!/usr/bin/env python
'''
    - retrieve data from the Tennholmen waverider
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf

!!! import to be implemented:
    - choose freely sdate and edate, currently only one month is possible
    - add model type to parsing options
    - add leadtime to parsing options
!!!

'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from stationmod import get_buoy, get_loc_idx
from modelmod import get_model, check_date
from collocmod import collocate
from buoy_specs import buoy_dict
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

from ncmod import dumptonc_coll_ts_Tennholmen
from ncmod import dumptonc_coll_ts_buoy
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter

import numpy as np

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieves data from the waverider Tennholmen and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./collocate_buoy.py -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-buoy", metavar='buoyname',
    help="name of buoy")

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

if args.buoy is None:
    buoy = 'Tennholmen'
else:
    buoy=args.buoy


# retrieve PID
grab_PID()
model='mwam4'
basetime = buoy_dict[buoy]['basetime']

tmpdate = deepcopy(sdate)
element = 0
buoy_lst = []
mod_lst = []
dist_lst = []
mod_lon_lst = []
mod_lat_lst = []
buoy_lat_lst = []
buoy_lon_lst = []
time_dt_lst = []
time_s_lst = []

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
        for i in range(len(leadtimes)):
            try:
                element = leadtimes[i]
                print('leadtimes: ', element)
                fc_date = ( tmpdate
                            - timedelta(hours=tdeltas[i])
                            )
                print('fc_date: ', fc_date)
                init_date = tmpdate-timedelta(hours=init_step)
                print('init_date: ', init_date)
                # get wave rider data
                time_s, time_dt, Hm0, Tm02, lons, lats = \
                            get_buoy(fc_date,fc_date,buoyname=buoy)
                # get wave model
                check_date(model,fc_date=fc_date,leadtime=element)
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc",model=model,fc_date=fc_date,
                    init_date=init_date,leadtime=element)
                # collocate with wave model
                idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
                                        model_lats,model_lons,\
                                        lats,lons,\
                                        mask=None)
                # append values
                mod_lst.append(model_Hs.squeeze()[idx,idy][0])
                mod_lat_lst.append(picked_lat[0])
                mod_lon_lst.append(picked_lon[0])
                buoy_lst.append(Hm0[0])
                buoy_lat_lst.append(lats[0])
                buoy_lon_lst.append(lons[0])
                dist_lst.append(distM[idx,idy][0])
                time_dt_lst.append(time_dt[0])
                time_s_lst.append(time_s[0])
            except SystemExit:
                print('error: --> leadtime is not available')
            except (ValueError,IOError,KeyError,IndexError) as e:
                print(e)
    else:
        for i in range(len(leadtimes)):
            try:
                element = leadtimes[i]
                print('leadtimes: ', element)
                fc_date = ( tmpdate
                            - timedelta(hours=tdeltas[i])
                            )
                print('fc_date: ', fc_date)
                init_date = tmpdate-timedelta(hours=init_step)
                print('init_date: ', init_date)
                # get wave rider data
                time_s, time_dt, Hm0, Tm02, lons, lats = \
                            get_buoy(fc_date,fc_date,buoyname=buoy)
                # get wave model
                check_date(model,fc_date=fc_date,leadtime=element)
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc",model=model,fc_date=fc_date,
                    init_date=init_date,leadtime=element)
                # collocate with wave model
                idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
                                        model_lats,model_lons,\
                                        lats,lons,\
                                        mask=None)
                # append values
                mod_lst.append(model_Hs.squeeze()[idx,idy][0])
                mod_lat_lst.append(picked_lat[0])
                mod_lon_lst.append(picked_lon[0])
                buoy_lst.append(Hm0[0])
                buoy_lat_lst.append(lats[0])
                buoy_lon_lst.append(lons[0])
                dist_lst.append(distM[idx,idy][0])
                time_dt_lst.append(time_dt[0])
                time_s_lst.append(time_s[0])
            except SystemExit:
                print('error: --> leadtime is not available')
            except (ValueError,IOError,KeyError,IndexError) as e:
                print(e)
    tmpdate = tmpdate + timedelta(hours=init_step)

# collocation dictionary
coll_dict = {'basetime':basetime,
             'time':time_s_lst,
             'Hm0_model':mod_lst,
             'lats_model':mod_lat_lst,
             'lons_model':mod_lon_lst,
             'Hm0_buoy':buoy_lst,
             'lats_buoy':buoy_lat_lst,
             'lons_buoy':buoy_lon_lst,
             'buoy':buoy
            }
# dump tp nc-file
outpath = fc_date.strftime('/lustre/storeB/project/fou/om/'
                            + 'waveverification/mwam4/buoys/'
                            + 'CollocationFiles/'
                            + '%Y/%m/')
os.system('mkdir -p ' + outpath)
filename_ts=fc_date.strftime(model + "_" + buoy + "_coll_ts_%Y%m_bestguess.nc")
title_ts=(model + ' vs ' + buoy)
#if buoy == 'Tennholmen':
#    dumptonc_coll_ts_Tennholmen(outpath,filename_ts,title_ts,
#                            basetime,coll_dict,model)
#else:
#    dumptonc_coll_ts_buoy(outpath,filename_ts,title_ts,
#                            basetime,coll_dict,model)
dumptonc_coll_ts_buoy(outpath,filename_ts,title_ts,
                        basetime,coll_dict,model)
