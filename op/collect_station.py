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
from datetime import datetime, timedelta
from ncmod import dumptonc_ts_station
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
./collect_station.py -sd 2019010100 -ed 2019020200 -station ekofiskL -sensor waverider
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
    sdate = datetime(now.year,now.month,now.day,now.hour)-timedelta(hours=1)
else:
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
if args.ed is None:
    edate = datetime(now.year,now.month,now.day,now.hour)
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))

# retrieve PID
grab_PID()
statname = args.station
sensorname = args.sensor
mode = 'd22'
deltat = 10

tmpdate = deepcopy(sdate)
while tmpdate < edate:
    tmpedate = deepcopy(tmpdate + timedelta(minutes=deltat))
    print('#################')
    print(tmpdate)
    print('#################')
    try:
        sc_obj = sc(statname,tmpdate,tmpedate,mode=mode,sensorname=sensorname,deltat=deltat)
        stat_dict = {'basetime':sc_obj.basedate,
             'time':sc_obj.time[0:-1],
             'Hs_stat_10min':sc_obj.hs_obs[0:-1],
             'lats_stat':sc_obj.lat,
             'lons_stat':sc_obj.lon,
            }

        # dump tp nc-file
        outpath = tmpdate.strftime('/lustre/storeB/project/fou/om/'
                            + 'waveverification/obs/stations/'
                            + '%Y/%m/')
        os.system('mkdir -p ' + outpath)
        filename_ts=tmpdate.strftime(statname + "_" + sensorname + "_%Y%m.nc")
        title_ts=('Observations from ' + statname)
        dumptonc_ts_station(outpath,filename_ts,title_ts,sc_obj.basedate,stat_dict,statname,sensorname)
    except ValueError as e:
        print(e)
        pass
    tmpdate = tmpdate + timedelta(minutes=deltat)
