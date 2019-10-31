#!/usr/bin/env python
'''
    - retrieve and aggregate data from the Tennholmen waverider
    - dump to netcdf
'''
import sys
sys.path.append(r'../wavy')

import os
from stationmod import get_buoy
from buoy_specs import buoy_dict
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

from ncmod import dumptonc_ts_Tennholmen
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter
import pytz

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieve data from the waverider Tennholmen and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./collect_buoy.py -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")

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

# Tennholmen buoy
tmpdate = deepcopy(sdate)
while tmpdate < edate:
    print('#################')
    print(tmpdate)
    print('#################')
    try:
        time_s, time_dt, Hm0, Tm02, lons, lats = get_buoy(tmpdate,tmpdate)
        obs_dict = {'time_s': time_s,
                    'time_dt': time_dt,
                    'Hm0': Hm0,
                    'Tm02': Tm02,
                    'lons': lons,
                    'lats': lats,
                    }
        outpath = tmpdate.strftime('/lustre/storeB/project/fou/om/'
                                +  'waveverification/obs/buoys/%Y/%m/')
        os.system('mkdir -p ' + outpath)
        basetime=buoy_dict['Tennholmen']['basetime']
        filename_ts=tmpdate.strftime("Tennholmen_"
                                    + "%Y%m.nc")
        title_ts=('Observations from the Tennholmen waverider')
        dumptonc_ts_Tennholmen(outpath,filename_ts,title_ts,basetime,obs_dict)
    except ValueError as e:
        print(e)
        pass
    tmpdate = tmpdate + timedelta(hours=1)
