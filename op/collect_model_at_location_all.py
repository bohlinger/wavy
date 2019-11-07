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
import time

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
./collect_model_at_location_all.py -sd 2019020100 -ed 2019020200
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-mod", metavar='model',
    help="model to be used for collocation")
parser.add_argument("-var", metavar='varname',
    help="varname")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")

args = parser.parse_args()

now = datetime.now()

init_times = np.array([0,6,12,18]).astype('float')
init_diffs = now.hour - init_times
init_diffs[init_diffs<0] = np.nan
h_idx = np.where(init_diffs==np.min(init_diffs[~np.isnan(init_diffs)]))
h = int(init_times[h_idx[0][0]])

if args.sd is None:
    #sdate = datetime(now.year,now.month,now.day,now.hour)-timedelta(hours=1)
    sdate = datetime(now.year,now.month,now.day,h)
    sdatestr = sdate.strftime("%Y%m%d%H")
else:
    sdatestr = args.sd
if args.ed is None:
    #edate = datetime(now.year,now.month,now.day,now.hour)
    edate = datetime(now.year,now.month,now.day,h) + timedelta(hours=6)
    edatestr = edate.strftime("%Y%m%d%H")
else:
    edatestr = args.ed

if args.mod is None:
    model = 'mwam4'
else:
    model = args.mod
if args.var is None:
    var = 'Hs'
else:
    var = args.var

# retrieve PID
grab_PID()
basetime = datetime(1970,1,1)

stationlst = station_dict.keys()
for station in (stationlst):
    cmd = ("python collect_model_at_location.py"
            + " -sd " + sdatestr
            + " -ed " + edatestr 
            + " -station " + station 
            + " -mod " + model
            + " -var " + var)
    tmp=os.system(cmd)
    del tmp
