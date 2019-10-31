#!/usr/bin/env python
'''
    - retrieve data from station
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf
!!!

'''
import sys
sys.path.append(r'../wavy')

import yaml
import os
from stationmod import station_class as sc
from stationmod import parse_d22, extract_d22, matchtime, get_loc_idx
from modelmod import get_model, check_date
from collocmod import collocate
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np
import time

from ncmod import dumptonc_coll_ts_station
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter

# read yaml config files:
with open("station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieves data from a station and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./collocate_station_all_bestguess.py -sd 2018110112 -ed 2018110118
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
    sdatestr = sdate.strftime("%Y%m%d%H")
else:
    sdatestr = args.sd
if args.ed is None:
    edate = datetime(now.year,now.month,now.day,now.hour)
    edatestr = edate.strftime("%Y%m%d%H")
else:
    edatestr = args.ed

# retrieve PID
grab_PID()
model='mwam4'
basetime = datetime(1970,1,1)

if model == 'mwam4':
    init_step = 6
    init_start = np.min([0,6,12,18])
if model == 'mwam8':
    init_step = 12
    init_start = 12
tdeltas = range(1,init_step+1)[::-1]
leadtimes = range(init_step)

stationlst = station_dict.keys()
for station in (stationlst):
    for sensor in (station_dict[station]['sensor'].keys()):
        cmd = ("python collocate_station_bestguess.py -sd " 
            + sdatestr
            + " -ed " + edatestr 
            + " -station " + station 
            + " -sensor " + sensor)
        tmp=os.system(cmd)
        del tmp
