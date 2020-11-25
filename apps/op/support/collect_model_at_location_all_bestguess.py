#!/usr/bin/env python3
'''
    - retrieve data from station
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf
!!!

'''
import sys
import os
from datetime import datetime, timedelta
import numpy as np
import time
from copy import deepcopy
import argparse
from argparse import RawTextHelpFormatter
import yaml

sys.path.append('../../wavy')

from utils import grab_PID

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieves data from wave model at given location 
and dumps to monthly nc-file. If file exists, 
data is appended.

Usage:
./collect_model_at_location_all_bestguess.py -sd 2019020100 -ed 2019020200
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
print(args)
#now = datetime(2020,3,25,17)
now = datetime.now()

if args.mod is None:
    args.mod = 'mwam4'
if args.var is None:
    args.var = 'Hs'

if args.mod == 'mwam4':
    init_times = np.array([0,6,12,18]).astype('float')
    hd = 6
elif args.mod == 'ecwam':
    init_times = np.array([0,12]).astype('float')
    hd = 12

init_diffs = now.hour - init_times
init_diffs[init_diffs<0] = np.nan
h_idx = np.where(init_diffs==np.min(init_diffs[~np.isnan(init_diffs)]))
h = int(init_times[h_idx[0][0]])

if args.sd is None:
    sdate = datetime(now.year,now.month,now.day,h)
    #sdate = datetime(now.year,now.month,now.day,now.hour) - timedelta(hour=1)
    sdatestr = sdate.strftime("%Y%m%d%H")
else:
    sdatestr = args.sd
if args.ed is None:
    edate = datetime(now.year,now.month,now.day,h) + timedelta(hours=hd)
    #edate = datetime(now.year,now.month,now.day,now.hour)
    edatestr = edate.strftime("%Y%m%d%H")
else:
    edatestr = args.ed

# retrieve PID
grab_PID()
basetime = datetime(1970,1,1)

with open("../../config/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

op_support_path = os.path.abspath(os.path.dirname( __file__ ))

for station in (station_dict):
    cmd = ("python " + op_support_path
            + "/collect_model_at_location_bestguess.py"
            + " -sd " + sdatestr
            + " -ed " + edatestr 
            + " -stat " + station
            + " -mod " + args.mod
            + " -var " + args.var)
    tmp=os.system(cmd)
    del tmp
