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
sys.path.append(r'/home/patrikb/wavy/wavy')
#sys.path.append(os.getenv("HOME") + "/met-ecflow-support/lib/python")
#import python
#python.module('load', 'compiler/intelPE2018')
#python.module('load', 'hdf5/1.10.5-intel2018')
#python.module('load', 'netcdf/4.7.0-intel2018')
#python.module('load', 'Python/3.7.2')

from datetime import datetime, timedelta
import numpy as np
import time

from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter
import yaml

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

#now = datetime(2019,11,20,11)
now = datetime.now()

init_times = np.array([0,6,12,18]).astype('float')
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
    edate = datetime(now.year,now.month,now.day,h) + timedelta(hours=6)
    #edate = datetime(now.year,now.month,now.day,now.hour)
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

with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

for station in (station_dict):
    cmd = ("python /home/patrikb/wavy/op/collect_model_at_location_bestguess.py"
            + " -sd " + sdatestr
            + " -ed " + edatestr 
            + " -station " + station 
            + " -mod " + model
            + " -var " + var)
    tmp=os.system(cmd)
    del tmp
