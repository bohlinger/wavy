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

with open("../../config/pathfinder.yaml", 'r') as stream:
    pathfinder=yaml.safe_load(stream)
wavy_dir = pathfinder['wavy_dir']
sys.path.append(wavy_dir + '/wavy')

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

#now = datetime(2019,11,20,11)
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

sdatestr = sdate.strftime("%Y%m%d%H")
edatestr = edate.strftime("%Y%m%d%H")

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

with open(wavy_dir + "/config/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

for station in (station_dict):
    cmd = ("python " + wavy_dir
            + "/op/collect_model_at_location.py"
            + " -sd " + sdatestr
            + " -ed " + edatestr
            + " -station " + station 
            + " -mod " + model
            + " -var " + var)
    tmp=os.system(cmd)
    del tmp
