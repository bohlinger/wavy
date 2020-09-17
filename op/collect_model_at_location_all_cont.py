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
./collect_model_at_location_all_cont.py -sd 2019020100 -ed 2019020200
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

sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))

sdatestr = sdate.strftime("%Y%m%d%H")
edatestr = edate.strftime("%Y%m%d%H")

# retrieve PID
grab_PID()

with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

for station in (station_dict):
    cmd = ("python /home/patrikb/wavy/op/collect_model_at_location_cont.py"
            + " -sd " + sdatestr
            + " -ed " + edatestr
            + " -stat " + station 
            + " -mod " + args.mod
            + " -var " + args.var)
    tmp=os.system(cmd)
    del tmp
