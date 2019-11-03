#!/usr/bin/env python
'''
    - retrieve and aggregate data from the Tennholmen waverider
    - dump to netcdf
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/op')

import os
from datetime import datetime, timedelta

import yaml
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter

with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_specs=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieve data from the waverider Tennholmen and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./collect_station_all.py -sd 2019100100 -ed 2019110100 -var varname
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-var", metavar='varname',
    help="variable name")

args = parser.parse_args()

now = datetime.now()
if args.sd is None:
    sdate = datetime(now.year,now.month,now.day,now.hour)-timedelta(hours=1)
    sdstr = sdate.strftime("%Y%m%d%H")
else:
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
    sdstr = sdate.strftime("%Y%m%d%H")
if args.ed is None:
    edate = datetime(now.year,now.month,now.day,now.hour)
    edstr = edate.strftime("%Y%m%d%H")
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))
    edstr = edate.strftime("%Y%m%d%H")

# retrieve PID
grab_PID()

# stations
for station in station_specs.keys():
    for sensor in station_specs[station]['sensor'].keys():
        cmdstr = ("python /home/patrikb/wavy/op/collect_station.py"
            + " -sd "
            + sdstr
            + " -ed "
            + edstr
            + " -station " + station
            + " -sensor " + sensor
            + " -var " + args.var)
        str('Executing command: ')
        str(cmdstr)
        cmd = os.system(cmdstr)
