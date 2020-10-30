#!/usr/bin/env python
'''
    - retrieve and aggregate data from the Tennholmen waverider
    - dump to netcdf
'''
import sys
import os
from datetime import datetime, timedelta
import yaml
import argparse
from argparse import RawTextHelpFormatter

sys.path.append(r'../../wavy')

with open("../../config/station_specs.yaml", 'r') as stream:
    station_specs=yaml.safe_load(stream)

from utils import grab_PID

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

op_support_path = os.path.abspath(os.path.dirname( __file__ ))

# stations
for station in station_specs.keys():
    for sensor in station_specs[station]['sensor'].keys():
        cmdstr = ("python " + op_support_path
            + "/collect_station.py"
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
