#!/usr/bin/env python
'''
    - retrieve and aggregate data from the Tennholmen waverider
    - dump to netcdf
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patri/home/patrikb/wavy/wavy/op')

import os
from datetime import datetime, timedelta

from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter

from station_specs import station_dict

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieve data from the waverider Tennholmen and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./op_collect.py -sd 2018110112 -ed 2018110118
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

# buoys
cmdstr = ("python /home/patri/home/patrikb/wavy/wavy/op/collect_buoy.py -sd " 
            + sdstr 
            + " -ed " 
            + edstr)
str('Executing command: ')
str(cmdstr)
#cmd = os.system(cmdstr)

# stations
stationlst = station_dict.keys()
for station in (stationlst):
    for sensor in (station_dict[station]['sensor'].keys()):
        cmdstr = ("python /home/patri/home/patrikb/wavy/wavy/op/collect_station.py"
            + " -sd "
            + sdstr
            + " -ed "
            + edstr
            + " -station " + station
            + " -sensor " + sensor)
        str('Executing command: ')
        str(cmdstr)
        cmd = os.system(cmdstr)
