#!/usr/bin/env python
'''
    - retrieve data from the Tennholmen waverider
    - collocate with wave model
    - aggregate collocated time series
    - dump to netcdf

!!! import to be implemented:
    - choose freely sdate and edate, currently only one month is possible
    - add model type to parsing options
    - add leadtime to parsing options
!!!

'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import yaml
import os
from stationmod import get_buoy, get_loc_idx
from modelmod import get_model, check_date
from collocmod import collocate
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

from ncmod import dumptonc_coll_ts_Tennholmen
from ncmod import dumptonc_coll_ts_buoy
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter

import numpy as np

# read yaml config files:
with open("/home/patrikb/wavy/wavy/buoy_specs.yaml", 'r') as stream:
    buoy_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Retrieves data from the waverider Tennholmen and dumps to monthly nc-file.
If file exists, data is appended.

Usage:
./collocate_buoy_all_bestguess.py -sd 2018110112 -ed 2018110118
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

if model == 'mwam4':
    init_step = 6
    init_start = np.min([0,6,12,18])
if model == 'mwam8':
    init_step = 12
    init_start = 12
tdeltas = range(1,init_step+1)[::-1]
leadtimes = range(init_step)

buoylst = buoy_dict.keys()
for buoy in (buoylst):
    cmd = ("python collocate_buoy_bestguess.py -sd "
            + sdatestr
            + " -ed " + edatestr
            + " -buoy " + buoy)
    tmp=os.system(cmd)
    del tmp
