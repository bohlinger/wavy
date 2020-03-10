#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import time
from datetime import datetime, timedelta
import calendar
from copy import deepcopy
import numpy as np
from dateutil.relativedelta import relativedelta
from calendar import monthrange
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from utils import grab_PID
from stationmod import matchtime
import os
import argparse
from argparse import RawTextHelpFormatter
from ncmod import get_arcmfc_stats

parser = argparse.ArgumentParser(
    description="""
Main program to run the monthly ARCMFC validation
with satellite.\n
Usage example in unix command line: 
./arcmfc_ncfile.py -d 201808\n
The argument consists of the year and month to be validated
If no date is given the last month is validated.
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument('-d',help="validation for given month",
                    type=int)
args = parser.parse_args()
#args, rest = parser.parse_known_args()

if len(sys.argv)<=2:
    now = datetime.now()-relativedelta(months=1)
else:
    nowstr = str(sys.argv[2])
    now = datetime(int(nowstr[0:4]),int(nowstr[4:6]),1)

# retrieve PID
grab_PID()

flatten = lambda l: [item for sublist in l for item in sublist]

pathstr=('/lustre/storeA/project/fou/om/waveverification/Arc-MFC/' 
        + 'monthly/' 
        + now.strftime('%Y') 
        + '_' 
        + now.strftime('%b') + '/')
filestr=('product_quality_stats_ARCTIC_ANALYSIS_FORECAST_WAV_002_006_' 
        + now.strftime('%Y') 
        + now.strftime('%m') 
        + '01-' 
        + now.strftime('%Y') + now.strftime('%m') 
        + str(monthrange(now.year, now.month)[1]) + '.nc')

nc = Dataset(pathstr + filestr, 'r+')
print ("\nAppending global attribute to netcdf validation file ...")
globalAttribs = {}
globalAttribs["product"] = "ARCTIC_ANALYSIS_FORECAST_WAV_002_014"
nc.setncatts(globalAttribs)
nc.sync()
print ("\nAdjust _FillValue")
plat = nc.variables['stats_VHM0_platform'][:]
platnan = (np.where(np.isnan(plat)))
if len(platnan[0])>0:
    print("NaNs found for platform data")
    for i in range(len(platnan[0])):
        plat[platnan[0][i],platnan[1][i],platnan[2][i],platnan[3][i],platnan[4][i]] = 9999.
alt = nc.variables['stats_VHM0_altimeter'][:]
altnan = (np.where(np.isnan(alt)))
if len(altnan[0])>0:
    print("NaNs found for altimeter data")
    for i in range(len(altnan[0])):
        alt[altnan[0][i],altnan[1][i],altnan[2][i],altnan[3][i],altnan[4][i]] = 9999.
nc.variables['stats_VHM0_platform'][:] = plat[:]
nc.variables['stats_VHM0_altimeter'][:] = alt[:]
nc.close()
print ("Attribute appended")
print ("Rename files to match naming of ARCMFC WAM 3km model product")
new_filestr=('product_quality_stats_ARCTIC_ANALYSIS_FORECAST_WAV_002_014_'
        + now.strftime('%Y')
        + now.strftime('%m')
        + '01-'
        + now.strftime('%Y') + now.strftime('%m')
        + str(monthrange(now.year, now.month)[1]) + '.nc')
cmd = "mv " + pathstr + filestr + " " + pathstr + new_filestr
h = os.system(cmd)
print ("### VALIDATION FINISHED ###")
