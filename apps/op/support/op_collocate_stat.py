# --- imports --- #

# standard
import yaml
import numpy as np
from datetime import datetime,timedelta
import os
import argparse
from argparse import RawTextHelpFormatter

# custom
from stationmod import station_class
from collocmod import collocation_class

# --- parser --- #
parser = argparse.ArgumentParser(
    description="""
Retrieve data from in-sistu stations, collocate with model,\n 
and dump data to monthly nc-file.
If file exists, data is appended.

Usage:
./op_collocate_stat.py -sd 2021010100 -ed 2021013123
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-var", metavar='varname',
    help="variable name")
parser.add_argument("-station", metavar='stationname',
    help="station name")
parser.add_argument("-sensor", metavar='sensorname',
    help="station name")
parser.add_argument("-model", metavar='modelname',
    help="model name")
parser.add_argument("-lt", metavar='leadtime',
    help="leadtime")
parser.add_argument("-dist", metavar='distance',
    help="distance limit for collocation")

args = parser.parse_args()

now = datetime.now()
if args.sd is None:
    args.sd = datetime(now.year,now.month,now.day,0)-timedelta(days=1)
else:
    args.sd = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
if args.ed is None:
    args.ed = datetime(now.year,now.month,now.day)-timedelta(minutes=1)
else:
    args.ed = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))-timedelta(minutes=1)

if args.var is None:
    args.var = 'Hs'

if args.model is None:
    args.model = 'mwam4'

if args.lt is None:
    args.lt = 0

if args.dist is None:
    args.dist = 6

print(args)

print( '# Start process of collecting platform'
        + ' data, collocate with model,\n'
        + ' and dump to nc-file #')

# --- prerequisites --- #

# get variable info
configfile = os.path.abspath(os.path.join(os.path.dirname( __file__ ), \
                            '../../../config/station_specs.yaml'))
with open(configfile,'r') as stream:
    station_dict=yaml.safe_load(stream)

# settings
if (args.station is None or args.station == 'all'):
    platformlst = station_dict['platform'].keys()
    sensorlst = None
else:
    platformlst = [args.station]
    sensorlst = [args.sensor]

date_incr = 1 # model time step

def get_and_store_data(platformlst,sensorlst,sd,ed,var,date_incr,model,lt,dist):
    for station in platformlst:
        if args.sensor is None:
            sensorlst = station_dict['platform'][station]['sensor'].keys()
        for sensor in sensorlst:
            try:
                print('station:',station,'; with sensor:',sensor)
                st_obj = station_class(station,sensor,sd,ed,
                                               varalias=var)
                col_obj = collocation_class(model=model,
                                            st_obj=st_obj,
                                            distlim=dist,
                                            leadtime=lt,
                                            date_incr=date_incr)
                # --- write to nc --- #
                col_obj.write_to_monthly_nc()
            except Exception as e:
                print(e)

# --- program body --- #
'''
Since station data comes in daily files I choose a daily loop
from sd to ed.
'''
get_and_store_data( platformlst,sensorlst,args.sd,args.ed,args.var,
                    date_incr,args.model,args.lt,args.dist)

print( '# Finished process of collecting platform'
        + ' data, collocate with model,\n'
        + ' and dump to nc-file #')
