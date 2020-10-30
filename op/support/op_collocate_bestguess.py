#!/usr/bin/env python
from pathlib import Path
home = str(Path.home())
wavy_script_dir = home + '/wavy/wavy'
import sys
sys.path.append(wavy_script_dir)

from datetime import datetime, timedelta
from satmod import satellite_altimeter as sa
from stationmod import station_class as sc
from stationmod import matchtime
from modelmod import get_model
from collocmod import collocate
from validationmod import validate
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter
from ncmod import get_nc_time, dumptonc_ts, check_vals_in_nc
import yaml
with open(home + "/wavy/config/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
                        '../..', 'config/variable_shortcuts.yaml'))
with open(moddir,'r') as stream:
    shortcuts_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Collocate wave model output and s3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./op_collocate_bestguess.py -mod mwam4 -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-mod", metavar='model',
    help="model to be used for collocation")
parser.add_argument("-sat", metavar='satellite',
    help="satellite mission to be used for collocation")
parser.add_argument("-reg", metavar='region',
    help="region of interest")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-twin", metavar='time window', type=int,
    help="time window for collocation")
parser.add_argument("-dist", metavar='distance limit', type=int,
    help="distance limit for collocation")
parser.add_argument("-path", metavar='outpath',
    help="path to where files are to be stored")

args = parser.parse_args()

now = datetime.now()

varlst = ['Hs']

if args.mod is None:
    args.mod = 'mwam4'

if args.sat is None:
    args.sat = ['s3a']
elif args.sat == "all":
    args.sat = ['s3a','s3b','c2','j3','al','cfo']
else: args.sat = [args.sat]

if args.reg is None:
    args.reg = args.mod

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

if args.path is None:
    args.path = '/lustre/storeB/project/fou/om/waveverification/'

# retrieve PID
grab_PID()

init_step = model_dict[args.mod]['init_step']

leadtimes = range(init_step)

# settings
if args.twin is None:
    args.twin = 30
if args.dist is None:
    args.dist = 6

for sat in args.sat:
    outpath = (args.path
           + args.mod + '/satellites/altimetry'
           + '/' + sat + '/'
           + 'CollocationFiles/')
    tmpdate = deepcopy(sdate)
    while tmpdate <= edate:
        for element in leadtimes:
            fc_date = tmpdate + timedelta(hours=element)
            init_date = deepcopy(tmpdate)
            print("leadtime: ", element, "h")
            print("fc_date: ", fc_date)
            print("init_date: ", init_date)
            # get sat values
            try:
                sa_obj = sa(fc_date,sat=sat,timewin=args.twin,polyreg=args.reg)
                if len(sa_obj.dtime)==0:
                    print("If possible proceed with another time step...")
                else:
                    d = model_dict[args.mod]['basetime']
                    # conversion from yaml-datetime.date to datetime.datetime
                    basetime = datetime(d.year,d.month,d.day)
                    filename_ts = fc_date.strftime(args.mod
                                    + "_vs_" + sat
                                    + "_for_" + args.reg
                                    + "_coll_ts_lt_best"
                                    + "_%Y%m.nc")
                    title_ts=('collocated time series for '
                        + ' model ' + args.mod
                        + ' vs ' + sat
                        + ' over region ' + args.reg
                        + ' with leadtime '
                        + ' bestguess')
                    # get_model
                    try:
                        model_Hs,\
                        model_lats,\
                        model_lons,\
                        model_time,\
                        model_time_dt = get_model(simmode="fc",
                                                model=args.mod,
                                                fc_date=fc_date,
                                                init_date=init_date,
                                                leadtime=element)
                        # collocation
                        results_dict = collocate(args.mod,model_Hs,
                                        model_lats,model_lons,
                                        model_time_dt,sa_obj,
                                        fc_date,distlim=args.dist)
                        dumptonc_ts(outpath + fc_date.strftime('%Y/%m/'), \
                            filename_ts,title_ts,basetime,results_dict)
                    except IOError as e:
                        print(e)
                    except ValueError as e:
                        print(e)
            except IndexError as e:
                print(e)
        tmpdate = tmpdate + timedelta(hours = init_step)
