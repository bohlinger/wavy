#!/usr/bin/env python

from datetime import datetime, timedelta
from copy import deepcopy
import argparse
from argparse import RawTextHelpFormatter
import os
import yaml
import sys

sys.path.append('../../wavy')

with open("../../config/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                        '../..', 'config/variable_shortcuts.yaml'))
with open(moddir,'r') as stream:
    shortcuts_dict=yaml.safe_load(stream)

from satmod import satellite_class as sa
from modelmod import get_model
from collocmod import collocate
from utils import grab_PID
from ncmod import dumptonc_ts

# parser
parser = argparse.ArgumentParser(
    description="""
Collocate wave model output and satellite data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./op_collocate.py -mod mwam4 -sd 2018110112 -ed 2018110118
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
    args.sat = ['s3a','s3b','c2','j3','al','cfo','h2b']
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

# define leadtimes
leadtimes = model_dict[args.mod]['leadtimes']
# define init_step
init_step = model_dict[args.mod]['init_step']

# settings
if args.twin is None:
    args.twin = 30
if args.dist is None:
    args.dist = 6

for sat in args.sat:
    outpath = (args.path + '/'
           + args.mod + '/satellites/altimetry'
           + '/' + sat + '/'
           + 'CollocationFiles/')
    tmpdate = deepcopy(sdate)
    while tmpdate <= edate:
        # get sat values
        if 'results_dict' in globals():
            del results_dict
        fc_date = deepcopy(tmpdate)
        sa_obj = sa(fc_date,sat=sat,timewin=args.twin,region=args.reg,
                    varlst=varlst)
        if ('vars' not in vars(sa_obj) or len(sa_obj.vars['time'])==0):
            print("If possible proceed with another time step...")
        else:
            # loop over all forecast lead times
            for element in leadtimes:
                print("leadtime: ", element, "h")
                print("fc_date: ", fc_date)
                filename_ts=fc_date.strftime(args.mod
                                            + "_vs_" + sat
                                            + "_for_" + args.reg
                                            + "_coll_ts_lt"
                                            + "{:0>3d}".format(element)
                                            + "h_%Y%m.nc")
                title_ts=('collocated time series for '
                        + ' model ' + args.mod
                        + ' vs ' + sat
                        + ' over region of ' + args.reg
                        + ' with leadtime '
                        + "{:0>3d}".format(element)
                        + ' h')
                init_date = fc_date - timedelta(hours=element)
                # get_model
                try:
                    model_var_dict = \
                        get_model(model=args.mod,fc_date=fc_date,
                        leadtime=element,init_date=init_date)
                    # collocation
                    if ('results_dict' in globals() 
                        and len(results_dict['idx_valid'])>0):
                        update_dict = collocate(args.mod,
                                            model_var_dict['model_var'],
                                            model_var_dict['model_lats'],
                                            model_var_dict['model_lons'],
                                            model_var_dict['model_time_dt'],
                                            sa_obj,shortcuts_dict[varlst[0]],
                                            fc_date,distlim=args.dist,
                                            idx_valid=results_dict['idx_valid'])
                        results_dict['model_matches']=\
                                            update_dict['model_matches']
                    else:
                        results_dict = collocate(args.mod,
                                            model_var_dict['model_var'],
                                            model_var_dict['model_lats'],
                                            model_var_dict['model_lons'],
                                            model_var_dict['model_time_dt'],
                                            sa_obj,shortcuts_dict[varlst[0]],
                                            fc_date,distlim=args.dist)
                    dumptonc_ts(outpath + fc_date.strftime('%Y/%m/'), \
                                filename_ts,title_ts,\
                                model_var_dict['model_time_unit'],\
                                results_dict)
                except IOError as e:
                    print(e)
                #    print('Model output not available')
                except ValueError as e:
                    print(e)
                #    print('Model wave field not available.')
                #    print('Continuing with next time step.')
                except IndexError as e:
                    print(e)
                #    print('time step not in model output file')
        tmpdate = tmpdate + timedelta(hours = init_step)
