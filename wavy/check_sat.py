#!/usr/bin/env python

# import libraries
import sys
import numpy as np
from datetime import datetime, timedelta
from satmod import satellite_class as sa
from validationmod import plot_sat
import argparse
from argparse import RawTextHelpFormatter
import os
from modelmod import get_model
from validationmod import comp_fig, validate, disp_validation
from collocmod import collocate
from ncmod import dumptonc_sat
import yaml

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_shortcuts.yaml'))
with open(moddir,'r') as stream:
    shortcuts_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Check availability of satellite SWH data. Example:
./check_sat.py -sat s3a -reg mwam4 -mod mwam4 -sd 2020100112 -lt 0 -twin 30 --col --show
./check_sat.py -sat s3a -reg mwam4 -sd 2020100100 -ed 2020101000 --show
./check_sat.py -sat s3a -reg mwam4 -sd 2020100100 -ed 2020101000 -dump outpath/
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-reg", metavar='region',
    help="region to check")
parser.add_argument("-sat", metavar='satellite',
    help="satellite mission, currently available: \
        \ns3a - Sentinel-3A\
        \ns3b - Sentinel-3B\
        \nj3 - Jason-3 (reference mission)\
        \nc2 - Cryosat-2\
        \nal - SARAL/AltiKa\
        \ncfo - CFOSAT\
        \nh2b - HaiYang-2B")
parser.add_argument('-l', metavar='satellite list', 
    help='delimited list input for sats', type=str)
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period to check")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period to check")
parser.add_argument("-mod", metavar='model',
    help="chosen wave model")
parser.add_argument("-lt", metavar='lead time', type=int,
    help="lead time from initialization")
parser.add_argument("-twin", metavar='time window', type=int,
    help="time window for collocation")
parser.add_argument("-dist", metavar='distance limit', type=int,
    help="distance limit for collocation")
parser.add_argument("--col",metavar="collocation",
    help="collocation",action='store_const',const=True)
parser.add_argument("--show",
    help="show figure",action='store_const',const=True)
parser.add_argument("-dump", metavar="outpath",
    help="dump data to .nc-file")

args = parser.parse_args()
print ("Parsed arguments: ",args)

flatten = lambda l: [item for sublist in l for item in sublist]

# setup
varlst = ['Hs']

sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))

if args.twin is None:
    timewin = 30
else:
    timewin = args.twin
if args.dist is None:
    dist = 10
else:
    dist = args.dist

if args.ed is None:
    edate = sdate
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                    int(args.ed[6:8]),int(args.ed[8:10]))
    timewin = 0
if args.twin is None:
    args.lt = 0

# get data
if args.sat == 'all':
    satlist = ['s3a','s3b','j3','c2','al','cfo','h2b']
    lats = []
    lons = []
    var = []
    time = []
    for sat in satlist:
        try:
            sa_obj_tmp = sa(sdate,sat=sat,edate=edate,timewin=timewin,
                            region=args.reg,varlst=varlst)
            if ('vars' in vars(sa_obj_tmp).keys() 
            and len(sa_obj_tmp.vars['time'])>0):
                sa_obj = sa_obj_tmp
                lats.append(sa_obj.vars['latitude'])
                lons.append(sa_obj.vars['longitude'])
                var.append(sa_obj.vars[shortcuts_dict[varlst[0]]])
                time.append(sa_obj.vars['time'])
        except:
            print(sat + ' not available')
            pass
    lats = flatten(lats)
    lons = flatten(lons)
    var = flatten(var)
    time = flatten(time)
    sa_obj.vars['latitude'] = np.array(lats)
    sa_obj.vars['longitude'] = np.array(lons)
    sa_obj.vars[shortcuts_dict[varlst[0]]] = np.array(var)
    sa_obj.vars['time'] = time
    sa_obj.region = args.reg
    sa_obj.sat = str(satlist)
elif args.sat == 'multi':
    satlist = args.l.split(',')
    lats = []
    lons = []
    var = []
    time = []
    for sat in satlist:
        try:
            sa_obj_tmp = sa(sdate,sat=sat,edate=edate,timewin=timewin,
                            region=args.reg,varlst=varlst)
            if ('vars' in vars(sa_obj_tmp).keys()
            and len(sa_obj_tmp.vars['time'])>0):
                sa_obj = sa_obj_tmp
                lats.append(sa_obj.vars['latitude'])
                lons.append(sa_obj.vars['longitude'])
                var.append(sa_obj.vars[shortcuts_dict[varlst[0]]])
                time.append(sa_obj.vars['time'])
        except:
            pass
    lats = flatten(lats)
    lons = flatten(lons)
    var = flatten(var)
    time = flatten(time)
    sa_obj.vars['time'] = time
    sa_obj.vars['latitude'] = np.array(lats)
    sa_obj.vars['longitude'] = np.array(lons)
    sa_obj.vars[shortcuts_dict[varlst[0]]] = np.array(var)
    sa_obj.region = args.reg
    sa_obj.sat = str(satlist)
else:
    sa_obj = sa(sdate,sat=args.sat,edate=edate,timewin=timewin,
                region=args.reg,varlst=varlst)

# plot
if bool(args.show)==True:
    if args.mod is None:
        plot_sat(sa_obj,shortcuts_dict[varlst[0]])
    elif (args.mod is not None and args.col is True):
        # get model collocated values
        # get_model
        init_date = edate - timedelta(hours=args.lt)
        model_var_dict = get_model(model=args.mod,fc_date=edate,
                                leadtime=args.lt,init_date=init_date)
        # collocation
        results_dict = collocate(args.mod,model_var_dict['model_var'],
            model_var_dict['model_lats'],model_var_dict['model_lons'],
            model_var_dict['model_time_dt'],sa_obj,
            shortcuts_dict[varlst[0]],edate,distlim=dist)
        valid_dict=validate(results_dict)
        disp_validation(valid_dict)
        comp_fig(args.mod,sa_obj,model_var_dict['model_var'],
                model_var_dict['model_lons'],model_var_dict['model_lats'],
                results_dict,shortcuts_dict[varlst[0]])
    else:
        # get model collocated values
        # get_model
        init_date = edate - timedelta(hours=args.lt)
        model_var_dict = get_model(model=args.mod,fc_date=edate,
                                leadtime=args.lt,init_date=init_date)
        results_dict = {'valid_date':[edate],
                        'date_matches':[edate-timedelta(minutes=timewin),
                                        edate+timedelta(minutes=timewin)],
                        'model_lons_matches':sa_obj.vars['longitude'],
                        'model_lats_matches':sa_obj.vars['latitude'],
                        'sat_matches':sa_obj.vars[shortcuts_dict[varlst[0]]]}
        comp_fig(args.mod,sa_obj,model_var_dict['model_var'],
                model_var_dict['model_lons'],model_var_dict['model_lats'],
                results_dict,shortcuts_dict[varlst[0]])

# dump to .ncfile
if args.dump is not None:
    dumptonc_sat(sa_obj,args.dump)
