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

# parser
parser = argparse.ArgumentParser(
    description="""
Check availability of satellite SWH data. Example:
./check_sat.py -sat s3a -reg mwam4 -mod mwam4 -sd 2019060112 -lt 0 -twin 30 --col --show
./check_sat.py -sat s3a -reg mwam4 -sd 2019060112 -ed 2019060318 --show
./check_sat.py -sat s3a -reg mwam4 -sd 2019060112 -ed 2019060318 -dump outpath/
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-reg", metavar='region',
    help="region to check")
parser.add_argument("-sat", metavar='satellite',
    help="source satellite mission")
parser.add_argument('-l', metavar='list of satellites', 
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
    dtime = []
    for sat in satlist:
        sa_obj = sa(sdate,sat=sat,edate=edate,timewin=timewin,
                    region=args.reg,varlst=varlst)
        lats.append(sa_obj.vardict['latitude'])
        lons.append(sa_obj.vardict['longitude'])
        Hs.append(sa_obj.Hs)
        time.append(sa_obj.time)
        dtime.append(sa_obj.dtime)
    loc0 = flatten(loc0)
    loc1 = flatten(loc1)
    Hs = flatten(Hs)
    time = flatten(time)
    dtime = flatten(dtime)
    loc = [loc0,loc1]
    sa_obj.loc = np.array(loc)
    sa_obj.Hs = np.array(Hs)
    sa_obj.time = time
    sa_obj.dtime = dtime
    sa_obj.region = args.reg
    sa_obj.sat = str(satlist)
elif args.sat == 'multi':
    #satlist = [int(item) for item in args.list.split(',')]
    satlist = args.l.split(',')
    loc0 = []
    loc1 = []
    Hs = []
    time = []
    dtime = []
    for sat in satlist:
        sa_obj = sa(sdate,sat=sat,edate=edate,
                    timewin=timewin,region=args.reg)
        loc0.append(sa_obj.loc[0])
        loc1.append(sa_obj.loc[1])
        Hs.append(sa_obj.Hs)
        time.append(sa_obj.time)
        dtime.append(sa_obj.dtime)
    loc0 = flatten(loc0)
    loc1 = flatten(loc1)
    Hs = flatten(Hs)
    sa_obj.time = time
    sa_obj.dtime = dtime
    loc = [loc0,loc1]
    sa_obj.loc = loc
    sa_obj.Hs = Hs
    sa_obj.region = args.reg
    sa_obj.sat = str(satlist)
else:
    sa_obj = sa(sdate,sat=args.sat,edate=edate,timewin=timewin,
                region=args.reg,varlst=varlst)

# plot
if bool(args.show)==True:
    if args.mod is None:
        plot_sat(sa_obj,'Hs')
    elif (args.mod is not None and args.col is True):
        # get model collocated values
        #get_model
        init_date = edate - timedelta(hours=args.lt)
        model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(model=args.mod,fc_date=edate,
            leadtime=args.lt,init_date=init_date)
        #collocation
        results_dict = collocate(args.mod,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,edate,distlim=dist)
        valid_dict=validate(results_dict)
        disp_validation(valid_dict)
        comp_fig(args.mod,sa_obj,model_Hs,model_lons,model_lats,results_dict)
    else:
        # get model collocated values
        #get_model
        init_date = edate - timedelta(hours=args.lt)
        model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(model=args.mod,fc_date=edate,
            leadtime=args.lt,init_date=init_date)
        results_dict = {'valid_date':[edate],
                        'date_matches':[edate-timedelta(minutes=timewin),
                                        edate+timedelta(minutes=timewin)],
                        'model_lons_matches':sa_obj.loc[1],
                        'model_lats_matches':sa_obj.loc[0],
                        'sat_Hs_matches':sa_obj.Hs}
        comp_fig(args.mod,sa_obj,model_Hs,model_lons,model_lats,results_dict)

# dump to .ncfile
if args.dump is not None:
    dumptonc_sat(sa_obj,args.dump)
