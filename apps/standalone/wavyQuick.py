#!/usr/bin/env python

# import standard libraries
import sys
sys.path.append('../../wavy')

import numpy as np
from datetime import datetime, timedelta
from satmod import satellite_class as sa
import argparse
from argparse import RawTextHelpFormatter
import os
import yaml

# own imports
from wavy.modelmod import model_class as mc
from wavy.validationmod import validate, disp_validation
from wavy.quicklookmod import comp_fig, plot_sat
from wavy.collocmod import collocate
from wavy.collocmod import collocation_class as coll
from wavy.ncmod import dumptonc_sat

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                        '../../', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
                        '../../', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Check availability of satellite SWH data. Example:
./wavyQuick.py -sat s3a -reg mwam4 -mod mwam4 -sd 2020100112 -lt 0 -twin 30 --col --show
./wavyQuick.py -sat s3a -reg mwam4 -sd 2020100100 -ed 2020101000 --show
./wavyQuick.py -sat s3a -reg mwam4 -sd 2020100100 -ed 2020101000 -dump outpath/
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
parser.add_argument("-var", metavar='varalias',
    help="alias for chosen variable")
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
parser.add_argument("-savep", metavar="savepath",
    help="save figure to path")
parser.add_argument("-dump", metavar="outpath",
    help="dump data to .nc-file")

args = parser.parse_args()
print ("Parsed arguments: ",args)

flatten = lambda l: [item for sublist in l for item in sublist]

# setup
if args.var is None:
    args.var = 'Hs'

sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))

if args.twin is None:
    twin = 30
else:
    twin = args.twin
if args.dist is None:
    args.dist = int(10)

if args.ed is None:
    edate = sdate
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                    int(args.ed[6:8]),int(args.ed[8:10]))
    twin = 0

# get data
if args.sat == 'all':
    satlist = ['s3a','s3b','j3','c2','al','cfo','h2b']
    lats = []
    lons = []
    var = []
    time = []
    sats = []
    satnamelst = []
    for sat in satlist:
        try:
            sa_obj_tmp = sa(sdate,sat=sat,edate=edate,twin=twin,
                            region=args.reg,varalias=args.var)
            if ('vars' in vars(sa_obj_tmp).keys()
            and len(sa_obj_tmp.vars['time'])>0):
                sa_obj = sa_obj_tmp
                lats.append(sa_obj.vars['latitude'])
                lons.append(sa_obj.vars['longitude'])
                var.append(sa_obj.vars[variable_info[args.var]\
                                            ['standard_name']])
                time.append(sa_obj.vars['time'])
                sats.append(sat)
                satnamelst.append([sat]*len(sa_obj.vars['time']))
        except:
            print(sat + ' not available')
            pass
    lats = flatten(lats)
    lons = flatten(lons)
    var = flatten(var)
    time = flatten(time)
    satnames = flatten(satnamelst)
    sa_obj.vars['latitude'] = np.array(lats)
    sa_obj.vars['longitude'] = np.array(lons)
    sa_obj.vars[variable_info[args.var]['standard_name']] = np.array(var)
    sa_obj.vars['time'] = time
    sa_obj.region = args.reg
    sa_obj.sat = str(sats)
    sa_obj.satname_ts = satnames
elif args.sat == 'multi':
    satlist = args.l.split(',')
    lats = []
    lons = []
    var = []
    time = []
    sats = []
    satnamelst = []
    for sat in satlist:
        try:
            sa_obj_tmp = sa(sdate,sat=sat,edate=edate,twin=twin,
                            region=args.reg,varalias=args.var)
            if ('vars' in vars(sa_obj_tmp).keys()
            and len(sa_obj_tmp.vars['time'])>0):
                sa_obj = sa_obj_tmp
                lats.append(sa_obj.vars['latitude'])
                lons.append(sa_obj.vars['longitude'])
                var.append(sa_obj.vars[variable_info[args.var]\
                                            ['standard_name']])
                time.append(sa_obj.vars['time'])
                sats.append(sat)
                satnamelst.append([sat]*len(sa_obj.vars['time']))
        except:
            pass
    lats = flatten(lats)
    lons = flatten(lons)
    var = flatten(var)
    time = flatten(time)
    sa_obj.vars['time'] = time
    sa_obj.vars['latitude'] = np.array(lats)
    sa_obj.vars['longitude'] = np.array(lons)
    sa_obj.vars[variable_info[args.var]['standard_name']] = np.array(var)
    sa_obj.region = args.reg
    sa_obj.sat = str(sats)
    sa_obj.satname_ts = satnames
else:
    sa_obj = sa(sdate,sat=args.sat,edate=edate,twin=twin,
                region=args.reg,varalias=args.var)

# plot
if (args.mod is None and sa_obj.region not in model_dict):
    #plot_sat(sa_obj,variable_info[args.var]['standard_name'],
    plot_sat(sa_obj,savepath=args.savep,showfig=args.show)
elif (args.mod is None and sa_obj.region in model_dict):
    print('Chosen region is a specified model domain')
    mc_obj = mc(model=sa_obj.region,
                fc_date=model_dict[sa_obj.region]['grid_date'],
                varalias=args.var)
    plot_sat(sa_obj,mc_obj=mc_obj,savepath=args.savep,showfig=args.show)
elif (args.mod is not None and args.col is True):
    # get model collocated values
    mc_obj = mc(model=args.mod,fc_date=edate,leadtime=args.lt,
                varalias=args.var)
    #collocation
    coll_obj = coll(mc_obj,
                    obs_obj=sa_obj,
                    distlim=args.dist )
    valid_dict=validate(coll_obj.vars)
    disp_validation(valid_dict)
    comp_fig(sa_obj=sa_obj,mc_obj=mc_obj,coll_obj=coll_obj,
            savepath=args.savep,showfig=args.show)
else:
    # get model collocated values
    mc_obj = mc(model=sa_obj.region,fc_date=edate,leadtime=args.lt,
                varalias=args.var)
    results_dict = {'valid_date':[edate],
                    'date_matches':[edate-timedelta(minutes=twin),
                                    edate+timedelta(minutes=twin)],
                    'model_lons_matches':sa_obj.vars['longitude'],
                    'model_lats_matches':sa_obj.vars['latitude'],
                    'sat_matches':sa_obj.vars[variable_info[args.var]\
                                                    ['standard_name']]}
    comp_fig(args.mod,sa_obj,
            mc_obj.vars[variable_info[args.var]['standard_name']],
            mc_obj.vars['longitude'],mc_obj.vars['latitude'],
            results_dict,variable_info[args.var]['standard_name'],
            savepath=args.savep,showfig=args.show)

# dump to .ncfile
if args.dump is not None:
    dumptonc_sat(sa_obj,args.dump)
