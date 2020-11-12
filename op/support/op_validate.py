#!/usr/bin/env python

import sys
import netCDF4
from datetime import datetime, timedelta
from copy import deepcopy
import argparse
from argparse import RawTextHelpFormatter
import os
import yaml
import numpy as np

sys.path.append('../../wavy')

with open("../../config/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)

from stationmod import matchtime
from validationmod import validate
from ncmod import dumptonc_stats
from utils import grab_PID

# parser
parser = argparse.ArgumentParser(
    description="""
Validate wave model output against satellite data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./op_validate.py -mod mwam4 -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-mod", metavar='model',
    help="model to be evaluated")
parser.add_argument("-sat", metavar='satellite',
    help="satellite mission to be used for collocation")
parser.add_argument("-reg", metavar='region',
    help="region of interest")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")
parser.add_argument("-path", metavar='outpath',
    help="path to where files are to be stored")

args = parser.parse_args()

now = datetime.now()

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

for sat in args.sat:
    inpath = (args.path
           + args.mod + '/satellites/altimetry/'
           + sat + '/' + 'CollocationFiles/')
    outpath = (args.path
           + args.mod + '/satellites/altimetry/' 
           + sat + '/' + 'ValidationFiles/')
    tmpdate = deepcopy(sdate)
    while tmpdate <= edate:
        print(tmpdate)
        for element in leadtimes:
            # settings
            fc_date = deepcopy(tmpdate)
            # get model collocated values
            from ncmod import get_nc_ts
            filename_ts = fc_date.strftime(args.mod
                                        + "_vs_" + sat
                                        + "_for_" + args.reg
                                        + "_coll_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
            if not os.path.exists(inpath 
                                + fc_date.strftime('%Y/%m/') 
                                + filename_ts):
                print(filename_ts + ' not found!')
            else:
                coll_dict = get_nc_ts(inpath 
                                        + fc_date.strftime('%Y/%m/') 
                                        + filename_ts,
                                        ['mHs','sHs']
                                        )
                del filename_ts
                # find collocations for given model time step and validate
                from stationmod import matchtime
                time_lst = list(coll_dict['time'].filled(np.nan))
                ctime, idx = matchtime(tmpdate,tmpdate,time_lst,
                                       coll_dict['time_unit'],
                                       timewin=30)
                if len(idx)==0:
                    pass
                else:
                    dtime = netCDF4.num2date(coll_dict['time'].filled(np.nan),
                                             coll_dict['time_unit'])
                    results_dict = {'date_matches':dtime[idx],
                        'model_matches':coll_dict['mHs'][idx].filled(np.nan),
                        'sat_matches':coll_dict['sHs'][idx].filled(np.nan)}
                    valid_dict=validate(results_dict)
                    print(valid_dict)
                    # dump to nc-file: validation
                    title_stat='validation file'
                    filename_stat=fc_date.strftime(args.mod
                                            + "_vs_" + sat
                                            + "_for_" + args.reg
                                            + "_val_ts_lt"
                                            + "{:0>3d}".format(element)
                                            + "h_%Y%m.nc")
                    time_dt = fc_date
                    dumptonc_stats(outpath + fc_date.strftime('%Y/%m/'), \
                        filename_stat,title_stat,time_dt,\
                        coll_dict['time_unit'],valid_dict)
        tmpdate = tmpdate + timedelta(hours=init_step)
