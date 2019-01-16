#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

from datetime import datetime, timedelta
from satmod import sentinel_altimeter as sa
from stationmod import station_class as sc
from stationmod import matchtime, get_model
from modelmod import get_model, collocate
from satmod import validate
from copy import deepcopy
from model_specs import model_dict
from custom_nc import dumptonc_stats
import argparse
from argparse import RawTextHelpFormatter
from utils import grab_PID

# parser
parser = argparse.ArgumentParser(
    description="""
Validate wave model output against S3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./op_validate.py -m mwam4 -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-m", metavar='model',
    help="model to be evaluated")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")

args = parser.parse_args()

now = datetime.now()

if args.m is None:
    model = 'mwam4'
else:
    model = args.m

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

# retrieve PID
grab_PID()

# Settings
region = args.m
inpath = ('/lustre/storeB/project/fou/om/waveverification/S3a/' 
        + model
        + '/CollocationFiles/')
outpath = ('/lustre/storeB/project/fou/om/waveverification/S3a/'
        + model 
        + '/ValidationFiles/')

tmpdate = deepcopy(sdate)

leadtimes = [0, 6, 12, 18, 24]

while tmpdate <= edate:
    print(tmpdate)
    for element in leadtimes:
        # settings
        fc_date = deepcopy(tmpdate)
        basetime=model_dict[model]['basetime']
        # get model collocated values
        from custom_nc import get_arcmfc_ts
        filename_ts=fc_date.strftime(model + "_coll_ts_lt"
                                            + "{:0>3d}".format(element)
                                            + "h_%Y%m.nc")
        dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_ts)
        del filename_ts
        # find collocations for given model time step and validate
        from stationmod import matchtime
        time_lst = []
        for dt in dtime:
            time_lst.append((dt-basetime).total_seconds())
        ctime,idx = matchtime(tmpdate, tmpdate, time_lst, basetime, timewin=30)
        if len(idx)==0:
            pass
        else:
            results_dict = {'date_matches':dtime[idx],
                        'model_Hs_matches':mHs[idx],
                        'sat_Hs_matches':sHs[idx]}
            valid_dict=validate(results_dict)
            print(valid_dict)
            # dump to nc-file: validation
            title_stat='validation file'
            filename_stat=fc_date.strftime(model + "_val_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
            time_dt = fc_date
            dumptonc_stats(outpath,filename_stat,title_stat,basetime,time_dt,valid_dict)
    tmpdate = tmpdate + timedelta(hours=6)
