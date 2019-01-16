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
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter
from custom_nc import get_nc_time, dumptonc_ts
from model_specs import model_dict

# parser
parser = argparse.ArgumentParser(
    description="""
Collocate wave model output and S3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./op_collocate.py -m mwam4 -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-m", metavar='model',
    help="model to be used for collocation")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")

args = parser.parse_args()

now = datetime.now()

if args.m is None:
    sdate = 'mwam4'
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

leadtimes = [0, 6, 12, 18, 24]

# settings
timewin = 30
distlim = 6
region = args.m
#model = args.m
outpath = ('/lustre/storeB/project/fou/om/waveverification/S3a/' 
           + args.m 
           + '/CollocationFiles/')

tmpdate = deepcopy(sdate)
while tmpdate <= edate:
    # get S3a values
    fc_date = deepcopy(tmpdate)
    sa_obj = sa(fc_date,timewin=timewin,region=region)
    if len(sa_obj.rtime)==0:
        print("If possible proceed with another time step...")
    else:
        # loop over all forecast lead times
        for element in leadtimes:
            print("leadtime: ", element, "h")
            print("fc_date: ", fc_date)
            basetime=model_dict[model]['basetime']
            filename_ts=fc_date.strftime(args.m + "_coll_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
            title_ts=('collocated time series for model '
                    + args.m
                    + ' with leadtime '
                    + "{:0>3d}".format(element)
                    + ' h')
            #dtime=get_nc_time(outpath+filename_ts)
            init_date = fc_date - timedelta(hours=element)
            #get_model
            try:
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc",model=model,fc_date=fc_date,
                    init_date=init_date,leadtime=element)
                #collocation
                results_dict = collocate(model,model_Hs,model_lats,
                    model_lons,model_time_dt,sa_obj,fc_date,distlim=distlim)
                dumptonc_ts(outpath,filename_ts,title_ts,basetime,results_dict)
            except IOError:
                print('Model output not available')
            except ValueError:
                print('Model wave field not available.')
                print('Continuing with next time step.')
    tmpdate = tmpdate + timedelta(hours=6)
