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

# retrieve PID
grab_PID()

# parser
parser = argparse.ArgumentParser(
    description="""
Collocate wave model output and S3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./arcmfc_collocate.py
./arcmfc_collocate.py -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period")

args = parser.parse_args()

now = datetime.now()
if args.sd is None:
    sdate = datetime(now.year,now.month,now.day)-timedelta(days=1)
else:
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))

if args.ed is None:
    edate = datetime(now.year,now.month,now.day,now.hour)
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))


forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]

# settings
timewin = 30
region = 'ARCMFC'
model = 'ARCMFC'
distlim = 6

# get S3a values
#sa_obj = sa(fc_date,timewin=timewin,region=region)

#date = (datetime(now.year,now.month,now.day)
#        - timedelta(days=1) 
#        + timedelta(hours=12))
#tmpdate = deepcopy(date)

sdate = datetime(2018,10,1)
edate = datetime(2018,12,10)
dummydate = deepcopy(startdate)
sa_obj = sa(sdate,edate=enddate,timewin=timewin,
            region="ARCMFC",mode="ARCMFC")


while dummydate<edate:
    fc_date = dummydate + timedelta(hours=12)
    for element in forecasts:
        print("leadtime: ", element, "h")
        init_date = fc_date - timedelta(hours=element)
        #get_model
        model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(simmode="fc",model=model,fc_date=fc_date,
            init_date=init_date)
        #collocation
        results_dict = collocate(model,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=distlim)
        try:
            from model_specs import model_dict
            # dump to nc-file: ts
            from custom_nc import dumptonc_ts
            basetime=model_dict[model]['basetime']
            outpath=('/lustre/storeB/project/fou/om/ARCMFC/'
                    + 'S3a/CollocationFiles/')
            filename_ts=dummydate.strftime("ARCMFC_coll_ts_lt" 
                                        + "{:0>3d}".format(element) 
                                        + "h_%Y%m.nc")
            title_ts=('collocated time series for ARCMFC with leadtime '
                    + "{:0>3d}".format(element)
                    + ' h')
            dumptonc_ts(outpath,filename_ts,title_ts,basetime,results_dict)
        except:
            print('no matches')
            pass
    dummydate = dummydate + timedelta(hours=24)
