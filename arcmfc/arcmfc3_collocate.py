#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

from datetime import datetime, timedelta
from satmod import satellite_altimeter as sa
from stationmod import station_class as sc
from stationmod import matchtime
from modelmod import get_model, check_date
from collocmod import collocate
from validationmod import validate
from copy import deepcopy
from utils import grab_PID
import argparse
from argparse import RawTextHelpFormatter
from ncmod import get_nc_time, dumptonc_ts
import yaml

with open("/home/patrikb/wavy/wavy/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Collocate wave model output and s3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./arcmfc3_collocate.py
./arcmfc3_collocate.py -sd 2019070100 -ed 2019070200
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

args = parser.parse_args()

now = datetime.now()

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

forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]

# settings
timewin = 30
region = args.reg
model = args.mod
distlim = 4

#satlist = ['s3a','s3b','j3','c2','al']
satlist = [args.sat]

tmpdate = deepcopy(sdate)
while tmpdate <= edate:
    # loop over all satellites
    for sat in satlist:
        # get s3a values
        fc_date = deepcopy(tmpdate)
        sa_obj = sa(fc_date,sat=sat,timewin=timewin,polyreg=region)
        if len(sa_obj.dtime)==0:
            print("If possible proceed with another time step...")
        else:
            # loop over all forecast lead times
            for element in forecasts:
                print("leadtime: ", element, "h")
                print("fc_date: ", fc_date)
                basetime=model_dict[model]['basetime']
                outpath=('/lustre/storeB/project/fou/om/' 
                    + 'waveverification/'
                    + model
                    + '/satellites/altimetry/'
                    + sat
                    + '/CollocationFiles/'
                    + fc_date.strftime("%Y/%m/"))
                filename_ts=fc_date.strftime(model + "_vs_"
                                        + sat 
                                        + '_for_'
                                        + region
                                        + "_coll_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
                title_ts=('collocated time series for ' 
                    + model + ' vs ' + sat
                    + ' for region ' + region
                    + ' with leadtime '
                    + "{:0>3d}".format(element)
                    + ' h')
                #dtime=get_nc_time(outpath+filename_ts)
                init_date = fc_date - timedelta(hours=element)
                #get_model
                try:
                    check_date(model,fc_date=fc_date,leadtime=element)
                    model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                        get_model(simmode="fc",model=model,fc_date=fc_date,
                        init_date=init_date,leadtime=element)
                    #collocation
                    results_dict = collocate(model,model_Hs,model_lats,
                        model_lons,model_time_dt,sa_obj,fc_date,distlim=distlim)
                    dumptonc_ts(outpath,filename_ts,title_ts,\
                                basetime,results_dict)
                except Exception as e:
                    print(e)
                    pass
    tmpdate = tmpdate + timedelta(hours=6)
