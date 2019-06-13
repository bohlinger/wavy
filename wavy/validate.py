#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

from datetime import datetime, timedelta
from satmod import satellite_altimeter as sa
from stationmod import station_class as sc
from stationmod import matchtime
from modelmod import get_model, check_date
from collocmod import collocate
#from satmod import validate
from validationmod import comp_fig, validate
from utils import identify_outliers
import os
import argparse
from argparse import RawTextHelpFormatter
from  utils import disp_validation

# parser
parser = argparse.ArgumentParser(
    description="""
Validate a wave model (mwam4, mwam8, ARCMFC) against observations 
(platform, satellite, buoys). Examples:
./validate.py -m ARCMFC -fc 2018080218 -sat s3a
./validate.py -m mwam4 -fc 2018080218 -sat s3a
./validate.py -m mwam4 -fc 2018080218 -lt 6 -sat s3a
./validate.py -m mwam4 -fc 2018120412 -lt 18 -sat s3a --show
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-m", metavar='model', type=str,
    help="model to check (mwam4,mwam8,ARCMFC)")
parser.add_argument("-r", metavar='region', type=str,
    help="region to check")
parser.add_argument("-lt", metavar='leadtime', type=int,
    help="leadtime in hours")
parser.add_argument("-fc", metavar='fcdate',
    help="forecasted date to check")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period to be evaluated")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period to be evaluated")
parser.add_argument("-plat", metavar='platform',
    help="name of platform")
parser.add_argument("-sat", metavar='satellite',
    help="name of satellite")
parser.add_argument("-buoy", metavar='buoy',
    help="name of buoy")
parser.add_argument("-dts", metavar='outpath',
    help="dump time series to nc files in folder outpath")
parser.add_argument("-dval", metavar='outpath',
    help="dump validation statistics to nc files in folder outpath")
parser.add_argument("-sfig", metavar='outpath',
    help="NOT AVAILABLE (save figure(s) in folder outpath)")
parser.add_argument("--diag", action='store_const', const=True,
    help="make diagnostics")
parser.add_argument("--show", action='store_const', const=True,
    help="show figures")

args = parser.parse_args()
print ("Parsed arguments: ",args)

# setup
if (args.lt is None and args.fc is not None):
    fc_date = datetime(int(args.fc[0:4]),int(args.fc[4:6]),
                int(args.fc[6:8]),int(args.fc[8:10]))
    timewin = 30
    init_date = fc_date
elif (args.lt is not None and args.fc is not None):
    fc_date = datetime(int(args.fc[0:4]),int(args.fc[4:6]),
                int(args.fc[6:8]),int(args.fc[8:10]))
    init_date = fc_date - timedelta(hours=args.lt)
    timewin = 30

if (args.ed is None and args.sd is not None):
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
    edate = sdate
elif (args.ed is not None and args.sd is not None):
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))
if (args.fc is None and args.sd is None):
    sys.exit("-> Error: A date or time period needs to be given!")
if (args.plat is None and args.sat is None and args.buoy is None):
    sys.exit("-> Error: A source of observations needs to be given!")
if (args.m is None):
    sys.exit("-> Error: A model to validate needs to be given!")
if (args.r is None):
    args.r = args.m
if (args.sat is None):
    sat = 's3a'
else: sat = args.sat

# Get sat data
#sa_obj = sa(fc_date,timewin=timewin,region=args.m)
#sa_obj = sa(fc_date,timewin=timewin,polyreg='BarentsSea')
sa_obj = sa(fc_date,sat=sat,timewin=timewin,polyreg=args.r)
if len(sa_obj.dtime)==0:
    print("If possible proceed with another time step...")
else:
    if (args.m != 'ARCMFC' and args.m != 'MoskNC'\
        and args.m != 'MoskWC'):
        # get model collocated values
        check_date(args.m,fc_date=fc_date,leadtime=args.lt)
        #get_model
        model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(simmode="fc",model=args.m,fc_date=fc_date,
            leadtime=args.lt)
        #collocation
        results_dict = collocate(args.m,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=6)
        valid_dict=validate(results_dict)
        #print(valid_dict)
        disp_validation(valid_dict)

    if (args.m == 'ARCMFC'):
        # get model collocated values
        check_date(args.m,fc_date=fc_date,leadtime=args.lt)
        model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(simmode="fc",model=args.m,fc_date=fc_date,
            init_date=init_date)
        #collocation
        results_dict = collocate(args.m,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=6)
        valid_dict=validate(results_dict)
        #print(valid_dict)
        disp_validation(valid_dict)

    if (args.m == 'MoskNC' or args.m == 'MoskWC'):
        region = "Mosk_dom"
        # get model collocated values
        model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(simmode="fc",model=args.m,fc_date=fc_date,
            init_date=init_date)
        #collocation
        results_dict = collocate(args.m,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=6)
        #print results_dict
        valid_dict=validate(results_dict)
        #print(valid_dict)
        disp_validation(valid_dict)

    if args.show is True:
        comp_fig(args.m,sa_obj,model_Hs,model_lons,model_lats,results_dict)

    if args.dts is not None:
        # dump to nc-file
        from model_specs import model_dict
        from ncmod import dumptonc_ts
        basetime=model_dict[args.m]['basetime']
        outpath=args.dts
        filename_ts= args.m + "_nc_ts.nc"
        title_ts='collocated time series'
        dumptonc_ts(outpath,filename_ts,title_ts,basetime,results_dict)

    if args.dval is not None:
        # dump to nc-file
        from model_specs import model_dict
        from ncmod import dumptonc_stats
        basetime=model_dict[args.m]['basetime']
        outpath=args.dval
        filename_stat= args.m + "_nc_val.nc"
        title_stat='validation file'
        time_dt = fc_date
        dumptonc_stats(outpath,filename_stat,title_stat,basetime,time_dt,valid_dict)

    if args.plat is not None:
        sc_obj = sc(args.plat,sdate,edate)
        ctime, cidx = matchtime(fc_date,fc_date,sc_obj.time,sc_obj.basedate)
