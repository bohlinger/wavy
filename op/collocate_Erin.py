#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

from datetime import datetime, timedelta
from satmod import satellite_altimeter as sa
from stationmod import station_class as sc
from stationmod import matchtime
from modelmod import get_model
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

args = parser.parse_args()

now = datetime.now()

if args.mod is None:
    model = 'mwam4'
else:
    model = args.mod

if args.sat is None:
    sat = 's3a'
else:
    sat = args.sat

if args.reg is None:
    region = model
else:
    region = args.reg

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

leadtimes = [0]

# settings
timewin = 30
distlim = 3
outpath = ('/lustre/storeB/project/fou/om/waveverification/'
           + model + '/satellites/altimetry'
           + '/' + sat + '/'
           + 'CollocationFiles/')

tmpdate = deepcopy(sdate)
while tmpdate <= edate:
    # get s3a values
    if 'results_dict' in globals():
        del results_dict
    fc_date = deepcopy(tmpdate)
    sa_obj = sa(fc_date,sat=sat,timewin=timewin,polyreg=region)
    if len(sa_obj.dtime)==0:
        print("If possible proceed with another time step...")
    else:
        # loop over all forecast lead times
        for element in leadtimes:
            print("leadtime: ", element, "h")
            print("fc_date: ", fc_date)
            basetime=model_dict[model]['basetime']
            filename_ts=fc_date.strftime(model 
                                        + "_vs_" + sat
                                        + "_" + region
                                        + "_coll_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
            title_ts=('collocated time series for '
                    + ' model ' + model
                    + ' vs ' + sat
                    + ' over region of ' + region
                    + ' with leadtime '
                    + "{:0>3d}".format(element)
                    + ' h')
            init_date = fc_date - timedelta(hours=element)
            #get_model
            try:
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc",model=model,fc_date=fc_date,
                    init_date=init_date,leadtime=element)
                #collocation
                if 'results_dict' in globals():
                    update_dict = collocate(model,model_Hs,model_lats,
                                        model_lons,model_time_dt,
                                        sa_obj,fc_date,distlim=distlim,
                                        idx_valid=results_dict['idx_valid'])
                    results_dict['mode_Hs_matches']=\
                                        update_dict['model_Hs_matches']
                else:
                    results_dict = collocate(model,model_Hs,model_lats,
                                        model_lons,model_time_dt,
                                        sa_obj,fc_date,distlim=distlim)
                dumptonc_ts(outpath + fc_date.strftime('%Y/%m/'), \
                            filename_ts,title_ts,basetime,results_dict)
            except IOError as e:
                print(e)
            #    print('Model output not available')
            except ValueError as e:
                print(e)
            #    print('Model wave field not available.')
            #    print('Continuing with next time step.')
    tmpdate = tmpdate + timedelta(hours=1)
