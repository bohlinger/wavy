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
from ncmod import dumptonc_stats
import argparse
from argparse import RawTextHelpFormatter
from utils import grab_PID
import yaml

with open("/home/patrikb/wavy/wavy/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)

# parser
parser = argparse.ArgumentParser(
    description="""
Validate wave model output against s3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./arcmfc3_validate.py
./arcmfc3_validate.py -sd 2020010112 -ed 2020020118
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

tmpdate = deepcopy(sdate)

forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]
region = args.reg
model = args.mod
sat = args.sat

while tmpdate <= edate:
    print(tmpdate)
    for element in forecasts:
        # settings
        fc_date = deepcopy(tmpdate)
        basetime=model_dict[model]['basetime']
        # get model collocated values
        from ncmod import get_arcmfc_ts

        inpath = ('/lustre/storeB/project/fou/om/waveverification/'
                + region + '/satellites/altimetry/'
                + sat + '/CollocationFiles/'
                + fc_date.strftime('%Y/%m/'))

        filename_ts=fc_date.strftime(   model + "_vs_" + sat 
                                        + "_for_" + region
                                        + "_coll_ts_lt"
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

            outpath=('/lustre/storeB/project/fou/om/waveverification/'
                    + model + '/satellites/altimetry/'
                    + sat + '/ValidationFiles/'
                    + fc_date.strftime('%Y/%m/'))

            title_stat=('validation ts for '
                    + model + ' vs ' + sat
                    + ' for region ' + region
                    + ' with leadtime '
                    + "{:0>3d}".format(element)
                    + ' h')

            filename_stat=fc_date.strftime(model + "_vs_"
                                        + sat
                                        + '_for_'
                                        + region
                                        + "_val_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
            time_dt = fc_date
            dumptonc_stats(outpath,filename_stat,title_stat,basetime,
                          time_dt,valid_dict)
    tmpdate = tmpdate + timedelta(hours=6)
