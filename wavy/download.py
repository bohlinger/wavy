#!/usr/bin/env python
import time
from datetime import datetime as dt
from datetime import datetime, timedelta
from satmod import get_remotefiles
from pathfinder import satpath_lustre, satpath_ftp_014_001
import argparse
from argparse import RawTextHelpFormatter

# parser
parser = argparse.ArgumentParser(
    description="""
Validate a wave model (mwam4, mwam8, ARCMFC) against observations 
(platform, satellite, buoys). Examples:
./download.py
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period to be downloaded")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period to be downloaded")

args = parser.parse_args()

now=dt.now()
if args.sd is None:
    sdate = dt(now.year,now.month,now.day)-timedelta(days=1)
else:
    sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))

if args.ed is None:
    edate = dt(now.year,now.month,now.day,now.hour)
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                int(args.ed[6:8]),int(args.ed[8:10]))

start_time = time.time()
sa_obj = get_remotefiles(satpath_ftp_014_001,satpath_lustre,
            sdate,edate,timewin=30,corenum=8,download=True)
time1 = time.time() - start_time
print("Time used for collecting data: ", time1, " seconds")
