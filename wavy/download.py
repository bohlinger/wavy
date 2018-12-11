#!/usr/bin/env python
import time
from datetime import datetime, timedelta
from satmod import get_remotefiles
from pathfinder import satpath_lustre, satpath_ftp_014_001
import argparse
from argparse import RawTextHelpFormatter

# parser
parser = argparse.ArgumentParser(
    description="""
Download S3a netcdf from Copernicus DU.

Usage:
./download.py
./download.py -sd 2018110112 -ed 2018110118
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period to be downloaded")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period to be downloaded")

args = parser.parse_args()

now=datetime.now()
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

start_time = time.time()
sa_obj = get_remotefiles(satpath_ftp_014_001,satpath_lustre,
            sdate,edate,timewin=30,corenum=8,download=True)
time1 = time.time() - start_time
print("Time used for collecting data: ", time1, " seconds")
