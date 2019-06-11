#!/usr/bin/env python

# import libraries
from datetime import datetime, timedelta
from satmod import sentinel_altimeter as sa
from validationmod import plot_S3a
import argparse
from argparse import RawTextHelpFormatter
import os

# parser
parser = argparse.ArgumentParser(
    description="""
Check Sentinel-3a data. Example:
./check_S3a.py -r ARCMFC -sd 2018080112 -ed 2018080718 -m -a --show -save outpath/ -dump outpath/
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-r", metavar='region',
    help="region to check")
parser.add_argument("-sd", metavar='startdate',
    help="start date of time period to check")
parser.add_argument("-ed", metavar='enddate',
    help="end date of time period to check")
parser.add_argument("-m",
    help="make map",action='store_const',const=True)
parser.add_argument("-a",
    help="compute availability",action='store_const',const=True)
parser.add_argument("--show",
    help="show figure",action='store_const',const=True)
parser.add_argument("-save",metavar='outpath',
    help="save figure(s)")
parser.add_argument("-dump", metavar="outpath",
    help="dump data to .nc-file")

args = parser.parse_args()
print ("Parsed arguments: ",args)

# setup
sdate = datetime(int(args.sd[0:4]),int(args.sd[4:6]),
                int(args.sd[6:8]),int(args.sd[8:10]))
if args.ed is None:
    timewin = 30
    edate = sdate
else:
    edate = datetime(int(args.ed[0:4]),int(args.ed[4:6]),
                    int(args.ed[6:8]),int(args.ed[8:10]))
    timewin = 0
print(timewin)

# get data
sa_obj = sa(sdate,edate=edate,timewin=timewin,polyreg=args.r)

# plot
if bool(args.m)==True:
    plot_S3a(sa_obj)

# check availability
if bool(args.a)==True:
    freqlst,datelst=sa_obj.bintime()
    sa_obj.plotavail(datelst,freqlst,show=bool(args.show),save=bool(args.save))

if args.save is not None:
    os.system('mkdir -p ' + args.save)
    os.system('mv altimeter*.pdf ' + args.save)

# dump to .ncfile
if args.dump is not None:
    sa_obj.dumptonc(args.dump)
