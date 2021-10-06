#!/usr/bin/env python3
"""
download satellite data from Copernicus
"""
# --- imports -------------------------------------------------------- #
# standard library imports
import os
import time
from datetime import datetime, timedelta
import yaml
import argparse
from argparse import RawTextHelpFormatter

# own import
from wavy.satmod import get_remote_files
from wavy.wconfig import load_or_default
# -------------------------------------------------------------------- #


def main():
    # read yaml config files:
    satellite_dict = load_or_default('satellite_specs.yaml')

# parser
    parser = argparse.ArgumentParser(description="""
    Download satellite nc-files

    Usage:
    ./wavyDownload.py -sat s3a -sd 2020100100 -ed 2020101000
        """,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument("-sd",
                        metavar='startdate',
                        help="start date of time period to be downloaded")
    parser.add_argument("-ed",
                        metavar='enddate',
                        help="end date of time period to be downloaded")
    parser.add_argument("-sat",
                        metavar='satellite',
                        help="satellite mission, currently available: \
            \ns3a - Sentinel-3A\
            \ns3b - Sentinel-3B\
            \nj3 - Jason-3 (reference mission)\
            \nc2 - Cryosat-2\
            \nal - SARAL/AltiKa\
            \ncfo - CFOSAT\
            \nh2b - HaiYang-2B\
            \nall - all availabe satellites")
    parser.add_argument("-path",
                        metavar='path',
                        help="destination for downloaded data")
    parser.add_argument("-provider",
                        metavar='provider',
                        help="institution providing the data")
    parser.add_argument("-api_url",
                        metavar='api_url',
                        help="source of eumetsat L2 data")
    parser.add_argument("-nproc",
                        metavar='nproc',
                        help="number of simultaneous processes",
                        type=int)

    args = parser.parse_args()

    # settings
    instr = 'altimeter'

    now = datetime.now()
    if args.sat is None:
        satlst = ['s3a']
    elif args.sat == 'all':
        satlst = satellite_dict['cmems']['mission'].keys()
    else:
        satlst = [args.sat]

    if args.sd is None:
        sdate = datetime(now.year, now.month, now.day,
                         now.hour) - timedelta(hours=24)
    else:
        sdate = datetime(int(args.sd[0:4]), int(args.sd[4:6]),
                         int(args.sd[6:8]), int(args.sd[8:10]))

    if args.ed is None:
        edate = datetime(now.year, now.month, now.day, now.hour, now.minute)
    else:
        edate = datetime(int(args.ed[0:4]), int(args.ed[4:6]),
                         int(args.ed[6:8]), int(args.ed[8:10]))

    if args.nproc is None:
        args.nproc = 1

    if args.provider is None:
        args.provider = 'cmems'

    print(args)

    twin = 30

    for sat in satlst:
        for i in range(1):
            #    try:
            print("Attempting to download data for:", sat)
            print("Time period:", str(sdate), "to", str(edate))
            start_time = time.time()
            dict_for_sub = {'mission':sat}
            sa_obj = get_remote_files(\
                                args.path, sdate, edate, twin,
                                args.nproc, args.provider,
                                args.api_url, sat, dict_for_sub)
            time1 = time.time() - start_time
            print("Time used for collecting data: ", time1, " seconds")


#    except Exception as e:
#        print('Experienced error when downloading data for',sat,
#            '\nwith the error:',e,
#            '\nSkip and continue ...')

if __name__ == "__main__":
    main()
