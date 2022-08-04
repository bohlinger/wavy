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
from wavy.sat_collectors import get_remote_files
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
                        help="satellite mission, currently available\n \
            \ncmems_L3_NRT:\
            \n s3a - Sentinel-3A\
            \n s3b - Sentinel-3B\
            \n j3 - Jason-3 (reference mission)\
            \n c2 - Cryosat-2\
            \n al - SARAL/AltiKa\
            \n cfo - CFOSAT\
            \n h2b - HaiYang-2B\
            \n\
            \ncmems_L3_s6a:\
            \n s6a - Sentinel-6A Michael Freilich\
            \n\
            \neumetsat_L2:\
            \n s3a - Sentinel-3A\
            \n s3b - Sentinel-3B\
            \n\
            \ncci_L2P:\
            \n j1 - Jason-1\
            \n j2 - Jason-2\
            \n j3 - Jason-3\
            \n c2 - Cryosat-2\
            \n envisat - Envisat\
            \n ers1 - European Remote-Sensing Satellite-1\
            \n ers2 - European Remote-Sensing Satellite-2\
            \n topex - TOPEX/Poseidon\
            \n al - SARAL/AltiKa\
            \n gfo - GEOSAT Follow-On\
            \n\
            \ncci_L3:\
            \n multi - multimission product 1991-2018\
            \n\
            \ncfo_swim_L2P:\
            \n cfo - CFOSAT\
            \n\
            \n")
    parser.add_argument("-path",
                        metavar='path',
                        help="destination for downloaded data")
    parser.add_argument("-product",
                        metavar='product',
                        help="product name as specified in *_specs.yaml")
    parser.add_argument("-api_url",
                        metavar='api_url',
                        help="source of eumetsat L2 data")
    parser.add_argument("-nproc",
                        metavar='nproc',
                        help="number of possible simultaneous downloads",
                        type=int)

    args = parser.parse_args()

    # settings
    now = datetime.now()
    if args.sat is None:
        satlst = ['s3a']
    elif args.sat == 'cmems_L3_NRT':
        satlst = satellite_dict['cmems_L3_NRT']['mission'].keys()
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

    if args.product is None:
        args.product = 'cmems_L3_NRT'

    print(args)

    twin = 30

    for sat in satlst:
        for i in range(1):
            #    try:
            print("Attempting to download data for:", sat)
            print("Time period:", str(sdate), "to", str(edate))
            start_time = time.time()
            dict_for_sub = {'mission':sat}
            #path_local,sdate,edate,twin,nproc,product,api_url,sat,dict_for_sub
            #sa_obj = get_remote_files(\
            #                    args.path, sdate, edate, twin,
            #                    args.nproc, args.product,
            #                    args.api_url, sat,
            #                    dict_for_sub)
            sa_obj = get_remote_files(\
                                path_local=args.path,
                                sdate=sdate,edate=edate,
                                twin=twin,nproc=args.nproc,
                                product=args.product,
                                api_url=args.api_url,
                                mission=sat,
                                dict_for_sub=dict_for_sub)
            time1 = time.time() - start_time
            print("Time used for collecting data: ", time1, " seconds")


#    except Exception as e:
#        print('Experienced error when downloading data for',sat,
#            '\nwith the error:',e,
#            '\nSkip and continue ...')

if __name__ == "__main__":
    main()
