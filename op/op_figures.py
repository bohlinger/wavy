#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from datetime import datetime, timedelta
from copy import deepcopy
from graphicsmod import make_val_ts_fig_op, make_val_scatter_fig_op
from ncmod import get_arcmfc_stats, get_arcmfc_ts
import argparse
from argparse import RawTextHelpFormatter

# parser
parser = argparse.ArgumentParser(
    description="""
Validate wave model output against s3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./op_figures.py -mod mwam4
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-mod", metavar='model',
    help="model to be evaluated")
parser.add_argument("-sat", metavar='satellite',
    help="satellite mission to be used for collocation")
parser.add_argument("-reg", metavar='region',
    help="region of interest")

args = parser.parse_args()

now = datetime.now()

if args.mod is None:
    args.mod = 'mwam4'

if args.sat is None:
    args.sat = 's3a'

if args.reg is None:
    args.reg = model

fc_date = datetime.now()
forecasts = [0]

for element in forecasts:
    # Get stats ts
    inpath = ('/lustre/storeB/project/fou/om/waveverification/'
                + args.mod + '/satellites/altimetry'
                + '/' + args.sat + '/'
                + 'ValidationFiles/'
                + fc_date.strftime('%Y/%m/'))
    filename_stat = fc_date.strftime(args.mod
                                    + "_vs_" + args.sat
                                    + "_for_" + args.reg
                                    + "_val_ts_lt"
                                    + "{:0>3d}".format(element)
                                    + "h_%Y%m.nc")
    valid_dict, dtime = get_arcmfc_stats(inpath + filename_stat)

    # Make ts-plots
    for val_name in valid_dict:
        filename_fig = fc_date.strftime(args.mod 
                                + "_vs_" + args.sat
                                + "_for_" + args.reg
                                + "_fig_val" 
                                + "_ts_" + val_name
                                + "_lt{:0>3d}".format(element)
                                + "h_%Y%m.png")
        ts = valid_dict[val_name]
        make_val_ts_fig_op(val_name,ts,dtime,filename_fig)

    # Get collocation ts
    inpath = ('/lustre/storeB/project/fou/om/waveverification/'
            + args.mod + '/satellites/altimetry'
            + '/' + args.sat + '/'
            + 'CollocationFiles/'
            + fc_date.strftime('%Y/%m/'))
    filename_coll = fc_date.strftime(args.mod
                                + "_vs_" + args.sat
                                + "_for_" + args.reg
                                + "_coll_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_coll)

    # Make scatter-plots
    filename_fig = fc_date.strftime(args.mod
                                + "_vs_" + args.sat
                                + "_for_" + args.reg
                                + "_fig_val_scatter"
                                + "_lt{:0>3d}".format(element)
                                + "h_%Y%m.png")
    make_val_scatter_fig_op(mHs,sHs,filename_fig)
    outpath = ('/lustre/storeB/project/fou/om/waveverification/'
           + args.mod + '/satellites/altimetry'
           + '/' + args.sat + '/'
           + 'ValidationFigures/'
           + fc_date.strftime('%Y/%m/'))
    cmd = 'mkdir -p ' + outpath
    os.system(cmd)
    cmd = 'mv ' + model + '_fig_val*.png ' + outpath
    os.system(cmd)
