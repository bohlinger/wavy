#!/usr/bin/env python
import sys
sys.path.append(r'../wavy')

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
./op_figures.py -m mwam4
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-m", metavar='model',
    help="model to be evaluated")

args = parser.parse_args()

now = datetime.now()

if args.m is None:
    model = 'mwam4'
else:
    model = args.m

fc_date = datetime.now()
forecasts = [0]

for element in forecasts:
    # Get stats ts
#    inpath = ('/lustre/storeB/project/fou/om/waveverification/s3a/'
#         + model
#         + '/ValidationFiles/'
#         + fc_date.strftime('%Y/%m/'))
    inpath = ('/lustre/storeB/project/fou/om/waveverification/'
                + model
                + '/s3a/'
                + 'ValidationFiles/'
                + fc_date.strftime('%Y/%m/'))
    filename_stats = fc_date.strftime(model + "_val_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    valid_dict, dtime = get_arcmfc_stats(inpath + filename_stats)

    # Make ts-plots
    for val_name in valid_dict:
        filename_fig = fc_date.strftime(model + "_fig_val" 
                                + "_ts_" + val_name
                                + "_lt{:0>3d}".format(element)
                                + "h_%Y%m.png")
        ts = valid_dict[val_name]
        make_val_ts_fig_op(val_name,ts,dtime,filename_fig)

    # Get collocation ts
#    inpath = ('/lustre/storeB/project/fou/om/waveverification/s3a/'
#         + model
#         + '/CollocationFiles/'
#         + fc_date.strftime('%Y/%m/'))
    inpath = ('/lustre/storeB/project/fou/om/waveverification/'
            + model
            + '/s3a/'
            + 'CollocationFiles/'
            + fc_date.strftime('%Y/%m/'))
    filename_coll = fc_date.strftime(model + "_coll_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_coll)

    # Make scatter-plots
    filename_fig = fc_date.strftime(model + "_fig_val_scatter"
                            + "_lt{:0>3d}".format(element)
                            + "h_%Y%m.png")
    make_val_scatter_fig_op(mHs,sHs,filename_fig)
#    outpath = ('/lustre/storeB/project/fou/om/waveverification/s3a/'
#         + model
#         + '/ValidationFigures/'
#         + fc_date.strftime('%Y/%m/'))
    outpath = ('/lustre/storeB/project/fou/om/waveverification/'
           + model
           + '/s3a/'
           + 'ValidationFigures/'
           + fc_date.strftime('%Y/%m/'))
    cmd = 'mkdir -p ' + outpath
    os.system(cmd)
    cmd = 'mv ' + model + '_fig_val*.png ' + outpath
    os.system(cmd)
