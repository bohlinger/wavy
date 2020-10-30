#!/usr/bin/env python

import sys
import os
from datetime import datetime, timedelta
from copy import deepcopy
import argparse
from argparse import RawTextHelpFormatter
import yaml

sys.path.append('../../wavy')

from graphicsmod import make_val_ts_fig_op, make_val_scatter_fig_op
from ncmod import get_arcmfc_stats, get_arcmfc_ts

# parser
parser = argparse.ArgumentParser(
    description="""
Create validation figures for validation files based on satellite altimetry.
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
parser.add_argument("-d", metavar='date',
    help="month to be plotted fmt: %Y%m")
parser.add_argument("-path", metavar='outpath',
    help="path to where files are to be stored")

args = parser.parse_args()

varlst = ['Hs']
now = datetime.now()

if args.d is None:
    fc_date = datetime.now()
else:
    fc_date = datetime(int(args.d[0:4]),int(args.d[4:6]),1)

leadtimes = [0,24,48]

if args.mod is None:
    args.mod = 'mwam4'

if args.sat is None:
    args.sat = 's3a'

if args.reg is None:
    args.reg = args.mod

if args.path is None:
    args.path = '/lustre/storeB/project/fou/om/waveverification/'

# make a list of validation metrics for various lead times
rmsd_lst = []
bias_lst = []
corr_lst = []
SI_lst = []
nov_lst = []
dtime_lst = []
for element in leadtimes:
    # Get stats ts
    inpath = (args.path
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
    rmsd_lst.append(valid_dict['rmsd'])
    bias_lst.append(valid_dict['bias'])
    corr_lst.append(valid_dict['corr'])
    SI_lst.append(valid_dict['SI'])
    nov_lst.append(valid_dict['nov'])
    dtime_lst.append(dtime)

valid_dict_lst = {'rmsd':rmsd_lst,
                  'bias':bias_lst,
                  'corr':corr_lst,
                  'SI':SI_lst,
                  'nov':nov_lst}

# Make ts-plots
for val_name in valid_dict_lst:
    filename_fig = fc_date.strftime(args.mod 
                            + "_vs_" + args.sat
                            + "_for_" + args.reg
                            + "_fig_val" 
                            + "_ts_" + val_name
                            + "_%Y%m.png")
    ts = valid_dict_lst[val_name]
    make_val_ts_fig_op(val_name,ts,dtime_lst,filename_fig,leadtimes)

# Get collocation ts
dtime_lst = []
sHs_lst = []
mHs_lst = []
for element in leadtimes:
    inpath = (args.path
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
    dtime_lst.append(dtime)
    sHs_lst.append(sHs)
    mHs_lst.append(mHs)

# Make scatter-plots
for i in range(len(leadtimes)):
    filename_fig = fc_date.strftime(args.mod
                                + "_vs_" + args.sat
                                + "_for_" + args.reg
                                + "_fig_val_scatter"
                                + "_lt{:0>3d}".format(leadtimes[i])
                                + "h_%Y%m.png")
    make_val_scatter_fig_op(mHs_lst[i],sHs_lst[i],filename_fig,leadtimes,i)

# clean up
outpath = (args.path
           + args.mod + '/satellites/altimetry'
           + '/' + args.sat + '/'
           + 'ValidationFigures/'
           + fc_date.strftime('%Y/%m/'))
cmd = 'mkdir -p ' + outpath
os.system(cmd)
cmd = 'mv ' + args.mod + '*_fig_val*.png ' + outpath
os.system(cmd)
