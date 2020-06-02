#/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from datetime import datetime, timedelta
from copy import deepcopy
from graphicsmod import make_val_ts_fig_arcmfc, make_val_scatter_fig_arcmfc
from ncmod import get_arcmfc_stats, get_arcmfc_ts
import argparse
from argparse import RawTextHelpFormatter

# parser
parser = argparse.ArgumentParser(
    description="""
Collocate wave model output and s3a data and dump to monthly nc-file.
If file exists, data is appended.

Usage:
./arcmfc3_figures.py -mod ARCMFC3 -sat s3a -reg ARCMFC3
./arcmfc3_figures.py -mod ARCMFC3 -sat s3a -reg NordicSeas -d 202001
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument("-mod", metavar='model',
    help="model to be used for collocation")
parser.add_argument("-sat", metavar='satellite',
    help="satellite mission to be used for collocation")
parser.add_argument("-reg", metavar='region',
    help="region of interest")
parser.add_argument("-d", metavar='date',
    help="month to be plotted fmt: %Y%m")

args = parser.parse_args()

now = datetime.now()

if args.d is None:
    fc_date = datetime.now()
else:
    fc_date = datetime(int(args.d[0:4]),int(args.d[4:6]),1)


# settings
sat = args.sat
model = args.mod
region = args.reg
#fc_date = datetime.now()
forecasts = [12, 36, 60]
val_names = ['rmsd','bias','corr','SI','nov']

# Get stats ts
dtime_lst = []
rmsd_lst = []
bias_lst = []
corr_lst = []
SI_lst = []
nov_lst = []
for element in forecasts:

    inpath=('/lustre/storeB/project/fou/om/waveverification/'
                    + model + '/satellites/altimetry/'
                    + sat + '/ValidationFiles/'
                    + fc_date.strftime('%Y/%m/'))

    filename_stats = fc_date.strftime(  model + "_vs_"
                                        + sat
                                        + '_for_'
                                        + region
                                        + "_val_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
    valid_dict, dtime = get_arcmfc_stats(inpath + filename_stats)
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
for val_name in val_names:
    filename_fig = fc_date.strftime(
                            "ARCMFC3_for_" 
                            + region
#                            "ARCMFC3"
                            + "_fig_val" 
                            + "_ts_" + val_name
                            + "_%Y%m.png")
    ts = valid_dict_lst[val_name]
    make_val_ts_fig_arcmfc(val_name,ts,dtime_lst,filename_fig,forecasts)

# Get collocation ts
dtime_lst = []
sHs_lst = []
mHs_lst = []
for element in forecasts:

    inpath=('/lustre/storeB/project/fou/om/'
                    + 'waveverification/'
                    + model
                    + '/satellites/altimetry/'
                    + sat
                    + '/CollocationFiles/'
                    + fc_date.strftime("%Y/%m/"))

    filename_coll = fc_date.strftime(model + "_vs_"
                                    + sat
                                    + '_for_'
                                    + region
                                    + "_coll_ts_lt"
                                    + "{:0>3d}".format(element)
                                    + "h_%Y%m.nc")

    dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_coll)
    dtime_lst.append(dtime)
    sHs_lst.append(sHs)
    mHs_lst.append(mHs)

# Make scatter-plots
for i in range(len(forecasts)):
    filename_fig = fc_date.strftime(
                            "ARCMFC3_for_" + region
#                            "ARCMFC3"
                            + "_fig_val_scatter_lt"
                            + "{:0>3d}".format(forecasts[i])
                            + "h_%Y%m.png")
    make_val_scatter_fig_arcmfc(mHs_lst[i],sHs_lst[i],filename_fig,forecasts,i)

# clean up
outpath=('/lustre/storeB/project/fou/om/waveverification/'
        + model + '/satellites/altimetry/' + sat + '/ValidationFigures/'
        + fc_date.strftime('%Y') + '/' + fc_date.strftime('%m') + '/')
cmd = 'mkdir -p ' + outpath
os.system(cmd)
#cmd = 'mv ARCMFC3_fig_val*.png ' + outpath
cmd = 'mv ARCMFC3*_fig_val*.png ' + outpath
os.system(cmd)
print(cmd)
