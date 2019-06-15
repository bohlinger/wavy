#/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from datetime import datetime, timedelta
from copy import deepcopy
from graphicsmod import make_val_ts_fig_arcmfc, make_val_scatter_fig_arcmfc
from ncmod import get_arcmfc_stats, get_arcmfc_ts

# settings
fc_date = datetime.now()
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
    inpath=('/lustre/storeB/project/fou/om/ARCMFC/s3a/ValidationFiles/'
            + fc_date.strftime('%Y/%m/'))
    filename_stats = fc_date.strftime("ARCMFC_val_ts_lt"
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
    filename_fig = fc_date.strftime("ARCMFC_fig_val" 
                            + "_ts_" + val_name
                            + "_%Y%m.png")
    ts = valid_dict_lst[val_name]
    make_val_ts_fig_arcmfc(val_name,ts,dtime_lst,filename_fig,forecasts)

# Get collocation ts
dtime_lst = []
sHs_lst = []
mHs_lst = []
for element in forecasts:
    inpath=('/lustre/storeB/project/fou/om/ARCMFC/s3a/CollocationFiles/'
            + fc_date.strftime('%Y/%m/'))
    filename_coll = fc_date.strftime("ARCMFC_coll_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_coll)
    dtime_lst.append(dtime)
    sHs_lst.append(sHs)
    mHs_lst.append(mHs)

# Make scatter-plots
for i in range(len(forecasts)):
    filename_fig = fc_date.strftime("ARCMFC_fig_val_scatter_lt"
                            + "{:0>3d}".format(forecasts[i])
                            + "h_%Y%m.png")
    make_val_scatter_fig_arcmfc(mHs_lst[i],sHs_lst[i],filename_fig,forecasts,i)

# clean up
outpath=('/lustre/storeB/project/fou/om/ARCMFC/s3a/ValidationFigures/'
        + fc_date.strftime('%Y') + '/' + fc_date.strftime('%m') + '/')
cmd = 'mkdir -p ' + outpath
os.system(cmd)
cmd = 'mv ARCMFC_fig_val*.png ' + outpath
os.system(cmd)
