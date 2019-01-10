#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from datetime import datetime, timedelta
from satmod import validate
from copy import deepcopy
from graphicsmod import make_val_ts_fig_arcmfc, make_val_scatter_fig_arcmfc
from custom_nc import get_arcmfc_stats, get_arcmfc_ts

#fc_date = datetime(2018,11,1,0)
#forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]
fc_date = datetime.now()
forecasts = [12]

for element in forecasts:
    # Get stats ts
    inpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFiles/'
    filename_stats = fc_date.strftime("ARCMFC_val_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    valid_dict, dtime = get_arcmfc_stats(inpath + filename_stats)

    # Make ts-plots
    for val_name in valid_dict:
        filename_fig = fc_date.strftime("ARCMFC_fig_val" 
                                + "_ts_" + val_name
                                + "_lt{:0>3d}".format(element)
                                + "h_%Y%m.png")
        ts = valid_dict[val_name]
        make_val_ts_fig_arcmfc(val_name,ts,dtime,filename_fig)

    # Get collocation ts
    inpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/CollocationFiles/'
    filename_coll = fc_date.strftime("ARCMFC_coll_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_coll)

    # Make scatter-plots
    filename_fig = fc_date.strftime("ARCMFC_fig_val_scatter"
                            + "_lt{:0>3d}".format(element)
                            + "h_%Y%m.png")
    make_val_scatter_fig_arcmfc(mHs,sHs,filename_fig)

    outpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFigures/'
    cmd = 'mv ARCMFC_fig_val*.png ' + outpath
    os.system(cmd)
