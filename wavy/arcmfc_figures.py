#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

import os
from datetime import datetime, timedelta
from satmod import validate
from copy import deepcopy
from graphicsmod import make_val_ts_fig_arcmfc
from custom_nc import get_arcmfc_stats

fc_date = datetime(2018,11,1,0)

#forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]
forecasts = [12]

for element in forecasts:
    # settings
    region = 'ARCMFC'
    model = 'ARCMFC'
    inpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFiles/'
    filename_stats = fc_date.strftime("ARCMFC_val_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
    valid_dict, dtime = get_arcmfc_stats(inpath + filename_stats)
    # Make ts-plots
    for val_name in valid_dict:
        filename_fig = fc_date.strftime("ARCMFC_fig_val" 
                                + "_" + val_name
                                + "_lt{:0>3d}".format(element)
                                + "h_%Y%m.png")
        ts = valid_dict[val_name]
        make_val_ts_fig_arcmfc(val_name,ts,dtime,filename_fig)
        outpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFigures/'
        cmd = 'cp ARCMFC_fig_val*.png ' + outpath
        os.system(cmd)
    # Make scatter-plots
