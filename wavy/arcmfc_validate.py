#!/usr/bin/env python
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')

from datetime import datetime, timedelta
from satmod import sentinel_altimeter as sa
from stationmod import station_class as sc
from stationmod import matchtime, get_model
from modelmod import get_model, collocate
from satmod import validate
from copy import deepcopy

startdate = datetime(2018,11,1,0)
tmpdate = deepcopy(startdate)
enddate = datetime(2018,11,2,0)

while tmpdate <= enddate:
    # settings
    fc_date = tmpdate
    region = 'ARCMFC'
    model = 'ARCMFC'
    # get model collocated values
    results_dict = 

    try:
        #print results_dict
        valid_dict=validate(results_dict)
        from model_specs import model_dict
        # dump to nc-file: validation
        from custom_nc import dumptonc_stats
        basetime=model_dict[model]['basetime']
        outpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFiles/'
        filename_ts=tmpdate.strftime("ARCMFC_val_%Y%m.nc")
        title_stat='validation file'
        time_dt = fc_date
        dumptonc_stats(outpath,filename_stat,title_stat,basetime,time_dt,valid_dict)
    except:
        print('no matches')
        pass
    tmpdate = tmpdate + timedelta(hours=1)
