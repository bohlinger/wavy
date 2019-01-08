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
from model_specs import model_dict
from custom_nc import dumptonc_stats

sdate = datetime(2018,11,1,0)
tmpdate = deepcopy(sdate)
edate = datetime(2018,11,15,0)

forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]

while tmpdate <= edate:
    print(tmpdate)
    for element in forecasts:
        # settings
        fc_date = deepcopy(tmpdate)
        region = 'ARCMFC'
        model = 'ARCMFC'
        basetime=model_dict[model]['basetime']
        # get model collocated values
        from custom_nc import get_arcmfc_ts
        inpath = '/lustre/storeB/project/fou/om/ARCMFC/S3a/CollocationFiles/'
        filename_ts=fc_date.strftime("ARCMFC_coll_ts_lt"
                                            + "{:0>3d}".format(element)
                                            + "h_%Y%m.nc")
        dtime, sHs, mHs = get_arcmfc_ts(inpath + filename_ts)
        del filename_ts
        # find collocations for given model time step and validate
        from stationmod import matchtime
        time_lst = []
        for dt in dtime:
            time_lst.append((dt-basetime).total_seconds())
        ctime,idx = matchtime(tmpdate, tmpdate, time_lst, basetime, timewin=30)
        if len(idx)==0:
            pass
        else:
            results_dict = {'date_matches':dtime[idx],
                        'model_Hs_matches':mHs[idx],
                        'sat_Hs_matches':sHs[idx]}
            valid_dict=validate(results_dict)
            print(valid_dict)
            # dump to nc-file: validation
            outpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFiles/'
            title_stat='validation file'
            filename_stat=fc_date.strftime("ARCMFC_val_ts_lt"
                                        + "{:0>3d}".format(element)
                                        + "h_%Y%m.nc")
            time_dt = fc_date
            dumptonc_stats(outpath,filename_stat,title_stat,basetime,time_dt,valid_dict)
    tmpdate = tmpdate + timedelta(hours=6)
