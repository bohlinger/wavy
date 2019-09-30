import numpy as np
import pandas as pd
from stationmod import matchtime, get_loc_idx
from modelmod import get_model, check_date
from datetime import datetime, timedelta
from calendar import monthrange
import netCDF4

df=pd.read_pickle("./locations_mwam4_clean.pkl")
idx = df['idx'].values.astype('int')
idy = df['idy'].values.astype('int')
stationid_tmp = df['station:id'].values.astype('str')
stationid = np.zeros([len(stationid_tmp)]).astype('str')
for i in range(len(stationid)):
    stationid[i]=stationid_tmp[i][8::]
stationid = np.array(stationid, dtype='object')

#ltlst = np.arange(0,144,24)
ltlst = np.arange(0,72,24)
varlst = ['Hs','Tp', 'Tz', 'thq', 'u10', 'v10']

basedate = datetime(1970,1,1)
sdate = datetime(2019,6,1)
tmpdate = datetime(2019,6,1)
#edate = datetime(2019,6,monthrange(2019,6)[1])
edate = sdate

model = 'mwam4'

HsA=np.zeros([1,len(ltlst),99])*np.nan
TpA=np.zeros([1,len(ltlst),99])*np.nan
TzA=np.zeros([1,len(ltlst),99])*np.nan
thqA=np.zeros([1,len(ltlst),99])*np.nan
u10A=np.zeros([1,len(ltlst),99])*np.nan
v10A=np.zeros([1,len(ltlst),99])*np.nan
time_lst = []
# get wave model
lt = 0
substs=np.ones([99])*99.99

while tmpdate <= edate:
    print(tmpdate)
    for i in range(len(ltlst)):
        lt = ltlst[i]
        print(lt)
        try:
            fc_date = tmpdate
            init_date = tmpdate - timedelta(hours=lt)
            check_date(model,fc_date=fc_date,leadtime=lt)
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt,varname='Hs')
            HsA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt,varname='Tp')
            TpA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt,varname='Tz')
            TzA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt,varname='thq')
            thqA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model='mwam4force',fc_date=fc_date,\
                init_date=init_date,leadtime=lt,varname='u10')
            u10A[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model='mwam4force',fc_date=fc_date,\
                init_date=init_date,leadtime=lt,varname='v10')
            v10A[0,i,:]=model_var.squeeze()[idx,idy]
        except IOError:
            fc_date = tmpdate
            init_date = tmpdate - timedelta(hours=lt+6)
            check_date(model,fc_date=fc_date,leadtime=lt)
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt+6,varname='Hs')
            HsA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt+6,varname='Tp')
            TpA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt+6,varname='Tz')
            TzA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=fc_date,\
                init_date=init_date,leadtime=lt+6,varname='thq')
            thqA[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model='mwam4force',fc_date=fc_date,\
                init_date=init_date,leadtime=lt+6,varname='u10')
            u10A[0,i,:]=model_var.squeeze()[idx,idy]
            model_var,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model='mwam4force',fc_date=fc_date,\
                init_date=init_date,leadtime=lt+6,varname='v10')
            v10A[0,i,:]=model_var.squeeze()[idx,idy]
    time_lst.append((tmpdate-basedate).total_seconds())
    tmpdate = tmpdate + timedelta(days=1)

vardict = {
            'basetime':basedate,
            'time':time_lst,
            'lt':ltlst,
            'stationid':stationid,
            'Hs':HsA,
            'Tp':TpA,
            'Tz':TzA,
            'thq':thqA,
            'u10':u10A,
            'v10':v10A,
}

from ncmod import dumptonc_LCWVF
title = 'Collocated variables for LC-WVF from the operational mywavewam 4km at the Norwegian Meteorological Institute'
filename = 'test_LCWVF_METNO.nc'
outpath = 'out/'
dumptonc_LCWVF(outpath,filename,title,basedate,vardict)
