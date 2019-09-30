import numpy as np
import pandas as pd
from stationmod import matchtime, get_loc_idx
from modelmod import get_model, check_date
from datetime import datetime, timedelta

model = 'mwam4'
#lt = np.arange(0,144,24)
init_date = datetime(2019,1,1)
fc_date = datetime(2019,1,1)

# get valid locations
df = pd.read_table("SWH_all_201906.txt",delimiter="\s+")

lats = df['latitude']
lons = df['longitude']
statid = df['station:id']

# get wave model
check_date(model,fc_date=fc_date,leadtime=0)
model_var,model_lats,model_lons,model_time,model_time_dt = \
            get_model(simmode="fc",model=model,fc_date=fc_date,\
            init_date=init_date,leadtime=0,varname='Hs')

# collocate with wave model
idx_lst = []
idy_lst = []
for i in range(len(lats)):
    print('Position ' + str(i) + ' of ' + str(len(lats)))
    idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
            model_lats,model_lons,\
            lats[i],\
            lons[i],\
            mask=None)
    print('Collocation distance: ' + str(distM[idx,idy]))
    if distM[idx,idy]<4:
        print('Add location')
        idx_lst.append(idx.item())
        idy_lst.append(idy.item())
    else: 
        print('Reject location')
        idx_lst.append(np.nan)
        idy_lst.append(np.nan)

df['idx'] = idx_lst 
df['idy'] = idy_lst 

df.to_pickle('location_pickle.pkl')
