import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
from datetime import datetime, timedelta
import numpy as np
from plotly_graphs import ts_comp_figs
from ncmod import get_nc_ts

# stations

station = 'trollb'
sensor = 'MKIIIradar'

fobs = ('/lustre/storeB/project/fou/om/waveverification/'
        + 'obs/stations/2019/11/Hs_10min_trollb_MKIIIradar_201911.nc')
fmod = ('/lustre/storeB/project/fou/om/waveverification/'
        + 'mwam4/stations/CollocationFiles/2019/11/'
        + 'mwam4_Hs_at_trollb_ts_201911_bestguess.nc')

ts_obs = get_nc_ts(fobs,['Hs_10min'])
ts_mod = get_nc_ts(fmod,['Hs'])

ts_comp_figs(ts_obs['dtime'],ts_obs['Hs_10min'],
            ts_mod['dtime'],ts_mod['Hs'],
            'Hs',station,sensor)
