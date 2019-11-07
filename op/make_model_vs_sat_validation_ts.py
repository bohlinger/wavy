import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
from datetime import datetime, timedelta
import numpy as np
from plotly_graphs import ts_comp_fig, ts_comp_figs
from ncmod import get_nc_ts
import yaml

model = 'mwam4'

now = datetime.now()
#now = datetime(2019,10,1)

with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

# stations
for station in station_dict.keys():
    print(station)
    xobs_lst = []
    yobs_lst = []
    sensor_lst = []
    fmod = ('/lustre/storeB/project/fou/om/waveverification/'
        + model + '/stations/CollocationFiles/'
        + now.strftime('%Y') + '/'
        + now.strftime('%m') + '/'
        + model + '_Hs_at_' + station + '_ts_'
        + now.strftime('%Y%m')
        + '_bestguess.nc')
    for sensor in station_dict[station]['sensor']:
        print(sensor)
        fobs = ('/lustre/storeB/project/fou/om/waveverification/'
            + 'obs/stations/'
            + now.strftime('%Y') + '/'
            + now.strftime('%m') + '/'
            + 'Hs_10min_' 
            + station + '_' + sensor + '_'
            + now.strftime('%Y%m') + '.nc')
        ts_obs = get_nc_ts(fobs,['Hs_10min'])
        xobs_lst.append(ts_obs['dtime'])
        yobs_lst.append(ts_obs['Hs_10min'])
        sensor_lst.append(sensor)
    ts_mod = get_nc_ts(fmod,['Hs'])
    ts_comp_figs(ts_mod['dtime'],ts_mod['Hs'],xobs_lst,yobs_lst,
            'Hs',station,sensor_lst,model)
