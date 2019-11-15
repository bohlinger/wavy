import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
from datetime import datetime, timedelta
import numpy as np
from plotly_graphs import ts_comp_fig, ts_comp_figs
from ncmod import get_nc_ts
import yaml
from modelmod import get_model, get_latest_output_init_date
from stationmod import matchtime, get_loc_idx

# define additional fct
def get_fc_vals(ts_mod,model):
    fc_date = ts_mod['dtime'][-1]
    init_date = get_latest_output_init_date(model)
    leadtime = int((fc_date-init_date).seconds/60/60)
    ts_fc_time = []
    ts_fc_var = []
    for i in range(1,13):
        try:
            model_var,model_lats,model_lons,model_time,model_time_dt = get_model(simmode="fc",model=model,fc_date=fc_date + timedelta(hours=i),init_date=init_date,leadtime=leadtime + i,varname='Hs')
        except:
            print("file not yet not available")
            print("using previous run")
            model_var,model_lats,model_lons,model_time,model_time_dt = get_model(simmode="fc",model=model,fc_date=fc_date + timedelta(hours=i),init_date=init_date-timedelta(hours=6),leadtime=leadtime + i,varname='Hs')
        # collocate with wave model
        idx, idy, distM, picked_lat, picked_lon = get_loc_idx(\
                model_lats,model_lons,\
                station_dict[station]['coords']['lat'],\
                station_dict[station]['coords']['lon'],\
                mask=None)
        ts_fc_time.append( fc_date + timedelta(hours=i) )
        ts_fc_var.append( model_var.squeeze()[idx[0],idy[0]] )
    return ts_fc_time, ts_fc_var

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
        try:
#        for k in range(1):
            ts_obs = get_nc_ts(fobs,['Hs_10min'])
            xobs_lst.append(ts_obs['dtime'])
            yobs_lst.append(ts_obs['Hs_10min'])
            sensor_lst.append(sensor)
            ts_mod = get_nc_ts(fmod,['Hs'])
            ts_mod_fc_time, ts_mod_fc_vals = get_fc_vals(ts_mod,model)
            ts_mod_time = np.array(list(ts_mod['dtime'])+ts_mod_fc_time)
            ts_mod_vals = np.array(list(ts_mod['Hs'])+ts_mod_fc_vals)
            ts_comp_figs(ts_mod_time,ts_mod_vals,xobs_lst,yobs_lst,
            'Hs',station,sensor_lst,model)
        except:
            print "Unexpected error:", sys.exc_info()[0]
