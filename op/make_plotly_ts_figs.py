import sys
sys.path.append(r'../wavy')
# buoys
from datetime import datetime
from stationmod import get_buoy
from buoy_specs import buoy_dict
from datetime import datetime
import numpy as np
basetime = buoy_dict['Tennholmen']['basetime']
sdate = datetime(2019,1,1)
edate = datetime.now()
time_s, time_dt, Hm0, Tm02, lons, lats = get_buoy(sdate,edate)
# plot
from plotly_graphs import ts_figs
# boys
ts_figs(time_dt,Hm0,'Hm0','buoys','Tennholmen','waverider')
ts_figs(time_dt,Tm02,'Tm02','buoys','Tennholmen','waverider')

# stations
from stationmod import station_class as sc
from datetime import datetime
statname = 'ekofiskL'
sensorname = 'waverider'
mode = 'd22'
deltat = 10
sdate = datetime(2019,1,1,0)
edate = datetime.now()
#edate = datetime(2019,6,18,0)
sc_obj = sc(statname,sdate,edate,mode=mode,sensorname=sensorname,deltat=deltat)
time_dt = []
Hm0 = []
for i in range(len(sc_obj.hs)):
    try:
        if sc_obj.hs[i]>30:
            Hm0.append(np.nan)
        else:
            Hm0.append(sc_obj.hs[i])
        time_dt.append(sc_obj.timedt[i])
    except:
        pass
# plot
from plotly_graphs import ts_figs
# stations
ts_figs(time_dt,Hm0,'Hm0','stations',statname,sensorname)
#ts_figs(time_dt,Tm02,'Tm02','stations',statname,sensorname)
