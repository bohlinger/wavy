#https://plotly.com/python/plotly-express/
#https://plotly.com/python/dropdowns/
#https://plotly.com/python/subplots/
#https://towardsdatascience.com/how-to-create-a-plotly-visualization-and-embed-it-on-websites-517c1a78568b

# --- imports --- #
# general
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import yaml
import pandas as pd
# plotly
from plotly.subplots import make_subplots
#import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

# --- settings --- #
varstr = 'Hs'
modstr = 'mwam4'

# --- get specs --- #
with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream) 
with open("/home/patrikb/wavy/wavy/station_specs.yaml", 'r') as stream:
    station_dict=yaml.safe_load(stream)

# --- gather data and make figures --- #
fig = make_subplots(rows=len(station_dict.keys()), cols=2)
statlst = list(station_dict.keys())
for i in range(len(statlst)):
    senslst = list(station_dict[statlst[i]]['sensor'].keys())
    for j in range(len(senslst)):
        # model
        pstr = ("/lustre/storeB/project/fou/om/waveverification/"
            + modstr + "/stations/CollocationFiles/2020/03/")
        fstr = (modstr + "_" + varstr + "_at_" 
            + statlst[i] + "_ts_lt_best_202003.nc")
        nc = netCDF4.Dataset(pstr+fstr)
        varmod = nc.variables[varstr][:]
        time = nc.variables['time'][:]
        unit = nc.variables['time'].units
        ptmod = netCDF4.num2date(time,units = unit)
        nc.close()
        dfmod = pd.DataFrame({'time':ptmod,varstr:varmod})

        # obs
        pstr = ("/lustre/storeB/project/fou/om/"
            + "waveverification/obs/stations/2020/03/")
        fstr = ("Hs_10min_" + statlst[i] + "_" + senslst[j] + "_202003.nc")
        nc = netCDF4.Dataset(pstr+fstr)
        varobs = nc.variables[varstr + '_10min'][:]
        time = nc.variables['time'][:]
        unit = nc.variables['time'].units
        ptobs = netCDF4.num2date(time,units = unit)
        nc.close()
        dfobs = pd.DataFrame({'time':ptobs,varstr:varobs})

        fig.add_trace(
            px.scatter(dfmod, x="time", y=varstr),
            row=i+1, col=1
        )
        fig.add_trace(
            px.scatter(dfobs, x="time", y=varstr),
            row=i+1, col=1
        )
        dfscat = pd.DataFrame({'obs':varmod,'mod':varmod})
        fig.add_trace(
            px.scatter(dfscat, x = 'obs', y = 'mod'),
            row=i+1, col=2
        )

        fig.update_layout(height=200, width=800, title_text="Subplots")

pio.write_html(fig, file='index.html', auto_open=True)
