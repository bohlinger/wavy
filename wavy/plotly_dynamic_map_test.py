import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
import numpy as np
from datetime import datetime, timedelta
from modelmod import get_model
import plotly
import plotly.plotly as py
from region_specs import poly_dict
from matplotlib.patches import Polygon
from copy import deepcopy
import os

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go

from stationlist_arcmfc import stationlist_arcmfc
from station_specs import station_dict
from buoy_specs import buoy_dict

# get data
from ncmod import get_nc_ts
stationvars = {}
vardict = get_nc_ts('/lustre/storeB/project/fou/om/waveverification/mwam4/buoys/CollocationFiles/2019/06/mwam4_Tennholmen_coll_ts_lt000h_201906.nc',['Hm0_model','Hm0_buoy'])
stationvars['Tennholmen']={}
stationvars['Tennholmen']['waverider']=vardict
#stationlst = station_dict.keys()
stationlst = ['draugen','ekofiskL','heimdal','snorrea','sleipner','norne','trolla','trollb','oseberg','gjoa','heidrun']
for station in (stationlst):
    stationvars[station]={}
    for sensor in (station_dict[station]['sensor']):

        vardict = get_nc_ts('/lustre/storeB/project/fou/om/waveverification/obs/stations/2019/01/' + station + '_' + sensor +'_201901.nc',['Hs'])
        stationvars[station][sensor]=vardict

locations = stationlst
locations.append('Tennholmen')

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# create figure

station_lats = []
station_lons= [] 
station_hovertxt = []
for element in stationlist_arcmfc:
    station_lats.append(station_dict[element]['coords']['lat'])
    station_lons.append(station_dict[element]['coords']['lon'])
    station_hovertxt.append(element)
buoy_lats = []
buoy_lons = []
buoy_hovertxt = []
for element in buoy_dict:
    buoy_lats.append(buoy_dict[element]['lat'])
    buoy_lons.append(buoy_dict[element]['lon'])
    buoy_hovertxt.append(element)
mapbox_access_token = "pk.eyJ1IjoicHJpeWF0aGFyc2FuIiwiYSI6ImNqbGRyMGQ5YTBhcmkzcXF6YWZldnVvZXoifQ.sN7gyyHTIq1BSfHQRBZdHA"
traces=[]
traces.append(go.Scattermapbox(
    lat=station_lats,
    lon=station_lons, mode='markers',
    hoverinfo='text', 
    marker={'symbol': "circle",
                'size': 8,
                'opacity':.8
                 },
    text=station_hovertxt,
    name = 'offshore'
    ))
traces.append(go.Scattermapbox(
    lat=buoy_lats, 
    lon=buoy_lons, mode='markers',
    hoverinfo='text', 
    marker={'symbol': "circle",
                'size': 8,
                'opacity':.8
                 },
    text=buoy_hovertxt,
    name = 'near-shore'
    ))
layout1 = go.Layout(title='Locations of available wave observations',           
                    autosize=True, hovermode='closest',
                    showlegend=True,
                    mapbox={'accesstoken':
                                mapbox_access_token,
                                'bearing': 0,
                                'center': {'lat': 61, 'lon': 7},
                                'pitch': 45, 'zoom': 4,
                                "style": 'mapbox://styles/mapbox/light-v9'
                        },
                    height=600,
                    )
app.layout = html.Div([
    html.Div([
        dcc.Graph(
            id='crossfilter-indicator-scatter',
            hoverData={'points': [{'customdata': 'Tennholmen'}]},
            figure={
                'data': traces,
                'layout':layout1
            }
        )
    ], style={'width': '49%','height':'99%', 'display': 'inline-block'}),
    html.Div([
        dcc.Graph(id='y-time-series'),
    ], style={'display': 'inline-block', 'width': '49%', 'height':'99%'}),
])

def create_time_series(stationvars,varname,stationstr):
    if (stationstr in station_dict and 
    len(station_dict[stationstr]['sensor'].keys())>1):
        traces=[]
        for sensor in station_dict[stationstr]['sensor']:
            x=stationvars[stationstr][sensor]['dtime']
            y=stationvars[stationstr][sensor][varname]
            traces.append(go.Scatter( x=x, y=y, name=sensor))
        data = traces
    else:
        for sensor in stationvars[stationstr].keys():
            pass
        x=stationvars[stationstr][sensor]['dtime']
        y=stationvars[stationstr][sensor][varname]
        trace = go.Scatter( x=x, y=y)
        data = [trace]
    layout = dict(
            title = varname + ' from ' + stationstr + ' ' + sensor,
            xaxis = dict(
                rangeselector=dict(
                    buttons=list([
                        dict(count=1,
                             label='1d',
                             step='day',
                             stepmode='backward'),
                        dict(count=1,
                             label='1m',
                             step='month',
                             stepmode='backward'),
                        dict(count=1,
                            label='YTD',
                            step='year',
                            stepmode='todate'),
                        dict(count=1,
                            label='1y',
                            step='year',
                            stepmode='backward'),
                        dict(step='all')
                    ])
                ),
                rangeslider=dict(
                    visible = True
                ),
                type='date',
                tickformat = '%Y-%m-%d UTC%H:%M',
            ),
            yaxis = dict(
                hoverformat = '.2f',
                range=[0, np.nanmax(y)]
            ),
            height=600
    )
    return {
        'data':data,
        'layout':layout
        }

@app.callback(
    dash.dependencies.Output('y-time-series', 'figure'),
    [dash.dependencies.Input('crossfilter-indicator-scatter', 'hoverData'),])

def update_timeseries(hoverData):
    #print(hoverData['points'][0])
    #{'pointNumber': 0, 'lon': 13.5618, 'curveNumber': 1, 'text': 'Tennholmen', 'lat': 67.3406, 'pointIndex': 0}
    try:
        stationstr = hoverData['points'][0]['customdata']
    except:
        stationstr = hoverData['points'][0]['text']
    if stationstr == 'Tennholmen':
        varname = 'Hm0_buoy'
    else:
        varname = 'Hs'
#        Hs = stationvars[stationstr]['waverider']['Hm0_buoy']
#        dtime = stationvars[stationstr]['waverider']['dtime']
#    else:
#        Hs = stationvars[stationstr]['Hs']
#        dtime = stationvars[stationstr]['dtime']
    return create_time_series(stationvars,varname,stationstr)


#app.css.config.serve_locally = True
#app.scripts.config.serve_locally = True

if __name__ == '__main__':
    #paste this into the browser: http://xvis-m3a:8050/
#    app.run_server(debug = True,host = '0.0.0.0')
    app.run_server(debug = False,host = '0.0.0.0')
