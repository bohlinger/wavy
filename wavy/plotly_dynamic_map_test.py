import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
import yaml
import numpy as np
from datetime import datetime, timedelta
from modelmod import get_model
import plotly
import plotly.plotly as py
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

from datetime import *; from dateutil.relativedelta import *
import calendar

flatten = lambda l: [item for sublist in l for item in sublist]

# get data
from ncmod import get_nc_ts
sdate = datetime(2019,1,1)
tmpdate = datetime(2019,1,1)
#edate = datetime(2019,1,31)
edate = datetime(2019,2,28)

stationlst = buoy_dict.keys()
#stationlst = ['Tennholmen','Fauskane','SulaA','SulaB','SulaD','SulaF','SulaG']
stationvars = {}
for station in stationlst:
    print(station)
    stationvars[station]={}
    dtime_tmp=[]
    Hm0_model_tmp=[]
    Hs_station_1h_tmp=[]
    tmpdate = datetime(2019,1,1)
    while tmpdate<edate:
        try:
            vardict = get_nc_ts('/lustre/storeB/project/fou/om/waveverification/mwam4/buoys/CollocationFiles/' + tmpdate.strftime('%Y') + '/' + tmpdate.strftime('%m') + '/mwam4_' + station + '_coll_ts_' + tmpdate.strftime('%Y') + tmpdate.strftime('%m') + '_bestguess.nc',['Hm0_model','Hm0_buoy'])
            dtime_tmp.append(list(vardict['dtime']))
            Hm0_model_tmp.append(list(vardict['Hm0_model']))
            Hs_station_1h_tmp.append(list(vardict['Hm0_buoy']))
            tmpdate = tmpdate + relativedelta(months=+1)
        except:
            tmpdate = tmpdate + relativedelta(months=+1)
            pass
    vardict={
            'dtime':np.array(flatten(dtime_tmp)),
            'Hm0_buoy':np.array(flatten(Hs_station_1h_tmp)),
            'Hm0_model':np.array(flatten(Hm0_model_tmp))
            }
    stationvars[station]['sensor']=vardict

#print(stationvars)

#stationlst = station_dict.keys()
stationlst = ['draugen','ekofiskL','heimdal','snorrea','sleipner','norne','trolla','trollb','oseberg','gjoa','heidrun','grane','stafjorda','visund']
#stationlst = ['draugen','ekofiskL','snorrea','sleipner','trollb']
for station in (stationlst):
    stationvars[station]={}
    for sensor in (station_dict[station]['sensor']):
        print(station,sensor)
        dtime_tmp=[]
        Hm0_model_tmp=[]
        Hs_station_1h_tmp=[]
        tmpdate = datetime(2019,1,1)
        while tmpdate<edate:
            try:
#            vardict = get_nc_ts('/lustre/storeB/project/fou/om/waveverification/obs/stations/' + tmpdate.strftime('%Y') + '/' + tmpdate.strftime('%m') +'/' + station + '_' + sensor +'_' + tmpdate.strftime('%Y%m') + '.nc',['Hs'])
                filestr = ('/lustre/storeB/project/fou/om/'
                                + 'waveverification/mwam4/stations/'
                                + 'CollocationFiles/'
                                + tmpdate.strftime('%Y') + '/'
                                + tmpdate.strftime('%m') +'/'
                                + 'mwam4_'
                                + station + '_'
                                + sensor +'_coll_ts_'
                                + tmpdate.strftime('%Y%m')
                                + '_bestguess.nc')
                vardict = get_nc_ts(filestr,['Hm0_model','Hs_station_1h'])
                dtime_tmp.append(list(vardict['dtime']))
                Hm0_model_tmp.append(list(vardict['Hm0_model']))
                Hs_station_1h_tmp.append(list(vardict['Hs_station_1h']))
                tmpdate = tmpdate + relativedelta(months=+1)
            except:
                tmpdate = tmpdate + relativedelta(months=+1)
                pass
        vardict={
                'dtime':np.array(flatten(dtime_tmp)),
                'Hm0_model':np.array(flatten(Hm0_model_tmp)),
                'Hs_station':np.array(flatten(Hs_station_1h_tmp))
                }
        stationvars[station][sensor]=vardict

#locations = stationlst
#locations.append('Tennholmen')

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
                                'pitch': 45, 'zoom': 3.6,
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
        dcc.Graph(id='scatter-plot'),
    ], style={'display': 'inline-block', 'width': '49%'}),
    html.Div([
        dcc.Graph(id='y-time-series'),
    ], style={'display': 'inline-block', 'width': '99%'}),
])

def create_scatter(stationvars,varname1,varname2,stationstr):
    if (stationstr in station_dict and
    len(station_dict[stationstr]['sensor'].keys())>1):
        traces=[]
        ymax=[]
        for sensor in station_dict[stationstr]['sensor']:
            x=stationvars[stationstr][sensor]['dtime']
            y1=stationvars[stationstr][sensor][varname1]
            y2=stationvars[stationstr][sensor][varname2]
            ymax.append(np.nanmax(stationvars[stationstr][sensor][varname1]))
            traces.append(go.Scatter( x=y1, y=y2, name=sensor, mode = 'markers'))
        data = traces
        sensor='all sensors'
    else:
        ymax=[]
        for sensor in stationvars[stationstr].keys():
            pass
        x=stationvars[stationstr][sensor]['dtime']
        y1=stationvars[stationstr][sensor][varname1]
        y2=stationvars[stationstr][sensor][varname2]
        ymax.append(np.nanmax(y1))
        trace = go.Scatter( x=y1, y=y2, mode = 'markers')
        data = [trace]
    layout = dict(
            title = 'scatter plot for ' + stationstr + ' ' + sensor,
            xaxis = dict(
                hoverformat = '.2f',
                range=[0, (np.nanmax(ymax)+.5)],title='Hs from model [m]'
            ),
            yaxis = dict(
                hoverformat = '.2f',
                range=[0, (np.nanmax(ymax)+.5)],title='Hs from obs [m]'
            ),
            height=600
    )
    return {
        'data':data,
        'layout':layout
        }

@app.callback(
    dash.dependencies.Output('scatter-plot', 'figure'),
    [dash.dependencies.Input('crossfilter-indicator-scatter', 'hoverData'),])

def update_scatter(hoverData):
    #print(hoverData['points'][0])
    #{'pointNumber': 0, 'lon': 13.5618, 'curveNumber': 1, 'text': 'Tennholmen', 'lat': 67.3406, 'pointIndex': 0}
    try:
        stationstr = hoverData['points'][0]['customdata']
    except:
        stationstr = hoverData['points'][0]['text']
    if stationstr in buoy_dict.keys():
        varname1 = 'Hm0_model'
        varname2 = 'Hm0_buoy'
    else:
        varname1 = 'Hm0_model'
        varname2 = 'Hs_station'
#        Hs = stationvars[stationstr]['waverider']['Hm0_buoy']
#        dtime = stationvars[stationstr]['waverider']['dtime']
#    else:
#        Hs = stationvars[stationstr]['Hs']
#        dtime = stationvars[stationstr]['dtime']
    return create_scatter(stationvars,varname1,varname2,stationstr)

def create_time_series(stationvars,varname1,varname2,stationstr):
    if (stationstr in station_dict and 
    len(station_dict[stationstr]['sensor'].keys())>1):
        traces=[]
        ymax=[]
        for sensor in station_dict[stationstr]['sensor']:
            x=stationvars[stationstr][sensor]['dtime']
            y1=stationvars[stationstr][sensor][varname1]
            y2=stationvars[stationstr][sensor][varname2]
            ymax.append(np.nanmax(stationvars[stationstr][sensor][varname1]))
            traces.append(go.Scatter( x=x, y=y2, name=sensor))
        traces.append(go.Scatter( x=x, y=y1, name='mwam4'))
        data = traces
        sensor='all sensors'
    else:
        traces=[]
        ymax=[]
        for sensor in stationvars[stationstr].keys():
            pass
        if stationstr in buoy_dict.keys():
            sensor='sensor'
        x=stationvars[stationstr][sensor]['dtime']
        y1=stationvars[stationstr][sensor][varname1]
        y2=stationvars[stationstr][sensor][varname2]
        ymax.append(np.nanmax(y1))
        traces.append(go.Scatter( x=x, y=y1, name='mwam4'))
        traces.append(go.Scatter( x=x, y=y2, name=sensor))
        data = traces
    layout = dict(
            title = 'Hs from ' + stationstr + ' ' + sensor,
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
                range=[0, (np.nanmax(ymax)+.5)],
                title='Hs [m]'
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
    if stationstr in buoy_dict.keys():
        varname1 = 'Hm0_model'
        varname2 = 'Hm0_buoy'
    else:
        varname1 = 'Hm0_model'
        varname2 = 'Hs_station'
#        Hs = stationvars[stationstr]['waverider']['Hm0_buoy']
#        dtime = stationvars[stationstr]['waverider']['dtime']
#    else:
#        Hs = stationvars[stationstr]['Hs']
#        dtime = stationvars[stationstr]['dtime']
    return create_time_series(stationvars,varname1,varname2,stationstr)


#app.css.config.serve_locally = True
#app.scripts.config.serve_locally = True

if __name__ == '__main__':
    #paste this into the browser: http://xvis-m3a:8050/
#    app.run_server(debug = True,host = '0.0.0.0')
    app.run_server(debug = False,host = '0.0.0.0')
