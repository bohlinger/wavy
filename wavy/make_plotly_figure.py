import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
import numpy as np
from datetime import datetime, timedelta
from modelmod import get_model
import plotly
import plotly.plotly as py
import pandas as pd
from region_specs import poly_dict
from matplotlib.patches import Polygon

def ggt(a, b):
    if a < b: a,b = b,a
    while a%b != 0:
        a,b = b,a%b
    return b

def make_plotly_figure(sa_obj=None,\
                       region=None,domain=None,proj=None,\
                       model=None,grid_date=None,tiling=None):
    if grid_date is None:
        if model == 'mwam8':
            grid_date = datetime(2019,2,1,6)
        elif (model == 'MoskNC' or model == 'MoskWC'):
            grid_date = datetime(2018,3,1)
        else:
            grid_date = datetime(2019,2,1)
    if region is None and sa_obj is None:
        region == 'NorwegianSea'
    if region is None and model is not None:
        region = model
    if domain is None:
        domain = 'Global'

    if model is not None:
        # get model grid
        if (model == 'ARCMFC' or model == 'MoskNC' or model=='MoskWC'):
            model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=grid_date,
                init_date=grid_date)
        else:
            model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                get_model(simmode="fc",model=model,fc_date=grid_date,
                leadtime=0)
        if tiling is not None:
            # tiling of model domain
            in1,in2 = model_lats.shape[0],model_lats.shape[1]
            if isinstance(tiling,int):
                g = tiling
            else:
                if ggt(in1,in2)==1:
                    g = ggt(in1-1,in2-1)
                else:
                    g = ggt(in1,in2)
            print('tile size is: ' + str(g))
            xtilesize = model_lats.shape[0]/g
            ytilesize = model_lats.shape[1]/g
            xidx = []
            for i in range(g):
                xidx.append(xtilesize*i)        
            xidx.append(model_lats.shape[0])
            yidx = []
            for i in range(g):
                yidx.append(ytilesize*i)
            yidx.append(model_lats.shape[1])
            print(xidx)
            print(yidx)
            tiles = []
            for i in range(g):
                for j in range(g):
                    tiles.append([model_lons[xidx[i]:xidx[i+1],yidx[j]:yidx[j+1]],model_lats[xidx[i]:xidx[i+1],yidx[j]:yidx[j+1]]])

    if(model is not None and region is not model):
        # get only edge indices
        lontmp1=list(model_lons[:,0])
        lontmp2=list(model_lons[:,-1])
        lontmp3=list(model_lons[0,:])
        lontmp4=list(model_lons[-1,:])
        mlons = lontmp1 + lontmp2 + lontmp3 + lontmp4
        lattmp1=list(model_lats[:,0])
        lattmp2=list(model_lats[:,-1])
        lattmp3=list(model_lats[0,:])
        lattmp4=list(model_lats[-1,:])
        mlats = lattmp1 + lattmp2 + lattmp3 + lattmp4
        pdmlons = pd.Series(mlons)
        pdmlats = pd.Series(mlats)

    if region == model:
        if tiling is not None:
            mlats_lst = []
            mlons_lst = []
            for i in range(len(tiles)):
                model_lons = tiles[i][0]
                model_lats = tiles[i][1]
                # create polygon of model domain
                idx = 1
                lontmp1=list(model_lons[:,0])
                lontmp2=list(model_lons[:,-1])
                lontmp3=list(model_lons[0,:])
                lontmp4=list(model_lons[-1,:])
                mlons = (lontmp3[::-1][::idx] + lontmp1[::idx] \
                    + lontmp4[::idx] + lontmp2[::-1][::idx] \
                    + [lontmp3[::-1][0]] + [lontmp3[::-1][0]])
                lattmp1=list(model_lats[:,0])
                lattmp2=list(model_lats[:,-1])
                lattmp3=list(model_lats[0,:])
                lattmp4=list(model_lats[-1,:])
                mlats = (lattmp3[::-1][::idx] + lattmp1[::idx] \
                    + lattmp4[::idx] + lattmp2[::-1][::idx] \
                    + [lattmp3[::-1][0]] + [lattmp3[::-1][0]])
                mlats_lst = mlats_lst + mlats
                mlons_lst = mlons_lst + mlons
            pdlons = pd.Series(mlons_lst)
            pdlats = pd.Series(mlats_lst)
            poly = Polygon(list(zip(pdlons,pdlats)), closed=False)
            pdlons = pd.Series(poly.xy[:,0])
            pdlats = pd.Series(poly.xy[:,1])
            # interactive text
            names = []
            for i in range(len(pdlons)):
                names.append('polygon node: E'
                           + '{:0.2f}'.format(poly.xy[i,0])
                           + ', N'
                           + '{:0.2f}'.format(poly.xy[i,1]))
            pdnames = pd.Series(names)
        else:
            # create polygon of model domain
            idx = 10
            lontmp1=list(model_lons[:,0])
            lontmp2=list(model_lons[:,-1])
            lontmp3=list(model_lons[0,:])
            lontmp4=list(model_lons[-1,:])
            mlons = (lontmp3[::-1][::idx] + lontmp1[::idx] \
                + lontmp4[::idx] + lontmp2[::-1][::idx] \
                + [lontmp3[::-1][0]] + [lontmp3[::-1][0]])
            lattmp1=list(model_lats[:,0])
            lattmp2=list(model_lats[:,-1])
            lattmp3=list(model_lats[0,:])
            lattmp4=list(model_lats[-1,:])
            mlats = (lattmp3[::-1][::idx] + lattmp1[::idx] \
                + lattmp4[::idx] + lattmp2[::-1][::idx] \
                + [lattmp3[::-1][0]] + [lattmp3[::-1][0]])
            pdlons = pd.Series(mlons)
            pdlats = pd.Series(mlats)
            poly = Polygon(list(zip(pdlons,pdlats)), closed=False)
            pdlons = pd.Series(poly.xy[:,0])
            pdlats = pd.Series(poly.xy[:,1])
            # interactive text
            names = []
            for i in range(len(pdlons)):
                names.append('polygon node: E' 
                       + '{:0.2f}'.format(poly.xy[i,0]) 
                       + ', N' 
                       + '{:0.2f}'.format(poly.xy[i,1]))
            pdnames = pd.Series(names)
    elif (region != model and model is not None):
        # get region
        poly = Polygon(list(zip(poly_dict[region]['lons'],\
                                poly_dict[region]['lats'])),\
                                closed=True)
        pdlons = pd.Series(poly.xy[:,0])
        pdlats = pd.Series(poly.xy[:,1])
        # interactive text for polygon
        names = []
        for i in range(len(pdlons)):
            names.append('polygon node: E' 
                       + '{:0.2f}'.format(poly.xy[i,0]) 
                       + ', N' 
                       + '{:0.2f}'.format(poly.xy[i,1]))
        pdnames = pd.Series(names)
        # for model
        mnames=[]
        for i in range(len(pdmlons)):
            mnames.append('model domain: E' 
                        + '{:0.2f}'.format(pdmlons[i]) 
                        + ', N' 
                        + '{:0.2f}'.format(pdmlats[i]))
            pdmnames = pd.Series(mnames)
    elif (region != model and model is None):
        # get region
        poly = Polygon(list(zip(poly_dict[region]['lons'],\
                                poly_dict[region]['lats'])),\
                                closed=True)
        pdlons = pd.Series(poly.xy[:,0])
        pdlats = pd.Series(poly.xy[:,1])
        # interactive text for polygon
        names = []
        for i in range(len(pdlons)):
            names.append('polygon node: E' 
                       + '{:0.2f}'.format(poly.xy[i,0]) 
                       + ', N' 
                       + '{:0.2f}'.format(poly.xy[i,1]))
        pdnames = pd.Series(names)

    # compute center of figure
    def angular_mean(angles_deg):
        N = len(angles_deg)
        mean_c = 1.0 / N * np.sum(np.exp(1j * angles_deg * np.pi/180.0))
        return np.angle(mean_c, deg=True)

    if (model is not None and model is not region):
        latmean = (angular_mean(pdmlats) + angular_mean(pdlats))/2
        lonmean = (angular_mean(pdmlons) + angular_mean(pdlons))/2
        # compute range of figure
        coefflon = 0.1
        coefflat = 0.04
        lonrange = [ np.min(list(mlons)+list(pdlons)) - \
                coefflon * np.min(list(mlons)+list(pdlons)), \
                np.max(list(mlons)+list(pdlons)) + \
                coefflon * np.max(list(mlons)+list(pdlons))]
        latrange = [ np.min(list(mlats)+list(pdlats)) - \
                coefflat * np.min(list(mlats)+list(pdlats)), \
                np.max(list(mlats)+list(pdlats)) + \
                coefflat * np.max(list(mlats)+list(pdlats))]
    else:
        latmean = angular_mean(pdlats)
        lonmean = angular_mean(pdlons)
        # compute range of figure
        coefflon = 0.1
        coefflat = 0.04
        lonrange = [ np.min(list(pdlons)) - \
                coefflon * np.min(list(pdlons)), \
                np.max(list(pdlons)) + \
                coefflon * np.max(list(pdlons))]
        latrange = [ np.min(list(pdlats)) - \
                coefflat * np.min(list(pdlats)), \
                np.max(list(pdlats)) + \
                coefflat * np.max(list(pdlats))]

    # region polygon
    print("adding polygon nodes to map")
    airports = [ dict(
            type = 'scattergeo',
            lon = pdlons,
            lat = pdlats,
            hoverinfo = 'text',
            text = pdnames,
            mode = 'markers',
            marker = dict(
                size=2,
                color='red',
                line = dict(
                    width=.5,
                    color='rgba(68, 68, 68, 0)'
                )
            ))]
    data = airports

    if model is not region:    
#    if model is not None:    
        # connecting polygon nodes
        print("connecting polygon nodes")
        flight_paths = []
        for i in range( len( poly.xy[:,0] ) - 1):
            flight_paths.append(
                dict(
                    type = 'scattergeo',
                    lon = [ poly.xy[:,0][i], poly.xy[:,0][i+1] ],
                    lat = [ poly.xy[:,1][i], poly.xy[:,1][i+1] ],
                    mode = 'lines',
                    line = dict(
                        width = 2,
                        color = 'red',
                    ),
                    opacity = float(0.8),
                )
            )
        data = data + flight_paths

    if (model is not None and region is not model):
        print("adding model boundary points to map")
        # model grid
        mgrid = [ dict(
            type = 'scattergeo',
            lon = pdmlons,
            lat = pdmlats,
            hoverinfo = 'text',
            text = pdmnames,
            mode = 'markers',
            marker = dict(
                size=2,
                color='black',
                line = dict(
                    width=3,
                    color='rgba(68, 68, 68, 0)'
                )
            ))]
        data = data + mgrid

    # adding S3a Hs data
    sanames = []
    if sa_obj is not None:
        print("adding S3a hovering legend")
        thin = 1
        if len(sa_obj.rHs) > 5000:
            thin = 2
        if len(sa_obj.rHs) > 15000:
            thin = 3
        if len(sa_obj.rHs) > 25000:
            thin = 5
        if len(sa_obj.rHs) > 60000:
            thin = 10
        for i in range(len(sa_obj.rHs[::thin])):
            sanames.append('Hs: ' 
                        + '{:0.2f}'.format(sa_obj.rHs[::thin][i])
                        + ' (' 
                        + sa_obj.rtime[::thin][i].strftime('%Y-%m-%d %H:%M:%S')
                        + ')')
            pdsanames = pd.Series(sanames)
        # model grid
        sa_grid = [ dict(
            type = 'scattergeo',
            lon = sa_obj.rloc[1][::thin],
            lat = sa_obj.rloc[0][::thin],
            hoverinfo = 'text',
            text = pdsanames,
            mode = 'markers',
#            marker = dict(
#                size=5,
#                color=sa_obj.rHs,
#                line = dict(
#                    width=3,
#                    color='rgba(68, 68, 68, 0)'
#                )
#            )
            marker = dict(
                color = sa_obj.rHs,
                colorscale = 'Portland',
                reversescale = False,
                opacity = 0.9,
                size = 5,
                colorbar = dict(
                    thickness = 10,
                    titleside = "right",
                    outlinecolor = "rgba(68, 68, 68, 0)",
                    ticks = "outside",
                    ticklen = 3,
                    showticksuffix = "last",
                    ticksuffix = " [m]",
                    dtick = 1
                )
            )
            )]
        data = data + sa_grid

    # make map
    print("make map")
    titlestr = (sa_obj.sdate.strftime('%Y-%m-%d %H:%M:%S UTC')
                + ' - '
                + sa_obj.edate.strftime('%Y-%m-%d %H:%M:%S UTC')
                + '<br>'
                + model 
                + ' domain (black) and chosen region ' 
                + region 
                + ' (red)'
                + '<br>'
                + 'displaying ' + str(len(sanames)) + ' footprints'
                + ' (applied thinning: ' + str(thin) + ')')
    if (proj=='ortho' and domain=='auto'):
        layout = dict(
                title = titlestr,
                showlegend = False,
                geo = dict(
                    resolution = 50,
                    showland = True,
                    showcountries = True,
                    showocean = True,
                    countrywidth = 0.5,
                    landcolor = 'rgb(243, 243, 243)',
                    oceancolor = 'lightblue',
                    projection = dict(
                        type = 'orthographic',
                        rotation = dict(
                            lon = lonmean,
                            lat = latmean,
                            roll = 0,
                        )
                    ),
                    lataxis = dict(
                        range = latrange,
                        showgrid = True,
                        dtick = 5
                    ),
                    lonaxis = dict(
                        range = lonrange,
                        showgrid = True,
                        dtick = 5
                    ),
                )
            )
    elif (proj==None and (domain == 'Global' or domain == None)):
        layout = dict(
            title = titlestr,
            showlegend = False,
            geo = dict(
                scope='world',
                projection=dict( type='kavrayskiy7' ),
                showland = True,
                showocean = True,
                landcolor = 'rgb(243, 243, 243)',
                countrycolor = 'rgb(204, 204, 204)',
                oceancolor = 'lightblue',
                lataxis = dict(
                    showgrid = True,
                    dtick = 5
                ),
                lonaxis = dict(
                    showgrid = True,
                    dtick = 5
                ),
            )
        )
    elif (proj=='ortho' and (domain == 'Global' or domain == None)):
        layout = dict(
            title = titlestr,
            showlegend = False,
            geo = dict(
                scope='world',
                projection=dict( 
                        type='orthographic',
                        rotation = dict(
                            lon = 0,
                            lat = 90,
                            roll = 0,
                        )
                    ),
                showland = True,
                showocean = True,
                landcolor = 'rgb(243, 243, 243)',
                countrycolor = 'rgb(204, 204, 204)',
                oceancolor = 'lightblue',
                lataxis = dict(
                    showgrid = True,
                    dtick = 5
                ),
                lonaxis = dict(
                    showgrid = True,
                    dtick = 5
                ),
            )
        )
    fig = dict( data=data, layout=layout )
    print("save figure")
    plotly.offline.plot( fig, filename=('/lustre/storeB/project/fou/om/'
                                       +'waveverification/d3-flight-paths')
                                       , auto_open=False)
