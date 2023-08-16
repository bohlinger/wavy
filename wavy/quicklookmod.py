"""
    Module for quicklook fct
"""
# imports
import numpy as np
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from abc import abstractmethod

# own imports
from wavy.wconfig import load_or_default

# read yaml config files:
region_dict = load_or_default('region_specs.yaml')
variable_info = load_or_default('variable_info.yaml')
model_dict = load_or_default('model_specs.yaml')
quicklook_dict = load_or_default('quicklook_specs.yaml')

# define global functions

def comp_fig(sa_obj=None,mc_obj=None,coll_obj=None,**kwargs):
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    from cartopy.mpl.ticker import LatitudeLocator, LongitudeLocator

    # sort out data/coordinates for plotting
    sat = "NA"
    model = "NA"

    """
    If sa_obj is not None get satellite_altimetry data for plotting
    """
    if sa_obj is not None:
        slons, slats = sa_obj.vars['longitude'],sa_obj.vars['latitude']
        svar = sa_obj.vars[sa_obj.stdvarname]
        stdvarname = sa_obj.stdvarname
        sat = sa_obj.mission
    """
    If mc_obj is not None get model data for plotting
    """
    if mc_obj is not None:
        mlons = mc_obj.vars['longitude']
        mlats = mc_obj.vars['latitude']
        mvar = mc_obj.vars[mc_obj.stdvarname]
        # inflate coords if regular lat/lon grid
        if (len(mlons.shape)==1):
            mlons, mlats = np.meshgrid(mlons, mlats)
        stdvarname = mc_obj.stdvarname
        model = mc_obj.model

    """
    If sa_obj is not None get satellite_altimetry data for plotting
    """
    if sa_obj is None:
        slons, slats = coll_obj.vars['obs_lons'],coll_obj.vars['obs_lats']
        svar = coll_obj.vars['obs_values']

    """
    Get all misc in **kwargs
    """

    """
    Prepare plotting
    """

    # determine region bounds
    if sa_obj is not None:
        if sa_obj.region in region_dict['rect']:
            latmin = region_dict['rect'][sa_obj.region]['llcrnrlat']
            latmax = region_dict['rect'][sa_obj.region]['urcrnrlat']
            lonmin = region_dict['rect'][sa_obj.region]['llcrnrlon']
            lonmax = region_dict['rect'][sa_obj.region]['urcrnrlon']
        elif sa_obj.region in region_dict['geojson']:
            latmin = np.min(sa_obj.vars['latitude']) - .5
            latmax = np.max(sa_obj.vars['latitude']) + .5
            lonmin = np.min(sa_obj.vars['longitude']) - .5
            lonmax = np.max(sa_obj.vars['longitude']) + .5
        elif sa_obj.region in region_dict['poly']:
            latmin = np.min(region_dict['poly'][sa_obj.region]['lats'])-.5
            latmax = np.max(region_dict['poly'][sa_obj.region]['lats'])+.5
            lonmin = np.min(region_dict['poly'][sa_obj.region]['lons'])-.5
            lonmax = np.max(region_dict['poly'][sa_obj.region]['lons'])+.5
        elif sa_obj.region in model_dict:
            # model bounds
            latmin = np.min(mlats)
            latmax = np.max(mlats)
            lonmin = np.min(mlons)
            lonmax = np.max(mlons)
        else: print("Error: Region not defined!")
    elif (sa_obj is None and mc_obj is not None):
            # model bounds
            latmin = np.min(mlats)
            latmax = np.max(mlats)
            lonmin = np.min(mlons)
            lonmax = np.max(mlons)

    # determine projection
    """
    Here, a routine is needed to determine a suitable projection.
    As for now, there is Mercator as default.
    """
    projection_default = ccrs.Mercator(
                         central_longitude=(lonmin+lonmax)/2.,
                         min_latitude=latmin, max_latitude=latmax,
                         globe=None,
                         latitude_true_scale=(latmin+latmax)/2.,
                         false_easting=0.0, false_northing=0.0,
                         scale_factor=None)
    projection = kwargs.get('projection',projection_default)

    land = cfeature.GSHHSFeature(scale='i', levels=[1],
                    facecolor=cfeature.COLORS['land'])

    # make figure
    fig, ax = plt.subplots(nrows=1, ncols=1,
                        subplot_kw=dict(projection=projection),
                        figsize=(9, 9))
    # plot domain extent
    ax.set_extent([lonmin, lonmax,latmin, latmax],crs = ccrs.PlateCarree())

    # plot model domain if model is available
    if mc_obj is not None:
        ax.plot(mlons[0,:], mlats[0,:], '-', transform= ccrs.PlateCarree(),
                color = 'gray', linewidth =2)
        ax.plot(mlons[-1,:], mlats[-1,:], '-', transform= ccrs.PlateCarree(),
                color = 'gray', linewidth =2)
        ax.plot(mlons[:,0], mlats[:,0], '-', transform= ccrs.PlateCarree(),
                color = 'gray', linewidth =2)
        ax.plot(mlons[:,-1], mlats[:,-1], '-', transform= ccrs.PlateCarree(),
                color = 'gray', linewidth =2)

    # plot polygon if defined
    if sa_obj.region in region_dict['poly']:
        ax.plot(region_dict['poly'][sa_obj.region]['lons'],
                region_dict['poly'][sa_obj.region]['lats'],
                '-', transform= ccrs.PlateCarree(),
                color = 'gray', linewidth =2)
    # plot polygon from geojson if defined
    if sa_obj.region in region_dict['geojson']:
        import geojson
        import matplotlib.patches as patches
        from matplotlib.patches import Polygon
        fstr = region_dict['geojson'][sa_obj.region]['fstr']
        with open(fstr) as f:
            gj = geojson.load(f)
        fidx = region_dict['geojson'][sa_obj.region].get('fidx')
        if fidx is not None:
            geo = {'type': 'Polygon',
                   'coordinates':\
                    gj['features'][fidx]['geometry']['coordinates'][0]}
            poly = Polygon([tuple(l)
                            for l in geo['coordinates'][0]],
                            closed=True)
            ax.add_patch(\
                patches.Polygon(\
                    poly.get_path().to_polygons()[0],
                    transform=ccrs.PlateCarree(),
                    facecolor='None',
                    edgecolor='red',
                    lw = 3,
                    alpha=0.2))
        else:
            for i in range(len(gj['features'])):
                geo = {'type': 'Polygon',
                        'coordinates':\
                        gj['features'][i]['geometry']['coordinates'][0]}
                poly = Polygon([tuple(l)
                                for l in geo['coordinates'][0]],
                                closed=True)
                ax.add_patch(\
                    patches.Polygon(\
                        poly.get_path().to_polygons()[0],
                        transform=ccrs.PlateCarree(),
                        facecolor='None',
                        edgecolor='red',
                        lw = 3,
                        alpha=0.2))

    # colors
    if stdvarname == 'sea_surface_wave_significant_height':
        cmap = cmocean.cm.amp
        levels = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                3,3.25,3.5,3.75,4,4.5,5,6,7,8,9,10,11,12,13,14,15,
                16,17,18]
    elif stdvarname == 'wind_speed':
        cmap = cmocean.cm.amp
        levels = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,
                36,40,44,48]
    if 'cmap' in kwargs.keys():
        cmap = kwargs['cmap']
    if 'levels' in kwargs.keys():
        levels = kwargs['levels']
    if 'scl' in kwargs.keys():
        scl = kwargs['scl']
    else: scl = 18
    if 'icl' in kwargs.keys():
        icl = kwargs['icl']
    else: icl = 1

    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    extend = 'neither'
    if 'extend' in kwargs.keys():
        extend = kwargs['extend']

    if sa_obj.region in quicklook_dict:
        if ('cm_levels' in \
        quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]\
        and quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]['cm_levels'] is not None):
            levels = quicklook_dict[sa_obj.region]['varspec']\
                    [sa_obj.varalias]['cm_levels']
        if ('scl' in \
        quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]\
        and quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]['scl'] is not None):
            scl = quicklook_dict[sa_obj.region]['varspec']\
                    [sa_obj.varalias]['scl']
        if ('icl' in \
        quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]\
        and quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]['icl'] is not None):
            scl = quicklook_dict[sa_obj.region]['varspec']\
                    [sa_obj.varalias]['icl']

    # draw figure features
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fs = 11

    # plot lats/lons
    gridcolor = 'gray'
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(),
                linewidth = 1,color = gridcolor, alpha = 0.4,
                linestyle = '-')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fs, 'color': gridcolor}
    gl.ylabel_style = {'size': fs, 'color': gridcolor}
    #gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.xlocator = LongitudeLocator()
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()

    # - model contours
    im = ax.contourf(mlons, mlats, mvar, levels = levels,
                    transform = ccrs.PlateCarree(),
                    cmap = cmocean.cm.amp, norm = norm, extend = extend)
    imc = ax.contour(mlons, mlats, mvar, levels = levels[scl::icl],
                    transform = ccrs.PlateCarree(),
                    colors='w', linewidths = 0.3)
    ax.clabel(imc, fmt='%2d', colors='w', fontsize=fs)

    # - add coastline
    ax.add_geometries(land.intersecting_geometries(
                    [lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)

    # - add land color
    ax.add_feature( land, facecolor = 'burlywood', alpha = 0.5 )

    # - add satellite
    if len(slats)>0:
        sc = ax.scatter(slons,slats,s=10,
                c=svar,
                marker='o', edgecolor = 'face',
                cmap=cmocean.cm.amp, norm = norm,
                transform=ccrs.PlateCarree())

    # - point of interests depending on region
    if ('quicklook_dict' in globals()
    and sa_obj.region in quicklook_dict
    and 'poi' in quicklook_dict[sa_obj.region]):
        for poi in quicklook_dict[sa_obj.region]['poi']:
            pname = quicklook_dict[sa_obj.region]['poi'][poi]['name']
            plat = quicklook_dict[sa_obj.region]['poi'][poi]['lat']
            plon = quicklook_dict[sa_obj.region]['poi'][poi]['lon']
            scp = ax.scatter(plon,plat,s=20, c='b',
                    marker = quicklook_dict[sa_obj.region]['poi'][poi]['marker'],
                    transform = ccrs.PlateCarree())
            ax.text(plon,plat,pname,transform = ccrs.PlateCarree())

    # - colorbar
    axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
    cbar = fig.colorbar(im, cax=axins)
    cbar.ax.set_ylabel( stdvarname + ' [' +
                        variable_info[sa_obj.varalias]['units']
                        + ']',size=fs)
    cbar.ax.tick_params(labelsize=fs)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

    # - title
    ax.set_title(model + ' model time step: '
        + coll_obj.vars['valid_date'][0].strftime("%Y-%m-%d %H%M:%S UTC")
        + '\n'
        + sat
        + ' coverage \n from '
        + coll_obj.vars['datetime'][0].strftime("%Y-%m-%d %H:%M:%S UTC" )
        + ' to '
        + coll_obj.vars['datetime'][-1].strftime("%Y-%m-%d %H:%M:%S UTC")
        ,fontsize=fs)

    # - save figure
    if ('savepath' in kwargs.keys() and kwargs['savepath'] != None):
        plt.savefig( kwargs['savepath'] + '/' + model
                + '_vs_satellite_'
                + coll_obj.vars['valid_date'][0].strftime("%Y%m%d")
                + 'T'
                + coll_obj.vars['valid_date'][0].strftime("%H")
                + 'Z.png', format = 'png', dpi=200)

    # - show figure
    if ('showfig' in kwargs.keys() and kwargs['showfig'] == True):
        plt.show()

def plot_sat(sa_obj,**kwargs):

    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker

    stdvarname = sa_obj.stdvarname
    # sort out data/coordinates for plotting
    slons, slats = sa_obj.vars['longitude'],sa_obj.vars['latitude']
    if 'col_obj' in kwargs.keys():
        model_lats = kwargs['col_obj'].vars['model_lats']
        model_lons = kwargs['col_obj'].vars['model_lons']
        model_vals = kwargs['col_obj'].vars['model_values']
    if (sa_obj.region in model_dict and 'mc_obj' in kwargs.keys()):
        grid_date = model_dict[sa_obj.region]['grid_date']
        model_var_dict = kwargs['mc_obj'].vars
    # check region and determine projection
    if (sa_obj.region == 'global'
        or (sa_obj.region in region_dict['rect']
            and 'boundinglat' in region_dict['rect'][sa_obj.region].keys())
        ):
        # Polar Stereographic Projection
        polarproj = ccrs.NorthPolarStereo(
                            central_longitude=0.0,
                            true_scale_latitude=66,
                            globe=None)
        projection = polarproj
        land = cfeature.GSHHSFeature(scale='i', levels=[1],
                        facecolor=cfeature.COLORS['land'])
    else:
        if sa_obj.region in region_dict['rect']:
            latmin = region_dict['rect'][sa_obj.region]['llcrnrlat']
            latmax = region_dict['rect'][sa_obj.region]['urcrnrlat']
            lonmin = region_dict['rect'][sa_obj.region]['llcrnrlon']
            lonmax = region_dict['rect'][sa_obj.region]['urcrnrlon']
        elif sa_obj.region in region_dict['geojson']:
            latmin = np.min(sa_obj.vars['latitude']) - .5
            latmax = np.max(sa_obj.vars['latitude']) + .5
            lonmin = np.min(sa_obj.vars['longitude']) -.5
            lonmax = np.max(sa_obj.vars['longitude']) + .5
        elif sa_obj.region in region_dict['poly']:
            latmin = np.min(region_dict['poly'][sa_obj.region]['lats'])-.5
            latmax = np.max(region_dict['poly'][sa_obj.region]['lats'])+.5
            lonmin = np.min(region_dict['poly'][sa_obj.region]['lons'])-.5
            lonmax = np.max(region_dict['poly'][sa_obj.region]['lons'])+.5
        elif sa_obj.region in model_dict:
            # model bounds
            latmin = np.min(model_var_dict['latitude'])
            latmax = np.max(model_var_dict['latitude'])
            lonmin = np.min(model_var_dict['longitude'])
            lonmax = np.max(model_var_dict['longitude'])
        else: print("Error: Region not defined!")
        projection = ccrs.Mercator(
                        central_longitude=(lonmin+lonmax)/2.,
                        min_latitude=latmin, max_latitude=latmax,
                        globe=None, latitude_true_scale=(latmin+latmax)/2.,
                        false_easting=0.0, false_northing=0.0,
                        scale_factor=None)
        land = cfeature.GSHHSFeature(scale='i', levels=[1],
                        facecolor=cfeature.COLORS['land'])

    # make figure
    fig, ax = plt.subplots(nrows=1, ncols=1,
                        subplot_kw=dict(projection=projection),
                        figsize=(9, 9))
    # plot domain extent
    if 'polarproj' not in locals():
        ax.set_extent([lonmin, lonmax,latmin, latmax],crs = ccrs.PlateCarree())
    else:
        ax.set_extent([-180, 180,40, 90],crs = ccrs.PlateCarree())

    # plot model domain if region is a model domain
    if (sa_obj.region in model_dict
    and len(model_var_dict['latitude'].shape)==1):
        lenlons = len(model_var_dict['longitude'][:])
        lenlats = len(model_var_dict['latitude'][:])
        ax.plot([model_var_dict['longitude'][0]]*lenlats,
                model_var_dict['latitude'][:], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['longitude'][:],
                [model_var_dict['latitude'][-1]]*lenlons, '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot([model_var_dict['longitude'][-1]]*lenlats,
                model_var_dict['latitude'][::-1], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['longitude'][::-1],
                [model_var_dict['latitude'][0]]*lenlons, '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
    if (sa_obj.region in model_dict
    and len(model_var_dict['latitude'].shape)==2):
        ax.plot(model_var_dict['longitude'][0,:],
                model_var_dict['latitude'][0,:], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['longitude'][-1,:],
                model_var_dict['latitude'][-1,:], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['longitude'][:,0],
                model_var_dict['latitude'][:,0], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['longitude'][:,-1],
                model_var_dict['latitude'][:,-1], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)

    # plot polygon if defined
    if sa_obj.region in region_dict['poly']:
        ax.plot(region_dict['poly'][sa_obj.region]['lons'],
                region_dict['poly'][sa_obj.region]['lats'],
                '-', transform= ccrs.PlateCarree(),
                color = 'gray', linewidth =2)
    # colors
    cmap = cmocean.cm.amp
    if stdvarname == 'sea_surface_wave_significant_height':
        cmap = cmocean.cm.amp
        levels = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                3,3.25,3.5,3.75,4,4.5,5,6,7,8,9,10,11,12,13,14,15,
                16,17,18]
    elif stdvarname == 'wind_speed':
        cmap = cmocean.cm.amp
        levels = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,
                36,40,44,48]
    if 'cmap' in kwargs.keys():
        cmap = kwargs['cmap']
    if 'levels' in kwargs.keys():
        levels = kwargs['levels']

    # draw figure features
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fs = 11

    if sa_obj.region in quicklook_dict:
        if ('cm_levels' in \
        quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]\
        and quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]['cm_levels'] is not None):
            levels = quicklook_dict[sa_obj.region]['varspec']\
                    [sa_obj.varalias]['cm_levels']
        if ('scl' in \
        quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]\
        and quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]['scl'] is not None):
            scl = quicklook_dict[sa_obj.region]['varspec']\
                    [sa_obj.varalias]['scl']
        if ('icl' in \
        quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]\
        and quicklook_dict[sa_obj.region]['varspec'][sa_obj.varalias]['icl'] is not None):
            scl = quicklook_dict[sa_obj.region]['varspec']\
                    [sa_obj.varalias]['icl']

    # plot lats/lons
    gridcolor = 'gray'
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(),
                linewidth = 1,color = gridcolor, alpha = 0.4,
                linestyle = '-')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': fs, 'color': gridcolor}
    gl.ylabel_style = {'size': fs, 'color': gridcolor}

    # - add coastline
    if 'polarproj' not in locals():
        ax.add_geometries(land.intersecting_geometries(
                    [lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)
    else:
        ax.add_geometries(land.intersecting_geometries(
                    [-180, 180, 0, 90]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)

    # - add land color
    ax.add_feature( land, facecolor = 'burlywood', alpha = 0.5 )

    # scale colormap
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    extend = 'neither'
    if 'extend' in kwargs.keys():
        extend = kwargs['extend']

    # - add satellite
    sc = ax.scatter(slons,slats,s=10,
                #c='k',#sa_obj.vars[sa_obj.stdvarname],
                c = sa_obj.vars[sa_obj.stdvarname],
                marker='o', edgecolor = 'face',
                cmap=cmocean.cm.amp, norm = norm,
                transform=ccrs.PlateCarree())
    if 'col_obj' in kwargs.keys():
        # - add satellite
        sc2 = ax.scatter(model_lons,model_lats,s=10,
                #c='r',#model_vals,
                c = model_vals,
                marker='o', edgecolor = 'face',
                cmap=cmocean.cm.amp, norm = norm,
                transform=ccrs.PlateCarree())

    # - point of interests depending on region
    if ('quicklook_dict' in globals()
    and sa_obj.region in quicklook_dict
    and 'poi' in quicklook_dict[sa_obj.region]):
        for poi in quicklook_dict[sa_obj.region]['poi']:
            pname = quicklook_dict[sa_obj.region]['poi'][poi]['name']
            plat = quicklook_dict[sa_obj.region]['poi'][poi]['lat']
            plon = quicklook_dict[sa_obj.region]['poi'][poi]['lon']
            scp = ax.scatter(plon,plat,s=20,
                    c = quicklook_dict[sa_obj.region]['poi'][poi].get('color','b'),
                    marker = quicklook_dict[sa_obj.region]['poi'][poi]['marker'],
                    transform = ccrs.PlateCarree())
            ax.text(plon,plat,pname,transform = ccrs.PlateCarree())

    # - plot polygon
    if sa_obj.region in region_dict['poly']:
        ax.plot(
            region_dict['poly'][sa_obj.region]['lons'],
            region_dict['poly'][sa_obj.region]['lats'],
            'k:',transform=ccrs.PlateCarree()
            )

    # plot polygon from geojson if defined
    if sa_obj.region in region_dict['geojson']:
        import geojson
        import matplotlib.patches as patches
        from matplotlib.patches import Polygon
        fstr = region_dict['geojson'][sa_obj.region]['fstr']
        with open(fstr) as f:
            gj = geojson.load(f)
        fidx = region_dict['geojson'][sa_obj.region].get('fidx')
        if fidx is not None:
            geo = {'type': 'Polygon',
                   'coordinates':\
                    gj['features'][fidx]['geometry']['coordinates'][0]}
            poly = Polygon([tuple(l)
                            for l in geo['coordinates'][0]],
                            closed=True)
            ax.add_patch(\
                patches.Polygon(\
                    poly.get_path().to_polygons()[0],
                    transform=ccrs.PlateCarree(),
                    facecolor='None',
                    edgecolor='red',
                    lw = 3,
                    alpha=0.2))
        else:
            for i in range(len(gj['features'])):
                geo = {'type': 'Polygon',
                        'coordinates':\
                        gj['features'][i]['geometry']['coordinates'][0]}
                poly = Polygon([tuple(l)
                                for l in geo['coordinates'][0]],
                                closed=True)
                ax.add_patch(\
                    patches.Polygon(\
                        poly.get_path().to_polygons()[0],
                        transform=ccrs.PlateCarree(),
                        facecolor='None',
                        edgecolor='red',
                        lw = 3,
                        alpha=0.2))

    # - colorbar
    axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
    cbar = fig.colorbar(sc, cax=axins)
    cbar.ax.set_ylabel(sa_obj.stdvarname + ' [' +
                        variable_info[sa_obj.varalias]['units']
                        + ']')
    cbar.ax.tick_params(labelsize=fs)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

    ax.set_title(sa_obj.mission
            + ' with '
            + str(len(sa_obj.vars[sa_obj.stdvarname]))
            + ' footprints: '
            + '\n'
            + sa_obj.sdate.strftime("%Y-%m-%d %H:%M:%S UTC" )
            + ' to '
            + sa_obj.edate.strftime("%Y-%m-%d %H:%M:%S UTC" )
            ,fontsize=fs)

    # - save figure
    if ('savepath' in kwargs.keys() and kwargs['savepath'] != None):
        plt.savefig( kwargs['savepath'] + '/'
            + 'satellite_coverage_for_'
            + sa_obj.region + '_from_'
            + sa_obj.sdate.strftime("%Y%m%d")
            + 'T'
            + sa_obj.sdate.strftime("%H")
            + 'Z'
            + '_to_'
            + sa_obj.edate.strftime("%Y%m%d")
            + 'T'
            + sa_obj.edate.strftime("%H")
            + 'Z'
            + '.png', format = 'png', dpi=300)

    # - show figure
    if ('showfig' in kwargs.keys() and kwargs['showfig'] == True):
        plt.show()

def ts_fig(results_dict):
    import numpy as np
    from datetime import datetime, timedelta
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    modelval = results_dict["model_matches"]
    satval = results_dict["sat_matches"]
    time = results_dict["date_matches"]
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    fs = 12
    plt.plot(time,satval,'ko',label='sHs')
    plt.plot(time,modelval,'ro',label='mHs')
    plt.legend(fontsize=fs,loc='best')
    plt.ylabel('Hs [m]',fontsize=fs)
    #plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
    plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=1))
    #plt.gca().xaxis.set_minor_locator(mdates.SecondLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d %H:%M:%S'))
    plt.gcf().autofmt_xdate()
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.show()

class quicklook_class_sat:

    def quicklook(self, a=False, projection=None, **kwargs):
        """
        Enables to explore the class object (and retrieved results)
        by plotting time series and map.

        param:
            m - map figure (True/False)
            ts - time series (True/False)
            a - all figures (True/False)
            projection - specified projection for cartopy

        return:
            figures
        """
        from wavy.consolidate import consolidate_class
        from wavy.collocmod import collocation_class
        # set plots
        m = kwargs.get('m', a)
        ts = kwargs.get('ts', a)
        mode = kwargs.get('mode', 'comb') # comb,indiv
        if ('ocos' in vars(self).keys()
        and isinstance(self, consolidate_class)
        and isinstance(self.ocos[0], collocation_class)):
            sc = kwargs.get('sc', a)
            if m:
                import cartopy.crs as ccrs
                import cmocean
                import matplotlib.pyplot as plt
                import matplotlib.cm as mplcm
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                lons = self.vars['obs_lons']
                lats = self.vars['obs_lats']
                var = self.vars['obs_values']
                if projection is None:
                    projection = ccrs.PlateCarree()
                vartype = variable_info[self.varalias].get('type','default')
                if kwargs.get('cmap') is None:
                    if vartype == 'cyclic':
                        cmap = mplcm.twilight
                    else:
                        cmap = cmocean.cm.amp
                else:
                    cmap = kwargs.get('cmap')
                lonmax,lonmin = np.nanmax(lons),np.nanmin(lons)
                latmax,latmin = np.nanmax(lats),np.nanmin(lats)
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1, projection=projection)
                ax.set_extent(  [lonmin, lonmax,latmin, latmax],
                                crs = projection )
                sc = ax.scatter(lons,lats,s=10,
                                c = var,
                                marker='o', edgecolor = 'face',
                                cmap=cmap,
                                transform=ccrs.PlateCarree())
                axins = inset_axes(ax,
                           width="5%",  # width = 5% of parent_bbox width
                           height="100%",  # height : 50%
                           loc='lower left',
                           bbox_to_anchor=(1.01, 0., 1, 1),
                           bbox_transform=ax.transAxes,
                           borderpad=0,
                           )
                fig.colorbar(sc, cax=axins, label=self.varalias
                                            + ' [' + self.units + ']')
                ax.coastlines()
                gl = ax.gridlines(draw_labels=True,crs=projection,
                                  linewidth=1, color='grey', alpha=0.4,
                                  linestyle='-')
                gl.top_labels = False
                gl.right_labels = False
                plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                ax.set_title(self.obsname + ' (' + self.obstype + ')\n'
                          + 'from ' + str(self.vars['datetime'][0])
                          + ' to ' + str(self.vars['datetime'][-1]))
                #fig.suptitle('', fontsize=16) # unused
                plt.show()
            if ts:
                import matplotlib.pyplot as plt
                import matplotlib.dates as mdates
                fig = plt.figure(figsize=(9,3.5))
                ax = fig.add_subplot(111)
                colors = ['k','orange']
                if self.obstype == 'insitu':
                    ax.plot(self.vars['datetime'],self.vars['obs_values'],
                        linestyle='None',color=colors[0],
                        label='obs ( ' + self.nID + ' - '
                                       + self.sensor + ' )',
                        marker='o',alpha=.5,ms=2)
                elif self.obstype == 'satellite_altimeter':
                    ax.plot(self.vars['datetime'],self.vars['obs_values'],
                        linestyle='None',color=colors[0],
                        label='obs (' + self.mission + ')',
                        marker='o',alpha=.5,ms=2)
                else:
                    ax.plot(self.vars['datetime'],self.vars['obs_values'],
                        linestyle='None',color=colors[0],
                        label='obs (' + self.obsname + ')',
                        marker='o',alpha=.5,ms=2)
                ax.plot(self.vars['datetime'],self.vars['model_values'],
                        linestyle='None',color=colors[1],
                        label='model (' + self.model + ')',
                        marker='o',alpha=.8,ms=2)
                plt.ylabel(self.varalias + '[' + self.units + ']')
                plt.legend(loc='best')
                plt.tight_layout()
                #ax.set_title()
                plt.show()
            if sc:
                from wavy.utils import compute_quantiles
                lq = np.arange(0.01,1.01,0.01)
                lq = kwargs.get('lq',lq)
                modq = compute_quantiles(np.array(self.vars['model_values']),lq)
                obsq = compute_quantiles(np.array(self.vars['obs_values']),lq)
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(4,4))
                ax = fig.add_subplot(111)
                colors = ['k']
                if self.obstype == 'insitu':
                    ax.plot(self.vars['obs_values'],self.vars['model_values'],
                        linestyle='None',color=colors[0],
                        marker='o',alpha=.5,ms=2)
                    plt.xlabel='obs ( ' + self.nID + ' - '\
                                         + self.sensor + ' )'
                elif self.obstype == 'satellite_altimeter':
                    ax.plot(self.vars['obs_values'],self.vars['model_values'],
                        linestyle='None',color=colors[0],
                        label='obs (' + self.mission + ')',
                        marker='o',alpha=.5,ms=2)
                    plt.xlabel('obs ( ' + self.mission + ' )')
                elif self.obstype == 'POI':
                    ax.plot(self.vars['obs_values'],self.vars['model_values'],
                        linestyle='None',color=colors[0],
                        marker='o',alpha=.5,ms=2)
                    plt.xlabel('obs ( ' + self.obstype + ' )')
                else:
                    ax.plot(self.vars['obs_values'],self.vars['model_values'],
                    linestyle='None',color=colors[0],
                    label='obs (' + self.obstype + ')',
                    marker='o',alpha=.5,ms=2)
                    plt.xlabel('obs ( ' + self.mission + ' )')
                # add quantiles
                ax.plot(obsq,modq,'r')
                # 45 degree line for orientation
                ax.axline((0, 0), (1, 1), lw=.5, color='grey',ls='--')
                plt.ylabel('models ( ' + self.model + ' )')
                vartype = variable_info[self.varalias].get('type','default')
                if vartype == 'cyclic':
                    plt.xlim([0,360])
                    plt.ylim([0,360])
                else:
                    maxv = np.nanmax([self.vars['model_values'],
                                      self.vars['obs_values']])
                    minv = 0
                    plt.xlim([minv,maxv+0.15*maxv])
                    plt.ylim([minv,maxv+0.15*maxv])
                ax.set_title(self.varalias + '[' + self.units + ']')
                plt.tight_layout()
                #ax.set_title()
            plt.show()
        else:
            if m:
                import cartopy.crs as ccrs
                import cmocean
                import matplotlib.pyplot as plt
                import matplotlib.cm as mplcm
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                lons = self.vars['longitude']
                lats = self.vars['latitude']
                var = self.vars[self.stdvarname]
                if projection is None:
                    projection = ccrs.PlateCarree()
                # parse kwargs
                vartype = variable_info[self.varalias].get('type', 'default')
                if kwargs.get('cmap') is None:
                    if vartype == 'cyclic':
                        cmap = mplcm.twilight
                    else:
                        cmap = cmocean.cm.amp
                else:
                    cmap = kwargs.get('cmap')
                lonmax,lonmin = np.max(lons),np.min(lons)
                latmax,latmin = np.max(lats),np.min(lats)
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1, projection=projection)
                # plot track if applicable
                lonmax,lonmin = np.max(lons),np.min(lons)
                latmax,latmin = np.max(lats),np.min(lats)
                if kwargs.get('poi') is not None:
                    plats = kwargs.get('poi').get('latitude')
                    platsmax,platsmin = np.max(plats), np.min(plats)
                    plons = kwargs.get('poi').get('longitude')
                    plonsmax,plonsmin = np.max(plons), np.min(plons)
                    tc = ax.plot(plons,plats,color='cornflowerblue',
                                 ls='-',lw=1,
                                 zorder=-1)
                    tc = ax.plot(plons,plats,color='cornflowerblue',
                                 ls='None',marker='o',ms=5,
                                 markeredgecolor='k',
                                 zorder=-1)
                    lonmax,lonmin = np.max([lonmax,plonsmax]),\
                                    np.min([lonmin,plonsmin])
                    latmax,latmin = np.max([latmax,platsmax]),\
                                    np.min([latmin,platsmin])
                # plot sat
                sc = ax.scatter(lons,lats,s=10,
                                c = var,
                                marker='o', edgecolor = 'face',
                                cmap=cmap,
                                transform=ccrs.PlateCarree())
                axins = inset_axes(ax,
                           width="5%",  # width = 5% of parent_bbox width
                           height="100%",  # height : 50%
                           loc='lower left',
                           bbox_to_anchor=(1.01, 0., 1, 1),
                           bbox_transform=ax.transAxes,
                           borderpad=0,
                           )
                fig.colorbar(sc, cax=axins, label=self.varalias
                                            + ' [' + self.units + ']')
                ax.set_extent(  [lonmin-1, lonmax+1,latmin-1, latmax+1],
                                 crs = projection )
                ax.coastlines(color='k')
                gl = ax.gridlines(draw_labels=True,crs=projection,
                                  linewidth=1, color='grey', alpha=0.4,
                                  linestyle='-')
                gl.top_labels = False
                gl.right_labels = False
                plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                ax.set_title(self.mission + ' (' + self.provider + ')\n'
                          + 'from ' 
                          + (self.vars['datetime'][0]).\
                                  strftime('%Y-%m-%d %H:%M:%S')
                          + ' to '
                          + (self.vars['datetime'][-1]).\
                                  strftime('%Y-%m-%d %H:%M:%S'))
                # plot from quickloop config file
                if ('region' in vars(self).keys()
                    and self.region in quicklook_dict
                    and 'poi' in quicklook_dict[self.region]):
                    for poi in quicklook_dict[self.region]['poi']:
                        pname = quicklook_dict[self.region]['poi'][poi]['name']
                        plat = quicklook_dict[self.region]['poi'][poi]['lat']
                        plon = quicklook_dict[self.region]['poi'][poi]['lon']
                        scp = ax.scatter(plon, plat, s = 20,
                                         c = quicklook_dict\
                                                [self.region]['poi'][poi]\
                                                .get('color','b'),
                                         marker = quicklook_dict[self.region]\
                                                ['poi'][poi]['marker'],
                                         transform = ccrs.PlateCarree())
                    ax.text(plon,plat,pname,transform = ccrs.PlateCarree())
                #fig.suptitle('', fontsize=16) # unused
                plt.show()
            if (ts and mode == 'comb'):
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(9,3.5))
                ax = fig.add_subplot(111)
                colors = ['k']
                ax.plot(self.vars['datetime'],
                        self.vars[self.stdvarname],
                        linestyle='None',color=colors[0],
                        label=self.label,
                        marker='o',alpha=.5,ms=2)
                plt.ylabel(self.varalias + ' [' + self.units + ']')
                plt.legend(loc='best')
                plt.tight_layout()
                #ax.set_title()
                plt.show()
            elif (ts and mode == 'indiv'):
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(9,3.5))
                ax = fig.add_subplot(111)
                for oco in self.ocos:
                    ax.plot(oco.vars['datetime'],
                            oco.vars[oco.stdvarname],
                            linestyle='None',
                            label=oco.label,
                            marker='o',alpha=.5,ms=2)
                plt.ylabel(self.varalias + ' [' + self.units + ']')
                plt.legend(loc='best')
                plt.tight_layout()
                #ax.set_title()
                plt.show()

"""
from abc import abstractmethod


class Quicklook:
    @property
    @abstractmethod
    def projection():
        pass

    @property
    @abstractmethod
    def vars():
        pass

    def quicklook(self, blah):
        proj = self.projection

        for var in self.vars:
            # do something
            pass

class Sat(Quicklook):
    vars = [ "hs", "foo" ]
    projection = pyproj.Proj('+proj=latlong')
"""
