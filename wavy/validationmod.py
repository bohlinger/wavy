"""
    Module to organize the validation procedure
    could include:
        1. Function to create two types of validation files
            One netcdf files comprising the observational and collocated 
            model time series together and one netcdf file including the 
            computed validation statistics for each considered model time 
            step
            May use functions from netcdf module ncmod
            Any new time step (e.g. continuing for a month) is appended
        2. Function to produce chosen validation figure
        3. 
"""
import yaml
import numpy as np

# read yaml config files:
with open("../config/region_specs.yaml", 'r') as stream:
    region_dict=yaml.safe_load(stream)

# define global functions
def rmsd(a,b):
    '''
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    '''
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    n = len(a1)
    diff2 = (a1-b1)**2
    msd = diff2.sum()/n
    rmsd = np.sqrt(msd)
    return msd, rmsd

def scatter_index(obs,model):
    msd,rmsd = rmsd(obs,model)
    stddiff = np.nanstd(obs-model)
    SIrmse = rmsd/np.nanmean(obs)*100.
    SIstd = stddiff/np.nanmean(obs)*100.
    return SIrmse,SIstd

def corr(a,b):
    '''
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    '''
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    corr = np.corrcoef(a1,b1)[1,0]
    return corr

def bias(a,b):
    """
    if nans exist the prinziple of marginalization is applied
    """
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    bias = np.sum(a1-b1)/N
    return bias

def mad(a,b):
    """
    mean absolute deviation
    if nans exist the prinziple of marginalization is applied
    """
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    mad = np.sum(np.abs(a1-b1))/N
    return mad

def disp_validation(valid_dict):
    print('\n')
    print('# ---')
    print('Validation stats')
    print('# ---')
    print('Correlation Coefficient: '
            + '{:0.2f}'.format(valid_dict['corr']))
    print('Root Mean Squared Difference: '
            + '{:0.2f}'.format(valid_dict['rmsd']))
    print('Mean Absolute Difference: ' + '{:0.2f}'.format(valid_dict['mad']))
    print('Bias: ' + '{:0.2f}'.format(valid_dict['bias']))
    print('Scatter Index: ' + '{:0.2f}'.format(valid_dict['SI'][1]))
    print('Mean of Model: ' + '{:0.2f}'.format(valid_dict['mop']))
    print('Mean of Observations: ' + '{:0.2f}'.format(valid_dict['mor']))
    print('Number of Collocated Values: ' + str(valid_dict['nov']))
    print('\n')
    pass


class validation_class():


    def __init__(self,date):
        print ('# ----- ')
        print (" ### Initializing validation_class instance ###")
        print ('# ----- ')

def validate(results_dict,boot=None):
    import numpy as np
    """
    produced metrics:
    mean of product --> mop
    mean of reference --> mor
    mean square difference --> msd
    number of data values --> nov
    scatter index --> SI
    """
    date_matches = np.array(results_dict['date_matches'])
    model_Hs_matches = np.array(results_dict['model_Hs_matches'])
    sat_Hs_matches = np.array(results_dict['sat_Hs_matches'])
    if (boot is None or boot ==  False):
        mop = np.nanmean(model_Hs_matches)
        mor = np.nanmean(sat_Hs_matches)
        msd, rmsd = rmsd(model_Hs_matches,sat_Hs_matches)
        nov = len(sat_Hs_matches)
        mad = mad(model_Hs_matches,sat_Hs_matches)
        corr = corr(model_Hs_matches,sat_Hs_matches)
        bias = bias(model_Hs_matches,sat_Hs_matches)
        SI = scatter_index(model_Hs_matches,sat_Hs_matches)
        arcmfc_validation_dict = {
            'mop':mop,
            'mor':mor,
            'msd':msd,
            'nov':nov,
            'rmsd':rmsd,
            'corr':corr,
            'mad':mad,
            'bias':bias,
            'SI':SI}
    elif boot is True:
        from utils import bootstr, marginalize
        reps=1000
        newmodel,newsat,newidx = marginalize(model_Hs_matches,sat_Hs_matches)
        sat_boot,boot_idx=bootstr(newsat,reps)
        print (len(sat_boot[np.isnan(sat_boot)]))
        RMSD=np.zeros(reps)*np.nan
        MSD=np.zeros(reps)*np.nan
        BIAS=np.zeros(reps)*np.nan
        CORR=np.zeros(reps)*np.nan
        for i in range(reps):
            results_dict = {'date_matches':date_matches[newidx[boot_idx[:,i]]],
                        'model_Hs_matches':newmodel[boot_idx[:,i]],
                        'sat_Hs_matches':newsat[boot_idx[:,i]]}
            try:
                RMSD[i]=validate(results_dict)['rmsd']
                MSD[i]=validate(results_dict)['mad']
                BIAS[i]=validate(results_dict)['bias']
                CORR[i]=validate(results_dict)['corr']
            except IndexError as e:
                print (e)
        arcmfc_validation_dict = {'rmsd':RMSD,'mad':MSD,'bias':BIAS,'corr':CORR}
    return arcmfc_validation_dict

def comp_fig(model,sa_obj,MHs,Mlons,Mlats,results_dict,mode=None,path=None):

    # imports
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    # sort out data/coordinates for plotting
    slons, slats = sa_obj.vars['longitude'],sa_obj.vars['latitude']
    clons, clats = results_dict["model_lons_matches"],\
                    results_dict["model_lats_matches"]
    mhs = MHs.squeeze() 
    mhs[mhs<0] = np.nan
    mhs[mhs>30] = np.nan

    # inflate coords if regular lat/lon grid
    if (len(Mlons.shape)==1):
        Mlons, Mlats = np.meshgrid(Mlons, Mlats)    
    # check region and determine projection
    polarproj=None

    if (sa_obj.region == 'ARCMFC' or sa_obj.region == 'mwam8'\
        or sa_obj.region == 'ARCMFC3' or sa_obj.region == 'ARCMFC3_hc'
        or sa_obj.region == 'CustomSvalbard'):
        # Polar Stereographic Projection
        polarproj = ccrs.NorthPolarStereo(
                            central_longitude=0.0,
                            true_scale_latitude=66,
                            globe=None)
        projection = polarproj
        land = cfeature.GSHHSFeature(scale='i', levels=[1],
                        facecolor=cfeature.COLORS['land'])
    else:
        # other projection for regional visualization
        if sa_obj.region in region_dict['rect']:
            latmin = region_dict['rect'][sa_obj.region]['llcrnrlat']
            latmax = region_dict['rect'][sa_obj.region]['urcrnrlat']
            lonmin = region_dict['rect'][sa_obj.region]['llcrnrlon']
            lonmax = region_dict['rect'][sa_obj.region]['urcrnrlon']
        elif sa_obj.region in region_dict['poly']:
            latmin = np.min(region_dict['poly'][sa_obj.region]['lats'])-.5
            latmax = np.max(region_dict['poly'][sa_obj.region]['lats'])+.5
            lonmin = np.min(region_dict['poly'][sa_obj.region]['lons'])-.5
            lonmax = np.max(region_dict['poly'][sa_obj.region]['lons'])+.5
        else: print("Error: Region not defined!")

        land = cfeature.GSHHSFeature(scale='i', levels=[1], 
                        facecolor=cfeature.COLORS['land'])
        projection = ccrs.LambertAzimuthalEqualArea(
                        central_longitude=0.0, 
                        central_latitude=60.0)

    # make figure
    fig, ax = plt.subplots(nrows=1, ncols=1, 
                        subplot_kw=dict(projection=projection),
                        figsize=(9, 9))
    # plot domain extent
    if projection != polarproj:
        #ax.set_extent([lonmin, 52,latmin, 81],crs = ccrs.PlateCarree())
        ax.set_extent([-26, 32.,45, 85.],crs = ccrs.PlateCarree())
        #ax.set_extent([lonmin, lonmax,latmin, latmax],crs = ccrs.PlateCarree())
        #ax.set_extent([-26, 32,latmin, latmax],crs = ccrs.PlateCarree())
        ax.plot(Mlons[0,:], Mlats[0,:], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
        ax.plot(Mlons[-1,:], Mlats[-1,:], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
        ax.plot(Mlons[:,0], Mlats[:,0], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
        ax.plot(Mlons[:,-1], Mlats[:,-1], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
    else:
        latmin = 40
        latmax = 90
        lonmin = -180
        lonmax = 180
        ax.set_extent([lonmin, lonmax, latmin, latmax],crs = ccrs.PlateCarree())

    # colors
    if mode == 'dir':
        cmap = cmocean.cm.phase
        levels = range(0,360,10)
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    else:
        #cmap = mplcm.GMT_haxby
        cmap = cmocean.cm.amp
        levels = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                3,3.25,3.5,3.75,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    # draw figure features
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fs = 12

    # - model contours
    im = ax.contourf(Mlons, Mlats, mhs, levels = levels, 
                    transform = ccrs.PlateCarree(), 
                    cmap = cmocean.cm.amp, norm = norm)
    #im = ax.contourf(Mlons, Mlats, mhs, transform = ccrs.PlateCarree())

    imc = ax.contour(Mlons, Mlats, mhs, levels = levels[18::1],
                    transform = ccrs.PlateCarree(), 
                    colors='w', linestyle = ':', linewidths = 0.3)
    #ax.clabel(imc, fmt='%.1f', colors='w', fontsize=fs)
    ax.clabel(imc, fmt='%2d', colors='w', fontsize=fs)

    if projection != polarproj:
        # - lons
        cs = ax.contour(Mlons, Mlats, Mlons, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = .6, alpha = 0.6, 
                    levels=range(-40,70,10))
        cs = ax.contour(Mlons, Mlats, Mlons, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = 2, alpha = 0.6, 
                    levels=range(0,1))
        # - lats
        cs = ax.contour(Mlons, Mlats, Mlats, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = .6, alpha = 0.6, 
                    levels=range(45,85,5))

        # - text for lats
        lat = np.arange (70, 90, 10)
        lon = np.repeat (40, len(lat))

        # - regular lat, lon projection
        for lon, lat in zip (lon, lat):
            plt.text (lon, lat, LATITUDE_FORMATTER.format_data(lat), 
                    transform = ccrs.PlateCarree(), fontsize=fs)
        # - text for lons
        lon = np.arange (30,50,10)
        lat = np.repeat (75,len(lon))

        # - regular lat, lon projection
        for lon, lat in zip (lon, lat):
            plt.text (lon, lat, LONGITUDE_FORMATTER.format_data(lon), 
                    transform = ccrs.PlateCarree(), fontsize=fs)

    # - add coastline
    #ax.add_geometries(land.intersecting_geometries(
    #                [lonmin, lonmax, latmin, latmax]),
    #                ccrs.PlateCarree(),
    #                facecolor=cfeature.COLORS['land'],
    #                edgecolor='black',linewidth=1)
    if projection != polarproj:
        ax.add_geometries(land.intersecting_geometries(
                    [lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)
    else:
        ax.add_geometries(land.intersecting_geometries(
                    [lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)

    # - add land color
    ax.add_feature( land, facecolor = 'burlywood', alpha = 0.5 )

    # - add satellite
    if len(clats)>0:
        sc = ax.scatter(clons,clats,s=10,
                c=results_dict['sat_Hs_matches'],
                marker='o',verts=levels, edgecolor = 'face',
                cmap=cmocean.cm.amp, norm = norm, 
                transform=ccrs.PlateCarree())

    # - colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='vertical',
                        fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Hs [m]',size=fs)
    cbar.ax.tick_params(labelsize=fs)

    # - title
    plt.title(model + ' model time step: '
            + results_dict['valid_date'][0].strftime("%Y-%m-%d %H:%M:%S UTC") 
            + '\n'
            + sa_obj.sat
            + ' coverage \n from ' 
            + results_dict['date_matches'][0].strftime("%Y-%m-%d %H:%M:%S UTC" )            + ' to '
            + results_dict['date_matches'][-1].strftime("%Y-%m-%d %H:%M:%S UTC")
            ,fontsize=12)
    if path != None:
        plt.savefig(path + model + '_test_' 
                + results_dict['valid_date'][0].strftime("%Y%m%d")
                + 'T' 
                + results_dict['valid_date'][0].strftime("%H") 
                + 'Z.png', format = 'png', dpi=200)
    plt.show()

def comp_wind(model,var,Mlons,Mlats,date,region,mode=None):

    # imports
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    # sort out data/coordinates for plotting
    var = var.squeeze()
    if model == 'ww3':
        var = (var - 180) % 360

    # check region and determine projection

    latmin = region_dict['rect'][region]['llcrnrlat']
    latmax = region_dict['rect'][region]['urcrnrlat']
    lonmin = region_dict['rect'][region]['llcrnrlon']
    lonmax = region_dict['rect'][region]['urcrnrlon']
    land = cfeature.GSHHSFeature(scale='i', levels=[1],
                    facecolor=cfeature.COLORS['land'])
    projection = ccrs.LambertAzimuthalEqualArea(
                    central_longitude=0.0,
                    central_latitude=60.0)

    # make figure
    fig, ax = plt.subplots(nrows=1, ncols=1,
                        subplot_kw=dict(projection=projection),
                        figsize=(9, 9))
    # plot domain extent
    ax.set_extent([-26, 32.,45, 85.],crs = ccrs.PlateCarree())
    #ax.set_extent([lonmin, lonmax,latmin, latmax],crs = ccrs.PlateCarree())
    #ax.set_extent([-26, 32,latmin, latmax],crs = ccrs.PlateCarree())
    ax.plot(Mlons[0,:], Mlats[0,:], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)
    ax.plot(Mlons[-1,:], Mlats[-1,:], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)
    ax.plot(Mlons[:,0], Mlats[:,0], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)
    ax.plot(Mlons[:,-1], Mlats[:,-1], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)

    # colors
    if mode == 'dir':
        cmap = cmocean.cm.phase
        levels = range(0,360,5)
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    else:
        #cmap = mplcm.GMT_haxby
        cmap = cmocean.cm.amp
        levels = [0,1,2,3,4,6,8,10,11,12,13,14,
             15,16,18,20,22,24,26,28,30,32,
             35,38,42]
        norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    # draw figure features
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fs = 12

    # - model contours
    im = ax.contourf(Mlons, Mlats, var, levels = levels,
                    transform = ccrs.PlateCarree(),
                    cmap = cmap, norm = norm)
    #im = ax.contourf(Mlons, Mlats, mhs, transform = ccrs.PlateCarree())
    if mode == 'dir':
        imc = ax.contour(Mlons, Mlats, var, levels = levels[::5],
                    transform = ccrs.PlateCarree(),
                    colors='w', linestyle = ':', linewidths = 0.3)
    else:
        imc = ax.contour(Mlons, Mlats, var, levels = levels[15::1],
                    transform = ccrs.PlateCarree(),
                    colors='w', linestyle = ':', linewidths = 0.3)

    ax.clabel(imc, fmt='%2d', colors='w', fontsize=fs)

    # - lons
    cs = ax.contour(Mlons, Mlats, Mlons, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = .6, alpha = 0.6,
                    levels=range(-40,70,10))
    cs = ax.contour(Mlons, Mlats, Mlons, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = 2, alpha = 0.6,
                    levels=range(0,1))
    # - lats
    cs = ax.contour(Mlons, Mlats, Mlats, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = .6, alpha = 0.6,
                    levels=range(45,85,5))

    # - text for lats
    lat = np.arange (70, 90, 10)
    lon = np.repeat (40, len(lat))

    # - regular lat, lon projection
    for lon, lat in zip (lon, lat):
        plt.text (lon, lat, LATITUDE_FORMATTER.format_data(lat),
                    transform = ccrs.PlateCarree(), fontsize=fs)
    # - text for lons
    lon = np.arange (30,50,10)
    lat = np.repeat (75,len(lon))

    # - regular lat, lon projection
    for lon, lat in zip (lon, lat):
        plt.text (lon, lat, LONGITUDE_FORMATTER.format_data(lon),
                    transform = ccrs.PlateCarree(), fontsize=fs)

    # - add coastline
    ax.add_geometries(land.intersecting_geometries(
                    [lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)

    # - add land color
    ax.add_feature( land, facecolor = 'burlywood', alpha = 0.5 )

    # - colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='vertical',
                        fraction=0.046, pad=0.04)
    if mode == 'dir':
        cbar.ax.set_ylabel('degree',size=fs)
    else:
        cbar.ax.set_ylabel('Wind speed [m/s]',size=fs)
    cbar.ax.tick_params(labelsize=fs)

    # - title
    plt.title(model + ' model time step: '
            + date.strftime("%Y-%m-%d %H:%M:%S UTC")
            ,fontsize=12)

    #plt.savefig(model + '_wind_test_' 
    #plt.savefig(model + '_dir_test_' 
    #            + date.strftime("%Y%m%d")
    #            + 'T' 
    #            + date.strftime("%H") 
    #            + 'Z.png', format = 'png', dpi=200)
    #plt.show()

def comp_wind_quiv(model,u,v,Mlons,Mlats,date,region):

    # imports
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    # sort out data/coordinates for plotting
    u = u.squeeze()
    v = v.squeeze()
    var = np.sqrt(u**2+v**2)

    # check region and determine projection

    latmin = region_dict['rect'][region]['llcrnrlat']
    latmax = region_dict['rect'][region]['urcrnrlat']
    lonmin = region_dict['rect'][region]['llcrnrlon']
    lonmax = region_dict['rect'][region]['urcrnrlon']
    land = cfeature.GSHHSFeature(scale='i', levels=[1],
                    facecolor=cfeature.COLORS['land'])
    projection = ccrs.LambertAzimuthalEqualArea(
                    central_longitude=0.0,
                    central_latitude=60.0)

    # make figure
    fig, ax = plt.subplots(nrows=1, ncols=1,
                        subplot_kw=dict(projection=projection),
                        figsize=(9, 9))
    # plot domain extent
    ax.set_extent([-26, 32.,45, 85.],crs = ccrs.PlateCarree())
    #ax.set_extent([lonmin, lonmax,latmin, latmax],crs = ccrs.PlateCarree())
    #ax.set_extent([-26, 32,latmin, latmax],crs = ccrs.PlateCarree())
    ax.plot(Mlons[0,:], Mlats[0,:], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)
    ax.plot(Mlons[-1,:], Mlats[-1,:], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)
    ax.plot(Mlons[:,0], Mlats[:,0], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)
    ax.plot(Mlons[:,-1], Mlats[:,-1], '-', transform= ccrs.PlateCarree(),
            color = 'gray', linewidth =2)

    # colors
    #cmap = mplcm.GMT_haxby
    cmap = cmocean.cm.amp
    levels = [0,1,2,3,4,6,8,10,11,12,13,14,
             15,16,18,20,22,24,26,28,30,32,
             35,38,42]
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    # draw figure features
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    fs = 12

    # - model contours
    im = ax.contourf(Mlons, Mlats, var, levels = levels,
                    transform = ccrs.PlateCarree(),
                    cmap = cmocean.cm.amp, norm = norm)
    #im = ax.contourf(Mlons, Mlats, mhs, transform = ccrs.PlateCarree())

    imc = ax.contour(Mlons, Mlats, var, levels = levels[15::1],
                    transform = ccrs.PlateCarree(),
                    colors='w', linestyle = ':', linewidths = 0.3)
    ax.clabel(imc, fmt='%2d', colors='w', fontsize=fs)
    # add quivers for wind
    qv = ax.quiver(Mlons, Mlats, u, v, color='k', transform=ccrs.PlateCarree(),scale=500)

    # - lons
    cs = ax.contour(Mlons, Mlats, Mlons, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = .6, alpha = 0.6,
                    levels=range(-40,70,10))
    cs = ax.contour(Mlons, Mlats, Mlons, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = 2, alpha = 0.6,
                    levels=range(0,1))
    # - lats
    cs = ax.contour(Mlons, Mlats, Mlats, transform = ccrs.PlateCarree(),
                    colors='k', linewidths = .6, alpha = 0.6,
                    levels=range(45,85,5))

    # - text for lats
    lat = np.arange (70, 90, 10)
    lon = np.repeat (40, len(lat))

    # - regular lat, lon projection
    for lon, lat in zip (lon, lat):
        plt.text (lon, lat, LATITUDE_FORMATTER.format_data(lat),
                    transform = ccrs.PlateCarree(), fontsize=fs)
    # - text for lons
    lon = np.arange (30,50,10)
    lat = np.repeat (75,len(lon))

    # - regular lat, lon projection
    for lon, lat in zip (lon, lat):
        plt.text (lon, lat, LONGITUDE_FORMATTER.format_data(lon),
                    transform = ccrs.PlateCarree(), fontsize=fs)

    # - add coastline
    ax.add_geometries(land.intersecting_geometries(
                    [lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black',linewidth=1)

    # - add land color
    ax.add_feature( land, facecolor = 'burlywood', alpha = 0.5 )

    # - colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='vertical',
                        fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Wind speed [m/s]',size=fs)
    cbar.ax.tick_params(labelsize=fs)

    # - title
    plt.title(model + ' model time step: '
            + date.strftime("%Y-%m-%d %H:%M:%S UTC")
            ,fontsize=12)

    #plt.savefig(model + '_wind_test_'
    #            + date.strftime("%Y%m%d")
    #            + 'T'
    #            + date.strftime("%H")
    #            + 'Z.png', format = 'png', dpi=200)
    #plt.show()


def plot_sat(sa_obj):
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    # sort out data/coordinates for plotting
    slons, slats = sa_obj.vars['longitude'],sa_obj.vars['latitude']

    # check region and determine projection

    polarproj=None

    if (sa_obj.region == 'ARCMFC' or sa_obj.region == 'mwam8'\
        or sa_obj.region == 'ARCMFC3' or sa_obj.region == 'ARCMFC3_hc'
        or sa_obj.region == 'ErinFix'):
        # Polar Stereographic Projection
        polarproj = ccrs.NorthPolarStereo(
                            central_longitude=0.0, 
                            true_scale_latitude=66, 
                            globe=None)
        projection = polarproj
        land = cfeature.GSHHSFeature(scale='i', levels=[1],
                        facecolor=cfeature.COLORS['land'])
        #pass
    else:
        if sa_obj.region in region_dict['rect']:
            latmin = region_dict['rect'][sa_obj.region]['llcrnrlat']
            latmax = region_dict['rect'][sa_obj.region]['urcrnrlat']
            lonmin = region_dict['rect'][sa_obj.region]['llcrnrlon']
            lonmax = region_dict['rect'][sa_obj.region]['urcrnrlon']
        elif sa_obj.region in region_dict['poly']:
            latmin = np.min(region_dict['poly'][sa_obj.region]['lats'])-.5
            latmax = np.max(region_dict['poly'][sa_obj.region]['lats'])+.5
            lonmin = np.min(region_dict['poly'][sa_obj.region]['lons'])-.5
            lonmax = np.max(region_dict['poly'][sa_obj.region]['lons'])+.5
        else: print("Error: Region not defined!")
        azimproj = ccrs.LambertAzimuthalEqualArea(
                        central_longitude=0.0,
                        central_latitude=60.0)
        projection = azimproj
        land = cfeature.GSHHSFeature(scale='i', levels=[1],
                        facecolor=cfeature.COLORS['land'])
        #projection = ccrs.Mercator()

    # make figure
    fig, ax = plt.subplots(nrows=1, ncols=1,
                        subplot_kw=dict(projection=projection),
                        figsize=(9, 9))
    # plot domain extent
    if projection != polarproj:
        ax.set_extent([lonmin+6, 45,latmin+6, 84],crs = ccrs.PlateCarree())
    else:
        ax.set_extent([-180, 180,40, 90],crs = ccrs.PlateCarree())


    # plot lats/lons
    #ax.plot(5, 60, transform = ccrs.PlateCarree())
    #ax.gridlines(draw_labels = True)
    #ax.gridlines(linewidth = 1.5,color = 'gray', alpha = 0.4, linestyle = '--')
    if (projection != polarproj and projection!=azimproj):
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    
    # colors
    #cmap = mplcm.GMT_haxby
    cmap = cmocean.cm.amp
    levels = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,
                    3,3.5,4,4.5,6,7,8,9,10,12,15,20]
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    # - add coastline
    if projection != polarproj:
        ax.add_geometries(land.intersecting_geometries(
                    #[lonmin, lonmax, latmin, latmax]),
                    [lonmin-20, lonmax+10, latmin-10, 90]),
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

    # - add satellite
    sc = ax.scatter(slons,slats,s=10,c=sa_obj.Hs,
                marker='o',verts=levels, edgecolor = 'face',
                cmap=cmocean.cm.amp, norm = norm,
                transform=ccrs.PlateCarree())
    
    # - colorbar
    cbar = fig.colorbar(sc, ax=ax, orientation='vertical',
                        fraction=0.03, pad=0.03)
                        #fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Hs [m]')
    plt.title(sa_obj.sat
            + ' with '
            + str(len(sa_obj.Hs))
            + ' footprints: '
            + '\n'
            + sa_obj.sdate.strftime("%Y-%m-%d %H:%M:%S UTC" )
            + ' to '
            + sa_obj.edate.strftime("%Y-%m-%d %H:%M:%S UTC" )
            ,fontsize=12)
    #plt.savefig(sa_obj.sat + '_coverage_from_'
    #        + sa_obj.sdate.strftime("%Y%m%d")
    #        + 'T'
    #        + sa_obj.sdate.strftime("%H")
    #        + 'Z'
    #        + '_to_'
    #        + sa_obj.edate.strftime("%Y%m%d")
    #        + 'T' 
    #        + sa_obj.edate.strftime("%H")
    #        + 'Z'
    #        + '.png', format = 'png', dpi=300)
    plt.show()

def ts_fig(results_dict):
    import numpy as np
    from datetime import datetime, timedelta
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    mHs = results_dict["model_Hs_matches"]
    sHs = results_dict["sat_Hs_matches"]
    time = results_dict["date_matches"]
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    fs = 12
    plt.plot(time,sHs,'ko',label='sHs')
    plt.plot(time,mHs,'ro',label='mHs')
    plt.legend(fontsize=fs,loc='best')
    plt.ylabel('Hs [m]',fontsize=fs)
    #plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
    plt.gca().xaxis.set_major_locator(mdates.SecondLocator(interval=1))
    #plt.gca().xaxis.set_minor_locator(mdates.SecondLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d %H:%M:%S'))
    plt.gcf().autofmt_xdate()
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.show()
