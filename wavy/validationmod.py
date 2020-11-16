"""
    Module to organize the validation procedure
"""
import yaml
import numpy as np
import os

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/region_specs.yaml'))
with open(moddir,'r') as stream:
    region_dict=yaml.safe_load(stream)
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_shortcuts.yaml'))
with open(moddir,'r') as stream:
    shortcuts_dict=yaml.safe_load(stream)
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/quicklook_specs.yaml'))
with open(moddir,'r') as stream:
    quicklook_dict=yaml.safe_load(stream)

# define global functions
def calc_rmsd(a,b):
    '''
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    '''
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    n = len(a1)
    diff2 = (a1-b1)**2
    msd = diff2.sum()/n
    rmsd = np.sqrt(msd)
    return msd, rmsd

def calc_scatter_index(obs,model):
    msd,rmsd = calc_rmsd(obs,model)
    stddiff = np.nanstd(obs-model)
    SIrmse = rmsd/np.nanmean(obs)*100.
    SIstd = stddiff/np.nanmean(obs)*100.
    return SIrmse,SIstd

def calc_corrcoef(a,b):
    '''
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    '''
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    corr = np.corrcoef(a1,b1)[1,0]
    return corr

def calc_bias(a,b):
    """
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    bias = np.sum(a1-b1)/N
    return bias

def calc_mad(a,b):
    """
    mean absolute deviation
    if nans exist the prinziple of marginalization is applied
    input: np.arrays with np.nan for invalids
    """
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
    vars in dict: np.arrays with np.nan for invalids
    produced metrics:
    mean of product --> mop
    mean of reference --> mor
    mean square difference --> msd
    number of data values --> nov
    scatter index --> SI
    """
    date_matches = results_dict['date_matches']
    model_matches = results_dict['model_matches']
    sat_matches = results_dict['sat_matches']
    if (boot is None or boot ==  False):
        mop = np.nanmean(model_matches)
        mor = np.nanmean(sat_matches)
        msd, rmsd = calc_rmsd(model_matches,sat_matches)
        nov = len(sat_matches)
        mad = calc_mad(model_matches,sat_matches)
        corr = calc_corrcoef(model_matches,sat_matches)
        bias = calc_bias(model_matches,sat_matches)
        SI = calc_scatter_index(model_matches,sat_matches)
        validation_dict = {
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
        newmodel,newsat,newidx = marginalize(model_matches,sat_matches)
        sat_boot,boot_idx=bootstr(newsat,reps)
        print (len(sat_boot[np.isnan(sat_boot)]))
        RMSD=np.zeros(reps)*np.nan
        MSD=np.zeros(reps)*np.nan
        BIAS=np.zeros(reps)*np.nan
        CORR=np.zeros(reps)*np.nan
        for i in range(reps):
            results_dict = {'date_matches':date_matches[newidx[boot_idx[:,i]]],
                        'model_matches':newmodel[boot_idx[:,i]],
                        'sat_matches':newsat[boot_idx[:,i]]}
            try:
                RMSD[i]=validate(results_dict)['rmsd']
                MSD[i]=validate(results_dict)['mad']
                BIAS[i]=validate(results_dict)['bias']
                CORR[i]=validate(results_dict)['corr']
            except IndexError as e:
                print (e)
        validation_dict = {'rmsd':RMSD,'mad':MSD,'bias':BIAS,'corr':CORR}
    return validation_dict

def comp_fig(model,sa_obj,MHs,Mlons,Mlats,results_dict,var,mode=None,**kwargs):
    # imports
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker

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
    elif sa_obj.region in model_dict:
        # model bounds
        latmin = np.min(Mlats)
        latmax = np.max(Mlats)
        lonmin = np.min(Mlons)
        lonmax = np.max(Mlons)
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
    ax.set_extent([lonmin, lonmax,latmin, latmax],crs = ccrs.PlateCarree())
    ax.plot(Mlons[0,:], Mlats[0,:], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
    ax.plot(Mlons[-1,:], Mlats[-1,:], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
    ax.plot(Mlons[:,0], Mlats[:,0], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)
    ax.plot(Mlons[:,-1], Mlats[:,-1], '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)

    # plot polygon if defined
    if sa_obj.region in region_dict['poly']:
        ax.plot(region_dict['poly'][sa_obj.region]['lons'],
                region_dict['poly'][sa_obj.region]['lats'], 
                '-', transform= ccrs.PlateCarree(), 
                color = 'gray', linewidth =2)

    # colors
    if mode == 'dir':
        cmap = cmocean.cm.phase
        levels = range(0,360,10)
    else:
        # cmap = mplcm.GMT_haxby
        cmap = cmocean.cm.amp
        levels = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                3,3.25,3.5,3.75,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18]

    if 'cmap' in kwargs.keys():
        cmap = kwargs['cmap']

    if (sa_obj.region in quicklook_dict
    and 'cm_levels' in quicklook_dict[sa_obj.region]
    and quicklook_dict[sa_obj.region]['cm_levels'] is not None):
        levels = quicklook_dict[sa_obj.region]['cm_levels']

    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

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
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': fs, 'color': gridcolor}
    gl.ylabel_style = {'size': fs, 'color': gridcolor}

    # - model contours
    im = ax.contourf(Mlons, Mlats, mhs, levels = levels, 
                    transform = ccrs.PlateCarree(), 
                    cmap = cmocean.cm.amp, norm = norm)
    scl = quicklook_dict[sa_obj.region]['scl']
    icl = quicklook_dict[sa_obj.region]['icl']
    imc = ax.contour(Mlons, Mlats, mhs, levels = levels[scl::icl],
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
    if len(clats)>0:
        sc = ax.scatter(clons,clats,s=10,
                c=results_dict['sat_matches'],
                marker = 'o', edgecolor = 'face',
                cmap = cmocean.cm.amp, norm = norm, 
                transform = ccrs.PlateCarree())

    # - point of interests depending on region
    if (sa_obj.region in quicklook_dict 
    and 'poi' in quicklook_dict[sa_obj.region]):
        for poi in quicklook_dict[sa_obj.region]['poi']:
            pname = quicklook_dict[sa_obj.region]['poi'][poi]['name']
            plat = quicklook_dict[sa_obj.region]['poi'][poi]['lat']
            plon = quicklook_dict[sa_obj.region]['poi'][poi]['lon']
            scp = ax.scatter(plon,plat,s=20, c='b',
                marker = quicklook_dict[sa_obj.region]['poi'][poi]['marker'],
                transform = ccrs.PlateCarree())
            ax.text(plon, plat, pname, transform = ccrs.PlateCarree())

    # - colorbar
    cbar = fig.colorbar(im, ax=ax, orientation='vertical',
                        fraction=0.035, pad=0.04)
    cbar.ax.set_ylabel(var + ' [m]',size=fs)
    cbar.ax.tick_params(labelsize=fs)

    # - title
    plt.title(model + ' model time step: '
            + results_dict['valid_date'][0].strftime("%Y-%m-%d %H:%M:%S UTC") 
            + '\n'
            + sa_obj.sat
            + ' coverage \n from ' 
            + results_dict['date_matches'][0].strftime("%Y-%m-%d %H:%M:%S UTC" )            + ' to '
            + results_dict['date_matches'][-1].strftime("%Y-%m-%d %H:%M:%S UTC")
            ,fontsize=fs)
    # - save figure
    if ('savepath' in kwargs.keys() and kwargs['savepath'] != None):
        plt.savefig( kwargs['savepath'] + '/' + model 
                + '_vs_satellite_'
                + results_dict['valid_date'][0].strftime("%Y%m%d")
                + 'T' 
                + results_dict['valid_date'][0].strftime("%H") 
                + 'Z.png', format = 'png', dpi=200)
    # - show figure
    if ('showfig' in kwargs.keys() and kwargs['showfig'] == True):
        plt.show()

def plot_sat(sa_obj,var,**kwargs):
    # imports
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cmocean
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker

    # sort out data/coordinates for plotting
    slons, slats = sa_obj.vars['longitude'],sa_obj.vars['latitude']
    if sa_obj.region in model_dict:
        from modelmod import get_model
        grid_date = model_dict[sa_obj.region]['grid_date']
        model_var_dict = \
                get_model(model=sa_obj.region,fc_date=grid_date)
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
        elif sa_obj.region in region_dict['poly']: 
            latmin = np.min(region_dict['poly'][sa_obj.region]['lats'])
            latmax = np.max(region_dict['poly'][sa_obj.region]['lats'])
            lonmin = np.min(region_dict['poly'][sa_obj.region]['lons'])
            lonmax = np.max(region_dict['poly'][sa_obj.region]['lons'])
        elif sa_obj.region in model_dict:
            # model bounds
            latmin = np.min(model_var_dict['model_lats'])
            latmax = np.max(model_var_dict['model_lats'])
            lonmin = np.min(model_var_dict['model_lons'])
            lonmax = np.max(model_var_dict['model_lons'])
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
        ax.set_extent([-180, 180,-90, 90],crs = ccrs.PlateCarree())

    # plot model domain if region is a model domain
    if (sa_obj.region in model_dict
    and len(model_var_dict['model_lats'].shape)==1):
        lenlons = len(model_var_dict['model_lons'][:])
        lenlats = len(model_var_dict['model_lats'][:])
        ax.plot([model_var_dict['model_lons'][0]]*lenlats,    
                model_var_dict['model_lats'][:], '-', 
                transform = ccrs.PlateCarree(), 
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['model_lons'][:], 
                [model_var_dict['model_lats'][-1]]*lenlons, '-',    
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot([model_var_dict['model_lons'][-1]]*lenlats,      
                model_var_dict['model_lats'][::-1], '-', 
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['model_lons'][::-1], 
                [model_var_dict['model_lats'][0]]*lenlons, '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
    if (sa_obj.region in model_dict 
    and len(model_var_dict['model_lats'].shape)==2):
        ax.plot(model_var_dict['model_lons'][0,:],
                model_var_dict['model_lats'][0,:], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['model_lons'][-1,:],
                model_var_dict['model_lats'][-1,:], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['model_lons'][:,0],
                model_var_dict['model_lats'][:,0], '-',
                transform = ccrs.PlateCarree(),
                color = 'gray', linewidth = 2)
        ax.plot(model_var_dict['model_lons'][:,-1],
                model_var_dict['model_lats'][:,-1], '-',
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
    levels = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                3,3.25,3.5,3.75,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
    if (sa_obj.region in quicklook_dict 
    and 'cm_levels' in quicklook_dict[sa_obj.region]
    and quicklook_dict[sa_obj.region]['cm_levels'] is not None):
        levels = quicklook_dict[sa_obj.region]['cm_levels']

    if 'cmap' in kwargs.keys():
        cmap = kwargs['cmap']

    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

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

    # - add satellite
    sc = ax.scatter(slons,slats,s=10,c=sa_obj.vars[var],
                marker='o', edgecolor = 'face',
                cmap = cmocean.cm.amp, norm = norm,
                transform = ccrs.PlateCarree())

    # - point of interests depending on region
    if (sa_obj.region in quicklook_dict
    and 'poi' in quicklook_dict[sa_obj.region]):
        for poi in quicklook_dict[sa_obj.region]['poi']:
            pname = quicklook_dict[sa_obj.region]['poi'][poi]['name']
            plat = quicklook_dict[sa_obj.region]['poi'][poi]['lat']
            plon = quicklook_dict[sa_obj.region]['poi'][poi]['lon']
            scp = ax.scatter(plon,plat,s=20, c='b',
                marker = quicklook_dict[sa_obj.region]['poi'][poi]['marker'],
                transform = ccrs.PlateCarree())
            ax.text(plon, plat, pname, transform = ccrs.PlateCarree())

    # - plot polygon
    if sa_obj.region in region_dict['poly']:
        ax.plot(
            region_dict['poly'][sa_obj.region]['lons'],
            region_dict['poly'][sa_obj.region]['lats'],
            'k:',transform=ccrs.PlateCarree()
            )
    # - colorbar
    cbar = fig.colorbar(sc, ax=ax, orientation='vertical',
                        fraction=0.035, pad=0.04)
    cbar.ax.set_ylabel(var + ' [m]')
    cbar.ax.tick_params(labelsize=fs)

    plt.title(sa_obj.sat
            + ' with '
            + str(len(sa_obj.vars[var]))
            + ' footprints: '
            + '\n'
            + sa_obj.sdate.strftime("%Y-%m-%d %H:%M:%S UTC" )
            + ' to '
            + sa_obj.edate.strftime("%Y-%m-%d %H:%M:%S UTC" )
            ,fontsize=fs)

    # - save figure
    if ('savepath' in kwargs.keys() and kwargs['savepath'] != None):
        plt.savefig( kwargs['savepath'] + '/' 
            + sa_obj.sat + '_coverage_for_' 
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
