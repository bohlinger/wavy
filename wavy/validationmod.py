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
class validation_class():


    def __init__(self,date):
        print ('# ----- ')
        print (" ### Initializing validation_class instance ###")
        print ('# ----- ')
        

def validate(results_dict,boot=None):
    import numpy as np
    from utils import rmsd, corr, mad, bias, scatter_index
    """
    produced metrics:
    mean of product --> mop
    mean of reference --> mor
    mean square difference --> msd
    number of data values --> nov
    scatter index --> SI
    """
    flatten = lambda l: [item for sublist in l for item in sublist]
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

def comp_fig(model,sa_obj,MHs,Mlons,Mlats,results_dict):
    from mpl_toolkits.basemap import Basemap, cm
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from region_specs import region_dict
    mHs = MHs.squeeze()
    mHs[np.where(mHs<0)[0],np.where(mHs<0)[1]]=np.nan
    if (sa_obj.region == 'MoskWC' or sa_obj.region == 'MoskNC'):
        clevs = np.arange(0,5,0.1)
    else:
        clevs = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,6,7,8,9,10,12,15,20]
    cmap=cm.GMT_haxby
    norm = mpl.colors.BoundaryNorm(clevs, cmap.N)
    if sa_obj.region == 'ARCMFC':
        # Polar Stereographic Projection
        m = Basemap(
            projection='npstere',
            boundinglat=region_dict[sa_obj.region]["boundinglat"],
            lon_0=0,
            resolution='l',area_thresh=1000
            )
    else:
        m = Basemap(width=4100000,
                height=4500000,
                resolution='l',
                projection='laea',
                lat_ts=66,lat_0=66,lon_0=1.)
    if (model == 'swan_karmoy250'):
        m.drawmeridians(np.arange(0,360,0.5))
        m.drawparallels(np.arange(-90,90,0.5))
    else:
        m.drawcoastlines()
        m.drawcountries()
        m.drawmeridians(np.arange(0,360,5))
        m.drawparallels(np.arange(-90,90,5))
    x, y = m(sa_obj.loc[1],sa_obj.loc[0])
    x2, y2 = m(results_dict["model_lons_matches"],results_dict["model_lats_matches"])
    if (model == 'swan_karmoy250'):
        Mlons, Mlats = np.meshgrid(Mlons,Mlats)
    lons, lats = m(Mlons,Mlats)
    cs = m.contourf(lons,lats,mHs,clevs,cmap=cmap,norm=norm)
    sc = m.scatter(x2,y2,s=30,c=results_dict['sat_Hs_matches'],marker='o',cmap=cmap,norm=norm,edgecolor='k',linewidths=0.05,verts=clevs)
    plt.title(model + ' model time step: '
            + results_dict['valid_date'][0].strftime("%Y-%m-%d %H:%M:%S UTC") 
            + '\n'
            + 'S3a coverage from ' 
            + results_dict['date_matches'][0].strftime("%Y-%m-%d %H:%M:%S UTC" )            + ' to '
            + results_dict['date_matches'][-1].strftime("%Y-%m-%d %H:%M:%S UTC" )
            ,fontsize=8)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Hs [m]')
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
