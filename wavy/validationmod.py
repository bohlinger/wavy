"""
    Module to organize the validation procedure
    could include:
        1. Function to create two types of validation files
            One netcdf files comprising the observational and collocated 
            model time series together and one netcdf file including the 
            computed validation statistics for each considered model time 
            step
            May use functions from netcdf module custom_nc
            Any new time step (e.g. continuing for a month) is appended
        2. Function to produce chosen validation figure
        3. 
"""
class validation_class():


    def __init__(self,date):
        print ('# ----- ')
        print (" ### Initializing validation_class instance ###")
        print ('# ----- ')
        

    #def valid_fig(self,figtype):

def comp_fig(model,sa_obj,MHs,Mlons,Mlats,results_dict):
    from mpl_toolkits.basemap import Basemap, cm
    import matplotlib.cm as mplcm
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    mHs = MHs.squeeze()
    mHs[np.where(mHs<0)[0],np.where(mHs<0)[1]]=np.nan
    clevs = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,6,7,8,9,10,12,15,20]
    cmap=cm.GMT_haxby
    norm = mpl.colors.BoundaryNorm(clevs, cmap.N)
    m = Basemap(width=4100000,
                height=4500000,
                resolution='l',
                projection='laea',
                lat_ts=66,lat_0=66,lon_0=1.)
    m.drawcoastlines()
    m.drawcountries()
    m.drawmeridians(np.arange(0,360,5))
    m.drawparallels(np.arange(-90,90,5))
    x, y = m(sa_obj.rloc[1],sa_obj.rloc[0])
    x2, y2 = m(results_dict["model_lons_matches"],results_dict["model_lats_matches"])
    lons, lats = m(Mlons,Mlats)
    cs = m.contourf(lons,lats,mHs,clevs,cmap=cmap,norm=norm)
    sc = m.scatter(x2,y2,s=30,c=results_dict['sat_Hs_matches'],marker='o',cmap=cmap,norm=norm,edgecolor='k',linewidths=0.05,verts=clevs)
    plt.title('Model time step: '
            + results_dict['valid_date'][0].strftime("%Y-%m-%d %H:%M:%S UTC") 
            + '\n'
            + 'S3a coverage from ' 
            + results_dict['date_matches'][0].strftime("%Y-%m-%d %H:%M:%S UTC" )            + ' to '
            + results_dict['date_matches'][-1].strftime("%Y-%m-%d %H:%M:%S UTC" )
            ,fontsize=8)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Hs [m]')
    plt.show()
