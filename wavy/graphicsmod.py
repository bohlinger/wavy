#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
Module to handle all visualization related functions.

I try to mostly follow the PEP convention for python code style. 
Constructive comments on style and effecient programming are most welcome!
'''

# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# control backend
import matplotlib
matplotlib.use('Agg')

# progress bar and other stuff
import sys

# ignore irrelevant warnings from matplotlib for stdout
import warnings
#warnings.filterwarnings("ignore")

# all class
import numpy as np
from datetime import datetime, timedelta
import argparse
from argparse import RawTextHelpFormatter
import os
import calendar
from dateutil.relativedelta import relativedelta
from copy import deepcopy
from mpl_toolkits import axes_grid1
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# get_altim
if sys.version_info <= (3, 0):
    from urllib import urlretrieve, urlcleanup # python2
else:
    from urllib.request import urlretrieve, urlcleanup # python3

# get necessary paths for module
import pathfinder

# --- global functions ------------------------------------------------#

def progress(count, total, status=''):
    "from: https://gist.github.com/vladignatyev/06860ec2040cb497f0f3"
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()

# flatten all lists before returning them
# define flatten function for lists
''' fct does the following:
flat_list = [item for sublist in TIME for item in sublist]
or:
for sublist in TIME:
for item in sublist:
flat_list.append(item)
'''
flatten = lambda l: [item for sublist in l for item in sublist]

def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    '''
    Add a vertical color bar to an image plot.
    '''
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, 
                                    aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)

val_fig_lim_dict =  {'SI':[0,60],
                    'msd':[0,1.5],
                    'corr':[-1,1],
                    'rmsd':[0,1.5],
                    'mor':[0,10],
                    'nov':[0,700],
                    'mop':[0,10],
                    'bias':[-1,1],
                    'mad':[0,1.5]
                    }
val_fig_varname_dict =  {'SI':'scatter index',
                         'msd':'mean square error',
                        'corr':'correlation coefficient',
                        'rmsd':'root mean square error',
                        'mor':'mean of observations',
                        'nov':'number of observations',
                        'mop':'mean of model',
                        'bias':'bias',
                        'mad':'mean absolute error'
                        }


def make_val_ts_fig_arcmfc(val_name,ts_lst,dtime_lst,filename_fig,forecasts):
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(111)
    fs = 14
    days_in_month = calendar.monthrange(dtime_lst[0][0].year,dtime_lst[0][0].month)[1]
    sdate = datetime(dtime_lst[0][0].year,dtime_lst[0][0].month,1)
    edate = datetime(dtime_lst[0][0].year,dtime_lst[0][0].month,days_in_month,23)
    ax.plot([sdate,edate],[np.zeros(len(dtime_lst[0])),np.zeros(len(dtime_lst[0]))],
        color='lightgray', linestyle='-',lw=1)
    pltcolors = ['k','skyblue','orange']
    pltlw = [2,2,2]
    pltms = [5,3,3]
    for i in range(len(forecasts)):
        ax.plot(dtime_lst[i],ts_lst[i],linestyle='-',
                color=pltcolors[i],lw=pltlw[i])
        ax.plot(dtime_lst[i],ts_lst[i],'o',markersize=pltms[i],
                color=pltcolors[i],label=str(forecasts[i]) + "h")
    plt.ylabel(val_fig_varname_dict[val_name],fontsize=fs)
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gcf().autofmt_xdate()
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.ylim([val_fig_lim_dict[val_name][0],val_fig_lim_dict[val_name][1]])
    plt.xlim([sdate,edate])
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(filename_fig,format='png',dpi=50)
    #plt.show()
    return

def make_val_ts_fig_op(val_name,ts,dtime,filename_fig):
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(111)
    fs = 14
    days_in_month = calendar.monthrange(dtime[0].year,dtime[0].month)[1]
    sdate = datetime(dtime[0].year,dtime[0].month,1)
    edate = datetime(dtime[0].year,dtime[0].month,days_in_month,23)
    ax.plot([sdate,edate],[np.zeros(len(dtime)),np.zeros(len(dtime))],
        color='lightgray', linestyle='-',lw=1)
    ax.plot(dtime,ts,'k-',lw=1)
    ax.plot(dtime,ts,'ko',markersize=3,label=val_name)
    plt.ylabel(val_fig_varname_dict[val_name],fontsize=fs)
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gcf().autofmt_xdate()
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.ylim([val_fig_lim_dict[val_name][0],val_fig_lim_dict[val_name][1]])
    plt.xlim([sdate,edate])
    plt.tight_layout()
    plt.savefig(filename_fig,format='png',dpi=50)
    #plt.show()
    return

def make_val_scatter_fig_arcmfc(ts_model_lst,ts_obs_lst,
                                filename_fig,forecasts,i):
    fig = plt.figure(figsize=(4*1.25,3*1.25))
    ax = fig.add_subplot(111)
    fs = 15
    pltcolors = ['k','skyblue','orange']
    plt.plot(ts_obs_lst,ts_model_lst,'o',markersize=5,
            color=pltcolors[i],alpha=.8,
            markeredgecolor=pltcolors[0],
            label=str(forecasts[i]) + 'h')
    lmin=0.
    #lmax=np.nanmax(list(mHs)+list(sHs))+.5
    lmax=14
    plt.plot([lmin, lmax], [lmin,lmax], ls="--", c=".3")
    plt.ylabel('model',fontsize=fs)
    plt.xlabel('observations',fontsize=fs)
    plt.ylim([0,lmax])
    plt.xlim([0,lmax])
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(filename_fig,format='png',dpi=50)
    #plt.show()
    return

def make_val_scatter_fig_op(ts_model,ts_obs,filename_fig):
    fig = plt.figure(figsize=(16*2/3.,9*2/3.))
    ax = fig.add_subplot(111)
    fs = 15
    plt.plot(ts_obs,ts_model,'ko',markersize=5,alpha=.8)
    lmin=0.
    #lmax=np.nanmax(list(mHs)+list(sHs))+.5
    lmax=14
    plt.plot([lmin, lmax], [lmin,lmax], ls="--", c=".3")
    plt.ylabel('model',fontsize=fs)
    plt.xlabel('observations',fontsize=fs)
    plt.ylim([0,lmax])
    plt.xlim([0,lmax])
    plt.tight_layout()
    plt.savefig(filename_fig,format='png',dpi=50)
    #plt.show()
    return

def polyonmap(poly,proj,mtype=None,\
    region=None,llclat=None,llclon=None,trclat=None,trclon=None):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon
    import numpy as np
    # Make the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # create the map object, m
    if proj=='cyl':
        m = Basemap(resolution='i', projection='cyl', \
            llcrnrlon=llclon, llcrnrlat=llclat, \
            urcrnrlon=trclon, urcrnrlat=trclat)
        # choose map illustration type
        # Drawing ArcGIS Basemap (only works with cylc projections??)
        # Examples of what each map looks like can be found here:
        # http://kbkb-wx-python.blogspot.com/2016/04/
        #       python-basemap-background-image-from.htm
        maps = ['ESRI_Imagery_World_2D',    # 0
            'ESRI_StreetMap_World_2D',  # 1
            'NatGeo_World_Map',         # 2
            'NGS_Topo_US_2D',           # 3
            'Ocean_Basemap',            # 4
            'USA_Topo_Maps',            # 5
            'World_Imagery',            # 6
            'World_Physical_Map',       # 7
            'World_Shaded_Relief',      # 8
            'World_Street_Map',         # 9
            'World_Terrain_Base',       # 10
            'World_Topo_Map'            # 11
            ]
        print "drawing image from arcGIS server...",
        m.arcgisimage(service=maps[mtype], xpixels=1000, verbose=False)
    elif (proj=='merc' or proj =='lcc'):
        if proj=='merc':
            m = Basemap(resolution='i',projection=proj,area_thresh=10000,\
                llcrnrlon=llclon, llcrnrlat=llclat,\
                urcrnrlon=trclon, urcrnrlat=trclat,\
                lat_ts=trclat-(trclat-llclat))
        elif proj=='lcc':
            m = Basemap(resolution='i',projection=proj,area_thresh=10000,\
                width=10000000,height=8000000,\
                rsphere=(6378137.00,6356752.3142),\
                #lat_0=trclat-(trclat-llclat),\
                lat_0=60,\
                lat_1=llclat,\
                lat_2=trclat,\
                lon_0=0)
        m.drawcoastlines()
        # convert polygon to points on map
        x,y = m(poly.xy[:,0],poly.xy[:,1])
        poly = Polygon(list(zip(x,y)), closed=True)
    elif (proj=='auto'):
        m=quim(region)
        # convert polygon to points on map
        x,y = m(poly.xy[:,0],poly.xy[:,1])
        poly = Polygon(list(zip(x,y)), closed=True)
    # Fill polygon shape
    patches = []
    patches.append(poly)
    ax.add_collection(PatchCollection(patches, facecolor='lightgreen',
        alpha=0.4, edgecolor='k', linewidths=1.5))
    plt.savefig('polytest.png')

def quip_ia(sa_obj, region=None, save=None, outpath=None, filetype=None):
    # ignore irrelevant warnings from matplotlib for stdout
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap, cm
    if region is None:
        region = 'Global'
    if outpath is None:
        outpath = 'outpath/'
    if save == True:
        os.system('mkdir -p ' + outpath)
        plt.savefig(
                    outpath
                    + 's3a_'
                    + region
                    + "_"
                    + sa_obj.sdate.strftime("%Y%m%d%H%M%S")
                    + "_"
                    + sa_obj.edate.strftime("%Y%m%d%H%M%S")
                    + '.pdf',
                    format='pdf'
                    )
# ---------------------------------------------------------------------#


class graphics_class():
    '''
    class to handle graphical applications
    '''
    from region_specs import region_dict

    def __init__(self,sdate,edate=None,timewin=None,download=None,region=None,
                corenum=None,mode=None):
        print ('# ----- ')
        print (" ### Initializing graphics_class instance ###")
        print ('# ----- ')
        if region is None:
            region='Global'
        self.region = region
        print ("graphics_class object initialized")

# --- help ------------------------------------------------------------#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
help text is coming
        """,
        formatter_class = RawTextHelpFormatter
        )
    args = parser.parse_args()
