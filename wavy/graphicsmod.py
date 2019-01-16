#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
Module to handle all visualization related functions.

I try to mostly follow the PEP convention for python code style. 
Constructive comments on style and effecient programming are most welcome!

'''
__version__ = "0.5.0"
__author__="Patrik Bohlinger, Norwegian Meteorological Institute"
__maintainer__ = "Patrik Bohlinger"
__email__ = "patrikb@met.no"
__status__ = "under development with operation ARCMFC branch"

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

val_fig_lim_dict =  {'SI':[0,50],
                    'msd':[0,1],
                    'corr':[-1,1],
                    'rmsd':[0,1],
                    'mor':[0,10],
                    'nov':[0,700],
                    'mop':[0,10],
                    'bias':[-.6,.6],
                    'mad':[0,1]
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


#def make_coll_ts_fig_arcmfc(name,ts,dtime):
#    import matplotlib.dates as mdates
#    fig = plt.figure(figsize=(20,5))
#    ax = fig.add_subplot(111)
#    fs = 14
#    ax.plot(all_dates,ts,'ko',markersize=6,label='Obs')
#    ax.plot(all_dates,ts,'ro',markersize=1.5,markeredgecolor='r',label='Model')
#    plt.legend(fontsize=fs,loc='upper left')
#    plt.ylabel('Hs [m]',fontsize=fs)
#    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=4))
#    #plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
#    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
#    plt.gcf().autofmt_xdate()
#    plt.tick_params(axis='both', which='major', labelsize=fs)
#    #plt.show()
#    return fig

def make_val_ts_fig_arcmfc(val_name,ts,dtime,filename_fig):
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

def make_val_scatter_fig_arcmfc(ts_model,ts_obs,filename_fig):
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
