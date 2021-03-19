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

import sys

# all class
import numpy as np
from datetime import datetime, timedelta
import argparse
from argparse import RawTextHelpFormatter
import os
import calendar
from copy import deepcopy
from mpl_toolkits import axes_grid1
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# --- global functions ------------------------------------------------#

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

def make_val_ts_fig_op(val_name,ts_lst,dtime_lst,filename_fig,forecasts):
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
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(111)
    fs = 14
    days_in_month = calendar.monthrange(dtime_lst[0][0].year,
                                        dtime_lst[0][0].month)[1]
    sdate = datetime(dtime_lst[0][0].year,dtime_lst[0][0].month,1)
    edate = datetime(dtime_lst[0][0].year,dtime_lst[0][0].month,
                    days_in_month,23)
    ax.plot([sdate,edate],[np.zeros(len(dtime_lst[0])),
                            np.zeros(len(dtime_lst[0]))],
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

def make_val_scatter_fig_op(ts_model_lst,ts_obs_lst,
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

# ---------------------------------------------------------------------#


class graphics_class():
    '''
    class to handle graphical applications
    '''
    def __init__(self,sdate,edate=None,twin=None,download=None,region=None,
                corenum=None,mode=None):
        print ('# ----- ')
        print (" ### Initializing graphics_class instance ###")
        print ('# ----- ')
        if region is None:
            region='global'
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
