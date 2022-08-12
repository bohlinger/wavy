#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to consolidate multiple sources of 
observations to perform a collective collocation/analysis
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import numpy as np
from datetime import datetime

# own imports
from wavy.utils import flatten
from wavy.satmod import satellite_class
from wavy.insitumod import insitu_class
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#

# read yaml config files:
variable_info = load_or_default('variable_info.yaml')

# --- global functions ------------------------------------------------#

def consolidate_obs():
    return

def consolidate_scos(scos):
    """
    consolidate sco.vars: 
        'sea_surface_wave_significant_height', 'time',
        'latitude', 'longitude', 'datetime'
    """
    varlst = []
    lonlst = []
    latlst = []
    timelst = []
    dtimelst = []
    for sco in scos:
        varlst.append(sco.vars[sco.stdvarname])
        lonlst.append(sco.vars['longitude'])
        latlst.append(sco.vars['latitude'])
        timelst.append(sco.vars['time'])
        dtimelst.append(sco.vars['datetime'])
    # flatten all, make arrays
    varlst = np.array(flatten(varlst))
    lonlst = np.array(flatten(lonlst))
    latlst = np.array(flatten(latlst))
    timelst = np.array(flatten(timelst))
    dtimelst = np.array(flatten(dtimelst))
    # sort according to time
    idx = np.argsort(timelst)
    # make dict
    vardict = {
            sco.stdvarname:varlst[idx],
            'longitude':lonlst[idx],
            'latitude':latlst[idx],
            'time':timelst[idx],
            'time_unit':sco.vars['time_unit'],
            'datetime':dtimelst[idx] }
    return vardict

def consolidate_icos():
    return

class consolidate_class():
    '''
    Class to handle multiple satellite_class objects
    '''
    def __init__(self,ocos):
        print('# ----- ')
        print(" ### Initializing consolidate_class object ###")
        print(" ")
        self.ocos = ocos
        self.varalias = ocos[0].varalias
        self.stdvarname = ocos[0].stdvarname
        self.varname = ocos[0].varname
        self.units = ocos[0].units
        self.sdate = ocos[0].sdate
        self.edate = ocos[0].edate
        if isinstance(ocos[0],satellite_class):
            self.vars = consolidate_scos(ocos)
        elif isinstance(ocos[0],insitu_class):
            self.vars = consolidate_icos(ocos)
        print(" ")
        print (" ### consolidate_class object initialized ###")
        print ('# ----- ')

    def rename_consolidate_object_parameters(self,**kwargs):
        # obsname, mission, obstype, nID, sensor
        self.obsname = kwargs.get('obsname','Consolidated observations')
        self.mission = kwargs.get('mission','')
        self.obstype = kwargs.get('obstype','')
        self.nID = kwargs.get('nID','')
        self.sensor = kwargs.get('sensor','')
        return

    def blend_obs_types(self):
        return

    def quicklook(self,a=False,projection=None,**kwargs):
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
        # set plots
        m = kwargs.get('m',a)
        ts = kwargs.get('ts',a)
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
            vartype = variable_info[self.varalias].get('type','default')
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
            ax.set_title( 'Consolidated obs\n'
                      + 'from ' + str(self.vars['datetime'][0])
                      + ' to ' + str(self.vars['datetime'][-1]))
            #fig.suptitle('', fontsize=16) # unused
            plt.show()
        if ts:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(9,3.5))
            ax = fig.add_subplot(111)
            colors = ['k']
            ax.plot(self.vars['datetime'],
                    self.vars[self.stdvarname],
                    linestyle='None',color=colors[0],
                    label='consolidated obs',
                    marker='o',alpha=.5,ms=2)
            plt.ylabel(self.varalias + ' [' + self.units + ']')
            plt.legend(loc='best')
            plt.tight_layout()
            #ax.set_title()
            plt.show()
