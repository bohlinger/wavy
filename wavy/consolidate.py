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
import xarray as xr

# own imports
from wavy.satellite_module import satellite_class as sco
from wavy.filtermod import filter_class as fc
from wavy.wconfig import load_or_default
from wavy.quicklookmod import quicklook_class_sat as qls
# ---------------------------------------------------------------------#

# read yaml config files:
variable_info = load_or_default('variable_def.yaml')

# --- global functions ------------------------------------------------#


def consolidate_ocos(ocos):
    """
    consolidate sco.vars:
        'sea_surface_wave_significant_height', 'time',
        'latitude', 'longitude', 'datetime'
    """
    ds_lst = [oco.vars for oco in ocos if 'vars' in list(vars(oco).keys())]
    ds = xr.concat(ds_lst, dim='time')
    return ds

def find_valid_oco(ocos):
    state = False
    for i, o in enumerate(ocos):
        if 'vars' in vars(o).keys():
            state = True
            break
    if state is True:
        return i
    else:
        print("")
        print("Caution: no valid observation object with values found!")
        print("")


# --------------------------------------------------------------------#

class consolidate_class(qls, fc):
    '''
    Class to consolidate multiple wavy objects
    '''
    def __init__(self, ocos, **kwargs):
        print('# ----- ')
        print(" ### Initializing consolidate_class object ###")
        print(" ")
        i = find_valid_oco(ocos)
        self.ocos = ocos
        self.varalias = ocos[i].varalias
        self.stdvarname = ocos[i].stdvarname
        self.units = ocos[i].units
        self.sd = ocos[i].sd
        self.ed = ocos[i].ed
        self.vars = consolidate_ocos(ocos)
        self.twin = kwargs.get('twin', 30)
        self.distlim = kwargs.get('distlim', 6)
        self.region = kwargs.get('region', 'global')
        self.obsname = 'consolidated-obs'
        self.obstype = 'consolidated-obs'
        self.label = 'consolidated-obs'
        self.name = 'consolidated-names'
        self.nID = 'consolidated-nIDs'
        self.sensor = 'consolidated-sensors'
        print(" ")
        print(" ### consolidate_class object initialized ###")
        print('# ----- ')

    def rename_consolidate_object_parameters(self, **kwargs):
        if kwargs.get('obsname') is not None:
            self.obsname = kwargs.get('obsname')
        if kwargs.get('mission') is not None:
            self.mission = kwargs.get('mission')
        if kwargs.get('obstype') is not None:
            self.obstype = kwargs.get('obstype')
        if kwargs.get('product') is not None:
            self.product = kwargs.get('product')
        if kwargs.get('nID') is not None:
            self.nID = kwargs.get('nID')
        if kwargs.get('sensor') is not None:
            self.sensor = kwargs.get('sensor')
        return
