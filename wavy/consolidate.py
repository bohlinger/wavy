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
from wavy.quicklookmod import quicklook_class_sat as qls
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

class consolidate_class(qls):
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
        self.obsname = 'Consolidated_Observations'
        self.obstype = 'Consolidated_Observations'
        self.mission = 'mission'
        self.product = 'product'
        self.provider = 'provider'
        self.nID = 'nID'
        self.sensor = 'sensor'
        print(" ")
        print (" ### consolidate_class object initialized ###")
        print ('# ----- ')

    def rename_consolidate_object_parameters(self,**kwargs):
        if kwargs.get('obsname') is not None:
            self.obsname = kwargs.get('obsname')
        if kwargs.get('mission') is not None:
            self.mission = kwargs.get('mission')
        if kwargs.get('obstype') is not None:
            self.obstype = kwargs.get('obstype')
        if kwargs.get('product') is not None:
            self.product = kwargs.get('product')
        if kwargs.get('provider') is not None:
            self.provider = kwargs.get('provider')
        if kwargs.get('nID') is not None:
            self.nID = kwargs.get('nID')
        if kwargs.get('sensor') is not None:
            self.sensor = kwargs.get('sensor')
        return

    def blend_obs_types(self):
        return
