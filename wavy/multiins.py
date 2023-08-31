#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from multiple insitu locations. 
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import numpy as np

# own imports
from wavy.utils import parse_date
from wavy.wconfig import load_or_default
from wavy.insitu_module import insitu_class as ic
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.consolidate import consolidate_class as cs
from wavy.utils import find_tagged_obs, expand_nID_for_sensors
# ---------------------------------------------------------------------#

# read yaml config files:
insitu_dict = load_or_default('insitu_specs.yaml')
# ---------------------------------------------------------------------#


class multiins_class(qls):
    '''
    Class to handle insitu based time series.
    '''

    def __init__(
    self, sdate, edate,
    nID = None, sensor = None,
    varalias = 'Hs', filterData = False, tags = None,
    **kwargs ):
        print('# ----- ')
        print(" ### Initializing multiins_class object ###")
        print(" ")
        self.obstype = 'insitu'
        if tags is None:
            # multiple nIDs
            nIDs = nID
            # products: either None, same as missions, or one product
            sensors = sensor
            # check if nIDs and sensors have same length
            assert len(nIDs) == len(sensors)
        else:
            nID = find_tagged_obs(tags,self.obstype)
            sensors = []
            nIDs = []
            for n in nID:
                nsensors = expand_nID_for_sensors(n,self.obstype)
                sensors += nsensors
                nIDs += [n]*len(nsensors)
            # check if nIDs and sensors have same length
            assert len(nIDs) == len(sensors)
        # retrieve
        icos = []
        for i,n in enumerate(nIDs):
            ico = ic( n,
                      sdate, edate,
                      varalias = varalias,
                      filterData = filterData,
                      sensor = sensors[i],
                      **kwargs )
            el = list(vars(ico).keys())
            if 'error' in el:
                print("Insitu location",n,"is not available and not appended")
                pass
            else:
                icos.append( ico )
        cso = cs(icos)
        # class object variables
        self.obsname = cso.obsname
        self.stdvarname = cso.stdvarname
        self.varalias = cso.varalias
        self.varname = cso.varname
        self.mission = cso.mission
        self.product = cso.product
        self.provider = cso.provider
        self.sdate = cso.sdate
        self.edate = cso.edate
        self.units = cso.units
        self.vars = cso.vars
        self.ocos = cso.ocos
        self.label = "multi-insitu-observations"

        print(" ")
        print (" ### multiins object initialized ###")
        print ('# ----- ')
