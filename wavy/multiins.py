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
from wavy.insitumod import insitu_class as ic
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.consolidate import consolidate_class as cs
# ---------------------------------------------------------------------#

# read yaml config files:
insitu_dict = load_or_default('insitu_specs.yaml')
# ---------------------------------------------------------------------#


class multiins_class(qls):
    '''
    Class to handle insitu based time series.
    '''

    def __init__(
    self, nID, sensor, sdate, edate,
    varalias = 'Hs', filterData = False,
    **kwargs ):
        print('# ----- ')
        print(" ### Initializing multiins_class object ###")
        print(" ")
        # multiple nIDs
        nIDs = nID
        # products: either None, same as missions, or one product
        sensors = sensor
        # check if nIDs and sensors have same length
        assert len(nIDs) == len(sensors)
        # retrieve
        icos = []
        for i,n in enumerate(nIDs):
            icos.append( ic( n, sensors[i],
                             sdate, edate,
                             varalias = varalias,
                             filterData = filterData,
                             **kwargs ) )
        cso = cs(icos)
        cso.rename_consolidate_object_parameters(obstype='insitu')
        # class variables
        self.obsname = cso.obsname
        self.stdvarname = cso.stdvarname
        self.varalias = cso.varalias
        self.varname = cso.varname
        self.obstype = cso.obstype
        self.mission = cso.mission
        self.product = cso.product
        self.provider = cso.provider
        self.sdate = cso.sdate
        self.edate = cso.edate
        self.units = cso.units
        self.vars = cso.vars
        self.ocos = cso.ocos

        print(" ")
        print (" ### multiins object initialized ###")
        print ('# ----- ')


