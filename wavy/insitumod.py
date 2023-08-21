#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from stations. I try to mostly follow the PEP
convention for python code style. Constructive comments on style and
effecient programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import numpy as np
from datetime import datetime, timedelta
import time
import os
from dateutil.relativedelta import relativedelta
import importlib.util
import dotenv
import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import dumptonc_ts_insitu
from wavy.ncmod import get_filevarname
from wavy.utils import make_pathtofile, get_pathtofile
from wavy.utils import finditem, make_subdict
from wavy.utils import parse_date
from wavy.utils import convert_meteorologic_oceanographic
from wavy.wconfig import load_or_default
from wavy.writermod import writer_class as wc
from wavy.filtermod import filter_class as fc
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.init_class_insitu import init_class
# ---------------------------------------------------------------------#

# read yaml config files:
insitu_dict = load_or_default('insitu_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')
# ---------------------------------------------------------------------#


class insitu_class(qls, wc, fc):
    '''
    Class to handle insitu based time series.
    '''
    def __init__(self, **kwargs):
        # parse and translate date input
        print('# ----- ')
        print(" ### Initializing insitu_class object ###")
        print(" ")
        print(" Given kwargs:")
        print(kwargs)
        # initializing useful attributes from config file
        dc = init_class('insitu', kwargs.get('nID'))
        # parse and translate date input
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        print('Chosen period: ' + str(self.sd) + ' - ' + str(self.ed))
        # add other class object variables
        self.nID = kwargs.get('nID')
        self.sensor = kwargs.get('sensor')
        self.varalias = kwargs.get('varalias', 'Hs')
        self.stdvarname = variable_def[self.varalias]['standard_name']
        self.units = variable_def[self.varalias].get('units')
        self.twin = int(kwargs.get('twin', 30))
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')
        self.cfg = dc
        print(" ")
        print(" ### insitu_class object initialized ### ")
        print('# ----- ')

    def list_input_files(self, show=False, **kwargs):
        print(" ## Find and list files ...")
        path = kwargs.get('path', None)
        wavy_path = kwargs.get('wavy_path', None)
        pathlst, _ = self._get_files(vars(self),
                                     path=path,
                                     wavy_path=wavy_path)
        print('source template:',
              insitu_dict[self.nID]['wavy_input']['src_tmplt'])
        if show is True:
            print(" ")
            print(pathlst)
            print(" ")
        return pathlst

    def _get_insitu_ts(self, **kwargs):
        """
        Main function to obtain data from insitu locations.
        reads files, apply region and temporal filter

        return: adjusted dictionary according to spatial and
                temporal constraints
        """

        # retrieve dataset
        ds = self.reader(**(vars(self)))
        self.vars = ds
        self.coords = list(self.vars.coords)
        return self

    @staticmethod
    def _enforce_longitude_format(ds):
        # adjust longitude -180/180
        attrs = ds.lons.attrs
        attrs['valid_min'] = -180
        attrs['valid_max'] = 180
        attrs['comments'] = 'forced to range: -180 to 180'
        ds.lons.values = ((ds.lons.values-180) % 360)-180
        return ds

    def _enforce_meteorologic_convention(self):
        ncvars = list(self.vars.variables)
        for ncvar in ncvars:
            if ('convention' in insitu_dict[self.nID].keys() and
            insitu_dict[self.nID]['convention'] == 'oceanographic'):
                print('Convert from oceanographic to meteorologic convention')
                self.vars[ncvar] =\
                    convert_meteorologic_oceanographic(self.vars[ncvar])
            elif 'to_direction' in self.vars[ncvar].attrs['standard_name']:
                print('Convert from oceanographic to meteorologic convention')
                self.vars[ncvar] =\
                    convert_meteorologic_oceanographic(self.vars[ncvar])

        return self

    def _change_varname_to_aliases(self):
        # variables
        ncvar = get_filevarname(self.varalias, variable_def,
                                insitu_dict[self.nID], self.meta)
        self.vars = self.vars.rename({ncvar: self.varalias})
        # coords
        coords = ['time', 'lons', 'lats']
        for c in coords:
            ncvar = get_filevarname(c, variable_def,
                                    insitu_dict[self.nID], self.meta)
            self.vars = self.vars.rename({ncvar: c}).set_index(time='time')
        return self

    def _change_stdvarname_to_cfname(self):
        # enforce standard_name for coordinate aliases
        self.vars['lons'].attrs['standard_name'] = \
            variable_def['lons'].get('standard_name')
        self.vars['lats'].attrs['standard_name'] = \
            variable_def['lats'].get('standard_name')
        self.vars['time'].attrs['standard_name'] = \
            variable_def['time'].get('standard_name')
        # enforce standard_name for variable alias
        self.vars[self.varalias].attrs['standard_name'] = \
            self.stdvarname
        return self

    def populate(self, **kwargs):
        print(" ### Read files and populate insitu_class object")

        # list input files or data origin if remote
        #lst = self.list_input_files(**kwargs)
        #self.pathlst = lst

        #print('')
        #print('Checking variables..')
        #self.meta = ncdumpMeta(self.pathlst[0])
        #ncvar = get_filevarname(self.varalias, variable_def,
        #                        insitu_dict[self.nID], self.meta)

        print('')
        print('Choosing reader..')
        # define reader
        dotenv.load_dotenv()
        WAVY_DIR = os.getenv('WAVY_DIR', None)
        if WAVY_DIR is None:
            print('###########')
            print('Environmental variable for WAVY_DIR needs to be defined!')
            print('###########')
        reader_str = kwargs.get('reader', self.cfg.reader)
        reader_mod_str = WAVY_DIR + '/wavy/insitu_readers.py'
        spec = importlib.util.spec_from_file_location(
                'insitu_readers.' + reader_str, reader_mod_str)

        # create reader module
        reader = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reader)

        # pick reader
        reader = getattr(reader, reader_str)
        self.reader = reader
        print('Chosen reader:', spec.name)
        print('')

        #if len(lst) > 0:
        if True:
            try:
                t0 = time.time()
                print('Reading..')
                self = self._get_insitu_ts(**kwargs)

                self = self._change_varname_to_aliases()
                self = self._change_stdvarname_to_cfname()
                self = self._enforce_meteorologic_convention()

                # convert longitude
                ds = self.vars
                ds_new = self._enforce_longitude_format(ds)
                self.vars = ds_new

                # adjust varalias if other return_var
                if kwargs.get('return_var') is not None:
                    newvaralias = kwargs.get('return_var')
                else:
                    newvaralias = self.varalias

                # define more class object variables
                if kwargs.get('return_var') is not None:
                    self.varalias = kwargs.get('return_var')
                    self.stdvarname = \
                        variable_def[newvaralias].get('standard_name')
                    self.units = variable_def[newvaralias].get('units')
                # create label for plotting
                self.label = self.mission
                t1 = time.time()
                print(" ")
                print(' ## Summary:')
                print(str(len(self.vars['time'])) + " footprints retrieved.")
                print("Time used for retrieving data:")
                print(round(t1-t0, 2), "seconds")
                print(" ")
                print(" ### insitu_class object populated ###")
                print('# ----- ')
            except Exception as e:
                logger.exception(e)
                #logger.debug(traceback.format_exc())
                print(e)
                print('Error encountered')
                print('insitu_class object not populated')
        else:
            print('No data found')
            print('insitu_class object not populated')
            print('# ----- ')
        return self

    def get_item_parent(self, item, attr):
        ncdict = self.vars['meta']
        lst = [i for i in ncdict.keys()
               if (attr in ncdict[i].keys()
               and item in ncdict[i][attr])
               ]
        if len(lst) >= 1:
            return lst
        else:
            return None

    def get_item_child(self, item):
        ncdict = self.vars['meta']
        parent = finditem(ncdict, item)
        return parent


#        if sensor is None:
#            sensor = list(insitu_dict[nID]['sensor'].keys())[0]
#            print('Sensor was not chosen')
#            print('Automatic choice:', sensor)
#        # for i in range(1):
#        try:
#            self.stdvarname = stdvarname
#            self.varalias = varalias
#            self.units = variable_info[varalias].get('units')
#            self.sensor = sensor
#            self.obstype = 'insitu'
#            if ('tags' in insitu_dict[nID].keys() and
#            len(insitu_dict[nID]['tags']) > 0):
#                self.tags = insitu_dict[nID]['tags']
#            print(" ")
#            print(" ## Read files ...")
#            t0 = time.time()
#            vardict, fifo, pathtofile = \
#                get_insitu_ts(nID=nID, sensor=sensor,
#                              sd=sd, ed=ed,
#                              varalias=varalias,
#                              basedate=self.basedate,
#                              dict_for_sub=vars(self),
#                              **kwargs)
#            self.vars = vardict
#            self.varname = varalias
#            if fifo == 'frost':
#                self.sensor = sensor
#            # create label for plotting
#            self.label = self.nID + '_' + self.sensor
#            t1 = time.time()
#            print(" ")
#            print('## Summary:')
#            print(str(len(self.vars['time'])) + " values retrieved.")
#            print("Time used for retrieving insitu data:",
#                   round(t1-t0, 2), "seconds")
#            print(" ")
#            print(" ### insitu_class object initialized ### ")
#        except Exception as e:
#            logger.exception(e)
#            print(e)
#            self.error = e
#            print("! No insitu_class object initialized !")
#        print('# ----- ')

    def get_item_parent(self, item, attr):
        ncdict = self.vars['meta']
        lst = [i for i in ncdict.keys() \
                if (attr in ncdict[i].keys() \
                and item in ncdict[i][attr]) \
                ]
        if len(lst) >= 1:
            return lst
        else:
            return None

    def get_item_child(self, item):
        ncdict = self.vars['meta']
        parent = finditem(ncdict, item)
        return parent
