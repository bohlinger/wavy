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
from copy import deepcopy
import glob
import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import get_filevarname
from wavy.ncmod import build_xr_ds_from_dict
from wavy.utils import make_pathtofile, get_pathtofile
from wavy.utils import finditem, make_subdict
from wavy.utils import parse_date
from wavy.utils import convert_meteorologic_oceanographic
from wavy.utils import date_dispatcher
from wavy.utils import flatten
from wavy.wconfig import load_or_default, load_dir
from wavy.filtermod import filter_class as fc
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.init_class_insitu import init_class
# ---------------------------------------------------------------------#

# read yaml config files:
variable_def = load_or_default('variable_def.yaml')
# ---------------------------------------------------------------------#

def check_date(filelst, date):
    '''
    Checks if str in lst according to desired date (sd, ed)

    return: idx for file
    '''
    idx = []
    for i in range(len(filelst)):
        element = filelst[i]
        tmp = element.find(date.strftime('%Y%m%d'))
        if tmp >= 0:
            idx.append(i)
    if len(idx) <= 0:
        idx = [0]
    return idx[0], idx[-1]


# ---------------------------------------------------------------------#

class insitu_class(qls, fc):
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
        # check for twinID
        nID = kwargs.get('nID')
        # initializing useful attributes from config file
        if kwargs.get('twinID') is None:
            dc = init_class('insitu', kwargs.get('nID'))
        else:
            dc = init_class('insitu', kwargs.get('twinID'))
            dc.nID = nID
        # parse and translate date input
        self.twin = int(kwargs.get('twin', 0))
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        if self.twin is not None:
            self.sd = self.sd - timedelta(minutes=self.twin)
            self.ed = self.ed + timedelta(minutes=self.twin)
        print('Chosen period: ' + str(self.sd) + ' - ' + str(self.ed))
        # add other class object variables
        self.nID = kwargs.get('nID')
        self.name = kwargs.get('name',
                               list(dc.name.keys())[0])
        self.varalias = kwargs.get('varalias', 'Hs')
        self.stdvarname = variable_def[self.varalias]['standard_name']
        self.units = variable_def[self.varalias].get('units')
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')
        self.cfg = dc
        print(" ")
        print(" ### insitu_class object initialized ### ")
        print('# ----- ')

    def download(self, path=None, nproc=1, **kwargs):
        print('')
        print('Choosing collector..')
        # define reader
        dotenv.load_dotenv()
        WAVY_DIR = os.getenv('WAVY_DIR', None)
        if WAVY_DIR is None:
            print('#')
            print('Environmental variable for WAVY_DIR not defined')
            print('Defaults are chosen')
            print('#')
            collector_mod_str = load_dir('insitu_collectors').name
        else:
            collector_mod_str = WAVY_DIR + '/wavy/insitu_collectors.py'

        collector_str = kwargs.get('collector', self.cfg.collector)
        spec = importlib.util.spec_from_file_location(
                'insitu_collectors.' + collector_str, collector_mod_str)

        # create collector module
        collector_tmp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(collector_tmp)

        # pick reader
        collector = getattr(collector_tmp, collector_str)
        self.collector = collector
        print('Chosen collector:', spec.name)
        print('')

        print("Downloading files ...")
        nproc = kwargs.get('nproc', 1)
        self.collector(nproc=nproc, path=path, **(vars(self)))

    def _create_pathlst(self, **kwargs):
        src_tmplt = self.cfg.wavy_input['src_tmplt']
        fl_tmplt = self.cfg.wavy_input['fl_tmplt']
        pth_tmplt = src_tmplt + '/' + fl_tmplt
        strsub = self.cfg.wavy_input['strsub']
        #dict_for_sub = vars(self.cfg)
        dict_for_sub = vars(self)
        file_date_incr = self.cfg.wavy_input['file_date_incr']
        file_date_incr_unit = self.cfg\
                              .wavy_input['file_date_incr_unit']
        subdict = make_subdict(strsub, class_object_dict=dict_for_sub)
        # loop from sdate to edate with dateincr
        tmpdate = deepcopy(self.sd)
        pathlst = []

        while (tmpdate <= self.ed):
            # get pathtofile
            pathtofile = get_pathtofile(pth_tmplt, strsub,
                                        subdict, tmpdate)
            pathlst.append(pathtofile)
            tmpdate = date_dispatcher(tmpdate,
                                      file_date_incr_unit,
                                      file_date_incr)
        return pathlst

    def list_input_files(self, show=False, **kwargs):
        try:
            if (kwargs.get('path') is None and kwargs.get('wavy_path') is None):
                pathlst = self._create_pathlst(**kwargs)
            else:
                # if defined path local
                print(" ## Find and list files based on given path...")
                path = kwargs.get('path', None)
                wavy_path = kwargs.get('wavy_path', None)
                pathlst = self._get_files(vars(self),
                                          path=path,
                                          wavy_path=wavy_path)

                # remove None values from pathlst
                pathlst = list(filter(lambda item: item is not None, pathlst))
                print(str(int(len(pathlst))) + " valid files found")

        except Exception as e:
            logger.exception(e)
            pathlst = None
            print(' no netcdf meta data retrieved')

        if show is True:
            print(" ")
            print("pathlst:")
            print(pathlst)
            print(" ")
        return pathlst

    def _get_files(self, dict_for_sub=None, path=None, wavy_path=None):
        """
        Function to retrieve list of files/paths for available
        locally stored data. This list is used for other functions
        to query and parsing.

        param:
            sd - start date (datetime object)
            ed - end date (datetime object)
            nID - nID as of model_cfg.yaml
            dict_for_sub - dictionary for substitution in templates
            path - a path if defined

        return:
            pathtotals - list of paths
        """
        filelst = []
        pathlst = []
        tmpdate = self.sd
        if wavy_path is not None:
            pathtotals = [wavy_path]
            filelst = [wavy_path]
        elif path is None:
            print('path is None -> checking config file')
            while (tmpdate <= date_dispatcher(self.ed,
            self.cfg.misc['date_incr_unit'],
            self.cfg.misc['date_incr'])):
                try:
                    # create local path for each time
                    path_template = \
                            self.cfg.nID['wavy_input'].get('src_tmplt')
                    strsublst = \
                            self.cfg.nID['wavy_input'].get('strsub')
                    subdict = \
                        make_subdict(strsublst,
                                     class_object_dict=dict_for_sub)
                    path = make_pathtofile(path_template,
                                           strsublst, subdict)
                    path = tmpdate.strftime(path)
                    if os.path.isdir(path):
                        tmplst = np.sort(os.listdir(path))
                        filelst.append(tmplst)
                        pathlst.append([os.path.join(path, e)
                                        for e in tmplst])
                    path = None
                except Exception as e:
                    logger.exception(e)

                tmpdate = date_dispatcher(tmpdate,
                            self.cfg.misc['date_incr_unit'],
                            self.cfg.misc['date_incr'])

            filelst = np.sort(flatten(filelst))
            pathlst = np.sort(flatten(pathlst))

            # check if type iterable
            try:
                iter(pathlst)
                pathtotals = pathlst
            except TypeError:
                pathtotals = [pathlst]
            else:
                print("Object is iterable")

            # limit to sd and ed based on file naming, see check_date
            idx_start, tmp = check_date(filelst, self.sd)
            tmp, idx_end = check_date(filelst, self.ed)
            if idx_end == 0:
                idx_end = len(pathlst)-1
            del tmp
            pathtotals = np.unique(pathtotals[idx_start:idx_end+1])
            filelst = np.unique(filelst[idx_start:idx_end+1])

        else:
            if os.path.isdir(path):
                pathlst = glob.glob(path+'/*')
            else:
                pathlst = glob.glob(path+'*')
            #
            # interesting other approach using pathlibs glob
            # https://stackoverflow.com/questions/3348753/\
            #        search-for-a-file-using-a-wildcard
            # from pathlib import Path
            # filelst = [p.name for p in Path(path).glob("*")]
            # pathlst = [str(p.parent) for p in Path(path).glob("*")]
            #
            # separate files from path
            filelst = [p.split('/')[-1] for p in pathlst]
            pathlst = [p[0:-len(f)] for p, f in zip(pathlst, filelst)]
            pathtotals = [os.path.join(p, f) for p, f in zip(pathlst, filelst)]

            # limit to sd and ed based on file naming, see check_date
            idx_start, tmp = check_date(filelst, self.sd)
            tmp, idx_end = check_date(filelst, self.ed)
            if idx_end == 0:
                idx_end = len(pathlst)-1
            del tmp
            pathtotals = np.unique(pathtotals[idx_start:idx_end+1])
            filelst = np.unique(filelst[idx_start:idx_end+1])

        return pathtotals


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
            if ('convention' in self.cfg.misc.keys() and
            self.cfg.misc['convention'] == 'oceanographic'):
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
                                vars(self.cfg), self.meta)
        self.vars = self.vars.rename({ncvar: self.varalias})
        # coords
        coords = ['time', 'lons', 'lats']
        for c in coords:
            try:  # because only if available
                ncvar = get_filevarname(c, variable_def,
                                        vars(self.cfg), self.meta)
                self.vars = self.vars.rename({ncvar: c})\
                            .set_index(time='time')
            except Exception as e:
                logger.exception(e)
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

    @staticmethod
    def _return_extension(fstr: str):
        from pathlib import Path
        # if m.endswith('.mp3'):
        # for case sensitive and multiple checks
        # m.lower().endswith(('.png', '.jpg', '.jpeg'))
        return Path(fstr).suffix

    def populate(self, **kwargs):
        print(" ### Read files and populate insitu_class object")

        try:
            self.pathlst = self.list_input_files(**kwargs)

            # only possible if netcdf
            if (self._return_extension(self.pathlst[0]) == '.nc'
                or self._return_extension(self.pathlst[0]) == '.ncml'):
                self.meta = ncdumpMeta(self.pathlst[0])
            else:
                self.meta = None
        except Exception as e:
            logger.exception(e)
            self.meta = None
            print(' no netcdf meta data retrieved')

        print('')
        print('Choosing reader..')
        # define reader
        dotenv.load_dotenv()
        WAVY_DIR = os.getenv('WAVY_DIR', None)
        if WAVY_DIR is None:
            print('#')
            print('Environmental variable for WAVY_DIR not defined')
            print('Defaults are chosen')
            print('#')
            reader_mod_str = load_dir('insitu_readers').name
        else:
            reader_mod_str = WAVY_DIR + '/wavy/insitu_readers.py'

        reader_str = kwargs.get('reader', self.cfg.reader)
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

                if self.meta is not None:
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
                t1 = time.time()
                print(" ")
                print(' ## Summary:')
                print(str(len(self.vars['time'])) + " values retrieved.")
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

class poi_class(qls, fc):
    '''
    Class to handle poi based time series.
    '''
    def __init__(self, poi: dict, **kwargs):
        # parse and translate date input
        print('# ----- ')
        print(" ### Initializing poi_class object ###")
        print(" ")
        print(" Given kwargs:")
        print(kwargs)

        # check for nID
        self.nID = poi.get('nID', kwargs.get('nID'))
        if self.nID is None:
            self.nID = 'NA'

        # check if all need variables are in poi
        # those are: lons, lats, time
        for v in ['lons', 'lats', 'time']:
            if v not in list(poi.keys()):
                raise Exception(v, 'is needed yet not in poi')

        # ancillary
        if 'var' not in list(poi.keys()):
            print(' ', 'No variable specified in poi')
            poi['var'] = np.zeros(len(poi['time']))*np.nan

        # check if varalias is specified, default is Hs
        self.varalias = poi.get('varalias', kwargs.get('varalias', 'Hs'))
        # rename var to varalias
        poi[self.varalias] = poi.pop('var')

        # parse all dates and return as datetime objects
        poi_times = [parse_date(d) for d in poi['time']]
        poi['time'] = poi_times

        self.sd = poi['time'][0]
        self.ed = poi['time'][-1]

        print('Chosen period: ' + str(self.sd) + ' - ' + str(self.ed))

        # add other class object variables
        self.stdvarname = variable_def[self.varalias]['standard_name']
        self.units = variable_def[self.varalias].get('units')
        self.twin = int(kwargs.get('twin', 30))
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')
        self.name = kwargs.get('name', self.nID)

        # build xarray dataset
        self.vars = self._build_xr_ds(poi)

        print(" ")
        print(" ### insitu_class object initialized ### ")
        print('# ----- ')

    def _build_xr_ds(self, poi):
        ds = build_xr_ds_from_dict(poi, 'time')
        return ds

    def crop_to_period(self, **kwargs):
        """
        Function to crop the variable dictionary to a given period
        """
        new = deepcopy(self)
        sd = parse_date(kwargs.get('sd', str(new.sd)))
        ed = parse_date(kwargs.get('ed', str(new.ed)))
        print('Crop to time period:', sd, 'to', ed)
        new.vars = new.vars.sel(time=slice(sd, ed))
        new.sd = sd
        new.ed = ed
        return new

