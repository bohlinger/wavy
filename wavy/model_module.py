#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to acquire, read, and prepare
geophysical variables from model output files for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import time
from functools import lru_cache
from tqdm import tqdm
import importlib.util
import dotenv
import os
import glob
from copy import deepcopy

import logging

# own imports
from wavy.utils import hour_rounder, make_fc_dates
from wavy.utils import finditem, parse_date
from wavy.utils import convert_meteorologic_oceanographic
from wavy.utils import date_dispatcher
from wavy.utils import find_direction_convention
from wavy.utils import flatten
from wavy.utils import make_pathtofile, make_subdict

from wavy.credentials import get_credentials

from wavy.ncmod import check_if_ncfile_accessible
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.ncmod import build_usr_pw_path

from wavy.wconfig import load_or_default, load_dir

from wavy.quicklookmod import quicklook_class_sat as qls

from wavy.init_class_mod import init_class

# ---------------------------------------------------------------------#

def crop_to_period(ds, sd, ed):
    """
    Function to crop the dataset to a given period
    """
    ds_sliced = ds.sel(time=slice(sd, ed))
    return ds_sliced

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

def generate_bestguess_leadtime(model, fc_date, lidx=None, **kwargs):
    """
    fct to return leadtimes for bestguess
    """
    logger = logging.getLogger(__name__)
    log_level = str(kwargs.get('logging', 'WARNING').upper())
    logger.setLevel(getattr(logging, log_level, logging.WARNING))

    if isinstance(fc_date, list):
        leadtime = \
            [generate_bestguess_leadtime(model, date, **kwargs)
                for date in fc_date]
    else:
        init_times = \
            np.array(model_dict[model]['misc']['init_times']).astype('float')
        diffs = fc_date.hour - np.array(init_times)

        # greater than zero
        gtz = diffs[diffs >= 0]

        if len(gtz) == 0:
            leadtime = int(np.abs(np.min(np.abs(diffs))
                           - model_dict[model]['misc']['init_step']))
        elif (len(gtz) > 0 and lidx is not None):
            leadtime = int(np.sort(diffs[diffs >= 0])[lidx])
        else:
            leadtime = int(np.min(diffs[diffs >= 0]))

        # multiple of date_incr
        if leadtime == 0:
            return leadtime
        elif leadtime % model_dict[model]['misc']['date_incr'] == 0:
            return leadtime
        elif (leadtime == 1 and model_dict[model]['misc']['date_incr'] == 1):
            return leadtime
        else:
            logger.warning('Lead time ' + str(leadtime) + 'h'
                           + ' not available for model ' + model)

# function used by satellite_module
# should probably be placed somewhere else
@lru_cache(maxsize=32)
def read_model_nc_output_lru(filestr, lonsname, latsname, timename):
    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    filestr = filestr.replace('\\', '')
    f = netCDF4.Dataset(filestr, 'r')
    # get coordinates and time
    model_lons = f.variables[lonsname][:]
    model_lats = f.variables[latsname][:]
    model_time = f.variables[timename]
    model_time_dt = list(
        netCDF4.num2date(model_time[:], units=model_time.units))
    f.close()
    return model_lons, model_lats, model_time_dt

# ---------------------------------------------------------------------#

# read yaml config files:
model_dict = load_or_default('model_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')


class model_class(qls):
    '''
    class to read and process model data
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and
    station classes.
    '''
    def __init__(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        logger.info('# ----- ')
        logger.info(" ### Initializing model_class object ###")
        logger.info(" ")
        logger.info(" Given kwargs:")
        logger.info(kwargs)

        # initializing useful attributes from config file
        dc = init_class('model', kwargs.get('nID'))
        # parse and translate date input
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        logger.info('Chosen period: ' + str(self.sd) + ' - ' + str(self.ed))

        # add other class object variables
        self.nID = kwargs.get('nID')
        self.model = kwargs.get('model', self.nID)
        self.varalias = kwargs.get('varalias', ['Hs'])
        if isinstance(self.varalias, str):
            self.varalias = [self.varalias]
        self.units = [variable_def[v].get('units') for v in self.varalias]
        self.stdvarname = [variable_def[v].get('standard_name') for v in\
                           self.varalias]                             
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')
        self.leadtime = kwargs.get('leadtime', 'best')
        self.cfg = dc

        logger.info(" ")
        logger.info(" ### model_class object initialized ### ")
        logger.info('# ----- ')


    def crop_to_period(self, **kwargs):
        """
        Function to crop the variable dictionary to a given period
        """
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        new = deepcopy(self)
        sd = parse_date(kwargs.get('sd', str(new.sd)))
        ed = parse_date(kwargs.get('ed', str(new.ed)))
        logger.info('Crop to time period:', sd, 'to', ed)
        new.vars = new.vars.sel(time=slice(sd, ed))
        new.sd = sd
        new.ed = ed
        return new

    def get_item_parent(self, item, attr):
        """
        Offers possibility to explore netcdf meta info.
        by specifying what you are looking for (item),
        e.g. part of a string, and in which attribute (attr),
        e.g. standard_name, this function returns the
        parent parameter name of the query string.

        param:
            item - (partial) string e.g. [m]
            attr - attribute e.g. units

        return: list of matching parameter strings

        e.g. for satellite_class object sco:

        .. code ::

            sco.get_item_parent('m','units')
        """

        lst = [i for i in self.meta.keys()
               if (attr in self.meta[i].keys()
               and item in self.meta[i][attr])
               ]
        if len(lst) >= 1:
            return lst
        else:
            return None

    def get_item_child(self, item):
        """
        Gets all attributes connected to given parameter name.

        param:
            item - (partial) string e.g. [m]

        return: matching parameter string

        e.g. for satellite_class object sco:

        .. code ::

            sco.get_item_child('time')
        """

        parent = finditem(self.meta, item)
        return parent

    def _get_model_filedate(self, fc_date, leadtime, **kwargs):
        '''
        get init_date for latest model output file and checks if available

        param:
            fc_date - datetime object
            leadtime - integer in hours

        return:
            suitable datetime to create model filename
        '''
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        if ('init_times' in vars(self.cfg)['misc'].keys()
                and vars(self.cfg)['misc']['init_times'] is not None):
            init_times = \
                np.array(vars(self.cfg)['misc']['init_times']).astype('float')
        else:
            logger.warning(
                'init_times for chosen model not specified in config file')
            logger.warning(
                'Assuming continuous simulation with hourly values')
            init_times = np.array(range(25)).astype('float')
        date = fc_date - timedelta(hours=leadtime)
        date_hour = hour_rounder(date).hour
        logger.info('date:', date)
        logger.info('date_hour:', date_hour)
        if date_hour in init_times:
            logger.info('Leadtime', leadtime, \
                        'available for date', fc_date)
            init_diffs = date_hour - init_times
            init_diffs[init_diffs < 0] = np.nan
            h_idx = np.where(init_diffs == \
                            np.min(init_diffs[~np.isnan(init_diffs)]))
            h = int(init_times[h_idx[0][0]])
            return datetime(date.year, date.month, date.day, h)
        else:
            logger.warning('leadtime for fc_date not available')
            logger.warning('-> returning None')
            return None

    def _make_model_filename(self, fc_date, leadtime, **kwargs):
        """
        creates/returns filename based on fc_date,leadtime

            param:
            fc_date - datetime object
            leadtime - integer in hours

        return:
            filename (consists of path + filename)

        comment:
                - special characters are escaped by adding "\\"
                - the escapes need to be removed for certain libraries
                  like xarray and netCDF4
        """
        if self.nID in model_dict:
            if 'xtra_h' in vars(self.cfg)['misc']:
                filedate = self._get_model_filedate(
                                fc_date, leadtime, **kwargs)
                pathdate = filedate + timedelta(hours=leadtime) \
                                    * vars(self.cfg)['misc']['lt_switch_p']
                tmpstr = vars(self.cfg)['wavy_input']['fl_tmplt']
                for i in range(vars(self.cfg)['misc']['nr_filedates']):
                    filedatestr = vars(self.cfg)['misc']['filedate_formats'][i]
                    replacestr = (filedate \
                                + timedelta(hours=leadtime\
                                    - (leadtime % \
                                        vars(self.cfg)['misc']['init_step']))
                                    * vars(self.cfg)['misc']['lt_switch_f'][i]
                                + timedelta(hours=\
                                      vars(self.cfg)['misc']['xtra_h'][i])).\
                                      strftime(filedatestr)
                    tmpstr = tmpstr.replace('filedate', replacestr, 1)
                filename = (
                            pathdate.strftime(\
                                vars(self.cfg)['wavy_input']['src_tmplt'])
                            + tmpstr)
            else:
                filedate = self._get_model_filedate(
                                fc_date, leadtime, **kwargs)
                if filedate is None:
                    filename = None
                else:
                    filename = (filedate.strftime(vars(self.cfg)
                                    ['wavy_input']['src_tmplt'])
                              + filedate.strftime(vars(self.cfg)
                                    ['wavy_input']['fl_tmplt']))
        else:
            raise ValueError("Chosen model is not specified in model_cfg.yaml")
        # replace/escape special characters
        if filename is not None:
            filename = filename.replace(" ", "\\ ")\
                               .replace("?", "\\?")\
                               .replace("&", "\\&")\
                               .replace("(", "\\(")\
                               .replace(")", "\\)")\
                               .replace("*", "\\*")\
                               .replace("<", "\\<")\
                               .replace(">", "\\>")
        return filename

    def _make_model_filename_wrapper(self, fc_date, leadtime, **kwargs):
        """
        Wrapper function of make_model_filename. Organizes various cases.

        param:
            model - modelname type(str)
            fc_date - datetime object
            leadtime - integer in hours

        return:
            filename
        """
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        remoteHostName = kwargs.get('remoteHostName',
                                    self.cfg.misc.get('remoteHostName'))

        if remoteHostName is not None:
            usr, pw = get_credentials(remoteHostName)

        if leadtime is None:
            leadtime = 'best'

        if (isinstance(fc_date, datetime) and leadtime != 'best'):
            filename = self._make_model_filename(fc_date, leadtime, **kwargs)
        elif (isinstance(fc_date, datetime) and leadtime == 'best'):
            switch = False
            leadtime = generate_bestguess_leadtime(self.nID, fc_date,
                                                   **kwargs)
            if leadtime is not None:
                while switch is False:
                    filename = self._make_model_filename(
                                    fc_date, leadtime, **kwargs)
                    # check if file is accessible
                    if remoteHostName is not None:
                        filename = build_usr_pw_path(filename,
                                                     remoteHostName,
                                                     usr, pw)
                    switch = check_if_ncfile_accessible(filename, **kwargs)
                    if (switch is False):
                        logger.warning(
                            "Desired file:", filename, " not accessible")
                        logger.warning(
                            "Continue to look for date" 
                            + " with extended leadtime")
                        leadtime = (leadtime
                                    + vars(self.cfg)['misc']['init_step'])
                    if (kwargs.get('max_lt') is not None
                        and leadtime > kwargs.get('max_lt')):
                        logger.warning("Leadtime:", leadtime,
                              "is greater as maximum allowed leadtime:",
                              str(kwargs.get('max_lt')))
                        break
            else:
                filename = None
        elif (isinstance(fc_date, list) and isinstance(leadtime, int)):
            filename = [self._make_model_filename(date, leadtime, **kwargs)
                        for date in fc_date]
        elif (isinstance(fc_date, list) and leadtime == 'best'):
            leadtime = generate_bestguess_leadtime(
                       self.nID, fc_date, **kwargs)
            filename = [self._make_model_filename(
                        fc_date[i], leadtime[i], **kwargs)
                        for i in range(len(fc_date))]
        elif leadtime is None:
            filename = None

        return filename

    def _make_list_of_model_filenames(self, fc_dates, lt, **kwargs):
        """
        return: flst - list of model files to be opened
                dlst - list of dates to be chosen within each file
        """
        flst = []
        for d in fc_dates:
            fn = self._make_model_filename_wrapper(d, lt, **kwargs)
            if fn is not None:
                flst.append(fn)
        return flst


    def _get_files(self, dict_for_sub=None, path=None,
    wavy_path=None, **kwargs):
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
            pathlst - list of paths
            filelst - list of files
        """
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        filelst = []
        pathlst = []
        tmpdate = self.sd
        if wavy_path is not None:
            pathtotals = [wavy_path]
            filelst = [wavy_path]
        elif path is None:
            logger.info('path is None -> checking config file')
            while (tmpdate <= date_dispatcher(self.ed,
            self.cfg.misc['date_incr_unit'], self.cfg.misc['date_incr'])):
                try:
                    # create local path for each time
                    path_template = \
                            vars(self.cfg)['wavy_input'].get('src_tmplt')
                    strsublst = vars(self.cfg)['wavy_input'].get('strsub')
                    subdict = \
                        make_subdict(strsublst,
                                     class_object_dict=dict_for_sub)
                    path = make_pathtofile(path_template,
                                           strsublst, subdict,
                                           **kwargs)
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
            pathtotals = [pathlst]

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
        print(str(int(len(pathtotals))) + " valid files found")
        return pathtotals, filelst

    def list_input_files(self, show=False, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        if (kwargs.get('path') is None and kwargs.get('wavy_path') is None):
            fc_dates = make_fc_dates(self.sd, self.ed,
                                     self.cfg.misc['date_incr_unit'],
                                     self.cfg.misc['date_incr'])

            pathlst = self._make_list_of_model_filenames(
                    fc_dates, self.leadtime, **kwargs)
        else:
            # if defined path local
            logger.info(" ## Find and list files ...")
            path = kwargs.get('path', None)
            wavy_path = kwargs.get('wavy_path', None)
            pathlst, _ = self._get_files(vars(self),
                                         path=path,
                                         wavy_path=wavy_path,
                                         **kwargs)

        if show is True:
            print(" ")
            print("Found files:")
            print(pathlst)
            print(" ")
        return pathlst

    def _get_model(self, **kwargs):
        """
        Main function to obtain data from satellite missions.
        reads files, apply region and temporal filter

        return: adjusted dictionary according to spatial and
                temporal constraints
        """

        # retrieve dataset
        ds = self.reader(**kwargs, **(vars(self)))

        self.vars = ds
        self.coords = list(self.vars.coords)
        return self

    def _enforce_longitude_format(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        new = deepcopy(self)
        # adjust longitude -180/180
        attrs = new.vars.lons.attrs
        logger.info(' enforcing lon max min = -180/180')
        attrs['valid_max'] = 180
        attrs['valid_min'] = -180
        attrs['comments'] = 'forced to range: -180 to 180'
        try:
            new.vars.lons.values = ((new.vars.lons.values + 180) % 360) - 180
        except Exception as e:
            logger.info('Exception in _enforce_longitude_format:')
            logger.info(e)
            new.vars.assign_coords({"lons":
                ((new.vars.lons.values + 180) % 360) - 180})
        return new

    def _enforce_meteorologic_convention(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        logger.info(' enforcing meteorological convention')
        for v in self.varalias:
            if ('convention' in vars(self.cfg)['misc'].keys() and
            vars(self.cfg)['misc']['convention'] == 'oceanographic'):
                logger.info(
                    'Convert from oceanographic to meteorologic convention')

                self.vars[v] = convert_meteorologic_oceanographic(self.vars[v])

            elif 'to_direction' in self.vars[v].attrs['standard_name']:
                logger.info(
                    'Convert from oceanographic to meteorologic convention')
                self.vars[v] = convert_meteorologic_oceanographic(self.vars[v])
        return self

    def _change_varname_to_aliases(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        logger.info(' changing variables to aliases')
        # variables
        for v in self.varalias:
            ncvar = get_filevarname(v, variable_def,
                                    vars(self.cfg), self.meta,
                                    **kwargs)
            if v in list(self.vars.keys()):
                logger.info('  ', ncvar, 'is alreade named correctly and'
                    + ' therefore not adjusted')
            else:
                self.vars = self.vars.rename({ncvar: v})
        # coords
        coords = ['time', 'lons', 'lats']
        for c in coords:
            ncvar = get_filevarname(c, variable_def,
                                    vars(self.cfg), self.meta,
                                    **kwargs)
            if c in list(self.vars.keys()):
                logger.info('  ', c, 'is alreade named correctly and'
                    + ' therefore not adjusted')
            else:
                self.vars = self.vars.rename({ncvar: c})\
                                        .set_index(time='time')
        return self

    def _change_stdvarname_to_cfname(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        logger.info(' complying to cf standard names')
        # enforce standard_name for coordinate aliases
        self.vars['lons'].attrs['standard_name'] = \
            variable_def['lons'].get('standard_name')
        self.vars['lats'].attrs['standard_name'] = \
            variable_def['lats'].get('standard_name')
        self.vars['time'].attrs['standard_name'] = \
            variable_def['time'].get('standard_name')
        # enforce standard_name for variable alias
        for i in range(len(self.varalias)):
            self.vars[self.varalias[i]].attrs['standard_name'] = \
                self.stdvarname[i]
        return self

    def populate(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        logger.info(" ### Read files and populate model_class object")

        fc_dates = make_fc_dates(self.sd, self.ed,
                                 self.cfg.misc['date_incr_unit'],
                                 self.cfg.misc['date_incr'])

        self.pathlst = self.list_input_files(**kwargs)

        if len(self.pathlst) > 0:
            logger.info('')
            logger.info('Checking variables..')
            self.meta = ncdumpMeta(self.pathlst[0])
            ncvar = [get_filevarname(v, variable_def,
                                    vars(self.cfg), 
                                    self.meta, **kwargs) for \
                                     v in self.varalias]
            logger.info('')
            logger.info('Choosing reader..')
            # define reader
            dotenv.load_dotenv()
            WAVY_DIR = os.getenv('WAVY_DIR', None)
            if WAVY_DIR is None:
                logger.debug('#')
                logger.debug('Environmental variable for WAVY_DIR not defined')
                logger.debug('Defaults are chosen')
                logger.debug('#')
                reader_mod_str = load_dir('model_readers').name
            else:
                reader_mod_str = WAVY_DIR + '/wavy/model_readers.py'

            reader_str = kwargs.get('reader', self.cfg.reader)
            spec = importlib.util.spec_from_file_location(
                    'model_readers.' + reader_str, reader_mod_str)

            # create reader module
            reader_tmp = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(reader_tmp)

            # pick reader
            reader = getattr(reader_tmp, reader_str)
            self.reader = reader
            logger.info('Chosen reader:', spec.name)
            logger.info('')

            # possible to select list of variables
            self.varname = ncvar

            # if credentials are needed
            remoteHostName = kwargs.get('remoteHostName',
                                        self.cfg.misc.get('remoteHostName'))

            kwargs['fc_dates'] = fc_dates

            try:
                t0 = time.time()
                logger.debug('Reading..')
                self = self._get_model(remoteHostName=remoteHostName,
                                       **kwargs)

                self = self._change_varname_to_aliases(**kwargs)
                self = self._change_stdvarname_to_cfname(**kwargs)
                self = self._enforce_meteorologic_convention(**kwargs)

                # convert longitude
                self = self._enforce_longitude_format(**kwargs)

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
                logger.info(" ")
                logger.info(' ## Summary:')
                logger.info(str(len(self.vars['time']))
                            + " time steps retrieved.")
                logger.info("Time used for retrieving data:")
                logger.info(round(t1-t0, 2), "seconds")
                logger.info(" ")
                logger.info(" ### model_class object populated ###")
                logger.info('# ----- ')
            except Exception as e:
                logger.exception(e)
                logger.error(e)
                logger.error('Error encountered')
                logger.error('model_class object not populated')
        else:
            logger.warning('No data data found')
            logger.warning('model_class object not populated')
            logger.warning('# ----- ')
        return self
