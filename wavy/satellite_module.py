#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to acquire, read, and prepare
geophysical variables related to wave height from satellite
altimetry files for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library igports
import sys
import numpy as np
from datetime import datetime, timedelta
import os
from copy import deepcopy
import time
from dateutil.relativedelta import relativedelta
import importlib.util
import pyproj
import dotenv
import glob
import xarray as xr
from tqdm import tqdm
#import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import get_filevarname
from wavy.ncmod import find_attr_in_nc

from wavy.utils import NoStdStreams
from wavy.utils import find_included_times
from wavy.utils import parse_date
from wavy.utils import make_pathtofile, make_subdict
from wavy.utils import finditem, haversineA
from wavy.utils import flatten
from wavy.utils import date_dispatcher
from wavy.utils import convert_meteorologic_oceanographic

from wavy.model_module import read_model_nc_output_lru
from wavy.model_module import model_class as mc

from wavy.insitu_module import poi_class as pc

from wavy.wconfig import load_or_default, load_dir

from wavy.filtermod import filter_class as fc

from wavy.quicklookmod import quicklook_class_sat as qls

from wavy.init_class_sat import init_class

from wavy.utils import footprint_pulse_limited_radius
# ---------------------------------------------------------------------#

# read yaml config files:
region_dict = load_or_default('region_cfg.yaml')
model_dict = load_or_default('model_cfg.yaml')
#satellite_dict = load_or_default('satellite_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')

# --- global functions ------------------------------------------------#

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

# ---------------------------------------------------------------------#


class satellite_class(qls, fc):
    '''
    Class to handle netcdf files containing satellite altimeter data
    e.g.: Hs[time], lat[time], lon[time], time[time]
    '''

    def __init__(self, **kwargs):
        print('# ----- ')
        print(" ### Initializing satellite_class object ###")
        print(" ")
        print(" Given kwargs:")
        print(kwargs)
        # initializing useful attributes from config file
        dc = init_class('satellite', kwargs.get('nID'))
        # parse and translate date input
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        # add other class object variables
        self.nID = kwargs.get('nID')
        self.name = dc.name[kwargs.get('name', 's3a')]
        self.varalias = kwargs.get('varalias', 'Hs')
        self.units = variable_def[self.varalias].get('units')
        self.stdvarname = variable_def[self.varalias].get('standard_name')
        self.twin = int(kwargs.get('twin', 30))
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        if isinstance(kwargs.get('region'), dict):
            self.region = kwargs.get('region')['name']
        else:
            self.region = kwargs.get('region', 'global')
        self.cfg = dc
        # poi should be poi_class
        self.poi = kwargs.get('poi', None)
        if self.poi is not None:
            self.sd = parse_date(str(self.poi.vars['time'].data[0]))
            self.ed = parse_date(str(self.poi.vars['time'].data[-1]))

        # super(config_class,self).__init__('satellite', kwargs.get('nID'))
        print(" ")
        print(" ### satellite_class object initialized ###")
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
            collector_mod_str = load_dir('satellite_collectors').name
        else:
            collector_mod_str = WAVY_DIR + '/wavy/satellite_collectors.py'

        collector_str = kwargs.get('collector', self.cfg.collector)
        spec = importlib.util.spec_from_file_location(
                'satellite_collectors.' + collector_str, collector_mod_str)

        # create collector module
        collector_tmp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(collector_tmp)

        # pick collector
        collector = getattr(collector_tmp, collector_str)
        self.collector = collector
        print('Chosen collector:', spec.name)
        print('')

        print("Downloading files ...")
        #kwargs_in = vars(self)
        #kwargs_in['path'] = path
        self.collector(nproc=nproc, path=path, **vars(self))


    def _get_files(self, dict_for_sub=None, path=None, wavy_path=None):
        """
        Function to retrieve list of files/paths for available
        locally stored satellite data. This list is used for
        other functions to query and parsing.

        param:
            sd - start date (datetime object)
            ed - end date (datetime object)
            twin - time window (temporal constraint) in minutes
            nID - nID as of satellite_cfg.yaml
            dict_for_sub - dictionary for substitution in templates
            path - a path if defined

        return:
            pathlst - list of paths
        """
        filelst = []
        pathlst = []
        tmpdate = self.sd-timedelta(minutes=self.twin)
        if wavy_path is not None:
            pathtotals = [wavy_path]
            filelst = [wavy_path]
        elif path is None:
            print('path is None -> checking config file')
            while (tmpdate <= date_dispatcher(self.ed,
            self.cfg.wavy_input['path_date_incr_unit'],
            self.cfg.wavy_input['path_date_incr'])):
                try:
                    # create local path for each time
                    path_template = \
                            self.cfg.wavy_input.get('src_tmplt')
                    strsublst = \
                        self.cfg.wavy_input.get('strsub')
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
                            self.cfg.wavy_input['path_date_incr_unit'],
                            self.cfg.wavy_input['path_date_incr'])

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
            idx_start, tmp = check_date(filelst,
                                        self.sd - timedelta(minutes=self.twin))
            tmp, idx_end = check_date(filelst,
                                      self.ed + timedelta(minutes=self.twin))
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
            idx_start, tmp = check_date(filelst,
                                        self.sd - timedelta(minutes=self.twin))
            tmp, idx_end = check_date(filelst,
                                      self.ed + timedelta(minutes=self.twin))
            if idx_end == 0:
                idx_end = len(pathlst)-1
            del tmp
            pathtotals = np.unique(pathtotals[idx_start:idx_end+1])
            filelst = np.unique(filelst[idx_start:idx_end+1])

        return pathtotals

    def list_input_files(self, show=False, **kwargs):
        print(" ## Find and list files ...")
        path = kwargs.get('path', None)
        wavy_path = kwargs.get('wavy_path', None)
        pathlst = self._get_files(vars(self),
                                     path=path,
                                     wavy_path=wavy_path)

        # remove None values from pathlst
        pathlst = list(filter(lambda item: item is not None, pathlst))
        print(str(int(len(pathlst))) + " valid files found")

        print('source template:',
              #satellite_dict[self.nID]['wavy_input']['src_tmplt'])
              self.cfg.wavy_input['src_tmplt'])
        if show is True:
            print(" ")
            print(pathlst)
            print(" ")
        return pathlst

    def crop_to_poi(self, **kwargs):
        if kwargs.get('poi') is not None:
            self.poi = kwargs.get('poi')
        if kwargs.get('distlim') is not None:
            self.distlim = kwargs.get('distlim')
        if kwargs.get('twin') is not None:
            self.twin = kwargs.get('twin')
        new = deepcopy(self)
        print('Crop to poi region')
        idx = new._match_poi(self.poi)
        new.vars = new.vars.sel(time=new.vars.time[idx])
        print('Region mask applied based on poi')
        print('For chosen poi region: ', len(new.vars['time']),
              'footprints found')
        return new

    def crop_to_region(self, region):
        new = deepcopy(self)
        print('Crop to region:', region)

        idx = new._match_region(new.vars['lats'].values,
                                new.vars['lons'].values,
                                region=region,
                                grid_date=new.sd)
        new.vars = new.vars.isel(time=idx)
        print('Region mask applied')
        print('For chosen region: ', len(new.vars['time']),
              'footprints found')
        return new

    def _get_sat_ts_blunt(self, **kwargs):
        """
        Main function to obtain data from satellite missions.
        reads files, apply region and temporal filter

        fct is insensitive to region

        return: adjusted dictionary according to spatial and
                temporal constraints
        """

        # retrieve dataset
        ds = self.reader(**(vars(self)))
        self.vars = ds
        self.coords = list(self.vars.coords)
        return self

    def _get_sat_ts(self, **kwargs):
        """
        Main function to obtain data from satellite missions.
        reads files, apply region and temporal filter

        fct will crop to region on the fly if given in init kwargs
        chunk_size determines number of files read before region filter

        return: adjusted dictionary according to spatial and
                temporal constraints
        """
        new = deepcopy(self)

        pathlst = self.pathlst
        chunk_size = kwargs.get('chunk_size', 1)
        region = kwargs.get('region', new.region)

        if isinstance(region, dict):
            new.region = region['name']

        ds_lst = []
        count = 0
        print('Reading', int((len(pathlst)+chunk_size)/chunk_size)+1,
              'chunks of files with chunk size', chunk_size)
        print('Total of', len(pathlst), 'files')

        for count in tqdm(range(0, len(pathlst)+chunk_size, chunk_size)):
            if count <= len(pathlst)-1:
                new.pathlst = pathlst[count:count+chunk_size]
                with NoStdStreams():
                    try:
                        # retrieve dataset
                        ds = new.reader(**(vars(new)))
                        new.vars = ds
                        new.coords = new.vars.coords
                        if (new.poi is None and
                        isinstance(region, str)):
                            ds_lst.append(new._change_varname_to_aliases()
                                          ._enforce_longitude_format()
                                          .crop_to_region(new.region).vars)
                        elif (new.poi is None and
                        isinstance(region, dict)):
                            ds_lst.append(new._change_varname_to_aliases()
                                          ._enforce_longitude_format()
                                          .crop_to_region(
                                              region['region']).vars)
                        else:
                            ds_lst.append(new._change_varname_to_aliases()
                                          ._enforce_longitude_format()
                                          .crop_to_region(region)
                                          .crop_to_poi().vars)
                    except Exception as e:
                        logger.exception(e)

        if len(ds_lst) > 1:
            combined = xr.concat(ds_lst, 'time',
                                 coords='minimal',
                                 data_vars='minimal',
                                 compat='override',
                                 combine_attrs='override',
                                 join='override')
        else:
            combined = ds_lst[0]

        new.vars = combined
        new.coords = list(new.vars.coords)

        return new

    def _enforce_longitude_format(self):
        new = deepcopy(self)
        # adjust longitude -180/180
        attrs = new.vars.lons.attrs
        print(' enforcing lon max min = -180/180')
        attrs['valid_max'] = 180
        attrs['valid_min'] = -180
        attrs['comments'] = 'forced to range: -180 to 180'
        new.vars.lons.values = ((new.vars.lons.values + 180) % 360) - 180
        return new

    def _enforce_meteorologic_convention(self):
        new = deepcopy(self)
        ncvars = list(new.vars.variables)
        #for ncvar in ncvars:
        #    if ('convention' in satellite_dict[new.nID].keys() and
        #    satellite_dict[new.nID]['convention'] == 'oceanographic'):

        for ncvar in ncvars:
            if ('convention' in vars(new.cfg)['misc'].keys() and
            vars(new.cfg)['misc']['convention'] == 'oceanographic'):
                print('Convert from oceanographic to meteorologic convention')
                new.vars[ncvar] =\
                    convert_meteorologic_oceanographic(new.vars[ncvar])
            elif 'to_direction' in new.vars[ncvar].attrs['standard_name']:
                print('Convert from oceanographic to meteorologic convention')
                new.vars[ncvar] =\
                    convert_meteorologic_oceanographic(new.vars[ncvar])

        return new


    def _change_varname_to_aliases(self):
        satellite_dict = load_or_default('satellite_cfg.yaml')
        print(' changing variables to aliases')
        new = deepcopy(self)
        # variables
        ncvar = get_filevarname(new.varalias, variable_def,
                                satellite_dict[new.nID], new.meta)
        if new.varalias in list(new.vars.keys()):
            print('  ', ncvar, 'is alreade named correctly and'
                + ' therefore not adjusted')
        else:
            new.vars = new.vars.rename({ncvar: new.varalias})
        # coords
        coords = ['time', 'lons', 'lats']
        for c in coords:
            try:
                ncvar = get_filevarname(c, variable_def,
                                        satellite_dict[new.nID], new.meta)
                if c in list(new.vars.keys()):
                    print('  ', c, 'is alreade named correctly and'
                        + ' therefore not adjusted')
                else:
                    new.vars = new.vars.rename({ncvar: c})#\
                                            #.set_index(time='time')
            except Exception as e:
                print(' ', ncvar, 'is not renamed')
                #logger.exception(e)
                #logger.debug(traceback.format_exc())
                print(e)

        return new


    def _change_stdvarname_to_cfname(self):
        # enforce standard_name for coordinate aliases
        new = deepcopy(self)
        new.vars['lons'].attrs['standard_name'] = \
            variable_def['lons'].get('standard_name')
        new.vars['lats'].attrs['standard_name'] = \
            variable_def['lats'].get('standard_name')
        new.vars['time'].attrs['standard_name'] = \
            variable_def['time'].get('standard_name')
        # enforce standard_name for variable alias
        new.vars[self.varalias].attrs['standard_name'] = \
            new.stdvarname
        return new

    def populate(self, **kwargs):
        print(" ### Read files and populate satellite_class object")

        satellite_dict = load_or_default('satellite_cfg.yaml')

        lst = self.list_input_files(**kwargs)
        self.pathlst = kwargs.get('pathlst', lst)

        print('')
        print('Checking variables..')
        self.meta = ncdumpMeta(self.pathlst[0])
        ncvar = get_filevarname(self.varalias, variable_def,
                                satellite_dict[self.nID], self.meta)
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
            reader_mod_str = load_dir('satellite_readers').name
        else:
            reader_mod_str = WAVY_DIR + '/wavy/satellite_readers.py'

        reader_str = kwargs.get('reader', self.cfg.reader)
        spec = importlib.util.spec_from_file_location(
                'satellite_readers.' + reader_str, reader_mod_str)

        # create reader module
        reader_tmp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reader_tmp)

        # pick reader
        reader = getattr(reader_tmp, reader_str)
        self.reader = reader
        print('Chosen reader:', spec.name)
        print('')

        # possible to select list of variables
        self.varname = ncvar

        if len(lst) > 0:
            try:
                t0 = time.time()
                print('Reading..')
                self = self._get_sat_ts(**kwargs)
                #self = self._get_sat_ts_blunt(**kwargs)

                self = self._change_varname_to_aliases()
                self = self._change_stdvarname_to_cfname()
                self = self._enforce_meteorologic_convention()

                # convert longitude
                self = self._enforce_longitude_format()

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
                print(str(len(self.vars['time'])) + " footprints retrieved.")
                print("Time used for retrieving data:")
                print(round(t1-t0, 2), "seconds")
                print(" ")
                print(" ### satellite_class object populated ###")
                print('# ----- ')

            except Exception as e:
                logger.exception(e)
                #logger.debug(traceback.format_exc())
                print(e)
                print('Error encountered')
                print('satellite_class object not populated')
        else:
            print('No data data found')
            print('satellite_class object not populated')
            print('# ----- ')
        return self

    def drop_duplicates(self, dim='time', keep='first'):
        print('Removing duplicates according to', dim)
        print('Keeping', keep, 'value for the duplicates')
        new = deepcopy(self)
        new.vars = self.vars.drop_duplicates(dim=dim, keep=keep)
        print(str(int(abs(len(self.vars[dim])-len(new.vars[dim])))),
              'values removed')
        print('New number of footprints is:', str(int(len(new.vars[dim]))))
        return new

    def _match_poi(self, poi):
        """
        return: idx that match to region
        """
        from tqdm import tqdm
        print('Match up poi locations')
        # increase with 1 degree in all directions to ensure
        # that distlim is still the determining factor
        region = {'llcrnrlat': poi.vars['lats'].min().data-1,
                  'urcrnrlat': poi.vars['lats'].max().data+1,
                  'llcrnrlon': poi.vars['lons'].min().data-1,
                  'urcrnrlon': poi.vars['lons'].max().data+1}

        ridx = match_region_rect(self.vars['lats'].data,
                                 self.vars['lons'].data,
                                 region=region)

        newdict = deepcopy(self.vars)
        idx = [self._poi_sat(newdict, self.twin, self.distlim, poi, ridx, i)
               for i in tqdm(range(len(poi.vars['time'])))]
        idx = list(np.array(ridx)[flatten(idx)])
        return idx

    @staticmethod
    def _poi_sat(ds, twin, distlim, poi, ridx, i):
        """
        return: indices for values matching the spatial and
                temporal constraints
        """
        target_t = parse_date(str(poi.vars['time'].data[i]))
        unfiltered_t = [parse_date(str(d)) for d in ds['time'][ridx].values]
        tidx = find_included_times(
                    unfiltered_t,
                    target_t=target_t,
                    twin=twin)
        slons = list(np.array(ds['lons'])[ridx][tidx])
        slats = list(np.array(ds['lats'])[ridx][tidx])
        plons = [poi.vars['lons'].data[i]]*len(slons)
        plats = [poi.vars['lats'].data[i]]*len(slats)
        dists_tmp = haversineA(slons, slats, plons, plats)
        sidx = np.argwhere(np.array(dists_tmp) <= distlim).flatten()
        #dists = list(np.array(dists_tmp)[sidx])
        return list(np.array(tidx)[sidx])

    @staticmethod
    def _match_region(LATS, LONS, region, grid_date):
        """
        Function to filter satellite data according to region      
        return:
            indices that match the region
        """
        # region in region_dict[poly]:
        # find values for given region
        if isinstance(region, dict) is True:
            print('Region is defined in custom dict')
            print(region)
            region = region
            ridx = match_region_rect(LATS, LONS, region=region)
        elif (region not in region_dict['poly'] and
        region not in model_dict and
        region not in region_dict['geojson']):
            if region is None:
                region = 'global'
            else:
                if (region not in region_dict['rect']
                and region not in region_dict['geojson']
                and isinstance(region, dict) is False):
                    sys.exit("Region is not defined")
                else:
                    print("Specified region: " + region + "\n"
                          + " --> Bounds: " +
                          str(region_dict['rect'][region]))
                    region = region_dict['rect'][region]
            ridx = match_region_rect(LATS, LONS, region=region)
        elif region in region_dict['geojson']:
            print("Region is defined as geojson")
            ridx = match_region_geojson(LATS, LONS, region=region)
        elif region in region_dict['poly']:
            ridx = match_region_poly(LATS, LONS, region=region,
                                     grid_date=grid_date)
        else:
            ridx = match_region_poly(LATS, LONS, region=region,
                                     grid_date=grid_date)
        return ridx

    def compute_pulse_limited_footprint_radius(self):
        """
        Compute pulse limited footprint size
        """
        new = deepcopy(self)
        Hs = new.vars['Hs'].values
        tau = new.cfg.misc['sat_specs'][new.name]['tau']*10**(-9)
        h = new.cfg.misc['sat_specs'][new.name]['h']*10**3
        fpr = footprint_pulse_limited_radius(Hs, h, tau)
        ds = new.vars
        ds = ds.assign({"fpr": (("time"), fpr)})
        ds["fpr"].attrs = variable_def["fpr"]
        new.vars = ds
        return new

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


def match_region_rect(LATS, LONS, region):
    """
    Takes care of regular grid regions
    """
    if (region is None or region == "global"):
        region = "global"
        ridx = range(len(LATS))
    else:
        llcrnrlat = region["llcrnrlat"]
        urcrnrlat = region["urcrnrlat"]
        llcrnrlon = region["llcrnrlon"]
        urcrnrlon = region["urcrnrlon"]
        ridx = np.where(
                    (np.array(LATS) > llcrnrlat)
                    & (np.array(LATS) < urcrnrlat)
                    & (np.array(LONS) > llcrnrlon)
                    & (np.array(LONS) < urcrnrlon)
                    )[0]
    print(len(ridx), " values found for chosen region and time frame.")
    return ridx


def match_region_geojson(LATS, LONS, region):
    """
    Takes care of regions defines as geojson
    """
    import geojson
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    points = np.c_[LONS, LATS]
    ridx = []
    fstr = region_dict['geojson'][region]['fstr']
    with open(fstr) as f:
        gj = geojson.load(f)
    fidx = region_dict['geojson'][region].get('fidx')
    if fidx is not None:
        geo = {'type': 'Polygon',
           'coordinates':\
                gj['features'][fidx]['geometry']['coordinates'][0]}
        poly = Polygon([tuple(l)
                        for l in geo['coordinates'][0]],
                       closed=True)
        hits = Path(poly.xy).contains_points(points, radius=1e-9)
        ridx_tmp = list(np.array(range(len(LONS)))[hits])
        ridx += ridx_tmp
    else:
        for i in range(len(gj['features'])):
            geo = {'type': 'Polygon',
               'coordinates':\
                    gj['features'][i]['geometry']['coordinates'][0]}
            poly = Polygon([tuple(l)
                            for l in geo['coordinates'][0]],
                           closed=True)
            hits = Path(poly.xy).contains_points(points, radius=1e-9)
            ridx_tmp = list(np.array(range(len(LONS)))[hits])
            ridx += ridx_tmp
    return ridx


def match_region_poly(LATS, LONS, region, grid_date):
    """
    Takes care of region defined as polygon
    """
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    import numpy as np
    if (region not in region_dict['poly'] and region not in model_dict):
        sys.exit("Region polygone is not defined")
    elif isinstance(region, dict) is True:
        print("Manuall specified region: \n"
              + " --> Bounds: " + str(region))
        poly = Polygon(list(zip(region['lons'],
                       region['lats'])), closed=True)
    elif (isinstance(region, str) is True and region in model_dict):
        # init model_class object
        mco = mc(nID=region)
        try:
            print('Use date for retrieving grid: ', grid_date)
            filestr = mco._make_model_filename_wrapper(grid_date, 'best')
            meta = ncdumpMeta(filestr)
            flon = get_filevarname('lons', variable_def,
                                   model_dict[region], meta)
            flat = get_filevarname('lats', variable_def,
                                   model_dict[region], meta)
            time = get_filevarname('time', variable_def,
                                   model_dict[region], meta)
            model_lons, model_lats, _ = \
                read_model_nc_output_lru(filestr, flon, flat, time)
        except (KeyError, IOError, ValueError) as e:
            print(e)
            if 'grid_date' in model_dict[region]:
                grid_date = model_dict[region]['grid_date']
                print('Trying default date ', grid_date)
            else:
                grid_date = datetime(
                                    datetime.now().year,
                                    datetime.now().month,
                                    datetime.now().day
                                    )
            filestr = mco._make_model_filename_wrapper(grid_date, 'best')
            meta = ncdumpMeta(filestr)
            flon = get_filevarname('lons', variable_def,
                                   model_dict[region], meta)
            flat = get_filevarname('lats', variable_def,
                                   model_dict[region], meta)
            time = get_filevarname('time', variable_def,
                                   model_dict[region], meta)
            model_lons, model_lats, _ = \
                read_model_nc_output_lru(filestr, flon, flat, time)
        if (len(model_lons.shape) == 1):
            model_lons, model_lats = np.meshgrid(
                                    model_lons,
                                    model_lats
                                    )
        print('Check if footprints fall within the chosen domain')
        ncdict = ncdumpMeta(filestr)
        try:
            proj4 = find_attr_in_nc('proj', ncdict=ncdict,
                                    subattrstr='proj4')
        except IndexError:
            print('proj4 not defined in netcdf-file')
            print('Using proj4 from model config file')
            proj4 = model_dict[region]['misc']['proj4']
        proj_model = pyproj.Proj(proj4)
        Mx, My = proj_model(model_lons, model_lats, inverse=False)
        Vx, Vy = proj_model(LONS, LATS, inverse=False)
        xmax, xmin = np.max(Mx), np.min(Mx)
        ymax, ymin = np.max(My), np.min(My)
        ridx = list(np.where((Vx > xmin) & (Vx < xmax) &
                             (Vy > ymin) & (Vy < ymax))[0])
    elif isinstance(region, str) is True:
        print("Specified region: " + region + "\n"
              + " --> Bounded by polygon: \n"
              + "lons: " + str(region_dict['poly'][region]['lons'])
              + "\n"
              + "lats: " + str(region_dict['poly'][region]['lats']))
        poly = Polygon(list(zip(region_dict['poly'][region]['lons'],
                       region_dict['poly'][region]['lats'])),
                       closed=True)
        # check if coords in region
        points = np.c_[LONS, LATS]
        # radius seems to be important to correctly define polygone
        # see discussion here:
        # https://github.com/matplotlib/matplotlib/issues/9704
        hits = Path(poly.xy).contains_points(points, radius=1e-9)
        ridx = list(np.array(range(len(LONS)))[hits])
    if (not ridx or len(ridx)<1):
        print("No values for chosen region and time frame!!!")
    return ridx
