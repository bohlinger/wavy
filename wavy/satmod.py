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
import netCDF4 as netCDF4
from dateutil.relativedelta import relativedelta
import pyproj
import zipfile
import tempfile
import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import get_filevarname
from wavy.ncmod import find_attr_in_nc, dumptonc_ts_sat
from wavy.utils import find_included_times, NoStdStreams
from wavy.utils import parse_date
from wavy.utils import make_pathtofile, make_subdict
from wavy.utils import finditem, haversineA
from wavy.utils import flatten
from wavy.utils import convert_meteorologic_oceanographic
from wavy.modelmod import make_model_filename_wrapper
from wavy.modelmod import read_model_nc_output_lru
from wavy.wconfig import load_or_default
from wavy.filtermod import filter_main, vardict_unique
from wavy.filtermod import rm_nan_from_vardict
from wavy.sat_collectors import get_remote_files
from wavy.sat_readers import read_local_files
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.writermod import writer_class as wc
from wavy.init_class import init_class
# ---------------------------------------------------------------------#

# read yaml config files:
region_dict = load_or_default('region_cfg.yaml')
model_dict = load_or_default('model_cfg.yaml')
satellite_dict = load_or_default('satellite_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')

# --- global functions ------------------------------------------------#

def crop_vardict_to_period(vardict,sdate,edate):
    """
    Function to crop the variable dictionary to a given period
    """
    for key in vardict:
        if (key != 'time_unit' and key != 'meta' and key != 'datetime'):
            vardict[key] =  list(np.array(vardict[key])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ])
        else:
            vardict[key] = vardict[key]
    vardict['datetime'] =   list(np.array(vardict['datetime'])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ])
    return vardict


def check_date(filelst,date):
    '''
    Checks if str in lst according to desired date (sdate,edate)

    return: idx for file
    '''
    idx = []
    for i in range(len(filelst)):
        element = filelst[i]
        tmp = element.find(date.strftime('%Y%m%d'))
        if tmp>=0:
            idx.append(i)
    if len(idx)<=0:
        idx=[0]
    return idx[0],idx[-1]

# ---------------------------------------------------------------------#


class satellite_class(qls, wc):
    '''
    Class to handle netcdf files containing satellite data e.g.
    Hs[time], lat[time], lon[time]
    '''

    def __init__(self, **kwargs):
        print('# ----- ')
        print(" ### Initializing satellite_class object ###")
        print(" ")
        dc = init_class('satellite', kwargs.get('nID'))
        # parse and translate date input
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        self.nID = kwargs.get('nID')
        self.mission = kwargs.get('mission', 's3a')
        self.varalias = kwargs.get('varalias', 'Hs')
        self.units = variable_info[self.varalias].get('units')
        self.stdvarname = variable_info[self.varalias].get('standard_name')
        self.twin = kwargs.get('twin', 0)
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')
        self.cfg = dc
        # self.poi = kwargs.get('poi',None)
        print(" ")
        print(" ### satellite_class object initialized ###")
        print('# ----- ')

    def download(self, path=None, nproc=1, **kwargs):
        print("Downloading files ...")
        get_remote_files(
                        path_local=path,
                        nproc=nproc,
                        twin=kwargs.get('twin', self.twin),
                        sd=kwargs.get('sd', self.sd),
                        ed=kwargs.get('ed', self.ed),
                        nID=self.nID,
                        mission=self.mission,
                        dict_for_sub=vars(self)
                        )

    def _get_files(self, dict_for_sub=None, path=None):
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
            local_path - a path if defined

        return:
            pathlst - list of paths
            filelst - list of files
        """
        filelst = []
        pathlst = []
        tmpdate = self.sd-timedelta(minutes=self.twin)
        if path is None:
            print('path is None -> checking config file')
            while (tmpdate <= self.ed + relativedelta(months=+1)):
                try:
                    # create local path for each time
                    path_template = \
                            satellite_dict[self.nID]['wavy_input'].get(
                                                  'src_tmplt')
                    strsublst = \
                        satellite_dict[self.nID]['wavy_input'].get('strsub')
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
                tmpdate = tmpdate + relativedelta(months=+1)
            filelst = np.sort(flatten(filelst))
            pathlst = np.sort(flatten(pathlst))
        else:
            filelst = np.sort(os.listdir(path))
            pathlst = [os.path.join(path, e) for e in filelst]
        idx_start, tmp = check_date(filelst,
                                    self.sd - timedelta(minutes=self.twin))
        tmp, idx_end = check_date(filelst,
                                  self.ed + timedelta(minutes=self.twin))
        if idx_end == 0:
            idx_end = len(pathlst)-1
        del tmp
        pathlst = np.unique(pathlst[idx_start:idx_end+1])
        filelst = np.unique(filelst[idx_start:idx_end+1])
        print(str(int(len(pathlst))) + " valid files found")
        return pathlst, filelst

    def list_files(self, **kwargs):
        print(" ## Find and list files ...")
        path = kwargs.get('path', None)
        pathlst, _ = self._get_files(vars(self), path=path)
        print('source template:',
              satellite_dict[self.nID]['wavy_input']['src_tmplt'])
        print(" ")
        print(pathlst)
        print(" ")
        return pathlst

    def _get_sat_ts(self, **kwargs):
        """
        Main function to obtain data from satellite missions.
        reads files, apply region and temporal filter

        return: adjusted dictionary according to spatial and
                temporal contarinst
        """
        cvardict = read_local_files(
                                    pathlst=self.pathlst,
                                    nID=self.nID,
                                    varalias=self.varalias,
                                    sd=self.sd,
                                    ed=self.ed,
                                    twin=self.twin,
                                    **kwargs
                                    )
        print('Total: ', len(cvardict['time']), ' footprints found')
        print('Apply region mask for region', self.region)
        ridx = self._match_region(cvardict['latitude'],
                                  cvardict['longitude'],
                                  region=self.region,
                                  grid_date=self.sd)
        print('Region mask applied')
        rvardict = {}
        for element in cvardict:
            if element != 'time_unit':
                rvardict[element] = list(np.array(
                                        cvardict[element]
                                        )[ridx])
            else:
                rvardict[element] = cvardict[element]
        del cvardict, ridx
        if len(rvardict['time']) > 0:
            rvardict['datetime'] = netCDF4.num2date(
                                        rvardict['time'],
                                        rvardict['time_unit'])
            print('For chosen region and time: ',
                    len(rvardict['time']), 'footprints found')
            # convert to datetime object
            timedt = rvardict['datetime']
            rvardict['datetime'] = [datetime(t.year, t.month, t.day,
                                             t.hour, t.minute, t.second,
                                             t.microsecond)
                                    for t in timedt]
        else:
            print('For chosen region and time: 0 footprints found!')
        if self.poi is not None:
            pvardict = {}
            #pidx = match_poi(rvardict, self.twin, self.distlim, self.poi)
            pidx = self._match_poi(rvardict, self.twin, self.distlim, self.poi)
            for element in rvardict:
                if element != 'time_unit':
                    pvardict[element] = list(np.array(
                                            rvardict[element]
                                            )[pidx])
                else:
                    pvardict[element] = rvardict[element]
            rvardict = pvardict
            print('For chosen poi: ', len(rvardict['time']),
                  'footprints found')
        # find variable name as defined in file
        ncdict = ncdumpMeta(self.pathlst[0])
        rvardict['meta'] = ncdict
        # adjust conventions
        if ('convention' in satellite_dict[self.nID].keys() and
        satellite_dict[self.nID]['convention'] == 'oceanographic'):
            print('Convert from oceanographic to meteorologic convention')
            rvardict[variable_info[self.varalias]['standard_name']] = \
                list(convert_meteorologic_oceanographic(np.array(
                    rvardict[variable_info[self.varalias]['standard_name']]
                    )))
        return rvardict

    def populate(self, **kwargs):
        print(" ### Read files and populate satellite_class object")
        lst = self.list_files(**kwargs)
        self.pathlst = lst
        self.poi = kwargs.get('poi', None)
        if len(lst) > 0:
            try:
                t0 = time.time()
                if self.filter is True:
                    print('Filter is activated,')
                    print('Make sure the chosen time window is appropriate')
                    # retrieve data
                    rvardict = self._get_sat_ts(**kwargs)
                    # adjust varalias if other return_var
                    if kwargs.get('return_var') is not None:
                        newvaralias = kwargs.get('return_var')
                    else:
                        newvaralias = self.varalias
                    # filter data
                    rvardict = filter_main(rvardict,
                                           varalias=newvaralias,
                                           **kwargs)
                else:
                    rvardict = self._get_sat_ts(**kwargs)
                    # adjust varalias if other return_var
                    if kwargs.get('return_var') is not None:
                        newvaralias = kwargs.get('return_var')
                    else:
                        newvaralias = self.varalias
                    # make ts in vardict unique
                    rvardict = vardict_unique(rvardict)
                    # rm NaNs
                    rvardict = rm_nan_from_vardict(newvaralias, rvardict)
                # find variable name as defined in file
                ncdict = ncdumpMeta(self.pathlst[0])
                with NoStdStreams():
                    filevarname = get_filevarname(
                                                  self.varalias,
                                                  variable_info,
                                                  satellite_dict[self.nID],
                                                  ncdict
                                                  )
                rvardict['meta'] = ncdict
                # define more class object variables
                self.vars = rvardict
                self.varname = filevarname
                if kwargs.get('return_var') is not None:
                    self.varalias = kwargs.get('return_var')
                    self.stdvarname = \
                        variable_info[newvaralias].get('standard_name')
                    self.units = variable_info[newvaralias].get('units')
                # create label for plotting
                self.label = self.mission
                t1 = time.time()
                print(" ")
                print(' ## Summary:')
                print(str(len(self.vars['time'])) + " footprints retrieved.")
                print("Time used for retrieving satellite data:")
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
            print('No satellite data found')
            print('satellite_class object not populated')
            print('# ----- ')
        return self

    def _match_poi(self, poi):
        """
        return: idx that match to region
        """
        from tqdm import tqdm
        print('Match up poi locations')
        # increase with 1 degree in all directions to ensure
        # that distlim is still teh determining factor
        region = {'llcrnrlat': np.min(poi['latitude']-1),
                  'urcrnrlat': np.max(poi['latitude']+1),
                  'llcrnrlon': np.min(poi['longitude']-1),
                  'urcrnrlon': np.max(poi['longitude']+1)}
        ridx = match_region_rect(self.vars['latitude'],
                                 self.vars['longitude'],
                                 region=region)
        newdict = deepcopy(self.vars)
        idx = [self._poi_sat(newdict, self.twin, self.distlim, poi, ridx, i)
               for i in tqdm(range(len(poi['datetime'])))]
        idx = list(np.array(ridx)[flatten(idx)])
        return idx

    @staticmethod
    def _poi_sat(indict, twin, distlim, poi, ridx, i):
        """
        return: indices for values matching the spatial and
                temporal constraints
        """
        tidx = find_included_times(
                    list(np.array(indict['datetime'])[ridx]),
                    target_t=poi['datetime'][i],
                    twin=twin)
        slons = list(np.array(indict['longitude'])[ridx][tidx])
        slats = list(np.array(indict['latitude'])[ridx][tidx])
        plons = [poi['longitude'][i]]*len(slons)
        plats = [poi['latitude'][i]]*len(slats)
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
        if (region not in region_dict['poly'] and
        region not in model_dict and
        region not in region_dict['geojson']):
            if region is None:
                region = 'global'
            else:
                if (region not in region_dict['rect']
                and region not in region_dict['geojson']
                and isinstance(region, dict) is False):
                    sys.exit("Region is not defined")
                elif isinstance(region, dict):
                    print("Region bound by min/max of poi coordinates")
                else:
                    print("Specified region: " + region + "\n"
                          + " --> Bounds: " +
                          str(region_dict['rect'][region]))
                    region = region_dict['rect'][region]
            ridx = match_region_rect(LATS, LONS, region=region)
        elif region in region_dict['geojson']:
            print("Region is defined as geojson")
            ridx = match_region_geojson(LATS, LONS, region=region)
        else:
            ridx = match_region_poly(LATS, LONS, region=region,
                                     grid_date=grid_date)
        return ridx


def poi_sat(indict, twin, distlim, poi, ridx, i):
    """
    return: indices for values matching the spatial and
            temporal constraints
    """
    tidx = find_included_times(
                list(np.array(indict['datetime'])[ridx]),
                target_t=poi['datetime'][i],
                twin=twin)
    slons = list(np.array(indict['longitude'])[ridx][tidx])
    slats = list(np.array(indict['latitude'])[ridx][tidx])
    plons = [poi['longitude'][i]]*len(slons)
    plats = [poi['latitude'][i]]*len(slats)
    dists_tmp = haversineA(slons, slats, plons, plats)
    sidx = np.argwhere(np.array(dists_tmp) <= distlim).flatten()
    #dists = list(np.array(dists_tmp)[sidx])
    return list(np.array(tidx)[sidx])

def match_poi(indict, twin, distlim, poi):
    """
    return: idx that match to region
    """
    from tqdm import tqdm
    print('Match up poi locations')
    region={'llcrnrlat':np.min(poi['latitude']),
            'urcrnrlat':np.max(poi['latitude']),
            'llcrnrlon':np.min(poi['longitude']),
            'urcrnrlon':np.max(poi['longitude'])}
    ridx = match_region_rect(indict['latitude'],
                             indict['longitude'],
                             region=region)
    sat_dict = deepcopy(indict)
    idx = [poi_sat(sat_dict,twin,distlim,poi,ridx,i) \
                  for i in tqdm(range(len(poi['datetime'])))]
    idx = list(np.array(ridx)[flatten(idx)])
    return idx

def match_region(LATS,LONS,region,grid_date):
    """
    Function to filter satellite data according to region

    return:
        indices that match the region
    """
    # region in region_dict[poly]:
    # find values for given region
    if (region not in region_dict['poly'] and \
    region not in model_dict and \
    region not in region_dict['geojson']):
        if region is None:
            region = 'global'
        else:
            if (region not in region_dict['rect'] \
            and region not in region_dict['geojson'] \
            and isinstance(region,dict)==False):
                sys.exit("Region is not defined")
            elif isinstance(region,dict):
                print("Region bound by min/max of poi coordinates")
            else:
                print("Specified region: " + region + "\n"
                + " --> Bounds: " + str(region_dict['rect'][region]))
                region = region_dict['rect'][region]
        ridx = match_region_rect(LATS,LONS,region=region)
    elif region in region_dict['geojson']:
        print("Region is defined as geojson")
        ridx = match_region_geojson(LATS,LONS,region=region)
    else:
        ridx = match_region_poly(LATS,LONS,region=region,
                                grid_date=grid_date)
    return ridx

def match_region_rect(LATS,LONS,region):
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
        ridx = np.where(\
                    (np.array(LATS)>llcrnrlat)\
                    &(np.array(LATS)<urcrnrlat)\
                    &(np.array(LONS)>llcrnrlon)\
                    &(np.array(LONS)<urcrnrlon)\
                    )[0]
    print (len(ridx), " values found for chosen region and time frame.")
    return ridx

def match_region_geojson(LATS,LONS,region):
    """
    Takes care of regions defines as geojson
    """
    import geojson
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    LATS = list(LATS)
    LONS = list(LONS)
    lats = np.array(LATS).ravel()
    lons = np.array(LONS).ravel()
    points = np.c_[lons,lats]
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
        hits = Path(poly.xy).contains_points(points,radius=1e-9)
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
            hits = Path(poly.xy).contains_points(points,radius=1e-9)
            ridx_tmp = list(np.array(range(len(LONS)))[hits])
            ridx += ridx_tmp
    return ridx

def match_region_poly(LATS,LONS,region,grid_date):
    """
    Takes care of region defined as polygon
    """
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    import numpy as np
    if (region not in region_dict['poly'] \
        and region not in model_dict):
        sys.exit("Region polygone is not defined")
    elif isinstance(region,dict)==True:
        print ("Manuall specified region: \n"
            + " --> Bounds: " + str(region))
        poly = Polygon(list(zip(region['lons'],
            region['lats'])), closed=True)
    elif (isinstance(region,str)==True and region in model_dict):
        try:
            print('Use date for retrieving grid: ', grid_date)
            filestr = make_model_filename_wrapper(\
                                    region,grid_date,'best')
            meta = ncdumpMeta(filestr)
            flon = get_filevarname('lons',variable_info,\
                                model_dict[region],meta)
            flat = get_filevarname('lats',variable_info,\
                                model_dict[region],meta)
            time = get_filevarname('time',variable_info,\
                                model_dict[region],meta)
            #M = xa.open_dataset(filestr, decode_cf=True)
            #model_lons = M[flon].data
            #model_lats = M[flat].data
            model_lons, model_lats, _ = \
                read_model_nc_output_lru(filestr,flon,flat,time)
        except (KeyError,IOError,ValueError) as e:
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
            filestr = make_model_filename_wrapper(\
                                    region,grid_date,'best')
            meta = ncdumpMeta(filestr)
            flon = get_filevarname('lons',variable_info,\
                                model_dict[region],meta)
            flat = get_filevarname('lats',variable_info,\
                                model_dict[region],meta)
            time = get_filevarname('time',variable_info,\
                                model_dict[region],meta)
            #M = xa.open_dataset(filestr, decode_cf=True)
            #model_lons = M[flon].data
            #model_lats = M[flat].data
            model_lons, model_lats, _ = \
                read_model_nc_output_lru(filestr,flon,flat,time)
        if (len(model_lons.shape)==1):
            model_lons, model_lats = np.meshgrid(
                                    model_lons,
                                    model_lats
                                    )
        print('Check if footprints fall within the chosen domain')
        ncdict = ncdumpMeta(filestr)
        try:
            proj4 = find_attr_in_nc('proj',ncdict=ncdict,
                                    subattrstr='proj4')
        except IndexError:
            print('proj4 not defined in netcdf-file')
            print('Using proj4 from model config file')
            proj4 = model_dict[region]['proj4']
        proj_model = pyproj.Proj(proj4)
        Mx, My = proj_model(model_lons,model_lats,inverse=False)
        Vx, Vy = proj_model(LONS,LATS,inverse=False)
        xmax, xmin = np.max(Mx), np.min(Mx)
        ymax, ymin = np.max(My), np.min(My)
        ridx = list( np.where((Vx>xmin) & (Vx<xmax) &
                              (Vy>ymin) & (Vy<ymax))[0] )
    elif isinstance(region,str)==True:
        print ("Specified region: " + region + "\n"
        + " --> Bounded by polygon: \n"
        + "lons: " + str(region_dict['poly'][region]['lons']) + "\n"
        + "lats: " + str(region_dict['poly'][region]['lats']))
        poly = Polygon(list(zip(region_dict['poly'][region]['lons'],
            region_dict['poly'][region]['lats'])), closed=True)
        # check if coords in region
        LATS = list(LATS)
        LONS = list(LONS)
        lats = np.array(LATS).ravel()
        lons = np.array(LONS).ravel()
        points = np.c_[lons,lats]
        # radius seems to be important to correctly define polygone
        # see discussion here:
        # https://github.com/matplotlib/matplotlib/issues/9704
        hits = Path(poly.xy).contains_points(points,radius=1e-9)
        ridx = list(np.array(range(len(LONS)))[hits])
    if not ridx:
        print ("No values for chosen region and time frame!!!")
    else:
        print ("Values found for chosen region and time frame.")
    return ridx
