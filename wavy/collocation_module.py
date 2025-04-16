#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
"""
- Module that should take care of collocation of points or swaths
- Needs input from modules that retrieve from observational platforms
  and models
"""
# --- import libraries ------------------------------------------------#
# standard library imports
import numpy as np
import netCDF4
from datetime import datetime, timedelta
import time
import pyresample
from tqdm import tqdm
from copy import deepcopy
import xarray as xr
import pandas as pd
import copy
#import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.utils import collocate_times
from wavy.utils import make_fc_dates
from wavy.utils import hour_rounder_pd, hour_rounder
from wavy.utils import NoStdStreams
from wavy.utils import parse_date
from wavy.utils import flatten
from wavy.utils import compute_quantiles
from wavy.utils import haversineA
from wavy.wconfig import load_or_default
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.model_module import model_class as mc
from wavy.gridder_module import gridder_class as gc
from wavy.grid_stats import apply_metric
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.validationmod import validate, disp_validation
# ---------------------------------------------------------------------#

# read yaml config files:
model_dict = load_or_default('model_cfg.yaml')
insitu_dict = load_or_default('insitu_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')


def collocate_observations(co1, co2, twin=5, dist_max=200):
    '''
     Collocates observations from two wavy objects, keeping only closest
     data from the second dataset within a given time range around each 
     observation from the first given dataset.

    Args:
        co1 (insitu_class or satellite_class object): first observation dataset
                                                    the data from the second
                                                    dataset will be collocated
                                                    around this one.
        co2 (insitu_class or satellite_class object): second observation 
                                                      dataset
        twin (int): time window length in minutes
        dist_max (int | float): Maximum distance
                                accepted to collocate
                                data, in km

    Returns:
        tuple (co1_filter (insitu_class or satellite_class object),
               co2_filter (insitu_class or satellite_class object)): first and
                                               second objects given as inputs, 
                                               with data filtered to keep only 
                                               collocated observations
    '''
    list_time_co2 = []
    list_time_co1 = []
    list_dist_min = []
    datetimes_co1 = [pd.to_datetime(t) for t in co1.vars['time'].values]

    for i, time_co1 in enumerate(datetimes_co1):

        # For each co1 measure, create a dataset
        # ds_tmp, of co2 data measured at times
        # between [time_co1 - twin, time_co1 + twin]
        time_sup = time_co1 + timedelta(minutes=twin)
        time_inf = time_co1 - timedelta(minutes=twin)
        ds_tmp = co2.vars.sel(time=slice(time_inf, time_sup))
        ds_tmp_size = ds_tmp.sizes['time']

        # if some co2 data are found
        if ds_tmp_size > 0:

            # Gets co1 corresponding datetime,
            # lat and lon and adds it to ds_tmp
            # as variables
            co1_vars_tmp = co1.vars.isel(time=i)
            lats_co1 = co1_vars_tmp['lats'].values
            lons_co1 = co1_vars_tmp['lons'].values

            ds_tmp = ds_tmp.assign(
                {
                    'time_co1': (
                        'time',
                        np.repeat(time_co1, ds_tmp_size)
                    ),
                    'lats_co1': (
                        'time',
                        np.repeat(lats_co1, ds_tmp_size)
                    ),
                    'lons_co1': (
                        'time',
                        np.repeat(lons_co1, ds_tmp_size)
                    )
                }
            )

            # Calculates distance between first and second
            # measurements and adds it as a variable to ds_tmp
            ds_tmp = ds_tmp.assign(
                {
                    'dist_co1_co2': (
                        'time',
                        haversineA(
                            ds_tmp['lons_co1'],
                            ds_tmp['lats_co1'],
                            ds_tmp['lons'],
                            ds_tmp['lats']
                        )[0]
                    )
                }
            )

            # Calculates the minimum value for distances
            min_dist_tmp = min(ds_tmp['dist_co1_co2'].values)

            # If the minimum value is lesser than defined max distance
            if min_dist_tmp <= dist_max:

                # Times for co1 and co2 measurements,
                # corresponding to the minimum distance are fetched
                dist = list(ds_tmp['dist_co1_co2'].values)
                ds_tmp = ds_tmp.isel(time=dist.index(min_dist_tmp))
                list_time_co2.append(ds_tmp['time'].data)
                list_dist_min.append(min_dist_tmp)
                list_time_co1.append(ds_tmp['time_co1'].data)

    co1_filter = copy.deepcopy(co1)
    co2_filter = copy.deepcopy(co2)

    # Original co1 and co2 object vars are filtered using
    # lists of times corresponding to minimum distance
    # between co1 and co2 observation
    co2_filter.vars = co2_filter.vars.sel(time=list_time_co2)
    co2_filter.vars = co2_filter.vars.assign(
        {
            'colloc_dist': (
                'time',
                list_dist_min
            )
        }
    )
    co1_filter.vars = co1_filter.vars.sel(time=list_time_co1)

    return co1_filter, co2_filter 


def collocation_fct(obs_lons, obs_lats, model_lons, model_lats):
    grid = pyresample.geometry.GridDefinition(
                                lats=model_lats,
                                lons=model_lons)
    # Define some sample points
    swath = pyresample.geometry.SwathDefinition(lons=obs_lons,
                                                lats=obs_lats)
    # Determine nearest (great circle distance) neighbour in the grid.
    # valid_input_index, valid_output_index, index_array, distance_array = \
    _, valid_output_index, index_array, distance_array = \
                            pyresample.kd_tree.get_neighbour_info(
                                source_geo_def=grid,
                                target_geo_def=swath,
                                radius_of_influence=1000000000,
                                neighbours=1)
    # get_neighbour_info() returns indices in the
    # flattened lat/lon grid. Compute the 2D grid indices:
    index_array_2d = np.unravel_index(index_array, grid.shape)
    return index_array_2d, distance_array, valid_output_index


def get_model_filename(nID, d, leadtime, **kwargs):
    mco = mc(nID=nID, sd=d, ed=d, leadtime=leadtime)
    return mco._make_model_filename_wrapper(parse_date(str(d)),
                                            leadtime, **kwargs)


def find_valid_fc_dates_for_model_and_leadtime(fc_dates, model,
                                               leadtime, colloc_time_method,
                                               **kwargs):
    '''
    Finds valid dates that are close to desired dates at a precision
    of complete hours
    '''
    #fc_dates_new = hour_rounder_pd(fc_dates)

    dt_tmp = [parse_date(str(d)) for d in fc_dates]
    fc_dates_new = [hour_rounder(d, method=colloc_time_method) for d in dt_tmp]

    fc_dates_new = np.unique(fc_dates_new)
    #if (leadtime is None or leadtime == 'best'):
    #    pass
    #else:
    fc_dates_new = [d for d in fc_dates_new
                    if get_model_filename(model, d, leadtime, **kwargs)
                    is not None]
    return fc_dates_new


def check_if_file_is_valid(fc_date, model, leadtime, **kwargs):
    fname = get_model_filename(model, fc_date, leadtime, **kwargs)
    print('Check if requested file:\n', fname, '\nis available and valid')
    try:
        nc = netCDF4.Dataset(fname, mode='r')
        time = nc.variables['time']
        dt = netCDF4.num2date(time[:], time.units)
        if fc_date in list(dt):
            print('File is available and contains requested date')
            return True
        else:
            print('Desired date ' + str(fc_date) + ' is not in', fname)
            return False
    except (FileNotFoundError, OSError) as e:
        print('File is not available or does not contain requested date')
        print(e)
        return False


def get_closest_date(overdetermined_lst, target_lst):
    idx = []
    for i in range(len(target_lst)):
        diffs = np.abs([(target_lst[i]
                        - overdetermined_lst[j]).total_seconds()
                        for j in range(len(overdetermined_lst))])
        mindiff = np.min(diffs)
        idx.append(list(diffs).index(mindiff))
    return idx


def adjust_dict_for_idx(indict, idx, excl_keys_lst):
    outdict = deepcopy(indict)
    for k in indict.keys():
        if k not in excl_keys_lst:
            outdict[k] = np.array(indict[k])[idx]
    return outdict


class collocation_class(qls):
    '''
    draft of envisioned collocation class object
    '''

    def __init__(self, oco=None, model=None, poi=None,
    distlim=None, leadtime=None, **kwargs):
        print('# ----- ')
        print(" ### Initializing collocation_class object ###")
        print(" ")
        # make clones to prevent overwriting
        self.varalias = oco.varalias
        self.varalias_obs = oco.varalias
        self.varalias_mod = kwargs.get('varalias', oco.varalias)
        self.model = model
        self.leadtime = leadtime
        self.oco = oco
        self.nID = oco.nID
        self.model = model
        self.obstype = str(type(oco))[8:-2]
        self.stdvarname = oco.stdvarname
        self.region = oco.region
        self.units = variable_def[self.varalias].get('units')
        self.sd = oco.sd
        self.ed = oco.ed
        self.twin = kwargs.get('twin', oco.twin)
        self.distlim = kwargs.get('distlim', 6)
        self.method = kwargs.get('method', 'closest')
        self.colloc_time_method = kwargs.get('colloc_time_method', 'nearest')
        print(" ")
        print(" ### Collocation_class object initialized ###")
        
    def populate(self, **kwargs):
    
        new = deepcopy(self)
        print(" ")
        print(" ## Collocate ... ")
        #for i in range(1):
        try:
            t0 = time.time()
            results_dict = new.collocate(**kwargs)
            new.model_time = results_dict['model_time']
            # build xarray dataset from results
            ds = new._build_xr_dataset(results_dict)
            ds = ds.assign_coords(time=ds.time.values)
            new.vars = ds
            new = new._drop_duplicates(**kwargs)
            t1 = time.time()
            print(" ")
            print(" ## Summary:")
            print(len(new.vars['time']), " values collocated.")
            print("Time used for collocation:", round(t1-t0, 2), "seconds")
            print(" ")

        except Exception as e:
            print(e)
            new.error = e
            new.vars = None
            print("! collocation_class object may be empty !")
        # add class variables
        print('# ----- ')
        
        return new

    def _build_xr_dataset(self, results_dict):
        ds = xr.Dataset({
                'time': xr.DataArray(
                        data=results_dict['obs_time'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['time'],
                        ),
                'dist': xr.DataArray(
                        data=results_dict['dist'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['dist'],
                        ),
                'obs_values': xr.DataArray(
                        data=results_dict['obs_values'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def[self.varalias],
                        ),
                'obs_lons': xr.DataArray(
                        data=results_dict['obs_lons'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['lons'],
                        ),
                'obs_lats': xr.DataArray(
                        data=results_dict['obs_lats'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['lats'],
                        ),
                'model_values': xr.DataArray(
                        data=results_dict['model_values'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def[self.varalias],
                        ),
                'model_lons': xr.DataArray(
                        data=results_dict['model_lons'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['lons'],
                        ),
                'model_lats': xr.DataArray(
                        data=results_dict['model_lats'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['lats'],
                        ),
                'colidx_x': xr.DataArray(
                        data=results_dict['collocation_idx_x'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['colidx_x'],
                        ),
                'colidx_y': xr.DataArray(
                        data=results_dict['collocation_idx_y'],
                        dims=['time'],
                        coords={'time': results_dict['obs_time']},
                        attrs=variable_def['colidx_y'],
                        ),
                    },
                attrs={'title': str(type(self))[8:-2] + ' dataset'}
                )
        return ds

    def _drop_duplicates(self, **kwargs):
        
        dim = kwargs.get('dim_duplicates', 'time')
        keep = kwargs.get('keep_duplicates', 'first')
        print('Removing duplicates according to', dim)
        print('Keeping', keep, 'value for the duplicates')
        new = deepcopy(self)
        new.vars = self.vars.drop_duplicates(dim=dim, keep=keep)
        print(str(int(abs(len(self.vars[dim])-len(new.vars[dim])))),
              'values removed')
        print('New number of footprints is:', str(int(len(new.vars[dim]))))
        return new

    def _collocate_field(self, mco, tmp_dict, **kwargs):
        """
        Some info
        """
        Mlons = mco.vars.lons.data
        Mlats = mco.vars.lats.data
        Mvars = mco.vars[mco.varalias].data
        if len(Mlons.shape) > 2:
            Mlons = Mlons[0, :].squeeze()
            Mlats = Mlats[0, :].squeeze()
            Mvars = Mvars[0, :].squeeze()
        elif len(Mlons.shape) == 2:
            Mlons = Mlons.squeeze()
            Mlats = Mlats.squeeze()
            Mvars = Mvars.squeeze()
        elif len(Mlons.shape) == 1:
            Mlons, Mlats = np.meshgrid(Mlons, Mlats)
            Mvars = Mvars.squeeze()
            assert len(Mlons.shape) == 2
        obs_vars = tmp_dict[self.varalias_obs]
        obs_lons = tmp_dict['lons']
        obs_lats = tmp_dict['lats']
        # Compare wave heights of satellite with model with
        # constraint on distance and time frame
        print("Perform collocation with distance limit\n",
              "distlim:", self.distlim)
        index_array_2d, distance_array, _ =\
                                    collocation_fct(
                                    obs_lons, obs_lats,
                                    Mlons, Mlats)
        # caution: index_array_2d is tuple
        # impose distlim
        dist_idx = np.where((distance_array < self.distlim*1000) &
                            (~np.isnan(Mvars[index_array_2d[0],
                             index_array_2d[1]])))[0]
        idx_x = index_array_2d[0][dist_idx]
        idx_y = index_array_2d[1][dist_idx]
        results_dict = {
                'dist': list(distance_array[dist_idx]),
                'model_values': list(Mvars[idx_x, idx_y]),
                'model_lons': list(Mlons[idx_x, idx_y]),
                'model_lats': list(Mlats[idx_x, idx_y]),
                'obs_values': list(obs_vars[dist_idx]),
                'obs_lons': list(obs_lons[dist_idx]),
                'obs_lats': list(obs_lats[dist_idx]),
                'collocation_idx_x': list(idx_x),
                'collocation_idx_y': list(idx_y),
                'time': tmp_dict['time'][dist_idx]
                }
        return results_dict

    def _collocate_track(self, **kwargs):
        """
        Some info
        """
        print("run: _collocate_track")

        # use only dates with a model time step closeby
        # given the time constrains
        print('Filtering valid dates ...')
        t1 = time.time()

        ndt = self.oco.vars.time.data
        ndt_datetime = [parse_date(str(d)) for d in ndt]

        ndt_valid = find_valid_fc_dates_for_model_and_leadtime(
                                    ndt, self.model, self.leadtime,
                                    self.colloc_time_method, **kwargs)

        ndt_valid = np.unique(ndt_valid)

        fc_date = ndt_valid
        t2 = time.time()

        print(f'... done, used {t2-t1:.2f} seconds')

        print("Start collocation ...")
        results_dict = {
                'model_time': [],
                'obs_time': [],
                'dist': [],
                'model_values': [],
                'model_lons': [],
                'model_lats': [],
                'obs_values': [],
                'obs_lons': [],
                'obs_lats': [],
                'collocation_idx_x': [],
                'collocation_idx_y': [],
                }

        for i in tqdm(range(len(fc_date))):
            print(fc_date[i])
        #for i in range(len(fc_date)):
            try:
                for j in range(1):
                #with NoStdStreams():
                    # filter needed obs within time period
                    target_date = [parse_date(str(fc_date[i]))]
                    
                    # if method is 'nearest', get the values that fall within
                    # a time window of +/- 30 minutes of model time by default 
                    if self.colloc_time_method=='nearest':
                        idx = collocate_times(ndt_datetime,
                                              target_t=target_date,
                                              twin=self.twin)
                    # if method is 'floor' get the values that fall between
                    # model time and model time + 1 hour
                    elif self.colloc_time_method=='floor':
                        sdate_colloc = target_date[0]
                        edate_colloc = target_date[0] + timedelta(hours=1)
                        idx = collocate_times(ndt_datetime,
                                              target_t=target_date,
                                              sdate=sdate_colloc,
                                              edate=edate_colloc,
                                              twin=0)
                    # if method is 'ceil' get the values that fall between
                    # model time - 1 hour and model time
                    elif self.colloc_time_method=='ceil':
                        sdate_colloc = target_date[0] - timedelta(hours=1)
                        edate_colloc = target_date[0] 
                        idx = collocate_times(ndt_datetime,
                                              target_t=target_date,
                                              sdate=sdate_colloc,
                                              edate=edate_colloc,
                                              twin=0)                                

                    print(len(idx), "footprints to be collocated")
                    # make tmp obs_obj with filtered data
                    tmp_dict = {}
                    tmp_dict['time'] = self.oco.vars['time'].values[idx]
                    tmp_dict['lats'] = self.oco.vars['lats'].values[idx]
                    tmp_dict['lons'] = self.oco.vars['lons'].values[idx]
                    tmp_dict[self.oco.varalias] = \
                        self.oco.vars[self.oco.varalias].values[idx]
                    mco = mc(sd=fc_date[i], ed=fc_date[i], nID=self.model,
                             leadtime=self.leadtime, varalias=self.varalias_mod,
                             **kwargs)
                    mco = mco.populate(**kwargs)
                    results_dict_tmp = self._collocate_field(
                                            mco, tmp_dict, **kwargs)
                    if (len(results_dict_tmp['model_values']) > 0):
                        # append to dict
                        results_dict['model_time'].append(fc_date[i])
                        results_dict['obs_time'].append(
                                results_dict_tmp['time'])
                        results_dict['dist'].append(
                                results_dict_tmp['dist'])
                        results_dict['model_values'].append(
                                results_dict_tmp['model_values'])
                        results_dict['model_lons'].append(
                                results_dict_tmp['model_lons'])
                        results_dict['model_lats'].append(
                                results_dict_tmp['model_lats'])
                        results_dict['obs_values'].append(
                                results_dict_tmp['obs_values'])
                        results_dict['obs_lats'].append(
                                results_dict_tmp['obs_lats'])
                        results_dict['obs_lons'].append(
                                results_dict_tmp['obs_lons'])
                        results_dict['collocation_idx_x'].append(
                                        results_dict_tmp['collocation_idx_x'])
                        results_dict['collocation_idx_y'].append(
                                        results_dict_tmp['collocation_idx_y'])
                    else:
                        pass
                    if 'results_dict_tmp' in locals():
                        del results_dict_tmp
            except (ValueError, FileNotFoundError, OSError) as e:
                # ValueError, pass if no collocation
                # FileNotFoundError, pass if file not accessible
                # OSError, pass if file not accessible from thredds
                logger.exception(e)
                print(e)
        # flatten all aggregated entries
        results_dict['model_time'] = results_dict['model_time']
        results_dict['obs_time'] = flatten(results_dict['obs_time'])
        results_dict['dist'] = flatten(results_dict['dist'])
        results_dict['model_values'] = flatten(results_dict['model_values'])
        results_dict['model_lons'] = flatten(results_dict['model_lons'])
        results_dict['model_lats'] = flatten(results_dict['model_lats'])
        results_dict['obs_values'] = flatten(results_dict['obs_values'])
        results_dict['obs_lats'] = flatten(results_dict['obs_lats'])
        results_dict['obs_lons'] = flatten(results_dict['obs_lons'])
        results_dict['collocation_idx_x'] = flatten(
                                    results_dict['collocation_idx_x'])
        results_dict['collocation_idx_y'] = flatten(
                                results_dict['collocation_idx_y'])
        return results_dict

    def _collocate_centered_model_value(self, time, lon, lat, **kwargs):

        #(time, lon, lat, nID_model, name_model, res):
   
        nID_model = self.model
        name_model = self.model 
        res = kwargs.get('res', (0.5,0.5))
        colloc_time_method = self.colloc_time_method
        
        print('Using resolution {}'.format(res))
        # ADD CHECK LIMITS FOR LAT AND LON
        res_dict = {}

        time = pd.to_datetime(time)
        time = hour_rounder(time, method=colloc_time_method)
        
        mco = mc(sd=time, ed=time,
                 nID=nID_model, name=name_model,
                 max_lt=12).populate(twin=5) # ADD AS PARAMETERS
   
        bb = (lon - res[0]/2, 
              lon + res[0]/2, 
              lat - res[1]/2, 
              lat + res[1]/2)

        gco = gc(lons=mco.vars.lons.squeeze().values.ravel(),
                 lats=mco.vars.lats.squeeze().values.ravel(),
                 values=mco.vars.Hs.squeeze().values.ravel(),
                 bb=bb, res=res,
                 varalias=mco.varalias,
                 units=mco.units,
                 sdate=mco.vars.time,
                 edate=mco.vars.time)
    
        gridvar, lon_grid, lat_grid = apply_metric(gco=gco)
    
        ts = gridvar['mor'].flatten()
        lon_flat = lon_grid.flatten()
        lat_flat = lat_grid.flatten()
    
        res_dict['hs'] = ts[0]
        res_dict['lon'] = lon_flat[0]
        res_dict['lat'] = lat_flat[0]
        res_dict['time'] = time

        return res_dict

    def _collocate_regridded_model(self, **kwargs):
    
        from joblib import Parallel, delayed
    
        hs_mod_list=[] 
        lon_mod_list=[]
        lat_mod_list=[]
        time_mod_list=[]
        
        nproc = kwargs.get('nproc', 16)
        
        oco_vars = self.oco.vars

        length = len(oco_vars.time.values)
        
        #Parallel should be optional, with nproc as parameter
        colloc_mod_list = Parallel(n_jobs=nproc)(
                               delayed(self._collocate_centered_model_value) (
                                           oco_vars.time.values[i],
                                           oco_vars.lons.values[i],
                                           oco_vars.lats.values[i],
                                           **kwargs) for i in range(length)
                                           )

        length_colloc_mod_list = len(colloc_mod_list)
        hs_mod_list = [colloc_mod_list[i]['hs'] for i in \
                       range(length_colloc_mod_list)]
        lon_mod_list = [colloc_mod_list[i]['lon'] for i in \
                       range(length_colloc_mod_list)]
        lat_mod_list = [colloc_mod_list[i]['lat'] for i in \
                       range(length_colloc_mod_list)]
        time_mod_list = [colloc_mod_list[i]['time'] for i in \
                       range(length_colloc_mod_list)]

        mod_colloc_vars = xr.Dataset(
            {
             "lats": (
                 ("time"),
                 lat_mod_list,
             ),
             "lons": (
                 ("time"), 
                 lon_mod_list
             ),
             "Hs": (
                 ("time"),
                 hs_mod_list
             )
          },
         coords={"time": time_mod_list},
        )

        results_dict = {
            'model_time': time_mod_list,
            'obs_time': oco_vars.time.values,
            'dist': [0]*length,
            'model_values': hs_mod_list,
            'model_lons': lon_mod_list,
            'model_lats': lat_mod_list,
            'obs_values': oco_vars.Hs.values,
            'obs_lons': oco_vars.lons.values,
            'obs_lats': oco_vars.lats.values,
            'collocation_idx_x': [0]*length,
            'collocation_idx_y': [0]*length,
            }
        
        return results_dict

    def collocate(self, **kwargs):
        """
        get obs value for model value for given
            temporal and spatial constraints
        """
        if (self.oco is None and len(self.oco.vars[self.oco.stdvarname]) < 1):
            raise Exception('\n###\n'
                    + 'Collocation not possible, '
                    + 'no observation values for collocation!'
                    + '\n###')
        if (self.model is None):
            raise Exception ('\n###\n'
                            +'Collocation not possible, '
                            +'no model available for collocation!'
                            +'\n###'
                            )
        if ((self.model is not None) and (self.oco is not None)):
            if self.method == 'closest':
                results_dict = self._collocate_track(**kwargs)
            elif self.method == 'regridded':
                results_dict = self._collocate_regridded_model(**kwargs)

        return results_dict


    def validate_collocated_values(self, **kwargs):
        times = self.vars['time']
        dtime = [parse_date(str(t.data)) for t in times]
        mods = self.vars['model_values']
        obs = self.vars['obs_values']
        sdate = dtime[0]
        edate = dtime[-1]
        validation_dict = validate_collocated_values(
                                dtime, obs, mods,
                                sdate=sdate, edate=edate,
                                **kwargs)
        return validation_dict


def validate_collocated_values(dtime, obs, mods, **kwargs):
    target_t, sdate, edate, twin = None, None, None, None
    if ('col_obj' in kwargs.keys() and kwargs['col_obj'] is not None):
        col_obj = kwargs['col_obj']
        mods = col_obj.vars['model_values']
        obs = col_obj.vars['obs_values']
        dtime = col_obj.vars['time']
    # get idx for date and twin
    if 'target_t' in kwargs.keys():
        target_t = kwargs['target_t']
    if 'sdate' in kwargs.keys():
        sdate = kwargs['sdate']
    if 'edate' in kwargs.keys():
        edate = kwargs['edate']
    if 'twin' in kwargs.keys():
        twin = kwargs['twin']
    idx = collocate_times(dtime,
                          target_t=target_t,
                          sdate=sdate,
                          edate=edate,
                          twin=twin)
    mods = np.array(mods)[idx]
    obs = np.array(obs)[idx]
    results_dict = {'model_values': mods, 'obs_values': obs}

    # validate
    validation_dict = validate(results_dict)
    disp_validation(validation_dict)

    return validation_dict

