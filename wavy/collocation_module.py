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
from datetime import datetime
import time
import pyresample
from tqdm import tqdm
from copy import deepcopy
import xarray as xr
#import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.utils import collocate_times
from wavy.utils import make_fc_dates
from wavy.utils import hour_rounder_pd
from wavy.utils import NoStdStreams
from wavy.utils import parse_date
from wavy.utils import flatten
from wavy.utils import compute_quantiles
from wavy.wconfig import load_or_default
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.model_module import model_class as mc
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.validationmod import validate, disp_validation
# ---------------------------------------------------------------------#

# read yaml config files:
model_dict = load_or_default('model_cfg.yaml')
insitu_dict = load_or_default('insitu_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')

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

def get_model_filename(nID, d, leadtime):
    mco = mc(nID, d, leadtime)
    return mco._make_model_filename_wrapper(d, leadtime)

def find_valid_fc_dates_for_model_and_leadtime(fc_dates, model, leadtime):
    '''
    Finds valid dates that are close to desired dates at a precision
    of complete hours
    '''
    fc_dates_new = hour_rounder_pd(fc_dates)
    if (leadtime is None or leadtime == 'best'):
        pass
    else:
        fc_dates_new = [d for d in fc_dates_new
                if get_model_filename(model, d, leadtime) != False]
    return fc_dates_new

def check_if_file_is_valid(fc_date, model, leadtime, max_lt=None):
    fname = get_model_filename(
                model, fc_date, leadtime, max_lt=max_lt)
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
    distlim=None, leadtime=None, max_lt=None, **kwargs):
        #twin=30, distlim=6):
        print('# ----- ')
        print(" ### Initializing collocation_class object ###")
        print(" ")
        # make clones to prevent overwriting
        self.varalias = oco.varalias
        self.model = model
        self.leadtime = leadtime
        self.oco = oco
        self.nID = oco.nID
        self.obstype = str(type(oco))[8:-2]
        self.stdvarname = oco.stdvarname
        self.region = oco.region
        self.units = variable_def[self.varalias].get('units')
        self.sd = oco.sd
        self.ed = oco.ed
        self.twin = kwargs.get('twin', oco.twin)
        self.distlim = kwargs.get('distlim', 6)
        print(" ")
        print(" ## Collocate ... ")
        try:
            t0 = time.time()
            results_dict = self.collocate(**kwargs)
            self.model_time = results_dict['model_time']
            # build xarray dataset from results
            ds = self._build_xr_dataset(results_dict)
            self.vars = ds
            t1 = time.time()
            print(" ")
            print(" ## Summary:")
            print(len(self.vars['time']), " values collocated.")
            print("Time used for collocation:", round(t1-t0, 2), "seconds")
            print(" ")
            print(" ### Collocation_class object initialized ###")
        except Exception as e:
            print(e)
            self.error = e
            print ("! No collocation_class object initialized !")
        # add class variables
        print('# ----- ')

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

    def _collocate_field(self, mco, tmp_dict, **kwargs):
        """
        Some info
        """
        Mlons = mco.vars.lons
        Mlats = mco.vars.lats
        Mvars = mco.vars[mco.varalias]
        if len(Mlons.shape) > 2:
            Mlons = Mlons[0, :].data.squeeze()
            Mlats = Mlats[0, :].data.squeeze()
            Mvars = Mvars[0, :].data.squeeze()
        else:
            Mlons = Mlons.data.squeeze()
            Mlats = Mlats.data.squeeze()
            Mvars = Mvars.data.squeeze()
        obs_vars = tmp_dict[self.varalias]
        obs_lons = tmp_dict['lons']
        obs_lats = tmp_dict['lats']
        # Compare wave heights of satellite with model with
        #   constraint on distance and time frame
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
                                    ndt, self.model, self.leadtime)

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
            try:
                with NoStdStreams():
                    # filter needed obs within time period
                    target_date = [parse_date(str(fc_date[i]))]
                    idx = collocate_times(ndt_datetime,
                                          target_t=target_date,
                                          twin=self.twin)

                    print(len(idx), "footprints to be collocated")
                    # make tmp obs_obj with filtered data
                    tmp_dict = {}
                    tmp_dict['time'] = self.oco.vars['time'].values[idx]
                    tmp_dict['lats'] = self.oco.vars['lats'].values[idx]
                    tmp_dict['lons'] = self.oco.vars['lons'].values[idx]
                    tmp_dict[self.oco.varalias] = \
                            self.oco.vars[self.oco.varalias].values[idx]
                    mco = mc(sd=fc_date[i], ed=fc_date[i], nID=self.model,
                             leadtime=self.leadtime, varalias=self.varalias,
                             **kwargs)
                    mco = mco.populate(**kwargs)
                    results_dict_tmp = self._collocate_field(
                                            mco, tmp_dict, **kwargs)
                    if (len(results_dict_tmp['model_values']) > 0):
                        # append to dict
                        results_dict['model_time'].append(fc_date[i])
                        results_dict['obs_time'].append(tmp_dict['time'])
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
                    else: pass
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
            results_dict = self._collocate_track(**kwargs)
        # same needed for insitu, multiins, multisat, consolidate
        return results_dict


    def write_to_nc(self, pathtofile=None, file_date_incr=None):
        if 'error' in vars(self):
            print('Erroneous collocation_class file detected')
            print('--> dump to netCDF not possible !')
        else:
            print('to be implemented ...')

    def write_to_pickle(self, pathtofile=None):
        import pickle
        # writing
        pickle.dump(self, open(pathtofile, "wb" ))
        print('collocation_class object written to:', pathtofile)
        # for reading
        # cco = pickle.load( open( pathtofile, "rb" ) )

    def validate_collocated_values(self,**kwargs):
        dtime = self.vars['time']
        mods = self.vars['model_values']
        obs = self.vars['obs_values']
        sdate = self.vars['time'][0]
        edate = self.vars['time'][-1]
        validation_dict = validate_collocated_values(
                                dtime,obs,mods,\
                                sdate=sdate,edate=edate,\
                                **kwargs)
        return validation_dict

def validate_collocated_values(dtime, obs, mods, **kwargs):
    target_t, sdate, edate, twin = None, None, None, None
    if ('col_obj' in kwargs.keys() and kwargs['col_obj'] is not None):
        col_obj = kwargs['col_obj']
        mods = col_obj.vars['model_values']
        obs = col_obj.vars['obs_values']
        dtime = col_obj.vars['datetime']
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

