import os
import numpy as np
from copy import deepcopy
import netCDF4
import pandas as pd
import roaring_landmask
from scipy.stats import circmean
from datetime import timedelta
from wavy.utils import flatten
import cartopy.io.shapereader as shpreader
from sklearn.neighbors import BallTree
from pyproj import Proj, Geod
from math import *
import logging

import pyresample as pr
from roaring_landmask import Shapes, LandmaskProvider
import shapely.wkb as wkb
from shapely.geometry import Polygon, mapping


# own imports
from wavy.utils import find_included_times, collocate_times
from wavy.wconfig import load_or_default

ROAR = None
variable_def = load_or_default('variable_def.yaml')


class filter_class:

    def apply_limits(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply limits (crude cleaning using valid range)')
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else:
            varalias = kwargs.get('varalias', new.varalias)
        llim = kwargs.get('llim',
                          variable_def[varalias]['valid_range'][0])
        ulim = kwargs.get('ulim',
                          variable_def[varalias]['valid_range'][1])
        ds = deepcopy(new.vars)
        y = ds[varalias]
        tmpdict = {'y': y}
        df = pd.DataFrame(data=tmpdict)
        dfmask = df['y'].between(llim, ulim, inclusive='both')
        idx = np.array(range(len(dfmask)))[dfmask]
        ds = new.vars.isel(time=idx)
        print(" Number of disregarded values:", len(dfmask[dfmask == False]))
        print(" Number of remaining values:", len(ds.time))
        new.vars = ds
        return new

    def filter_landMask(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply land mask')
        new = deepcopy(self)
        longitudes = np.array(new.vars['lons'])
        latitudes = np.array(new.vars['lats'])
        # apply land mask
        sea_mask = apply_land_mask(longitudes, latitudes, **kwargs)
        # impose on dataset
        idx = np.array(range(len(sea_mask)))
        idx = idx[sea_mask]
        ds = new.vars.isel(time=idx)

        indices = start_stop(sea_mask, True)
        indices_tmp = deepcopy(indices)
        len_of_indices_tmp = len(list(indices_tmp))

        del indices_tmp

        if len_of_indices_tmp == 0:
            pass

        else:
            no_chunks = 0
            for start_idx, stop_idx in indices:
                logger.info(' start_idx: ' +
                            str(start_idx) +
                            'stop_idx: ' + str(stop_idx))
                no_chunks += 1
                lenofchunk = len(list(range(start_idx, stop_idx)))
                logger.info(' -> Length of chunk: ' + str(lenofchunk))
        # add chunks to object for further use
        new.land_sea_chunks = indices
        # Assign back to class object
        new.vars = ds
        print(' Number of registered intersections with land:', no_chunks)
        print(' Number of disregarded values due to land intersections:',
              len(sea_mask[sea_mask is False]))
        print(' Number of remaining values:', len(new.vars['time']))
        print(' land_sea_chunks added to self')
        return new

    def filter_distance_to_coast(self, llim=0, ulim=100000000, **kwargs):
        """
        discards all values closer to shoreline than threshold
        """
        print("Apply distance_to_coast_mask")
        new = deepcopy(self)
        longitudes = np.array(new.vars['lons'])
        latitudes = np.array(new.vars['lats'])
        w = Shapes.wkb(LandmaskProvider.Gshhg)
        polys = wkb.loads(w)
        mapped = mapping(polys)
        c = mapped['coordinates']
        cA = np.vstack([np.flip(x[0][:]) for x in c])
        points_sdef = pr.geometry.SwathDefinition(longitudes, latitudes)
        coast_sdef = pr.geometry.SwathDefinition(cA[:, 1], cA[:, 0])
        _, _, _, distance_array = pr.kd_tree.get_neighbour_info(
            coast_sdef, points_sdef, 10000000, neighbours=1)
        # get rid of infs
        mask = np.where((distance_array > llim) & (distance_array < ulim))[0]
        # new.dist_to_coast = distance_array[mask]
        # impose on dataset
        ds = new.vars.isel(time=mask)
        # add to dataset
        ds = ds.assign({"dist_to_coast": (("time"), distance_array[mask])})
        ds["dist_to_coast"].attrs = variable_def["dist_to_coast"]
        # Assign back to class object
        new.vars = ds
        print(" Number of disregarded values:",
              (len(distance_array)-len(mask)))
        print(" Number of remaining values:", len(new.vars['time']))
        return new

    #def filter_blockMean(self, **kwargs):
    #    print('Apply blockMean')
    #    return self

    def filter_lanczos(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply lanczos filter')
        from wavy.utils import runmean
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else: 
            varalias = kwargs.get('varalias', new.varalias)
        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        ts_lst = []
        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j+1)
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    window = kwargs.get('window')
                    cutoff = kwargs.get('cutoff')
                    weights = lanczos_weights(window, cutoff)
                    ts, _ = runmean(y, window,
                                    mode='centered',
                                    weights=weights)
                    ts_lst.append(ts)
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))
        new.vars[varalias].values = flatten(ts_lst)
        return new

    def filter_runmean(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply running mean filter')
        from wavy.utils import runmean
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg="Variable {} selected by default. ".format(varalias)\
                + "If you wish to apply the filter to another variable, "+\
                      "please specify it using varalias."
                logger.warning(msg)
        else: 
            varalias = kwargs.get('varalias', new.varalias)

        print("Applying filter to {}".format(varalias))

        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        mode = kwargs.get('mode', 'centered')
        indices = new.slider_chunks(slider=win, overlap=ol)

        ts_lst = []
        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j)
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    window = kwargs.get('window')
                    ts, _ = runmean(y, window,
                                    mode=mode)
                    ts_lst.append(ts)
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))
        new.vars[varalias].values = flatten(ts_lst)
        return new

    def filter_GP(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply GPR filter')
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else:
            varalias = kwargs.get('varalias', new.varalias)
        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        ts_lst = []
        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j+1)
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    x = tmp_ds['time'].values[tmp_tgc_idx].astype(float)
                    X = x  # points for prediction
                    ts = smoother_GP(x, y, X, **kwargs)
                    ts_lst.append(ts)
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))
        new.vars[varalias].values = flatten(ts_lst)

        return new

    #def filter_NIGP(self, **kwargs):
    #    return self

    def filter_linearGAM(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply LinearGAM filter')
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else:
            varalias = kwargs.get('varalias', new.varalias)
        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        ts_lst = []
        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j+1)
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    x = tmp_ds['time'].values[tmp_tgc_idx].astype(float)
                    X = x  # points for prediction
                    ts = smoother_linearGAM(x, y, X, **kwargs)
                    ts_lst.append(ts)
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx])
                else:
                    print("Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))
        new.vars[varalias].values = flatten(ts_lst)

        return new

    @staticmethod
    def cleaner_blockStd(y, **kwargs):
        meanval = np.nanmean(y)
        stdval = np.nanstd(y)
        sigma_multiplyer = kwargs.get("sigma", 2)
        uplim = meanval + (sigma_multiplyer*stdval)
        lowlim = meanval - (sigma_multiplyer*stdval)
        idx = [i for i in range(len(y))
               if (y[i] < uplim and y[i] > lowlim)]
        return idx

    def despike_blockStd(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply blockStd despiking')
        """
        Uses slider blocks as basis
        """
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else: 
            varalias = kwargs.get('varalias', new.varalias)

        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j)
            logger.info('tmp_idx ' + str(tmp_idx))
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                logger.info('tmp_tgc_idx ' + str(tmp_idx))
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    idx = new.cleaner_blockStd(y, **kwargs)
                    logger.info('idx ' + str(idx))
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx][idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))

        print(" Number of disregarded values:",
              (len(self.vars.time)-len(new.vars.time)))
        print(" Number of remaining values:", len(new.vars['time']))

        return new

    @staticmethod
    def cleaner_blockQ(y, **kwargs):
        llim_pct = kwargs.get("llim_pct", .05)
        ulim_pct = kwargs.get("ulim_pct", .95)
        uplim = np.quantile(y, ulim_pct)
        lowlim = np.quantile(y, llim_pct)
        idx = [i for i in range(len(y))
               if (y[i] < uplim and y[i] > lowlim)]
        return idx

    def despike_blockQ(self, **kwargs):
        print('Apply blockStd despiking')
        """
        Uses slider blocks as basis
        """
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                    + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else:
            varalias = kwargs.get('varalias', new.varalias)

        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j+1)
            logger.info('tmp_idx ' + str(tmp_idx))
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                logger.info('tmp_tgc_idx ' + str(tmp_idx))
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    idx = new.cleaner_blockQ(y, **kwargs)
                    logger.info('idx ' + str(idx))
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx][idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))

        print(" Number of disregarded values:",
              (len(self.vars.time)-len(new.vars.time)))
        print(" Number of remaining values:", len(new.vars['time']))

        return new

    def despike_GP(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply GPR despiking')
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                    + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else:
            varalias = kwargs.get('varalias', new.varalias)

        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j+1)
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    x = tmp_ds['time'].values[tmp_tgc_idx].astype(float)
                    idx = cleaner_GP(x, y, **kwargs)
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx][idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))

        print(" Number of disregarded values:",
              (len(self.vars.time)-len(new.vars.time)))
        print(" Number of remaining values:", len(new.vars['time']))

        return new

    def despike_NIGP(self, **kwargs):
        return self

    def despike_linearGAM(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        print('Apply GAM despiking')
        new = deepcopy(self)
        if isinstance(new.varalias, list):
            varalias = kwargs.get('varalias', new.varalias[0])

            if (len(new.varalias) > 1) and\
             (kwargs.get('varalias', None) is None):
                msg = "Variable {} selected by default. ".format(varalias)\
                    + "If you wish to apply the filter to another variable, " +\
                      "please specify it using varalias."
                logger.warning(msg)
        else:
            varalias = kwargs.get('varalias', new.varalias)

        # apply slider if needed
        win = kwargs.get('slider', len(new.vars.time))
        ol = kwargs.get('overlap', 0)
        indices = new.slider_chunks(slider=win, overlap=ol)

        tgc_idx_lst = []
        for i, j in indices:
            tmp_idx = range(i, j+1)
            # create tmp dataset reduced to i:j
            tmp_ds = new.vars.isel(time=tmp_idx)
            # apply gap chunks if needed
            pdtimes = tmp_ds.time.to_pandas()
            tgc_indices = new.time_gap_chunks(pdtimes, **kwargs)
            for k, l in tgc_indices:
                tmp_tgc_idx = range(k, l+1)
                # apply min chunk size
                if len(tmp_tgc_idx) > kwargs.get("chunk_min", 5):
                    y = tmp_ds[varalias].values[tmp_tgc_idx]
                    x = tmp_ds['time'].values[tmp_tgc_idx].astype(float)
                    X = x  # points for prediction
                    idx = cleaner_linearGAM(X, y, **kwargs)
                    tgc_idx_lst.append(np.array(tmp_idx)[tmp_tgc_idx][idx])
                else:
                    logger.warning(
                        "Chunk size to small -> not filtered and rejected")
                    pass

        new.vars = new.vars.isel(time=flatten(tgc_idx_lst))

        return new

    def slider_chunks(self, **kwargs):
        """
        Purpose: chunk data to ease computational load
        """
        new = deepcopy(self)
        slider = kwargs['slider']
        overlap = kwargs.get('overlap', 0)

        print('Using slider')
        print(' Splitting up dataset in chunks of', slider)
        print(' with overlap', overlap)

        start_idx_lst = []
        stop_idx_lst = []

        # create slider chunks
        for i in range(0, len(new.vars.time), slider):
            start_idx = np.max([i-overlap, 0])
            stop_idx = np.min([i+slider+overlap,
                               len(new.vars.time)-1])
            if start_idx != stop_idx:
                start_idx_lst.append(start_idx)
                stop_idx_lst.append(stop_idx)
        indices = zip(start_idx_lst, stop_idx_lst)
        print(' Number of created chunks:', len(start_idx_lst))
        return indices

    @staticmethod
    def time_gap_chunks(pdtime, **kwargs):
        """
        Purpose: chunk data according to sampling gaps to make
                 neighbour points match up and make filtering
                 meaningful.
        """
        print("Using time_gap_chunks")

        sr = kwargs.get('sampling_rate_Hz', 20)
        mask = (pdtime.diff() > pd.to_timedelta((1./sr)*2,
                'seconds')).to_numpy().copy()

        start_idx_lst = []
        stop_idx_lst = []

        # create gap chunks
        start_idx_lst.append(0)
        while True in mask:
            idx = list(mask).index(True)
            stop_idx_lst.append(idx-1)
            start_idx_lst.append(idx)
            mask[idx] = False

        stop_idx_lst.append(len(mask)-1)

        assert len(start_idx_lst) == len(stop_idx_lst)

        indices = zip(start_idx_lst, stop_idx_lst)

        print(" Number of created chunks:", len(start_idx_lst,))

        return indices

    def filter_footprint_radius(self, llim=None, ulim=None):
        """
        Filters all data according to given limits (llim, ulim)
        of footprint size
        """
        print("Apply filter_footprint_radius")
        new = deepcopy(self)
        if "fpr" in list(new.vars.keys()):
            pass
        else:
            new = new.compute_pulse_limited_footprint_radius()
        new.vars = new.vars.where(new.vars.fpr < ulim, drop=True)
        new.vars = new.vars.where(new.vars.fpr > llim, drop=True)
        print(" Number of disregarded values:",
              (len(self.vars.time)-len(new.vars.time)))
        print(" Number of remaining values:", len(new.vars['time']))
        return new

    def filter_footprint_land_interaction(self, **kwargs):
        """
        Checks if footprint interacts with land based on footprint size.
        Filters away the ones that do interact and returns a clean data set.
        """
        print("Apply filter_footprint_land_interaction")
        new = deepcopy(self)
        if "fpr" in list(new.vars.keys()):
            pass
        else:
            new = new.compute_pulse_limited_footprint_radius()
        _, _, _, _, ls_idx_lst = new._generate_xtrack_footprints(**kwargs)
        # apply indices to dataset
        new.vars = new.vars.isel(time=ls_idx_lst)
        print(" Number of disregarded values:",
              (len(self.vars.lons)-len(ls_idx_lst)))
        print(" Number of remaining values:", len(new.vars['time']))
        return new

    def _generate_xtrack_footprints(self, **kwargs):
        logger = logging.getLogger(__name__)
        log_level = str(kwargs.get('logging', 'WARNING').upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        domain = kwargs.get('domain', 'lonlat')
        n = kwargs.get('number_of_seeds', 250) + 1
        new = deepcopy(self)
        lons = new.vars.lons.values
        lats = new.vars.lats.values
        if domain == 'cartesian':
            pass
        elif domain == 'lonlat':
            lats_perp_lst = []
            lons_perp_lst = []
            ls_idx_lst = []
            mask = []
            for i in range(len(lons)):
                ls_idx = []
                if i < (len(lons)-1):
                    P1 = (lons[i], lats[i])
                    P2 = (lons[i+1], lats[i+1])
                else:
                    P1 = (lons[i], lats[i])
                    P2 = (lons[i-1], lats[i-1])
                lons_perp_lst_tmp = []
                lats_perp_lst_tmp = []
                for s in range(n):
                    P_perp_minus, P_perp_plus = \
                        new._generate_xtrack_footprints_in_lonlat(
                            P1, P2, s)
                    # check if within pulse limited footprint
                    dist = new._distance(P1[0], P1[1],
                                         P_perp_minus[0], P_perp_plus[1])
                    if dist > new.vars.fpr.values[i]:
                        pass
                    else:
                        lons_perp = \
                            np.transpose(
                                np.array([P_perp_minus, P_perp_plus]))[0]
                        lats_perp = \
                            np.transpose(
                                np.array([P_perp_minus, P_perp_plus]))[1]
                    # check if footprints intersect with land
                    sea_mask = apply_land_mask(lons_perp, lats_perp)
                    if False in sea_mask:
                        logger.info(
                            'Polution by land is detected for index ' + str(i))
                        logger.info(' -> Footprint not included!')
                        ls_idx.append(False)
                    else:
                        ls_idx.append(True)
                    # gather perpendicular footprints in lonlat
                    lons_perp_lst_tmp.append(lons_perp)
                    lats_perp_lst_tmp.append(lats_perp)
                if False in ls_idx:
                    #ls_idx_lst.append(False)
                    ls_idx_lst.append(i)
                    mask.append(False)
                else:
                    ls_idx_lst.append(i)
                    mask.append(True)
                lons_perp_lst.append(lons_perp_lst_tmp)
                lats_perp_lst.append(lats_perp_lst_tmp)

            # mask bad neighbours if desired
            if (kwargs.get('rm_bad_neighbours', False) is True):
                print("removing bad neighbours")
                k = kwargs.get('kneigh', 1)

                invalid = np.array(ls_idx_lst)[~np.array(mask)]
                expanded = np.any(
                        np.abs(
                        np.array(ls_idx_lst)[:, None] - invalid)\
                        <= k, axis=1)
                mask = mask & ~expanded

            ls_idx_lst = np.array(ls_idx_lst)[mask]
            lons_perp = flatten([lons_perp_lst[i] for i in ls_idx_lst])
            lats_perp = flatten([lats_perp_lst[i] for i in ls_idx_lst])
        return lons_perp, lats_perp, lons_perp_lst, lats_perp_lst, ls_idx_lst

    @staticmethod
    def _generate_xtrack_footprints_in_lonlat(
            P1: tuple, P2: tuple, n=None):
        """
        Input are tuples (lon, lat) for points P1, P2
        """
        # create vector
        V = np.array([P1[0]-P2[0], P1[1]-P2[1]])
        # rotate 90 degree
        theta = np.deg2rad(90)
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        Vrot = np.dot(R, V)
        # produce footprints to either side
        n = n*0.1
        P_perp_minus = (P1[0] - n*Vrot[0], P1[1] - n*Vrot[1])
        P_perp_plus = (P1[0] + n*Vrot[0], P1[1] + n*Vrot[1])
        return P_perp_minus, P_perp_plus

    @staticmethod
    def _generate_xtrack_footprints_in_cartesian():
        import utm
        return

    @staticmethod
    def _lonlat_to_xy(lon, lat, utmzone):
        P = Proj(proj='utm', zone=utmzone,
                 ellps='WGS84', preserve_units=True)
        return P(lon, lat)

    @staticmethod
    def _xy_to_lonlat(x, y, utmzone):
        P = Proj(proj='utm', zone=utmzone,
                 ellps='WGS84', preserve_units=True)
        return P(x, y, inverse=True)

    @staticmethod
    def _distance(lon1, lat1, lon2, lat2):
        G = Geod(ellps='WGS84')
        return G.inv(lon1, lat1, lon2, lat2)[2]


def start_stop(a, trigger_val):
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False, np.equal(a, trigger_val), False]
    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])
    # Get the start and end indices with slicing along the shifting ones
    return zip(idx[::2], idx[1::2]-1)

def apply_land_mask(longitudes: np.ndarray, latitudes: np.ndarray):
    """ Mask out parts covering land

    Args:
        longitudes, latitudes

    Returns:
        vardict, sea_mask

    """
    global ROAR

    if ROAR is None:
        ROAR = roaring_landmask.RoaringLandmask.new()

    # ensure float64 type on input to ROAR
    longitudes = longitudes.astype(np.float64)
    latitudes = latitudes.astype(np.float64)

    land_mask = ROAR.contains_many(longitudes, latitudes)
    sea_mask = np.invert(land_mask)

    return sea_mask

def lanczos_weights(window, cutoff):
    """ Calculate weights for a low pass Lanczos filter

    args:
        window: (integer) the length of the filter window
        cutoff: (float) the cutoff frequency in inverse time steps

    returns: weights

    example: https://scitools.org.uk/iris/docs/v1.2/examples/
             graphics/SOI_filtering.html
    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]

def smoother_GP(x, y, X, **kwargs):
    from sklearn import gaussian_process
    from sklearn.gaussian_process.kernels import RBF
    from sklearn.gaussian_process.kernels import WhiteKernel
    from sklearn.gaussian_process.kernels import RationalQuadratic
    if isinstance(x, list):
        x = np.array(x)
    if isinstance(y, list):
        y = np.array(y)
    if isinstance(X, list):
        X = np.array(X)
    if X is None:
        X = x.reshape(-1, 1)
    else:
        X = X.reshape(-1, 1)
    x = x.reshape(-1, 1)
    # create a zero mean process
    Y = y.reshape(-1, 1) - np.nanmean(y)
    # define the kernel based on kwargs
    if 'kernel' in kwargs.keys():
        print('kernel is defined by user')
        kernel = kwargs['kernel']
    elif 'kernel_lst' in kwargs.keys():
        print('kernel constituents given by user')
        kernel = WhiteKernel(noise_level=1)
        if 'RBF' in kwargs['kernel_lst']:
            kernel += 1 * RBF(length_scale=1)
        if 'RationalQuadratic' in kwargs['kernel_lst']:
            kernel += 1 * RationalQuadratic(alpha=1, length_scale=1)
    else:
        print('default kernel')
        kernel = WhiteKernel(noise_level=1,
                             noise_level_bounds=(.5, 1000))\
                        + 1 * RBF(length_scale=1,
                                  length_scale_bounds=(2, 10.0))
    gp = gaussian_process.GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=10)
    gp.fit(x, Y)
    print(gp.kernel_)
    y_pred, _ = gp.predict(X, return_std=True)
    y_pred = y_pred + np.nanmean(y)
    return y_pred

def smoother_linearGAM(x, y, X, **kwargs):
    from pygam import LinearGAM, l, s
    if isinstance(x, list):
        x = np.array(x)
    x = x.reshape(len(x), 1)
    if isinstance(y, list):
        y = np.array(y)
    if isinstance(X, list):
        X = np.array(X)
    if X is None:
        X = x.reshape(len(x), 1)
    else:
        X = X.reshape(len(X), 1)
    #if 'n_splines' in kwargs.keys():
    #    n_splines = kwargs['n_splines']
    #else:
    #    # This is because the automatic approach is too smooth
    #    n_splines = int(len(y)/5)
    #gam = LinearGAM(n_splines=n_splines,\
    #                terms=s(0,basis='ps')\
    #                ).gridsearch(x, y)
    gam = LinearGAM(terms=s(0, basis='ps')).gridsearch(x, y)
    gam.summary()
    # sample on the input grid
    means = gam.predict(X)
    return means

def cleaner_GP(x, y, **kwargs):
    from sklearn import gaussian_process
    from sklearn.gaussian_process.kernels import RBF
    from sklearn.gaussian_process.kernels import WhiteKernel
    from sklearn.gaussian_process.kernels import RationalQuadratic
    if isinstance(x, list):
        x = np.array(x)
    if isinstance(y, list):
        y = np.array(y)
    X = x.reshape(-1, 1)
    # create a zero mean process
    ymean = np.nanmean(y)
    Y = y.reshape(-1, 1) - ymean
    # define the kernel based on kwargs
    if 'kernel' in kwargs.keys():
        print('kernel is defined by user')
        kernel = kwargs['kernel']
    elif 'kernel_lst' in kwargs.keys():
        print('kernel constituents given by user')
        kernel = WhiteKernel(noise_level=1)
        if 'RBF' in kwargs['kernel_lst']:
            kernel += RBF(length_scale=1)
        if 'RationalQuadratic' in kwargs['kernel_lst']:
            kernel += RationalQuadratic(alpha=1,\
                                        length_scale=1)
    else:
        print('default kernel')
        kernel = WhiteKernel(noise_level=1,
                             noise_level_bounds=(.5, 1000))\
                    + 1 * RBF(length_scale=1,
                              length_scale_bounds=(2, 10.0))
    gp = gaussian_process.GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=10)
    gp.fit(X, Y)
    print(gp.kernel_)
    y_pred, sigma = gp.predict(X, return_std=True)
    sigma_multiplyer = kwargs.get('sigma', 2)
    uplim = y_pred + (sigma_multiplyer*sigma)
    lowlim = y_pred - (sigma_multiplyer*sigma)
    idx = [i for i in range(len(Y))
           if (Y[i] < uplim[i] and Y[i] > lowlim[i])]
    return idx

def cleaner_linearGAM(x, y, **kwargs):
    from pygam import LinearGAM, l, s
    if isinstance(x, list):
        x = np.array(x)
    if isinstance(y, list):
        y = np.array(y)
    X = x.reshape(len(x), 1)
    #if 'n_splines' in kwargs.keys():
    #    n_splines = kwargs['n_splines']
    #else:
    #    # This is because the automatic approach is too smooth
    #    #n_splines = int(len(y)/5)
    #gam = LinearGAM(n_splines=n_splines,\
    #                terms=s(0,basis='ps')\
    #                ).gridsearch(X, y)
    gam = LinearGAM(terms=s(0, basis='ps')).gridsearch(X, y)
    gam.summary()
    #gam = LinearGAM(n_splines=n_splines,terms=s(0)).gridsearch(X, y)
    # sample on the input grid
    means = gam.predict(X)
    bounds = gam.prediction_intervals(X, width=kwargs.get('pct', .95))
    idx = [i for i in range(len(y)) \
           if (y[i] < bounds[i, 1] and y[i] > bounds[i, 0])]
    return idx
