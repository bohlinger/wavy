# Module to organize gridding data

# imports
import numpy as np
import tqdm
from wavy.grid_stats import apply_metric
from wavy.wconfig import load_or_default

validation_metric_abbreviations = load_or_default('validation_metrics.yaml')

class gridder_class():

    def __init__(
    self, oco=None, mco=None, cco=None, bb=None, grid='lonlat', res=(1, 1),
    **kwargs):
        """
        setup the gridder
        grid: lonlat or in m
        bb tupel: (lonmin, lonmax, latmin, latmax)
        res: resolution, tupel e.g. (.5,.5)
             in degree where tupel is (lon,lat) dimension
        """
        print('# ----- ')
        print(" ### Initializing gridder_class object ###")
        print(" ")
        self.mvals = None
        if oco is not None:
            self.olons = np.array(oco.vars['lons'].squeeze().values.ravel())
            self.olats = np.array(oco.vars['lats'].squeeze().values.ravel())
            self.ovals = np.array(oco.vars[oco.varalias].squeeze().values.ravel())
            self.stdvarname = oco.stdvarname
            self.varalias = oco.varalias
            self.units = oco.units
            self.sdate = oco.vars['time'][0]
            self.edate = oco.vars['time'][-1]
        elif cco is not None:
            self.olons = np.array(cco.vars['obs_lons'])
            self.olats = np.array(cco.vars['obs_lats'])
            self.ovals = np.array(cco.vars['obs_values'])
            self.mvals = np.array(cco.vars['model_values'])
            self.stdvarname = cco.stdvarname
            self.varalias = cco.varalias
            self.units = cco.units
            self.sdate = cco.vars['time'][0]
            self.edate = cco.vars['time'][-1]
        elif mco is not None:
            self.olons = np.array(mco.vars.lons.squeeze().values.flatten())
            self.olats = np.array(mco.vars.lats.squeeze().values.flatten())
            self.ovals = np.array(
                    mco.vars[mco.varalias].squeeze().values.flatten())
            self.stdvarname = mco.stdvarname
            self.varalias = mco.varalias
            self.units = mco.units
            self.sdate = mco.vars['time'][0]
            self.edate = mco.vars['time'][-1]
        else:
            self.olons = kwargs.get('lons')
            self.olats = kwargs.get('lats')
            self.ovals = kwargs.get('values')
            self.stdvarname = kwargs.get('stdvarname', None)
            self.varalias = kwargs.get('varalias', None)
            self.units = kwargs.get('units', None)
            self.sdate = kwargs.get('sdate', None)
            self.edate = kwargs.get('edate', None)

        self.bb = bb
        self.res = res
        self.grid = grid
        self.glons, self.glats = self.create_grid_coords()
        ovals, mvals, Midx = self.get_obs_grid_idx()
        self.ovals_clean = ovals
        self.mvals_clean = mvals
        self.Midx_clean = Midx
        print(" ")
        print(" ### gridder_class object initialized ###")
        print('# ----- ')

    def create_grid_coords(self):
        """
        returns grid coordinates
        """
        lons = np.arange(self.bb[0], self.bb[1]+self.res[0]/2, self.res[0])
        lats = np.arange(self.bb[2], self.bb[3]+self.res[1]/2, self.res[1])
        return np.array(lons), np.array(lats)

    def get_obs_grid_idx(self):
        Midx = self.assign_obs_to_grid(
                self.glons, self.glats,
                self.olons, self.olats,
                self.res)
        ovals, mvals, Midx = self.clean_Midx(
                Midx, self.ovals, self.mvals, self.glons, self.glats)
        return ovals, mvals, Midx

    @staticmethod
    def assign_obs_to_grid(glons, glats, olons, olats, res):
        """
        assigns observation coordinates to grid indices
        """
        lonidx = ((olons-np.min(glons))/res[0]).astype(int)
        latidx = ((olats-np.min(glats))/res[1]).astype(int)
        Midx = np.array([lonidx, latidx], dtype=object)
        return Midx

    @staticmethod
    def clean_Midx(Midx, ovals, mvals, glons, glats):
        """
        cleans Midx and observations (ovals) for grid cells outside bb
        """
        # clean x-dim (tlons)
        glons_idx = np.where((Midx[0] >= 0) & (Midx[0] < len(glons)))[0]
        Midx = Midx[:, glons_idx]
        ovals = ovals[glons_idx]
        if mvals is not None:
            mvals = mvals[glons_idx]
        # clean y-dim (tlats)
        glats_idx = np.where((Midx[1] >= 0) & (Midx[1] < len(glats)))[0]
        Midx = Midx[:, glats_idx]
        ovals = ovals[glats_idx]
        if mvals is not None:
            mvals = mvals[glats_idx]
        return ovals, mvals, Midx

    @staticmethod
    def region_filter():
        # filter grid cells for region of interest
        # region could be rectangular, polygon, ...
        # -> return grid
        return

    @staticmethod
    def get_grid_idx(Midx):
        return np.where((Midx[0] == Midx[0][0]) & (Midx[1] == Midx[1][0]))[0]

    @staticmethod
    def calc_mean(gidx, ovals):
        return np.mean(np.array(ovals)[gidx])

    @staticmethod
    def rm_used_idx_from_Midx(gidx, Midx):
        d1 = np.delete(Midx[0, :], gidx)
        d2 = np.delete(Midx[1, :], gidx)
        return np.array([d1, d2], dtype=object)

    @staticmethod
    def get_exteriors(glons, glats, res):
        xb = np.array([glons+res[0], glons+res[0],
                       glons, glons,
                       glons+res[0]])
        yb = np.array([glats,
                       glats+res[1], glats+res[1],
                       glats, glats])
        return xb, yb

    def grid_view(self, metric, mask_metric_llim, mask_metric, **kwargs):
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        import cmocean
        import matplotlib.pyplot as plt
        import matplotlib.cm as mplcm
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        import math
        from copy import deepcopy

        # shift coords for plotting
        lon_grid = kwargs.get('lon_grid') + self.res[0]/2.
        lat_grid = kwargs.get('lat_grid') + self.res[1]/2.

        # backup values
        all_grid = deepcopy(kwargs.get('val_grid'))
        mask_grid = all_grid[mask_metric]
        val_grid = all_grid[metric]

        # apply mask
        mask_llim_idx = np.where(mask_grid < mask_metric_llim)
        val_grid[mask_llim_idx[0], mask_llim_idx[1]] = np.nan

        if kwargs.get('projection') is None:
            projection = ccrs.PlateCarree()
        # parse kwargs
        if kwargs.get('cmap') is None:
            cmap = cmocean.cm.amp
        else:
            cmap = kwargs.get('cmap')

        # max/min for colorbar
        vmax = kwargs.get('vmax')
        vmin = kwargs.get('vmin')

        # plot track if applicable
        if kwargs.get('lonmax') is not None:
            lonmax = kwargs.get('lonmax')
        else:
            lonmax = np.max(lon_grid)
        if kwargs.get('latmax') is not None:
            latmax = kwargs.get('latmax')
        else:
            latmax = np.max(lat_grid)
        if kwargs.get('lonmin') is not None:
            lonmin = kwargs.get('lonmin')
        else:
            lonmin = np.min(lon_grid)
        if kwargs.get('latmin') is not None:
            latmin = kwargs.get('latmin')
        else:
            latmin = np.min(lat_grid)

        # land
        land = cfeature.GSHHSFeature(
                    scale=kwargs.get('land_mask_resolution', 'i'),
                    levels=[1],
                    facecolor=cfeature.COLORS['land'])

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=projection)
        # add land
        ax.add_geometries(land.intersecting_geometries(
                    [-180, 180, 0, 90]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black', linewidth=1)

        ax.set_extent([lonmin, lonmax, latmin, latmax], crs=projection)
        pc = ax.pcolormesh(
                lon_grid, lat_grid, val_grid,
                transform=projection, cmap=cmap,
                vmax=vmax, vmin=vmin)

        axins = inset_axes(ax,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.01, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )

        metric_name = validation_metric_abbreviations[metric].get('name')
        metric_units =\
            validation_metric_abbreviations[metric].get('units', self.units)
        if metric_units is None:
            cbar = fig.colorbar(pc, cax=axins, label=metric_name)
        else:
            cbar = fig.colorbar(pc, cax=axins,
                                label=metric_name
                                + ' [' + metric_units + ']')

        # ax.coastlines()
        gl = ax.gridlines(draw_labels=True, crs=projection,
                          linewidth=1, color='grey', alpha=0.4,
                          linestyle='-')
        gl.top_labels = False
        gl.right_labels = False
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        autotitle = ('Base variable: ' + self.varalias + '\n'
                     + 'from ' + str(self.sdate.data)
                     + ' to ' + str(self.edate.data))
        if kwargs.get('title') is None:
            ax.set_title(autotitle)
        else:
            ax.set_title(kwargs.get('title'))
        ax.title.set_size(11)
        # todo: add info on observation and model source for figure
        plt.show()

    def quicklook(self, metric='mor',
                  mask_metric_llim=10,
                  mask_metric='nov',
                  **kwargs):

        if metric == 'all':
            for key in kwargs['val_grid'].keys():
                self.grid_view(key, mask_metric_llim, mask_metric, **kwargs)
        else:
            self.grid_view(metric, mask_metric_llim, mask_metric, **kwargs)
