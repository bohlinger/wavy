"""
    Module for quicklook fct
"""
# imports
import numpy as np
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from abc import abstractmethod
import matplotlib.pyplot as plt

# own imports
from wavy.wconfig import load_or_default
from wavy.utils import parse_date
from wavy.utils import compute_quantiles
#from wavy.collocation_module import collocation_class as cc
#from wavy.insitu_module import poi_class as pc
#from wavy.insitu_module import insitu_class as ic
#from wavy.satellite_module import satellite_class as sc
#from wavy.model_module import insitu_class as mc
#from wavy.gridder_module import gridder_class as gc

# read yaml config files:
region_dict = load_or_default('region_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')
model_dict = load_or_default('model_cfg.yaml')
quicklook_dict = load_or_default('quicklook_cfg.yaml')


class quicklook_class_sat:

    def quicklook(self, a=False, projection=None, **kwargs):
        """
        Enables to explore the class object (and retrieved results)
        by plotting time series and map.

        param:
            m - map figure (True/False)
            ts - time series (True/False)
            a - all figures (True/False)
            projection - specified projection for cartopy

        return:
            figures
        """

        import matplotlib as mpl
        import cmocean

        # settings
        m = kwargs.get('m', a)
        ts = kwargs.get('ts', a)
        sc = kwargs.get('sc', False)
        hst = kwargs.get('hist', False)
        mode = kwargs.get('mode', 'comb')  # comb, indiv

        # set variables
        try:
            plot_var = self.vars[self.varalias]
            plot_lons = self.vars.lons
            plot_lats = self.vars.lats
        except Exception as e:
            plot_var = self.vars.obs_values
            plot_lons = self.vars.obs_lons
            plot_lats = self.vars.obs_lats
            plot_var_obs = self.vars.obs_values
            plot_lons_obs = self.vars.obs_lons
            plot_lats_obs = self.vars.obs_lats
            plot_var_model = self.vars.model_values
            plot_lons_model = self.vars.model_lons
            plot_lats_model = self.vars.model_lats

        fs = kwargs.get('fs', 12)

        vmin = kwargs.get('vmin', 0)
        vmax = kwargs.get('vmax', np.nanmax(plot_var))

        levels = kwargs.get('levels',
                            np.arange(0, vmax+.5, .5))
                            #np.arange(vmin, vmax, .5))

        cflevels = kwargs.get('cflevels', levels)
        clevels = kwargs.get('clevels', levels)

        # parse kwargs
        vartype = variable_info[self.varalias].get('type', 'default')
        if kwargs.get('cmap') is None:
            if vartype == 'cyclic':
                cmap = mpl.cm.twilight
            else:
                cmap = cmocean.cm.amp
        else:
            cmap = kwargs.get('cmap')

        norm = mpl.colors.BoundaryNorm(levels, cmap.N)

        if m is True:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            import matplotlib.cm as mplcm
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes

            zorder_land = kwargs.get('zorder_land', 10)

            # land
            land = cfeature.GSHHSFeature(
                    scale=kwargs.get('land_mask_resolution', 'i'),
                    levels=[1],
                    facecolor=cfeature.COLORS['land'])
            if projection is None:
                projection = ccrs.PlateCarree()

            lonmax, lonmin = np.max(plot_lons), np.min(plot_lons)
            latmax, latmin = np.max(plot_lats), np.min(plot_lats)

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection=projection)

            # add land
            ax.add_geometries(land.intersecting_geometries(
                    [-180, 180, 0, 90]),
                    ccrs.PlateCarree(),
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black', linewidth=1,
                    zorder=zorder_land)

            # - add land color
            ax.add_feature(land, facecolor='burlywood', alpha=0.5)

            # plot track if applicable
            lonmax, lonmin = np.max(plot_lons), np.min(plot_lons)
            latmax, latmin = np.max(plot_lats), np.min(plot_lats)
            if kwargs.get('poi') is not None:
                poi = kwargs.get('poi')
                plats = poi.vars.lats.data
                platsmax, platsmin = np.max(plats), np.min(plats)
                plons = poi.vars.lons.data
                plonsmax, plonsmin = np.max(plons), np.min(plons)
                tc = ax.plot(plons, plats, color='cornflowerblue',
                             ls='-', lw=1,
                             zorder=-1)
                tc = ax.plot(plons, plats, color='cornflowerblue',
                             ls='None', marker='o', ms=5,
                             markeredgecolor='k',
                             zorder=-1)
                lonmax, lonmin = np.max([lonmax, plonsmax]),\
                                 np.min([lonmin, plonsmin])
                latmax, latmin = np.max([latmax, platsmax]),\
                                 np.min([latmin, platsmin])
            # plot sat
            if kwargs.get("plot_xtrack_pulse_limited_fpr") is not None:
                domain = kwargs.get('domain', 'lonlat')
                number_of_seeds = kwargs.get('number_of_seeds', 100)
                lons_perp, lats_perp, _, _, ls_idx_lst = \
                    self._generate_xtrack_footprints(
                            domain=domain,
                            number_of_seeds=number_of_seeds)
                sc2 = ax.scatter(lons_perp, lats_perp,
                                 s=.2, c='b', marker='.',
                                 edgecolor='face',
                                 transform=ccrs.PlateCarree())
            if len(plot_var.shape) > 1:
                sc = ax.contourf(plot_lons.squeeze(),
                                 plot_lats.squeeze(),
                                 plot_var.squeeze(),
                                 cmap=cmap, levels=cflevels,
                                 vmin=vmin, vmax=vmax, norm=norm,
                                 transform=ccrs.PlateCarree())
                c = ax.contour(plot_lons.squeeze(),
                               plot_lats.squeeze(),
                               plot_var.squeeze(),
                               levels=clevels,
                               colors='w', linewidths=0.3,
                               transform=ccrs.PlateCarree())
            else:
                sc = ax.scatter(plot_lons, plot_lats, s=15,
                            c=plot_var,
                            marker='o',  # edgecolor='face',
                            edgecolors='k',
                            linewidths=0.1,
                            cmap=cmap, norm=norm,
                            transform=ccrs.PlateCarree())
            # axes for colorbar
            axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )

            # - colorbar
            cbar = fig.colorbar(sc, cax=axins,
                    label=self.varalias + ' [' + self.units + ']',
                    ticks=levels)
            cbar.ax.set_ylabel(self.units, size=fs)
            cbar.ax.tick_params(labelsize=fs)

            lon_range = (lonmax - lonmin)
            lat_range = (latmax - latmin)
            map_extent_multiplicator = kwargs.get(
                "map_extent_multiplicator", 0.1)
            map_extent_multiplicator_lon = kwargs.get(
                "map_extent_multiplicator_lon", map_extent_multiplicator)
            map_extent_multiplicator_lat = kwargs.get(
                "map_extent_multiplicator_lat", map_extent_multiplicator)
            ax.set_extent([lonmin-lon_range*map_extent_multiplicator_lon,
                           lonmax+lon_range*map_extent_multiplicator_lon,
                           latmin-lat_range*map_extent_multiplicator_lat,
                           latmax+lat_range*map_extent_multiplicator_lat],
                           crs=projection)

            #ax.coastlines(color='k')

            gl = ax.gridlines(draw_labels=True, crs=projection,
                              linewidth=1, color='grey', alpha=0.4,
                              linestyle='-')
            gl.top_labels = False
            gl.right_labels = False
            plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
            auto_title = (self.nID + '\n'
                          + 'from ' 
                          + (parse_date(str(self.vars['time'][0].values))).\
                              strftime('%Y-%m-%d %H:%M:%S')
                          + ' to '
                          + (parse_date(str(self.vars['time'][-1].values))).\
                              strftime('%Y-%m-%d %H:%M:%S'))
            title = kwargs.get('title',auto_title)
            ax.set_title(title)
            # plot from quickloop config file
            if ('region' in vars(self).keys()
                and self.region in quicklook_dict
                and 'poi' in quicklook_dict[self.region]):
                for poi in quicklook_dict[self.region]['poi']:
                    pname = quicklook_dict[self.region]['poi'][poi]['name']
                    plat = quicklook_dict[self.region]['poi'][poi]['lat']
                    plon = quicklook_dict[self.region]['poi'][poi]['lon']
                    scp = ax.scatter(plon, plat, s = 20,
                                     c=quicklook_dict\
                                           [self.region]['poi'][poi]\
                                           .get('color','b'),
                                     marker=quicklook_dict[self.region]\
                                            ['poi'][poi]['marker'],
                                     transform=ccrs.PlateCarree())
                ax.text(plon, plat, pname, transform=ccrs.PlateCarree())
            #fig.suptitle('', fontsize=16) # unused
            plt.show()

        if (ts is True and mode == 'comb'):
            fig = plt.figure(figsize=(9, 3.5))
            ax = fig.add_subplot(111)
            colors = ['k', 'r']
            ax.plot(self.vars['time'],
                    plot_var,
                    color=colors[0],
                    linestyle=kwargs.get('linestyle',''),
                    label=self.nID,
                    marker='o',alpha=.5,ms=2)
            try:
                if 'model' in vars(self):
                    label_scdplot = self.model
                else:
                    label_scdplot = self.nID
                ax.plot(self.vars['time'],
                        plot_var_model,
                        color=colors[1],
                        linestyle=kwargs.get('linestyle',''),
                        label=label_scdplot,
                        marker='o',alpha=.5,ms=2)
            except Exception as e:
                pass
            plt.ylabel(self.varalias + ' [' + self.units + ']')
            plt.legend(loc='best')
            plt.tight_layout()
            #ax.set_title()
            plt.show()

        elif (ts is True and mode == 'indiv'):
            fig = plt.figure(figsize=(9, 3.5))
            ax = fig.add_subplot(111)
            for oco in self.ocos:
                label = oco.name
                ax.plot(oco.vars['time'],
                        oco.vars[oco.varalias],
                        linestyle=kwargs.get('linestyle',''),
                        label=label,
                        marker='o',alpha=.5, ms=2)
            try:
                label = self.model
                ax.plot(self.vars['time'],
                        plot_var_model,
                        color=colors[1],
                        linestyle=kwargs.get('linestyle',''),
                        label=label,
                        marker='o',alpha=.5,ms=2)
            except Exception as e:
                pass
            plt.ylabel(self.varalias + ' [' + self.units + ']')
            plt.legend(loc='best')
            plt.tight_layout()
            #ax.set_title()
            plt.show()

        if sc is True:
            lq = np.arange(0.01, 1.01, 0.01)
            lq = kwargs.get('lq', lq)
            modq = compute_quantiles(plot_var_model, lq)
            obsq = compute_quantiles(plot_var_obs, lq)

            fig = plt.figure(figsize=(4, 4))
            ax = fig.add_subplot(111)
            colors = ['k']

            ax.plot(plot_var_obs, plot_var_model,
                    linestyle='None', color=colors[0],
                    marker='o', alpha=.5, ms=2)

            # add quantiles
            ax.plot(obsq, modq, 'r')

            # 45 degree line for orientation
            ax.axline((0, 0), (1, 1), lw=.5, color='grey',ls='--')

            # add axis labels
            plt.xlabel('obs (' + self.nID + ')')
            plt.ylabel('models (' + self.model + ')')

            vartype = variable_info[self.varalias].get('type','default')
            if vartype == 'cyclic':
                plt.xlim([0, 360])
                plt.ylim([0, 360])
            else:
                maxv = np.nanmax([self.vars['model_values'],
                                  self.vars['obs_values']])
                minv = 0
                plt.xlim([minv, maxv+0.15*maxv])
                plt.ylim([minv, maxv+0.15*maxv])
            ax.set_title(self.varalias + '[' + self.units + ']')
            plt.tight_layout()
            #ax.set_title()
            plt.show()

    def quick_anim():
        pass


"""
from abc import abstractmethod


class Quicklook:
    @property
    @abstractmethod
    def projection():
        pass

    @property
    @abstractmethod
    def vars():
        pass

    def quicklook(self, blah):
        proj = self.projection

        for var in self.vars:
            # do something
            pass

class Sat(Quicklook):
    vars = [ "hs", "foo" ]
    projection = pyproj.Proj('+proj=latlong')
"""
