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
from wavy.validationmod import linreg_evm, linreg_std

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
            ms - map figure + scatter (True/False)
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
        scat = kwargs.get('sc', False)
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
            plot_var_model = self.vars.model_values

        if str(type(self)) == "<class 'wavy.model_module.model_class'>":
            if len(plot_lons.shape) < 2:
                plot_lons, plot_lats = np.meshgrid(plot_lons, plot_lats)

        fs = kwargs.get('fs', 12)

        vartype = variable_info[self.varalias].get('type', 'default')
        if kwargs.get('cmap') is None:
            if vartype == 'cyclic':
                cmap = mpl.cm.twilight
            else:
                cmap = cmocean.cm.amp
        else:
            cmap = kwargs.get('cmap')


        if m is True:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            import matplotlib.cm as mplcm
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes

            vmin = kwargs.get('vmin', 0)
            vmax = kwargs.get('vmax', np.nanmax(plot_var))

            levels_incr = kwargs.get('levels_incr', .5)

            levels = kwargs.get('levels',
                                np.arange(vmin, vmax, levels_incr))

            cflevels = kwargs.get('cflevels', levels)
            clevels = kwargs.get('clevels', levels)

            norm = mpl.colors.BoundaryNorm(levels, cmap.N)

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
            #ax = fig.add_subplot(1, 1, 1, projection=ccrs.epsg(3857))

            # add land
            ax.add_geometries(land.intersecting_geometries(
                    [-180, 180, 0, 90]),
                    projection,
                    facecolor=cfeature.COLORS['land'],
                    edgecolor='black', linewidth=1,
                    zorder=zorder_land)

            # add sea map
            # ax.add_wmts("https://cache.kartverket.no/v1/wmts", 'sjokartraster')

            # - add land color
            ax.add_feature(land, facecolor='burlywood', alpha=0.5)

            # plot track if applicable
            if kwargs.get('lonmax') is not None:
                lonmax = kwargs.get('lonmax')
            else:
                lonmax = np.max(plot_lons)
            if kwargs.get('latmax') is not None:
                latmax = kwargs.get('latmax')
            else:
                latmax = np.max(plot_lats)
            if kwargs.get('lonmin') is not None:
                lonmin = kwargs.get('lonmin')
            else:
                lonmin = np.min(plot_lons)
            if kwargs.get('latmin') is not None:
                latmin = kwargs.get('latmin')
            else:
                latmin = np.min(plot_lats)

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
                                 transform=projection)
            if len(plot_var.shape) > 1:
                sc = ax.contourf(plot_lons.squeeze(),
                                 plot_lats.squeeze(),
                                 plot_var.squeeze(),
                                 cmap=cmap, levels=cflevels,
                                 vmin=vmin, vmax=vmax, norm=norm,
                                 transform=projection,
                                 transform_first=\
                                 kwargs.get('transform_first', False))
                c = ax.contour(plot_lons.squeeze(),
                               plot_lats.squeeze(),
                               plot_var.squeeze(),
                               levels=clevels,
                               colors='w', linewidths=0.3,
                               transform=projection,
                               transform_first=\
                               kwargs.get('transform_first', False))
            else:
                sc = ax.scatter(plot_lons, plot_lats, s=15,
                            c=plot_var,
                            marker='o',  # edgecolor='face',
                            edgecolors='k',
                            linewidths=0.1,
                            cmap=cmap, norm=norm,
                            transform=projection)

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
            if kwargs.get("cbar", True) is True:
                cbar = fig.colorbar(sc, cax=axins,
                        label=self.varalias + ' [' + self.units + ']',
                        ticks=levels)
                cbar.ax.set_ylabel(self.units, size=fs)
                cbar.ax.tick_params(labelsize=fs)

            # - add extend
            if kwargs.get("map_extent_llon") is None:
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
            elif kwargs.get('map_extend') is False:
                pass
            else:
                ax.set_extent([kwargs.get("map_extent_llon"),
                               kwargs.get("map_extent_ulon"),
                               kwargs.get("map_extent_llat"),
                               kwargs.get("map_extent_ulat")],
                              crs=projection)

            #ax.coastlines(color='k')
            if projection == ccrs.PlateCarree():
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
            title = kwargs.get('title', auto_title)
            ax.set_title(title)
            # plot from quickloop config file
            if ('region' in vars(self).keys()
                and self.region in quicklook_dict
                and 'poi' in quicklook_dict[self.region]):
                for poi in quicklook_dict[self.region]['poi']:
                    pname = quicklook_dict[self.region]['poi'][poi]['name']
                    plat = quicklook_dict[self.region]['poi'][poi]['lat']
                    plon = quicklook_dict[self.region]['poi'][poi]['lon']
                    scp = ax.scatter(plon, plat, s=20,
                                     c=quicklook_dict\
                                           [self.region]['poi'][poi]\
                                           .get('color','b'),
                                     marker=quicklook_dict[self.region]\
                                            ['poi'][poi]['marker'],
                                     transform=projection)
                ax.text(plon, plat, pname, transform=projection,
                        zorder=100)
            #fig.suptitle('', fontsize=16) # unused
            if kwargs.get("show", True) is True:
                plt.show()
            else:
                return fig, ax

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
            if kwargs.get("show", True) is True:
                plt.show()
            else:
                return fig, ax


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
            if kwargs.get("show", True) is True:
                plt.show()
            else:
                return fig, ax

        if scat is True:
            lq = np.arange(0.01, 1.01, 0.01)
            lq = kwargs.get('lq', lq)
            modq = compute_quantiles(plot_var_model, lq)
            obsq = compute_quantiles(plot_var_obs, lq)

            fig = plt.figure(figsize=(4, 4))
            ax = fig.add_subplot(111)
            colors = ['k']

            ax.plot(plot_var_obs, plot_var_model,
                    linestyle='None', color=colors[0],
                    marker='o', alpha=.5, ms=2,
                    label="data")

            # add quantiles
            ax.plot(obsq, modq, 'r', label='QQ')

            # 45 degree line for orientation
            ax.axline((0, 0), (1, 1), lw=.5, color='grey',
                      ls='--', label="45 deg")

            # linreg_evm line
            if kwargs.get('evm_regression_line') is True:
                rl = linreg_evm(plot_var_obs, plot_var_model, **kwargs)
                self.EVMreg = dict({'intercept': rl[1], "slope": rl[0]})
                ax.axline(xy1=(0, rl[1]), slope=rl[0],
                          color=kwargs.get('evm_regression_col', 'lightblue'),
                          lw=kwargs.get('evm_regression_lw', 1),
                          ls=kwargs.get('evm_regression_ls', '-'),
                          label="EVM-regr")

            # std linreg line
            if kwargs.get('std_regression_line') is True:
                rl = linreg_std(plot_var_obs, plot_var_model, **kwargs)
                self.linreg = dict({'intercept': rl['intercept'],
                                    'slope': rl['slope']})
                ax.axline(xy1=(0, rl['intercept']), slope=rl['slope'],
                          color=kwargs.get('std_regression_col', 'lightblue'),
                          lw=kwargs.get('std_regression_lw', 1),
                          ls=kwargs.get('std_regression_ls', '-'),
                          label="linregr")

            # add axis labels
            plt.xlabel('obs (' + self.nID + ')')
            plt.ylabel('models (' + self.model + ')')

            maxv = np.nanmax([self.vars['model_values'],
                              self.vars['obs_values']])
            minv = 0
            plt.xlim([minv, maxv*1.05])
            plt.ylim([minv, maxv*1.05])

            ax.set_title(self.varalias + '[' + self.units + ']')
            plt.legend()

            plt.tight_layout()

            #ax.set_title()
            if kwargs.get("show", True) is True:
                plt.show()
            else:
                return fig, ax

        if hst is True:
            lq = np.arange(0.01, 1.01, 0.01)
            lq = kwargs.get('lq', lq)
            modq = compute_quantiles(plot_var_model, lq)
            obsq = compute_quantiles(plot_var_obs, lq)

            fig = plt.figure(figsize=(5, 4))
            ax = fig.add_subplot(111)

            # 2d histogram
            lmin = 0
            lmax = np.nanmax([plot_var_obs, plot_var_model])*1.05
            plt.hist2d(plot_var_obs, plot_var_model,
                       bins=kwargs.get('bins', 100),
                       range=[[lmin, lmax], [lmin, lmax]],
                       norm=kwargs.get('norm', mpl.colors.LogNorm()),
                       cmap=kwargs.get('cmap', mpl.cm.gray),
                       cmin=kwargs.get('cmin', 1))
            cbar = plt.colorbar()
            cbar.ax.tick_params(labelsize=fs)
            cbar.set_label('Frequency', size=fs)
            plt.xlim([lmin, lmax])
            plt.ylim([lmin, lmax])

            # add quantiles
            ax.plot(obsq, modq, 'r', label='QQ')

            # 45 degree line for orientation
            ax.axline((0, 0), (1, 1), lw=.5, color='grey',
                      ls='--', label="45 deg")

            # linreg_evm line
            if kwargs.get('evm_regression_line') is True:
                rl = linreg_evm(plot_var_obs, plot_var_model, **kwargs)
                self.EVMreg = dict({'intercept': rl[1], "slope": rl[0]})
                ax.axline(xy1=(0, rl[1]), slope=rl[0],
                          color=kwargs.get('evm_regression_col', 'lightblue'),
                          lw=kwargs.get('evm_regression_lw', 1),
                          ls=kwargs.get('evm_regression_ls', '-'),
                          label="EVM-regr")

            # std linreg line
            if kwargs.get('std_regression_line') is True:
                rl = linreg_std(plot_var_obs, plot_var_model, **kwargs)
                self.linreg = dict({'intercept': rl['intercept'],
                                    'slope': rl['slope']})
                ax.axline(xy1=(0, rl['intercept']), slope=rl['slope'],
                          color=kwargs.get('std_regression_col', 'lightblue'),
                          lw=kwargs.get('std_regression_lw', 1),
                          ls=kwargs.get('std_regression_ls', '-'),
                          label="linregr")

            # add axis labels
            plt.xlabel('obs (' + self.nID + ')')
            plt.ylabel('models (' + self.model + ')')

            ax.set_title(self.varalias + '[' + self.units + ']')
            plt.legend()

            plt.tight_layout()

            if kwargs.get("show", True) is True:
                plt.show()
            else:
                return fig, ax


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
