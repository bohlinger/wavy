"""
Module for quicklook fct
"""

# imports
import numpy as np
from dataclasses import dataclass
from typing import Any, Optional, Tuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import cmocean
import logging

logger = logging.getLogger(__name__)

# own imports
from wavy.wconfig import load_or_default
from wavy.utils import parse_date
from wavy.utils import compute_quantiles
from wavy.validationmod import linreg_evm, linreg_std

# read yaml config files:
region_dict = load_or_default("region_cfg.yaml")
variable_info = load_or_default("variable_def.yaml")
model_dict = load_or_default("model_cfg.yaml")
quicklook_dict = load_or_default("quicklook_cfg.yaml")


@dataclass
class PlotContext:
    """Container for all data and parameters shared across plot methods.

    Build one via ``quicklook_class_sat._build_plot_context()`` before
    dispatching to a ``plot_*`` method.
    """

    varalias: str
    units: str
    plot_var: Any
    plot_lons: Any
    plot_lats: Any
    plot_var_obs: Any  # None when data is not collocated
    plot_var_model: Any  # None when data is not collocated
    fs: int
    cmap: Any
    projection: Any
    mode: str = "comb"


class quicklook_class_sat:
    """Mixin class providing quicklook plotting capability.

    Adding a new plot type
    ----------------------
    1. Define a method with the standard signature::

           def plot_<name>(self, ctx: PlotContext, **kwargs) -> Tuple[Figure, Axes]:
               ...

    2. Add one entry to ``_PLOT_REGISTRY`` at class level::

           "<kwarg>": {
               "method": "plot_<name>",
               "default": False,    # True → also triggered by quicklook(a=True)
               "needs_colloc": False,  # True → skipped when obs/model are absent
           }

    That is all.  ``quicklook()`` discovers and dispatches the method
    automatically.
    """

    # Registry drives quicklook() dispatching.
    # Each entry maps the enabling kwarg to plot metadata.
    _PLOT_REGISTRY: dict = {
        "m": {"method": "plot_map", "default": True, "needs_colloc": False},
        "ts": {"method": "plot_timeseries", "default": True, "needs_colloc": False},
        "sc": {"method": "plot_scat", "default": False, "needs_colloc": True},
        "hist": {"method": "plot_hist", "default": False, "needs_colloc": True},
    }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _check_varalias(self, **kwargs):
        """Resolve the active varalias and its display units from kwargs."""
        if isinstance(self.varalias, list):
            varalias = kwargs.get("varalias", self.varalias[0])
            assert varalias in self.varalias, "varalias must be one of {}".format(
                self.varalias
            )
            assert isinstance(varalias, str), "varalias argument should be a string"
            idx_units = np.argwhere(np.array(self.varalias) == varalias)[0][0]
            units_to_plot = self.units[idx_units]
        else:
            varalias = self.varalias
            units_to_plot = self.units

        return varalias, units_to_plot

    def _build_plot_context(self, projection=None, **kwargs) -> PlotContext:
        """Build and return a :class:`PlotContext` from the current object state."""
        varalias, units = self._check_varalias(**kwargs)

        try:
            plot_var = self.vars[varalias]
            plot_lons = self.vars.lons
            plot_lats = self.vars.lats
            plot_var_obs = None
            plot_var_model = None
        except Exception:
            list_vars = list(self.vars.variables)
            assert "model_" + varalias in list_vars, (
                f"model_{varalias} is missing in the dataset; "
                "specify varalias to validate a different variable."
            )
            assert "obs_" + varalias in list_vars, (
                f"obs_{varalias} is missing in the dataset; "
                "specify varalias to validate a different variable."
            )
            plot_var = self.vars["obs_" + varalias]
            plot_lons = self.vars.obs_lons
            plot_lats = self.vars.obs_lats
            plot_var_obs = self.vars["obs_" + varalias]
            plot_var_model = self.vars["model_" + varalias]

        if str(type(self)) == "<class 'wavy.model_module.model_class'>":
            if len(plot_lons.shape) < 2:
                plot_lons, plot_lats = np.meshgrid(plot_lons, plot_lats)

        fs = kwargs.get("fs", 12)
        vartype = variable_info[varalias].get("type", "default")
        if kwargs.get("cmap") is None:
            cmap = mpl.cm.twilight if vartype == "cyclic" else cmocean.cm.amp
        else:
            cmap = kwargs.get("cmap")

        return PlotContext(
            varalias=varalias,
            units=units,
            plot_var=plot_var,
            plot_lons=plot_lons,
            plot_lats=plot_lats,
            plot_var_obs=plot_var_obs,
            plot_var_model=plot_var_model,
            fs=fs,
            cmap=cmap,
            projection=_check_projection(projection),
            mode=kwargs.get("mode", "comb"),
        )

    def _add_regression_lines(self, ax, ctx: PlotContext, **kwargs):
        """Overlay optional EVM and/or standard regression lines on *ax*."""
        if kwargs.get("evm_regression_line") is True:
            rl = linreg_evm(ctx.plot_var_obs, ctx.plot_var_model, **kwargs)
            self.EVMreg = {"intercept": rl[1], "slope": rl[0]}
            ax.axline(
                xy1=(0, rl[1]),
                slope=rl[0],
                color=kwargs.get("evm_regression_col", "lightblue"),
                lw=kwargs.get("evm_regression_lw", 1),
                ls=kwargs.get("evm_regression_ls", "-"),
                label="EVM-regr",
            )
        if kwargs.get("std_regression_line") is True:
            rl = linreg_std(ctx.plot_var_obs, ctx.plot_var_model, **kwargs)
            self.linreg = {"intercept": rl["intercept"], "slope": rl["slope"]}
            ax.axline(
                xy1=(0, rl["intercept"]),
                slope=rl["slope"],
                color=kwargs.get("std_regression_col", "lightblue"),
                lw=kwargs.get("std_regression_lw", 1),
                ls=kwargs.get("std_regression_ls", "-"),
                label="linregr",
            )

    # ------------------------------------------------------------------
    # Plot methods  (uniform signature: plot_X(self, ctx, **kwargs))
    # ------------------------------------------------------------------

    def plot_sat(self, ax, **kwargs):
        """Hook: overlay additional satellite geometry on a map axis.

        Override in a subclass or extend here to render e.g. cross-track
        footprints.  The base implementation adds the pulse-limited
        footprint when ``plot_xtrack_pulse_limited_fpr=True`` is passed.
        """
        if kwargs.get("plot_xtrack_pulse_limited_fpr"):
            domain = kwargs.get("domain", "lonlat")
            number_of_seeds = kwargs.get("number_of_seeds", 100)
            lons_perp, lats_perp, _, _, _ = self._generate_xtrack_footprints(
                domain=domain, number_of_seeds=number_of_seeds
            )
            ax.scatter(
                lons_perp,
                lats_perp,
                s=0.2,
                c="b",
                marker=".",
                edgecolor="face",
                transform=ccrs.PlateCarree(),
            )

    def plot_map(self, ctx: PlotContext, **kwargs):
        """Render a geographic map of ``ctx.plot_var``.

        Keyword arguments:
        vmin, vmax           -- colour-scale limits (default: 0, nanmax)
        levels_incr          -- contour level increment (default: 0.5)
        levels               -- explicit contour levels
        cflevels             -- filled-contour levels (default: levels)
        clevels              -- contour-line levels (default: levels)
        zorder_land          -- z-order for land (default: 10)
        land_mask_resolution -- GSHHS scale letter (default: "i")
        lonmax/min, latmax/min -- explicit map bounds
        poi                  -- wavy object to overlay as a track
        cbar                 -- draw colorbar (default: True)
        title                -- axis title (default: auto-generated)
        map_extent_*         -- see :func:`_set_extent`
        """
        vmin = kwargs.get("vmin", 0)
        vmax = kwargs.get("vmax", np.nanmax(ctx.plot_var))
        levels_incr = kwargs.get("levels_incr", 0.5)
        levels = kwargs.get("levels", np.arange(vmin, vmax, levels_incr))
        cflevels = kwargs.get("cflevels", levels)
        clevels = kwargs.get("clevels", levels)
        norm = mpl.colors.BoundaryNorm(levels, ctx.cmap.N)
        zorder_land = kwargs.get("zorder_land", 10)

        land = cfeature.GSHHSFeature(
            scale=kwargs.get("land_mask_resolution", "i"),
            levels=[1],
            facecolor=cfeature.COLORS["land"],
        )

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=ctx.projection)

        ax.add_geometries(
            land.intersecting_geometries([-180, 180, 0, 90]),
            crs=ccrs.PlateCarree(),
            facecolor=cfeature.COLORS["land"],
            edgecolor="black",
            linewidth=1,
            zorder=zorder_land,
        )
        ax.add_feature(land, facecolor="burlywood", alpha=0.5)

        lonmax, lonmin, latmax, latmin = _set_lonlat_minmax(
            ctx.plot_lons, ctx.plot_lats, **kwargs
        )

        if poi := kwargs.get("poi"):
            lonmax, lonmin, latmax, latmin = plot_poi(
                ax, poi, lonmax, lonmin, latmax, latmin, **kwargs
            )

        self.plot_sat(ax, **kwargs)

        sc, _ = plot_var_field(
            ax,
            ctx.plot_var,
            ctx.plot_lons,
            ctx.plot_lats,
            vmin,
            vmax,
            ctx.cmap,
            cflevels,
            clevels,
            norm,
            **kwargs,
        )

        axins = inset_axes(
            ax,
            width="5%",
            height="100%",
            loc="lower left",
            bbox_to_anchor=(1.01, 0.0, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )

        if kwargs.get("cbar", True) is True:
            cbar = fig.colorbar(
                sc,
                cax=axins,
                label=f"{ctx.varalias} [{ctx.units}]",
                ticks=levels,
            )
            cbar.ax.set_ylabel(ctx.units, size=ctx.fs)
            cbar.ax.tick_params(labelsize=ctx.fs)

        _set_extent(ax, lonmax, lonmin, latmax, latmin, ctx.projection, **kwargs)

        if ctx.projection == ccrs.PlateCarree():
            gl = ax.gridlines(
                draw_labels=True,
                crs=ctx.projection,
                linewidth=1,
                color="grey",
                alpha=0.4,
                linestyle="-",
            )
            gl.top_labels = False
            gl.right_labels = False

        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)

        auto_title = (
            self.nID
            + "\nfrom "
            + parse_date(str(self.vars["time"][0].values)).strftime("%Y-%m-%d %H:%M:%S")
            + " to "
            + parse_date(str(self.vars["time"][-1].values)).strftime(
                "%Y-%m-%d %H:%M:%S"
            )
        )
        ax.set_title(kwargs.get("title", auto_title))

        if (
            "region" in vars(self)
            and self.region in quicklook_dict
            and "poi" in quicklook_dict[self.region]
        ):
            for poi_key, poi_cfg in quicklook_dict[self.region]["poi"].items():
                ax.scatter(
                    poi_cfg["lon"],
                    poi_cfg["lat"],
                    s=20,
                    c=poi_cfg.get("color", "b"),
                    marker=poi_cfg["marker"],
                    transform=ccrs.PlateCarree(),
                )
                ax.text(
                    poi_cfg["lon"],
                    poi_cfg["lat"],
                    poi_cfg["name"],
                    transform=ccrs.PlateCarree(),
                    zorder=100,
                )

        if kwargs.get("show", True) is True:
            plt.show()
        return fig, ax

    def plot_timeseries(self, ctx: PlotContext, **kwargs):
        """Render a time series of ``ctx.plot_var`` (and optionally the model)."""
        fig, ax = plt.subplots(figsize=(9, 3.5))

        if ctx.mode == "comb":
            ax.plot(
                self.vars["time"],
                ctx.plot_var,
                color="k",
                linestyle=kwargs.get("linestyle", ""),
                label=self.nID,
                marker="o",
                alpha=0.5,
                ms=2,
            )
            try:
                label = self.model if "model" in vars(self) else self.nID
                ax.plot(
                    self.vars["time"],
                    ctx.plot_var_model,
                    color="r",
                    linestyle=kwargs.get("linestyle", ""),
                    label=label,
                    marker="o",
                    alpha=0.5,
                    ms=2,
                )
            except Exception:
                pass

        elif ctx.mode == "indiv":
            for oco in self.ocos:
                ax.plot(
                    oco.vars["time"],
                    oco.vars[ctx.varalias],
                    linestyle=kwargs.get("linestyle", ""),
                    label=oco.name,
                    marker="o",
                    alpha=0.5,
                    ms=2,
                )
            try:
                ax.plot(
                    self.vars["time"],
                    ctx.plot_var_model,
                    color="r",
                    linestyle=kwargs.get("linestyle", ""),
                    label=self.model,
                    marker="o",
                    alpha=0.5,
                    ms=2,
                )
            except Exception:
                pass
        else:
            raise ValueError(f"mode must be 'comb' or 'indiv', got {ctx.mode!r}")

        ax.set_ylabel(f"{ctx.varalias} [{ctx.units}]")
        ax.legend(loc="best")
        fig.tight_layout()

        if kwargs.get("show", True) is True:
            plt.show()
        return fig, ax

    def plot_scat(self, ctx: PlotContext, **kwargs):
        """Render an obs-vs-model scatter plot with QQ overlay."""
        lq = kwargs.get("lq", np.arange(0.01, 1.01, 0.01))
        modq = compute_quantiles(ctx.plot_var_model, lq)
        obsq = compute_quantiles(ctx.plot_var_obs, lq)

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.plot(
            ctx.plot_var_obs,
            ctx.plot_var_model,
            linestyle="None",
            color="k",
            marker="o",
            alpha=0.5,
            ms=2,
            label="data",
        )
        ax.plot(obsq, modq, "r", label="QQ")
        ax.axline((0, 0), (1, 1), lw=0.5, color="grey", ls="--", label="45 deg")

        self._add_regression_lines(ax, ctx, **kwargs)

        ax.set_xlabel(f"obs ({self.nID})")
        ax.set_ylabel(f"model ({self.model})")
        maxv = np.nanmax(
            [self.vars["model_" + ctx.varalias], self.vars["obs_" + ctx.varalias]]
        )
        ax.set_xlim([0, maxv * 1.05])
        ax.set_ylim([0, maxv * 1.05])
        ax.set_title(f"{ctx.varalias} [{ctx.units}]")
        ax.legend()
        fig.tight_layout()

        if kwargs.get("show", True) is True:
            plt.show()
        return fig, ax

    def plot_hist(self, ctx: PlotContext, **kwargs):
        """Render a 2-D density histogram of obs vs model."""
        lq = kwargs.get("lq", np.arange(0.01, 1.01, 0.01))
        modq = compute_quantiles(ctx.plot_var_model, lq)
        obsq = compute_quantiles(ctx.plot_var_obs, lq)

        fig, ax = plt.subplots(figsize=(5, 4))
        lmin = 0
        lmax = np.nanmax([ctx.plot_var_obs, ctx.plot_var_model]) * 1.05
        _, _, _, im = ax.hist2d(
            ctx.plot_var_obs,
            ctx.plot_var_model,
            bins=kwargs.get("bins", 100),
            range=[[lmin, lmax], [lmin, lmax]],
            norm=kwargs.get("norm", mpl.colors.LogNorm()),
            cmap=kwargs.get("cmap", mpl.cm.gray),
            cmin=kwargs.get("cmin", 1),
        )
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.tick_params(labelsize=ctx.fs)
        cbar.set_label("Frequency", size=ctx.fs)
        ax.set_xlim([lmin, lmax])
        ax.set_ylim([lmin, lmax])

        ax.plot(obsq, modq, "r", label="QQ")
        ax.axline((0, 0), (1, 1), lw=0.5, color="grey", ls="--", label="45 deg")

        self._add_regression_lines(ax, ctx, **kwargs)

        ax.set_xlabel(f"obs ({self.nID})")
        ax.set_ylabel(f"model ({self.model})")
        ax.set_title(f"{ctx.varalias} [{ctx.units}]")
        ax.legend()
        fig.tight_layout()

        if kwargs.get("show", True) is True:
            plt.show()
        return fig, ax

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def quicklook(self, a=False, projection=None, **kwargs):
        """Explore the object by plotting one or more views.

        Parameters
        ----------
        a : bool
            Enable all "default" plots (map + time series).
        projection : cartopy projection, optional
            Map projection; defaults to PlateCarree.
        m : bool
            Map figure.
        ts : bool
            Time series.
        sc : bool
            Obs-vs-model scatter.
        hist : bool
            Obs-vs-model 2-D histogram.
        mode : str
            ``"comb"`` (all data on one plot) or ``"indiv"`` (per-source).
        show : bool
            Call ``plt.show()`` after each figure (default True).
            Pass ``show=False`` to suppress display and receive ``(fig, ax)``.

        Returns
        -------
        (fig, ax) when *show* is False, otherwise None.
        """
        log_level = str(kwargs.get("logging", "WARNING").upper())
        logger.setLevel(getattr(logging, log_level, logging.WARNING))

        ctx = self._build_plot_context(projection=projection, **kwargs)

        fig, ax = None, None
        for kwarg, spec in self._PLOT_REGISTRY.items():
            enabled = kwargs.get(kwarg, a if spec["default"] else False)
            if not enabled:
                continue
            if spec["needs_colloc"] and ctx.plot_var_obs is None:
                logger.warning(
                    "Skipping %r plot: no collocated obs/model data available.", kwarg
                )
                continue
            fig, ax = getattr(self, spec["method"])(ctx, **kwargs)

        if kwargs.get("show") is False:
            return fig, ax

    def quick_anim(self):
        pass


def _check_projection(projection):
    """
    Checks if a projection is specified. If not, defaults to PlateCarree projection."""
    if projection is None:
        projection = ccrs.PlateCarree()
        logger.debug("projection not specified, using default: PlateCarree")
    logger.info(f"projection: {projection}")
    return projection


def _set_lonlat_minmax(lons, lats, **kwargs):
    """
    Sets the min and max values for longitude and latitude based on provided keyword arguments or defaults to the min/max of the provided arrays.
    """
    if kwargs.get("lonmax") is not None:
        lonmax = kwargs.get("lonmax")
    else:
        lonmax = np.max(lons)
    if kwargs.get("latmax") is not None:
        latmax = kwargs.get("latmax")
    else:
        latmax = np.max(lats)
    if kwargs.get("lonmin") is not None:
        lonmin = kwargs.get("lonmin")
    else:
        lonmin = np.min(lons)
    if kwargs.get("latmin") is not None:
        latmin = kwargs.get("latmin")
    else:
        latmin = np.min(lats)
    return lonmax, lonmin, latmax, latmin


def _set_polar_extent(ax, lonmax, lonmin, latmax, latmin, projection, **kwargs):
    """
    Sets the extent of the map for polar stereographic projections based on the provided longitude and latitude ranges and keyword arguments.
    """
    map_extent_multiplicator = kwargs.get("map_extent_multiplicator", 0.1)
    map_extent_multiplicator_lon = kwargs.get(
        "map_extent_multiplicator_lon", map_extent_multiplicator
    )
    map_extent_multiplicator_lat = kwargs.get(
        "map_extent_multiplicator_lat", map_extent_multiplicator
    )
    lon_range = lonmax - lonmin
    lat_range = latmax - latmin
    lat0 = projection._proj4_params.get("lat_0")

    if lat0 > 60:
        # Northern hemisphere polar view
        extent = [
            max(-180, lonmin - lon_range * map_extent_multiplicator_lon),
            min(180, lonmax + lon_range * map_extent_multiplicator_lon),
            max(30, latmin - lat_range * map_extent_multiplicator_lat),
            90,
        ]
        logger.debug("Northern hemisphere projection detected")

    elif lat0 < -60:
        # Southern hemisphere polar view
        extent = [
            max(-180, lonmin - lon_range * map_extent_multiplicator_lon),
            min(180, lonmax + lon_range * map_extent_multiplicator_lon),
            -90,
            min(-30, latmax + lat_range * map_extent_multiplicator_lat),
        ]
        logger.debug("Southern hemisphere projection detected")
    else:
        extent = [
            lonmin - lon_range * map_extent_multiplicator_lon,
            lonmax + lon_range * map_extent_multiplicator_lon,
            latmin - lat_range * map_extent_multiplicator_lat,
            latmax + lat_range * map_extent_multiplicator_lat,
        ]
    return extent


def _set_extent(ax, lonmax, lonmin, latmax, latmin, projection, **kwargs):
    if kwargs.get("map_extent_llon") is None:

        map_extent_multiplicator = kwargs.get("map_extent_multiplicator", 0.1)
        map_extent_multiplicator_lon = kwargs.get(
            "map_extent_multiplicator_lon", map_extent_multiplicator
        )
        map_extent_multiplicator_lat = kwargs.get(
            "map_extent_multiplicator_lat", map_extent_multiplicator
        )

        logger.debug(
            f"map_extent_multiplicator_lon: {map_extent_multiplicator_lon}, map_extent_multiplicator_lat: {map_extent_multiplicator_lat}"
        )

        lon_range = lonmax - lonmin
        lat_range = latmax - latmin

        logger.debug(f"lon_range: {lon_range}, lat_range: {lat_range}")

        # handle polar stereographic projections differently
        if isinstance(
            projection,
            (ccrs.Stereographic, ccrs.NorthPolarStereo, ccrs.SouthPolarStereo),
        ):
            extent = _set_polar_extent(
                ax, lonmax, lonmin, latmax, latmin, projection, **kwargs
            )

        else:
            extent = [
                lonmin - lon_range * map_extent_multiplicator_lon,
                lonmax + lon_range * map_extent_multiplicator_lon,
                latmin - lat_range * map_extent_multiplicator_lat,
                latmax + lat_range * map_extent_multiplicator_lat,
            ]

        logger.debug(
            f"Final extent: {extent}"
        )  # BUG when latmin=40, and get Datarray with lon and lat range

        ax.set_extent(extent, crs=ccrs.PlateCarree())
    elif kwargs.get("map_extent_llon") is False:
        pass
    else:
        ax.set_extent(
            [
                kwargs.get("map_extent_llon"),
                kwargs.get("map_extent_ulon"),
                kwargs.get("map_extent_llat"),
                kwargs.get("map_extent_ulat"),
            ],
            crs=ccrs.PlateCarree(),
        )


def plot_poi(ax, poi, lonmax, lonmin, latmax, latmin, **kwargs):
    """Overlay a POI track on *ax* and return expanded lon/lat bounds.

    Parameters
    ----------
    ax : cartopy GeoAxes
    poi : wavy object with ``poi.vars.lons`` / ``poi.vars.lats``
    lonmax, lonmin, latmax, latmin : float
        Current map bounds, expanded to encompass the POI track.

    Returns
    -------
    lonmax, lonmin, latmax, latmin : float
    """
    plats = poi.vars.lats.data
    plons = poi.vars.lons.data
    platsmax, platsmin = np.max(plats), np.min(plats)
    plonsmax, plonsmin = np.max(plons), np.min(plons)
    # Draw track as a line and as individual points
    ax.plot(
        plons,
        plats,
        color="cornflowerblue",
        ls="-",
        lw=1,
        zorder=-1,
        transform=ccrs.PlateCarree(),
    )
    ax.plot(
        plons,
        plats,
        color="cornflowerblue",
        ls="None",
        marker="o",
        ms=5,
        markeredgecolor="k",
        zorder=-1,
        transform=ccrs.PlateCarree(),
    )
    lonmax = np.max([lonmax, plonsmax])
    lonmin = np.min([lonmin, plonsmin])
    latmax = np.max([latmax, platsmax])
    latmin = np.min([latmin, platsmin])
    return lonmax, lonmin, latmax, latmin


def plot_var_field(
    ax, plot_var, lons, lats, vmin, vmax, cmap, cflevels, clevels, norm, **kwargs
):
    if len(plot_var.shape) > 1:
        sc = ax.contourf(
            lons.squeeze(),
            lats.squeeze(),
            plot_var.squeeze(),
            cmap=cmap,
            levels=cflevels,
            vmin=vmin,
            vmax=vmax,
            norm=norm,
            transform=ccrs.PlateCarree(),
            transform_first=kwargs.get("transform_first", False),
        )
        c = ax.contour(
            lons.squeeze(),
            lats.squeeze(),
            plot_var.squeeze(),
            levels=clevels,
            colors="w",
            linewidths=0.3,
            transform=ccrs.PlateCarree(),
            transform_first=kwargs.get("transform_first", False),
        )
    else:
        sc = ax.scatter(
            lons,
            lats,
            s=15,
            c=plot_var,
            marker="o",  # edgecolor='face',
            edgecolors="k",
            linewidths=0.1,
            cmap=cmap,
            norm=norm,
            transform=ccrs.PlateCarree(),
        )
        c = None
    return sc, c
