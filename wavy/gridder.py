# Module to organize gridding data

# imports
import numpy as np
import tqdm

class gridder_class():

    def __init__(
        self, oco=None, cco=None, bb=None, grid='lonlat', res=(1,1)
        ):
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
        if oco is not None:
            self.olons = np.array(oco.vars['longitude'])
            self.olats = np.array(oco.vars['latitude'])
            self.ovals = np.array(oco.vars[oco.stdvarname])
        elif cco is not None:
            self.olons = np.array(cco.vars['longitude'])
            self.olats = np.array(cco.vars['latitude'])
            self.ovals = np.array(cco.vars['obs_values'])
            self.mvals = np.array(cco.vars['model_values'])
        self.bb = bb
        self.res = res
        self.grid = grid
        self.glons, self.glats = self.create_grid_coords()

    def create_grid_coords(self):
        """
        returns grid coordinates
        """
        lons = range(self.bb[0],self.bb[1],self.res[0])
        lats = range(self.bb[2],self.bb[3],self.res[1])
        return np.array(lons), np.array(lats)
    
    def get_obs_grid_idx(self):
        Midx = self.assign_obs_to_grid(
                self.glons,self.glats,
                self.olons,self.olats,
                self.ovals,self.res)
        ovals, Midx = self.clean_Midx(Midx,self.ovals,self.glons,self.glats)
        return ovals, Midx

    @staticmethod
    def assign_obs_to_grid(glons,glats,olons,olats,ovals,res):
        """
        assigns observation coordinates to grid indices
        """
        lonidx = ((olons-np.min(glons))/res[0]).astype(int)
        latidx = ((olats-np.min(glats))/res[1]).astype(int)
        Midx = np.array([lonidx,latidx],dtype=object)
        return Midx

    @staticmethod
    def clean_Midx(Midx,ovals,glons,glats):
        """
        cleans Midx and observations (ovals) for grid cells outside bb
        """
        # clean x-dim (tlons)
        glons_idx = np.where((Midx[0]>=0)&(Midx[0]<len(glons)))[0]
        Midx = Midx[:,glons_idx]
        ovals = ovals[glons_idx]
        # clean y-dim (tlats)
        glats_idx = np.where((Midx[1]>=0)&(Midx[1]<len(glats)))[0]
        Midx = Midx[:,glats_idx]
        ovals = ovals[glats_idx]
        return ovals, Midx

    @staticmethod
    def region_filter():
        # filter grid cells for region of interest
        # region could be rectangular, polygon, ...
        # -> return grid
        return

    @staticmethod
    def create_grid_polygons(bb,res):
        #def get_geometries(top_left, bottom_right, spacing=1):
        # bb tupel: (lonmin, lonmax, latmin, latmax)
        polygons = []
        points = []
        xmin = bb[0]
        xmax = bb[1]
        ymax = bb[3]
        y = bb[2]
        i = -1
        while True:
            if y > ymax:
                break
            x = xmin
            while True:
                if x > xmax:
                    break
                #components for polygon grid
                polygon = box(x, y, x+spacing, y+spacing)
                polygons.append(polygon)
                #components for point grid
                point = Point(x, y)
                points.append(point)
                i = i + 1
                x = x + spacing
            y = y + spacing
        return polygons, points

    @staticmethod
    def get_grid_idx(Midx):
        return np.where((Midx[0]==Midx[0][0]) & (Midx[1]==Midx[1][0]))[0]

    @staticmethod
    def calc_mean(gidx,ovals):
        return np.mean(np.array(ovals)[gidx])

    @staticmethod
    def rm_used_idx_from_Midx(gidx,Midx):
        d1 = np.delete(Midx[0,:],gidx)
        d2 = np.delete(Midx[1,:],gidx)
        return np.array([d1,d2],dtype=object)

    def grid_mean(self,Midx,ovals):
        # initialize grid
        val_grid = np.ones([len(self.glons),len(self.glats)])*np.nan
        lon_grid = np.ones([len(self.glons),len(self.glats)])*np.nan
        lat_grid = np.ones([len(self.glons),len(self.glats)])*np.nan
        pbarlen = len(ovals)
        novals = ovals.copy()
        # fill grid
        nMidx = Midx.copy()
        pbar = tqdm.tqdm(total=pbarlen)
        while len(novals)>0:
            l1 = len(novals)
            idx1,idx2 = nMidx[0][0], nMidx[1][0]
            gidx = self.get_grid_idx(nMidx)
            mval = self.calc_mean(gidx,novals)
            val_grid[idx1,idx2] = mval
            lon_grid[idx1,idx2] = self.glons[idx1]
            lat_grid[idx1,idx2] = self.glats[idx2]
            nMidx = self.rm_used_idx_from_Midx(gidx,nMidx)
            novals = np.delete(novals,gidx)
            l2 = len(novals)
            incr = l1-l2
            pbar.update(incr)
        return val_grid, lon_grid, lat_grid

    def apply_metric(self,Midx,ovals,**kwargs):
        '''
        dispatch table for various validation metrics
        '''
        dispatch_reader = {
                'mean':self.grid_mean,
                }
        metric = kwargs.get('metric')
        var_gridded = dispatch_reader[metric](Midx,ovals)
        return var_gridded

    @staticmethod
    def get_exteriors(glons,glats,res):
        xb = np.array([ glons+res[0],glons+res[0],
                        glons,glons,
                        glons+res[0] ])
        yb = np.array([ glats,
                        glats+res[1],glats+res[1],
                        glats,glats ])
        return xb,yb

    def quicklook(self,a=True,projection=None,**kwargs):
        lon_grid = kwargs.get('lon_grid')
        lat_grid = kwargs.get('lat_grid')
        val_grid = kwargs.get('val_grid')
        g = kwargs.get('g',a)
        m = kwargs.get('m',a)
        if m:
            import cartopy.crs as ccrs
            import cmocean
            import matplotlib.pyplot as plt
            import matplotlib.cm as mplcm
            import matplotlib as mpl
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            import math

            lons = lon_grid.ravel()
            lats = lat_grid.ravel()
            vals = val_grid.ravel()

            if projection is None:
                projection = ccrs.PlateCarree()
            # parse kwargs
            if kwargs.get('cmap') is None:
                cmap = cmocean.cm.amp
            else:
                cmap = kwargs.get('cmap')

            normalize = mpl.colors.Normalize(
                            vmin=np.nanmin(val_grid),
                            vmax=np.nanmax(val_grid) )

            lonmax,lonmin = np.nanmax(lons),np.nanmin(lons)
            latmax,latmin = np.nanmax(lats),np.nanmin(lats)

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_extent(  [lonmin, lonmax,latmin, latmax],
                            crs = projection )
            for i in range(len(val_grid.ravel())):
                xs,ys = self.get_exteriors(lons[i],lats[i],self.res)
                if math.isnan(np.mean(xs)):
                    pass
                else:
                    ax.fill(
                        xs, ys,transform=ccrs.Geodetic(),
                        color=cmap(normalize(vals[i])) )

            cbax = fig.add_axes([0.85, 0.12, 0.05, 0.78])
            cb = mpl.colorbar.ColorbarBase(cbax, 
                    cmap=cmap, norm=normalize, 
                    orientation='vertical')

            ax.coastlines()
            gl = ax.gridlines(draw_labels=True,crs=projection,
                              linewidth=1, color='grey', alpha=0.4,
                              linestyle='-')
            gl.top_labels = False
            gl.right_labels = False
            plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
            plt.show()

