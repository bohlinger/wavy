# Module to organize gridding data

class gridder_class():

    def __init__(
        obs_obj=None, col_obj=None,
        grid=None, grid_type='lonlat', resolution = 1
        ):
        # setup the gridder
        # grid_type: lonlat or in m
        # grid tupel: lonmin, lonmax, latmin, latmax
        # resolution: resolution, either integer or tupel e.g. (.5,.5)
        #             in degree where tupel is (lon,lat)

    def make_grid(self):
        # create grid
        # filter for land/sea
        # -> return grid
        return

    def regionalize_grid(self):
        # filter grid cells for region of interest
        # -> return grid
        return

    def assign_vals(self):
        # assign obs to grid cells
        # -> return index mask
        return

    def operate(self,operation="mean"):
        return

