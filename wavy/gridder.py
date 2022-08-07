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
        return

    def assign_vals(self):
        return

    def aggregate(self):
        return

