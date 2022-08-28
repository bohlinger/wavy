# wavy imports
from wavy.satmod import satellite_class as sc

class satellite_class():
    '''
    Class to handle netcdf files containing satellite data e.g.
    Hs[time], lat[time], lon[time]

    This class offers the following added functionality:
     - get swaths of desired days and read
     - get the closest time stamp(s)
     - get the location (lon, lat) for this time stamp
     - get Hs or 10m wind value for this time
     - region mask
    '''

    def __init__(
        self,sdate=None,mission=['s3a'],product='cmems_L3_NRT',
        edate=None,twin=30,download=False,path_local=None,
        region='mwam4',nproc=1,varalias='Hs',api_url=None,
        filterData=False,poi=None,distlim=None,**kwargs):
        print('# ----- ')
        print(" ### Initializing satellite_class object ###")
        print(" ")
        # parse and translate date input
        sdate_str = kwargs.get(sdate)
        sdate = parse_date(sdate)
        edate_str = kwargs.get(edate)
        edate = parse_date(edate)
        missions = kwargs.get('mission')
        # products: either None, same as missions, or one product
        products = kwargs.get('products') 
        twin = kwargs.get('twin',30)
        for m in missions:
            pass
