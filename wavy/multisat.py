# imports
import numpy as np
# wavy imports
from wavy.satmod import satellite_class as sc
from wavy.consolidate import consolidate_class as cs
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.utils import parse_date
from wavy.wconfig import load_or_default

# read yaml config files:
satellite_dict = load_or_default('satellite_specs.yaml')

class multisat_class(qls):
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
        self,
        sdate = None, edate = None, twin = 30,
        varalias = 'Hs', region = 'global',
        mission = ['s3a'], product = ['cmems_L3_NRT'],
        poi = None, distlim = None, filterData = False,
        download = False, path_local = None,
        nproc = 1, api_url = None,
        **kwargs):
        print('# ----- ')
        print(" ### Initializing multisat_class object ###")
        print(" ")
        missions = mission
        # products: either None, same as missions, or one product
        products = product
        providers = [satellite_dict[p].get('provider') for p in products]
        if len(products) != len(missions):
            if len(products) == 1:
                products = products * len(missions)
            else:
                print("products and missions need to correspond")
                assert len(products) == len(missions)
        print(missions)
        print(products)
        scos = []
        for i,m in enumerate(missions):
            sco = sc( sdate = sdate, edate = edate,
                         twin = twin, distlim = distlim,
                         mission = m, products = products[i],
                         region = region, varalias = varalias,
                         filterData = filterData, poi = poi,
                         nproc = nproc, api_url = api_url,
                         path_local = path_local,
                         **kwargs )
            if 'vars' in vars(sco):
                scos.append(sco)
        cso = cs(scos)
        missions = find_valid_missions(scos)
        cso.rename_consolidate_object_parameters(obstype='satellite_altimeter')
        cso.rename_consolidate_object_parameters(mission='-'.join(missions))
        if len(np.unique(products)) == 1:
            cso.rename_consolidate_object_parameters(product=products[0])
        else:
            cso.rename_consolidate_object_parameters(\
                            product='-'.join(products))
        if len(np.unique(providers)) == 1:
            cso.rename_consolidate_object_parameters(provider=providers[0])
        else:
            cso.rename_consolidate_object_parameters(\
                            provider='-'.join(providers))
        # class variables
        self.obsname = cso.obsname
        self.stdvarname = cso.stdvarname
        self.varalias = cso.varalias
        self.varname = cso.varname
        self.obstype = cso.obstype
        self.label = 'multi-mission-obs'
        self.mission = cso.mission
        self.product = cso.product
        self.provider = cso.provider
        self.sdate = cso.sdate
        self.edate = cso.edate
        self.units = cso.units
        self.region = region
        self.vars = cso.vars
        self.ocos = cso.ocos

        print(" ")
        print (" ### multisat object initialized ###")
        print ('# ----- ')

def find_valid_missions(scos):
    missions = [scos[0].mission]
    for i in range(1,len(scos)):
        if len(scos[i].vars['time']) > 0:
            missions.append(scos[i].mission)
    return missions
