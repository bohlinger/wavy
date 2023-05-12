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
    Class to combine multiple satellite datasets
    '''

    def __init__(self, **kwargs):
#        self,
#        sdate = None, edate = None, twin = 30,
#        varalias = 'Hs', region = 'global',
#        mission = ['s3a'], product = ['cmems_L3_NRT'],
#        poi = None, distlim = None, filterData = False,
#        download = False, path_local = None,
#        nproc = 1, api_url = None,
#        **kwargs):
        print('# ----- ')
        print(" ### Initializing multisat_class object ###")
        print(" ")
        # parse and translate date input
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        # add other class object variables
        self.nID = kwargs.get('nID') # list of nID
        self.mission = kwargs.get('mission', ['s3a']) # list of missions
        self.varalias = kwargs.get('varalias', 'Hs')
        self.units = variable_info[self.varalias].get('units')
        self.stdvarname = variable_info[self.varalias].get('standard_name')
        self.twin = int(kwargs.get('twin', 0))
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')

        missions = self.mission
        # products: either None, same as missions, or one product
        nIDs = self.nID
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
            sco = sc(sd=sd, ed=ed, nID=nID[i], mission=m,
                     twin=twin, distlim=distlim,
                     region=region, varalias=varalias,
                     **kwargs)
            sco = sco.populate()
            if 'vars' in vars(sco):
                scos.append(sco)
        # consolidate scos
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
