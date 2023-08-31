# imports
import numpy as np
from copy import deepcopy
# wavy imports
from wavy.satellite_module import satellite_class as sc
from wavy.consolidate import consolidate_class as cs
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.utils import parse_date
from wavy.wconfig import load_or_default

# read yaml config files:
satellite_dict = load_or_default('satellite_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')

class multisat_class(qls):
    '''
    Class to combine multiple satellite datasets
    '''

    def __init__(self, **kwargs):
        print('# ----- ')
        print(" ### Initializing multisat_class object ###")
        print(" ")
        # parse and translate date input
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        # add other class object variables
        self.nID = kwargs.get('nID')  # list of nID
        self.mission = kwargs.get('mission', ['s3a'])  # list of missions
        self.varalias = kwargs.get('varalias', 'Hs')
        self.units = variable_def[self.varalias].get('units')
        self.stdvarname = variable_def[self.varalias].get('standard_name')
        self.twin = int(kwargs.get('twin', 30))
        self.distlim = kwargs.get('distlim', 6)
        self.filter = kwargs.get('filter', False)
        self.region = kwargs.get('region', 'global')
        self.reader = kwargs.get('reader', ['read_local_ncfiles'])
        self.path = kwargs.get('path')

        missions = self.mission
        # products: either None, same as missions, or one product
        nIDs = self.nID
        if len(nIDs) != len(missions):
            if len(nIDs) == 1:
                products = nIDs * len(missions)
            else:
                print("products and missions need to correspond")
                assert len(products) == len(missions)
        print(missions)
        print(nIDs)
        scos = []
        for i, m in enumerate(missions):
            sco = sc(sd=self.sd, ed=self.ed,
                     nID=nIDs[i], mission=m,
                     twin=self.twin, distlim=self.distlim,
                     region=self.region, varalias=self.varalias,
                     reader=self.reader[i])
            #sco = sco.populate(reader=self.reader[i])
            sco = sco.populate(reader=self.reader[i],
                               path=self.path[i])
            if 'vars' in vars(sco):
                scos.append(deepcopy(sco))
            del sco
        # consolidate scos
        cso = cs(scos)
        missions = find_valid_missions(scos)
        #cso.rename_consolidate_object_parameters(obstype='satellite_altimeter')
        #cso.rename_consolidate_object_parameters(mission='-'.join(missions))
        # class variables
        self.label = 'multi-mission-obs'

        print(" ")
        print(" ### multisat object initialized ###")
        print('# ----- ')


def find_valid_missions(scos):
    missions = [scos[0].mission]
    for i in range(1, len(scos)):
        if len(scos[i].vars['time']) > 0:
            missions.append(scos[i].mission)
    return missions
