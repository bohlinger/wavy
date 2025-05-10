# imports
import numpy as np
from copy import deepcopy
import time 

# wavy imports
from wavy.satellite_module import satellite_class as sc
from wavy.consolidate import consolidate_class as cs
from wavy.quicklookmod import quicklook_class_sat as qls
from wavy.filtermod import filter_class as fc
from wavy.utils import parse_date
from wavy.wconfig import load_or_default

# read yaml config files:
satellite_dict = load_or_default('satellite_cfg.yaml')
variable_def = load_or_default('variable_def.yaml')

class multisat_class(qls, fc):
    '''
    Class to combine multiple satellite datasets
    '''

    def __init__(self, **kwargs):
        print('# ----- ')
        print(" ### Initializing multisat_class object ###")
        print(" ")
        # parse and translate date input
        self.nID = kwargs.get('nID', ['cmems_L3_NRT'])
        self.name = kwargs.get('name', ['s3a'])
        self.varalias = kwargs.get('varalias', 'Hs')
        self.stdvarname = variable_def[self.varalias].get('standard_name')
        self.sd = parse_date(kwargs.get('sd'))
        self.ed = parse_date(kwargs.get('ed', self.sd))
        self.varalias = kwargs.get('varalias', 'Hs')
        self.units = variable_def[self.varalias].get('units')
        self.twin = kwargs.get('twin', 30)
        self.distlim = kwargs.get('distlim', 6)
        self.region = kwargs.get('region', 'global')
        t0 = time.time()

        # products: either None, same as names, or one product
        if len(self.nID) != len(self.name):
            if len(self.nID) == 1:
                self.nID = self.nID * len(self.name)
            else:
                print("nIDs and names need to correspond")
                assert len(self.nID) == len(self.name)
        scos = []
        for i, n in enumerate(self.name):
            try:
                sco = sc(sd=self.sd, ed=self.ed,
                         nID=self.nID[i], name=n,
                         twin=self.twin, distlim=self.distlim,
                         region=self.region, varalias=self.varalias)
                sco = sco.populate()
                if 'vars' in list(vars(sco)):
                    scos.append(deepcopy(sco))
                del sco
            except Exception as e:
                print(e)
                print('no data found for', n)

        # consolidate scos
        cso = cs(scos)
        self.vars = cso.vars
        self.ocos = cso.ocos
        #cso.rename_consolidate_object_parameters(obstype='satellite_altimeter')
        #cso.rename_consolidate_object_parameters(mission='-'.join(missions))
        self.name = str(self.name)
        self.nID = str(np.unique(self.nID))
        self.obsname = 'consolidated-obs'
        self.obstype = 'consolidated-obs'
        self.label = 'consolidated-obs'
        t1 = time.time()

        print(" ")
        print(' ## Summary:')
        print(str(len(self.vars['time'])) + " footprints retrieved.")
        print("Time used for retrieving data:")
        print(round(t1-t0, 2), "seconds")

        print(" ")
        print(" ### multisat object initialized ###")
        print('# ----- ')


def find_valid_names(scos):
    names = [scos[0].name]
    for i in range(1, len(scos)):
        if len(scos[i].vars['time']) > 0:
            names.append(scos[i].name)
    return names
