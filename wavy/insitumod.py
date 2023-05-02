#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from stations. I try to mostly follow the PEP
convention for python code style. Constructive comments on style and
effecient programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import numpy as np
from datetime import datetime, timedelta
import time
import os
from dateutil.relativedelta import relativedelta
import traceback
import logging
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import dumptonc_ts_insitu
from wavy.ncmod import get_filevarname
from wavy.utils import make_pathtofile, get_pathtofile
from wavy.utils import finditem, make_subdict
from wavy.utils import parse_date
from wavy.wconfig import load_or_default
from wavy.insitu_readers import insitu_reader
from wavy.filtermod import filter_class as fc
# ---------------------------------------------------------------------#

# read yaml config files:
insitu_dict = load_or_default('insitu_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')
d22_dict = load_or_default('d22_var_dicts.yaml')
# ---------------------------------------------------------------------#


class insitu_class(fc):
    '''
    Class to handle insitu based time series.
    '''
    basedate = datetime(1970, 1, 1)

    def __init__(self, nID, sd, ed, varalias='Hs',
    sensor=None, **kwargs):
        # parse and translate date input
        sd = parse_date(sd)
        ed = parse_date(ed)
        print('# ----- ')
        print(" ### Initializing insitu_class object ###")
        print(" ")
        print('Chosen period: ' + str(sd) + ' - ' + str(ed))
        stdvarname = variable_info[varalias]['standard_name']
        if sensor is None:
            sensor = list(insitu_dict[nID]['sensor'].keys())[0]
            print('Sensor was not chosen')
            print('Automatic choice:', sensor)
        # for i in range(1):
        try:
            self.stdvarname = stdvarname
            self.varalias = varalias
            self.units = variable_info[varalias].get('units')
            self.sd = sd
            self.ed = ed
            self.nID = nID
            self.sensor = sensor
            self.obstype = 'insitu'
            if ('tags' in insitu_dict[nID].keys() and
            len(insitu_dict[nID]['tags']) > 0):
                self.tags = insitu_dict[nID]['tags']
            print(" ")
            print(" ## Read files ...")
            t0 = time.time()
            vardict, fifo, pathtofile = \
                get_insitu_ts(nID=nID, sensor=sensor,
                              sd=sd, ed=ed,
                              varalias=varalias,
                              basedate=self.basedate,
                              dict_for_sub=vars(self),
                              **kwargs)
            self.vars = vardict
            self.lat = np.nanmean(vardict['latitude'])
            self.lon = np.nanmean(vardict['longitude'])
            if fifo == 'nc':
                print(pathtofile)
                meta = ncdumpMeta(sd.strftime(pathtofile))
                self.vars['meta'] = meta
                varname = get_filevarname(varalias,
                                          variable_info,
                                          insitu_dict[nID],
                                          meta)
                self.varname = varname
            else:
                self.varname = varalias
            if fifo == 'frost':
                self.sensor = sensor
            # create label for plotting
            self.label = self.nID + '_' + self.sensor
            t1 = time.time()
            print(" ")
            print('## Summary:')
            print(str(len(self.vars['time'])) + " values retrieved.")
            print("Time used for retrieving insitu data:",
                   round(t1-t0, 2), "seconds")
            print(" ")
            print(" ### insitu_class object initialized ### ")
        except Exception as e:
            logger.exception(e)
            print(e)
            self.error = e
            print("! No insitu_class object initialized !")
        print('# ----- ')

    def get_item_parent(self, item, attr):
        ncdict = self.vars['meta']
        lst = [i for i in ncdict.keys() \
                if (attr in ncdict[i].keys() \
                and item in ncdict[i][attr]) \
                ]
        if len(lst) >= 1:
            return lst
        else:
            return None

    def get_item_child(self, item):
        ncdict = self.vars['meta']
        parent = finditem(ncdict, item)
        return parent

    def quicklook(self, a=False, projection=None, **kwargs):
        m = kwargs.get('m', a)
        ts = kwargs.get('ts', a)
        if ts is True:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(9, 3.5))
            ax = fig.add_subplot(111)
            colors = ['k']
            ax.plot(self.vars['datetime'],
                    self.vars[self.stdvarname],
                    linestyle='None', color=colors[0],
                    label=self.nID + ' ( ' + self.sensor + ' )',
                    marker='o', alpha=.5, ms=2)
            plt.ylabel(self.varalias + '[' + self.units + ']')
            plt.legend(loc='best')
            plt.tight_layout()
            # ax.set_title()
            plt.show()

    def write_to_nc(self,pathtofile=None,file_date_incr=None):
        # divide time into months by loop over months from sdate to edate
        if 'error' in vars(self):
            print('Erroneous insitu_class file detected')
            print('--> dump to netCDF not possible !')
        else:
            tmpdate = self.sd
            edate = self.ed
            while tmpdate <= edate:
                if pathtofile is None:
                    path_template = insitu_dict[self.nID]['dst']\
                                               ['path_template'][0]
                    file_template = insitu_dict[self.nID]['dst']\
                                                ['file_template']
                    strsublst = insitu_dict[self.nID]['dst']['strsub']
                    if 'filterData' in vars(self).keys():
                        file_template = 'filtered_' + file_template
                    tmppath = os.path.join(path_template,file_template)
                    pathtofile = make_pathtofile(tmppath,strsublst,
                                                 vars(self),
                                                 date=tmpdate)
                title = ( self.varalias + ' observations from '
                        + self.nID + ' ' + self.sensor )
                dumptonc_ts_insitu(self,pathtofile,title)
                # determine date increment
                if file_date_incr is None:
                    file_date_incr = insitu_dict[self.nID]\
                                    ['src'].get('file_date_incr','m')
                if file_date_incr == 'm':
                    tmpdate += relativedelta(months = +1)
                elif file_date_incr == 'Y':
                    tmpdate += relativedelta(years = +1)
                elif file_date_incr == 'd':
                    tmpdate += timedelta(days = +1)
        return

    def write_to_pickle(self, pathtofile=None):
        import pickle
        # writing
        pickle.dump( self, open( pathtofile, "wb" ) )
        print('insitu_class object written to:',pathtofile)
        # for reading
        # ico = pickle.load( open( pathtofile, "rb" ) )


def get_insitu_ts(nID, sensor, sd, ed, varalias, basedate,
dict_for_sub, **kwargs):
    # determine fifo
    fifo = kwargs.get('fifo', insitu_dict[nID]['fifo'])
    kwargs['fifo'] = fifo
    path_template = insitu_dict[nID]['src']['path_template']
    file_template = insitu_dict[nID]['src']['file_template']
    pathlst = [p + ('/' + file_template) for p in path_template]
    strsublst = insitu_dict[nID]['src']['strsub']
    if 'path_local' in kwargs.keys():
        pathlst = [kwargs['path_local'] + '/' + file_template]
    if fifo == 'frost':
        pathtofile = 'frost.api.no'
    else:
        subdict = make_subdict(strsublst, class_object_dict=dict_for_sub)
        pathtofile = get_pathtofile(pathlst, strsublst, subdict, sd)
    vardict = insitu_reader(nID=nID, sensor=sensor,
                            sd=sd, ed=ed,
                            varalias=varalias,
                            basedate=basedate,
                            pathlst=pathlst,
                            strsublst=strsublst,
                            dict_for_sub=dict_for_sub,
                            **kwargs)
    return vardict, fifo, pathtofile
