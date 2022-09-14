import os
from dateutil.relativedelta import relativedelta
#from datetime import datetime, timedelta

from wavy.ncmod import find_attr_in_nc, dumptonc_ts_sat
from wavy.utils import make_pathtofile, make_subdict
from wavy.wconfig import load_or_default

satellite_dict = load_or_default('satellite_specs.yaml')

class writer_class:

    def write_to_pickle(self,pathtofile=None):
        """
        Writes class object to pickle

        param:
            pathtofile - string
        """
        import pickle
        pickle.dump( self, open( pathtofile, "wb" ) )
        print( type(self), ' written to:', pathtofile )

    def write_to_nc(self,pathtofile=None,file_date_incr=None):
        tmpdate = self.sdate
        edate = self.edate
        while tmpdate <= edate:
            if pathtofile is None:
                path_template = satellite_dict[self.product]\
                                              ['dst']\
                                              ['path_template']
                file_template = satellite_dict[self.product]\
                                              ['dst']\
                                              ['file_template']
                strsublst = satellite_dict[self.product]\
                                          ['dst']['strsub']
                if 'filterData' in vars(self).keys():
                    file_template = 'filtered_' + file_template
                tmppath = os.path.join(path_template,file_template)
                subdict = make_subdict(strsublst,
                                       class_object_dict=vars(self))
                pathtofile = make_pathtofile(tmppath,strsublst,
                                             subdict,
                                             date=tmpdate)
            title = (self.obstype
                   + ' observations from '
                   + self.mission)
            dumptonc_ts_sat(self,pathtofile,title)
            # determine date increment
            if file_date_incr is None:
                file_date_incr = satellite_dict[self.product]\
                                ['dst'].get('file_date_incr','m')
            if file_date_incr == 'm':
                tmpdate += relativedelta(months = +1)
            elif file_date_incr == 'Y':
                tmpdate += relativedelta(years = +1)
            elif file_date_incr == 'd':
                tmpdate += timedelta(days = +1)

