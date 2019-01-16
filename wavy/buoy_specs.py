"""
file to specify model specifications such that models can be added 
and data can be imported easily
"""

from datetime import datetime, timedelta

buoy_dict={'Tennholmen':
            {'Hs':3,
            'Tm0':4,
            'time':2, # supposingly UTC
            'lon':13.5618,
            'lat':67.3406,
            'url_template':'157.249.178.22/datawell/waved/Tennholmen/%Y/%m/',
            'file_template':'Tennholmen{0x324}%Y-%m.csv',
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'delta_t':'0000-00-00 (01:00:00)'
            }
        }
