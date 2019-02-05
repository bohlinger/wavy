"""
file to specify model specifications such that models can be added 
and data can be imported easily
"""

from datetime import datetime, timedelta

buoy_dict={'Tennholmen':
            {'Hm0':3,
            'Tm02':4,
            'time':2, # supposingly UTC
            'lon':13.5618,
            'lat':67.3406,
            'lats':3,
            'lons':4,
            'Hs':3,
            'TI':4,
            'TE':5,
            'T1':6,
            'TZ':7,
            'T3':8,
            'Tc':9,
            'Tdw':10,
            'Tp':11,
            'Qp':12,
            'url_template':'157.249.178.22/datawell/waved/Tennholmen/%Y/%m/',
            # file 303: containing  - compressed spectral parameters
            #                         Hs, 
            #                         TI, TE, T1, TZ, T3, Tc, Tdw, Tp, 
            #                         Qp
            'file_template_spec_ext':'Tennholmen{0x303}%Y-%m.csv',
            # file 324: containing  - the significant wave height, 
            #                       - the mean period (Tz), 
            #                       - and the value of the spectral peak
            'file_template':'Tennholmen{0x324}%Y-%m.csv',
            # file 320: containing  - the power spectral density as a 
            #           function of frequency, i.e. the 1D spectrum
            'file_template_heave_spec':'Tennholmen{0x320}%Y-%m.csv',
            # file 321: containing  - the direction from which the waves 
            #           arrive, as well as the directional spread.
            #           Both are given in radians for each frequency 
            #           component. When converted to degrees we have
            #           0 = North, 90 = East, etc.
            'file_template_prim_dir_spec':'Tennholmen{0x321}%Y-%m.csv',
            # file 328: containing  -  additional spectral parameters 
            #           to determine the 2D spectrum, namely m2, n2, K.
            'file_template_sec_dir_spec':'Tennholmen{0x328}%Y-%m.csv',
            # file 380: containing  - the buoy position given in radians.
            'file_template_pos':'Tennholmen{0x380}%Y-%m.csv',
            # file 3C1: containing  - the battery status and some 
            #           additional info on accelerations (max? variance?)
            'file_template_tec':'Tennholmen{0x3C1}%Y-%m.csv',
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'delta_t':'0000-00-00 (01:00:00)'
            }
        }
