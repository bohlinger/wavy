"""
file to specify model specifications such that models can be added 
and data can be imported easily
"""

from datetime import datetime, timedelta

buoy_dict={
        'Tennholmen':
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
            #'Tdw':10,
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
            },
        'Fauskane':
            {'Hm0':'Significant_Wave_Height_Hm0',
            'Tm02':'Wave_Mean_Period_Tm02',
            'time':'time',
            'lon':5.725221,
            'lat':62.56689,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeB/project/met-obs/'
                            +'kystverket/buoy/%Y/%m/'),
            'file_template':('%Y%m_Kystverket-Smartbuoy-Fauskane'
                            +'_AanderaaMotusSensor.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00'
        },
        'SulaA':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.426267625,
            'lon':6.04578018,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_A_Sulafjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaB':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.40252876,
            'lon':6.079721455,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_B_Sulafjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaB1':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.405,
            'lon':6.0585,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_B1_Sulafjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaC':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.39201164,
            'lon':6.050519945,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_C_Sulafjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaC1':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.3969,
            'lon':6.047,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_C1_Sulafjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaD':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.440265655,
            'lon':5.93398094,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_D_Breisundet_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaF':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat':62.22093,
            'lon':5.90045,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_F_Vartalsfjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        },
        'SulaG':
            {'Hm0':'Hm0',
            'Tm02':'tm02',
            'time':'time',
            'lat': 63.08537674,
            'lon': 8.15519333,
            'lats':'longitude',
            'lons':'latitude',
            'path_template':('/lustre/storeA/project/SVV/E39/buoy/%Y/%m/'),
            'file_template':('%Y%m_E39_G_Vartalsfjorden_wave.nc'),
            'basetime':datetime(1970,1,1),
            'units_time':'seconds since 1970-01-01 00:00:00',
            'time_comment':'End of the 10 minute sampling period'
        }
    }
