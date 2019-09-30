"""
Additional information on variable retrieved from wave model
"""
var_dict = {
                'time':
                    {
                    'units':'seconds since',
                    'standard_name':'time',
                    },
                'Hs':
                    {
                    'units':'m',
                    'standard_name':'significant_height_from_wave_model',
                    'valid_range':[0., 30.]
                    },
                'Tp':
                    {
                    'units':'s',
                    'standard_name':'peak_period_from_wave_model',
                    'valid_range':[0., 30.]
                    },
                'Tm':
                    {
                    'units':'s',
                    'standard_name':'mean_period_from_wave_model',
                    'valid_range':[0., 30.]
                    },
                }
