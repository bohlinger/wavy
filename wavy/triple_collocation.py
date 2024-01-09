from wavy.satellite_module import satellite_class as sc
from wavy.insitu_module import insitu_class as ic
from wavy.utils import parse_date, haversineA
from wavy.utils import find_included_times
from wavy.validationmod import validate, disp_validation
import numpy as np
import pandas as pd
import copy
from datetime import datetime, timedelta


def collocate_sat_and_insitu(sco, ico, twin=5, dist_max=200):
    '''
     Collocates satellite and in-situ data, keeping only closest
     satellite data within a given time range around each in-situ
     measurements.

    Args:
        sco (satellite_class object): satellite data
        ico (insitu_class object): in-situ data
        twin (int): time window length in minutes
        dist_max (int | float): Maximum distance
                                accepted to collocate
                                data

    Returns:
        tuple (sco_filter (satellite_class object),
               ico_filter (insitu_class object)):
            Satellite and in-situ collocated data

    '''
    list_time_sat = []
    list_time_insitu = []
    datetimes_ico = [pd.to_datetime(t) for t in ico.vars['time'].values]

    for i, time_insitu in enumerate(datetimes_ico):

        # For each in-situ measure, create a dataset
        # ds_tmp, of satellite data measured at times
        # between [time_insitu-twin, time_insitu+twin]
        time_sup = time_insitu + timedelta(minutes=twin)
        time_inf = time_insitu - timedelta(minutes=twin)
        ds_tmp = sco.vars.sel(time=slice(time_inf, time_sup))
        ds_tmp_size = ds_tmp.sizes['time']

        # if some satellite data are found
        if ds_tmp_size > 0:

            # Gets in-situ corresponding datetime,
            # lat and lon and adds it to ds_tmp
            # as variables
            ico_vars_tmp = ico.vars.isel(time=i)
            lats_insitu = ico_vars_tmp['lats'].values
            lons_insitu = ico_vars_tmp['lons'].values

            ds_tmp = ds_tmp.assign(
                {
                    'insitu_time': (
                        'time',
                        np.repeat(time_insitu, ds_tmp_size)
                    ),
                    'lats_insitu': (
                        'time',
                        np.repeat(lats_insitu, ds_tmp_size)
                    ),
                    'lons_insitu': (
                        'time',
                        np.repeat(lons_insitu, ds_tmp_size)
                    )
                }
            )

            # Calculates distance between satellite and in-situ
            # measurements and adds it as a variable to ds_tmp
            ds_tmp = ds_tmp.assign(
                {
                    'dist_is_sat': (
                        'time',
                        haversineA(
                            ds_tmp['lons_insitu'],
                            ds_tmp['lats_insitu'],
                            ds_tmp['lons'],
                            ds_tmp['lats']
                        )[0]
                    )
                }
            )

            # Calculates the minimum value for distances
            min_dist_tmp = min(ds_tmp['dist_is_sat'].values)

            # If the minimum value is lesser than defined max distance
            if min_dist_tmp <= dist_max:

                # Times for in-situ and satellite measurements,
                # corresponding to the minimum distance are fetched
                dist = list(ds_tmp['dist_is_sat'].values)
                ds_tmp = ds_tmp.isel(time=dist.index(min_dist_tmp))
                list_time_sat.append(ds_tmp['time'].data)
                list_dist_min.append(min_dist_tmp)
                list_time_insitu.append(ds_tmp['insitu_time'].data)

    sco_filter = copy.deepcopy(sco)
    ico_filter = copy.deepcopy(ico)

    # Original sco and ico object vars are filtered using
    # lists of times corresponding to minimum distance
    # between satellite and in-situ observation
    sco_filter.vars = sco_filter.vars.sel(time=list_time_sat)
    sco_filter.vars = sco_filter.vars.assign(
        {
            'colloc_dist': (
                'time',
                list_dist_min
            )
        }
    )
    ico_filter.vars = ico_filter.vars.sel(time=list_time_insitu)

    return sco_filter, ico_filter


def remove_nan(A, B, C):
    '''
    Find indexes of nan values in each of three
    lists, and returns the filtered lists
    '''
    n = len(A)  # add test len are equal

    list_nan_indexes = []

    for i in range(n):
        if np.isnan(A[i]) | np.isnan(B[i]) | np.isnan(C[i]):
            list_nan_indexes.append(i)

    A = [A[i] for i in range(n) if i not in list_nan_indexes]
    B = [B[i] for i in range(n) if i not in list_nan_indexes]
    C = [C[i] for i in range(n) if i not in list_nan_indexes]

    return A, B, C


def mixed_second_moment(A, B):
    '''
    Returns the mixed second moment
    of given lists A and B
    '''
    N = len(A)
    mixed_serie_tmp = [A[i]*B[i] for i in range(N)]

    return np.sum(mixed_serie_tmp)/N


def covariance(A, B):
    '''
    Returns the covariance of lists
    A and B
    '''
    return mixed_second_moment(A, B) - np.mean(A)*np.mean(B)


def variance_estimates(A, B, C):
    '''
    Returns the estimate of vairance of random
    error from list A, according to triple
    collocation method with the covariance notation
    as in Gruber et al., 2016
    '''
    return covariance(A, A) -\
        (covariance(A, B)*covariance(A, C))/covariance(B, C)


def RMSE(A, B, C):
    '''
    Returns the Root Mean Square Error of
    error from A
    '''
    return np.sqrt(variance_estimates(A, B, C))


def SI(A, B, C, ref):
    '''
    Returns the scatter index from A, with
    reference to ref measurements
    '''
    return 100 * RMSE(A, B, C) / np.mean(ref)


def sensitivity_estimates(A, B, C):
    '''
    Returns sensitivity estimates of measurements
    A according to Gruber et al., 2016 definition
    '''
    return (covariance(A, B)*covariance(A, C)) / covariance(B, C)


def SNR(A, B, C, var_err_A):
    '''
    Returns signal-to-noise ratio of measurements A
    according to Gruber et al., 2016 definition
    '''
    return sensitivity_estimates(A, B, C)/var_err_A


def fMSE(A, B, C, var_err_A):
    '''
    Returns fractional Mean Square Error of measurements A
    according to Gruber et al., 2016 definition
    '''
    return 1/(1+SNR(A, B, C, var_err_A))


def SNR_dB(A, B, C, var_err_A):
    '''
    Returns signal-to-noise ratio in dB of measurements A
    according to Gruber et al., 2016 definition
    '''
    return 10*np.log10(sensitivity_estimates(A, B, C)/var_err_A)


def triple_collocation_validate(result_dict, ref=None):
    '''
    Runs the triple collocation given a dictionary
    containing three measurements, returns results
    in a dictionary.

    results_dict: {'name of measurement':list of values}
    ref: Name of one of the measurements, must correspond
    to one key of results_dict

    returns: dict of dict of the metrics for each measurement
    {'name of measurement': {'metric name':metric}}
    '''
    measure_names = list(result_dict.keys())
    measures = list(result_dict.values())

    tc_validate = {'data_sources':{key:{} for key in measure_names},
                     'reference_dataset':ref}

    mean_ref = np.mean(result_dict[ref])

    A = measures[0]
    B = measures[1]
    C = measures[2]

    cov_ab = np.cov(A,B)
    cov_bc = np.cov(B,C)
    cov_ac = np.cov(A,C)

    sens = [(cov_ab[0][1]*cov_ac[0][1])/cov_bc[0][1],
            (cov_ab[0][1]*cov_bc[0][1])/cov_ac[0][1],
            (cov_bc[0][1]*cov_ac[0][1])/cov_ab[0][1]]

    var_est = [cov_ab[0][0] - sens[0],
               cov_ab[1][1] - sens[1],
               cov_bc[1][1] - sens[2]]

    rmse = [np.sqrt(var_est[0]),
            np.sqrt(var_est[1]),
            np.sqrt(var_est[2])]

    si = [rmse[0]/mean_ref*100,
          rmse[1]/mean_ref*100,
          rmse[2]/mean_ref*100]

    snr = [sens[0]/var_est[0],
           sens[1]/var_est[1],
           sens[2]/var_est[2]]

    fmse = [1/(1+snr[0]),
            1/(1+snr[1]),
            1/(1+snr[2])]

    snr_db = [10*np.log10(snr[0]),
              10*np.log10(snr[1]),
              10*np.log10(snr[2])]

    for i,k in enumerate(measure_names):

        tc_validate['data_sources'][k]['var_est'] = var_est[i]
        tc_validate['data_sources'][k]['RMSE'] = rmse[i]
        tc_validate['data_sources'][k]['SI'] = si[i]
        tc_validate['data_sources'][k]['sensitivity'] = sens[i]
        tc_validate['data_sources'][k]['fMSE'] = fmse[i]
        tc_validate['data_sources'][k]['SNR_dB'] = snr_db[i]
        tc_validate['data_sources'][k]['mean'] = np.mean(measures[i])
        tc_validate['data_sources'][k]['std'] = np.std(measures[i])

    return tc_validate


def disp_tc_validation(tc_validate, dec=3):
    '''
    Displays results of triple collocation into a table
    '''

    ref = tc_validate["reference_dataset"]

    data_sources_res = tc_validate["data_sources"]
    data_sources_names = data_sources_res.keys()

    metrics_names = data_sources_res[ref].keys()

    data = np.array([np.round(list(data_sources_res[k].values()), dec)
                    for k in data_sources_names])
    data = np.transpose(data)

    format_row = "{:>12}" * (len(data_sources_names) + 1)
    print(format_row.format("", *data_sources_names))
    for ds, row in zip(metrics_names, data):
        print(format_row.format(ds, *row))

    print("\n The reference for the SI is:", ref)
