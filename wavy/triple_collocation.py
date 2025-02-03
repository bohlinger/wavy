from wavy.satellite_module import satellite_class as sc
from wavy.insitu_module import insitu_class as ic
from wavy.utils import parse_date, haversineA
from wavy.utils import find_included_times
from wavy.validationmod import validate, disp_validation
import numpy as np
import pandas as pd
import copy
from datetime import datetime, timedelta
import random


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


def triple_collocation_validate(result_dict,
                                metric_list=['var_est', 'rmse', 'si',
                                             'rho', 'mean', 'std'],
                                ref=None):
    '''
    Runs the triple collocation given a dictionary
    containing three measurements, returns results
    in a dictionary.

    results_dict: {'name of measurement':list of values}
    metric_list: Str "all" or List of the metrics to return, among 'var_est',
    'rmse', 'si', 'rho', 'sens', 'snr', 'snr_db', 'fmse', 'mean',
    'std'
    ref: Name of one of the measurements, must correspond
    to one key of results_dict

    returns: dict of dict of the metrics for each measurement
    {'name of measurement': {'metric name':metric}}
    '''
    measure_names = list(result_dict.keys())
    measures = list(result_dict.values())

    if ref is None:
        ref = measures_names[0]
        print('No reference was specified.',
              ref, 'was set as the reference.')

    tc_validate = {'data_sources': {key: {} for key in measure_names},
                   'reference_dataset': ref}

    mean_ref = np.mean(result_dict[ref])

    A = measures[0]
    B = measures[1]
    C = measures[2]

    cov_ab = np.cov(A, B)
    cov_bc = np.cov(B, C)
    cov_ac = np.cov(A, C)

    # Sensitivity
    sens = [(cov_ab[0][1]*cov_ac[0][1])/cov_bc[0][1],
            (cov_ab[0][1]*cov_bc[0][1])/cov_ac[0][1],
            (cov_bc[0][1]*cov_ac[0][1])/cov_ab[0][1]]

    # Estimate of the variance of random error
    var_est = [cov_ab[0][0] - sens[0],
               cov_ab[1][1] - sens[1],
               cov_bc[1][1] - sens[2]]

    # Root Mean Square Error
    rmse = [np.sqrt(var_est[0]),
            np.sqrt(var_est[1]),
            np.sqrt(var_est[2])]

    # Scatter Index
    if 'si' in metric_list or metric_list == 'all':
        si = [rmse[0]/mean_ref*100,
              rmse[1]/mean_ref*100,
              rmse[2]/mean_ref*100]

    # Signal to Noise Ratio
    snr = [sens[0]/var_est[0],
           sens[1]/var_est[1],
           sens[2]/var_est[2]]

    # Fractional Mean Squared Error
    if 'fmse' in metric_list or metric_list == 'all':
        fmse = [1/(1+snr[0]),
                1/(1+snr[1]),
                1/(1+snr[2])]

    # Signal to Noise Ratio (dB)
    if 'snr_db' in metric_list or metric_list == 'all':
        snr_db = [10*np.log10(snr[0]),
                  10*np.log10(snr[1]),
                  10*np.log10(snr[2])]

    # Data truth correlation
    if 'rho' in metric_list or metric_list == 'all':
        rho = [sens[0]/cov_ab[0][0],
               sens[1]/cov_ab[1][1],
               sens[2]/cov_bc[1][1]]

    for i, k in enumerate(measure_names):
        if 'var_est' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['var_est'] = var_est[i]
        if 'rmse' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['RMSE'] = rmse[i]
        if 'si' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['SI'] = si[i]
        if 'sens' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['sensitivity'] = sens[i]
        if 'rho' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['rho'] = rho[i]
        if 'snr' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['SNR'] = snr[i]
        if 'fmse' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['fMSE'] = fmse[i]
        if 'snr_db' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['SNR_dB'] = snr_db[i]
        if 'mean' in metric_list or metric_list == 'all':
            tc_validate['data_sources'][k]['mean'] = np.mean(measures[i])
        if 'std' in metric_list or metric_list == 'all':
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


def tc_analysis_wavy_objects(data, ref=None,
                             metric_list=['var_est', 'rmse', 'si',
                                          'rho', 'mean', 'std'],
                             calibrate=True, calibration_ref=None):
    '''
    Runs the triple collocation given a dictionary
    containing three wavy objects, returns results
    in a dictionary.

    data: {'name of measurement':wavy object of satellite,
                                 model or insitu class}
    metric_list: Str "all" or List of the metrics to return, among 'var_est',
    'rmse', 'si', 'rho', 'sens', 'snr', 'snr_db', 'fmse', 'mean',
    'std'
    ref: Name of one of the measurements, must correspond
    to one key of results_dict. Data source to use for metrics
    that depend on one reference data source.
    calibrate: True or False. Choose if calibration should be applied
    to the data or not.
    calibration_ref: Name of the measurement which will be used to calibrate
    the two other data sources.

    returns: dict of dict of the metrics for each measurement
    {'name of measurement': {'metric name':metric}}
    '''
    list_keys = list(data.keys())

    if ref is None:
        ref = list_keys[0]
        print("Since no reference was given,",
              ref, "is taken as a reference.")

    # Fetch the 1D-arrays for each data source
    data_bis = {}
    for k in list_keys:
        data_bis[k] = data[k].vars.Hs.values

    # Recreate dict of data, removing the nan values
    list_data = list(data_bis.values())
    list_data_na_rm = remove_nan(list_data[0],
                                 list_data[1],
                                 list_data[2])

    for i, k in enumerate(list_keys):
        data_bis[k] = list_data_na_rm[i]

    # Recreate dict of data, after calibrating
    if calibrate:
        if calibration_ref is None:
            calibration_ref = list_keys[0]
            print("Since no reference was given for the calibration,",
                  calibration_ref, "is taken as a reference.")

        # Get the two keys of data sources that are not the calibration
        # reference in order to give the right arguments for calibration
        keys_left = [k for k in list_keys if k != calibration_ref]
        # Apply calibration
        calibrated_data = calibration(data_bis[calibration_ref],
                                      data_bis[keys_left[0]],
                                      data_bis[keys_left[1]])
        # Overwrite the calibrated data in the data dictionary
        data_bis[keys_left[0]] = calibrated_data[0]
        data_bis[keys_left[1]] = calibrated_data[1]

    # Apply triple collocation analysis
    tc_results = triple_collocation_validate(data_bis,
                                             metric_list=metric_list,
                                             ref=ref)

    return tc_results


def calibration(R, A, B):
    '''
    Calibrate A and B relatively to R using triple collocation calibration
    constant estimates, following Gruber et al., 2016 method.

    R (list of floats): Reference data to use for calibration.
    A, B (lists of floats): Data series to calibrate relatively to the
    reference.

    returns:
    A_R, B_R (lists of floats): Calibrated data series
    '''
    c_AB = np.cov(A, B)[0, 1]
    c_RA = np.cov(R, A)[0, 1]
    c_RB = np.cov(R, B)[0, 1]

    a_A = c_AB/c_RB
    a_B = c_AB/c_RA
    a_R = 1

    A_R = (a_R/a_A)*(A - np.mean(A)) + np.mean(R)
    B_R = (a_R/a_B)*(B - np.mean(B)) + np.mean(R)

    return A_R, B_R


def bootstrap_ci(result_dict,
                 conf=95.,
                 sample_size=None,
                 metric='var_est',
                 n_bootstrap=100,
                 ref=None):
    '''
    Calculate 95% confidence intervals for a given estimate
    returned by triple_collocation_validate function, using
    bootstrap sampling.

    results_dict: {'name of measurement':list of values}
    sample_size: size of the bootstrap samples that will
    be drawn from the data.
    metric: metric for which the confidence interval is
    required. Must be one of the list that can be given
    to triple_collocation_validate metric_list argument.
    n_bootstrap: number of generated bootstrap samples.
    ref: Name of one of the measurements, must correspond
    to one key of results_dict

    returns: dict containing mean, upper and lower
    value of the 95% confidence interval for each
    data source.
    {'name of measurement': {'metric name':metric}}
    '''
    measure_names = list(result_dict.keys())
    measures = list(result_dict.values())

    if ref is None:
        ref = measure_names[0]
    if sample_size is None:
        sample_size = len(measures[0])

    indexes = np.arange(0, len(measures[0]))

    dict_res_tc = {m: [] for m in measure_names}

    for i in range(n_bootstrap):

        rand_ind_tmp = random.choices(indexes, k=sample_size)
        result_dict_tmp = {key: [result_dict[key][j] for j in rand_ind_tmp] for
                           key in measure_names}
        tc_validate_tmp = triple_collocation_validate(result_dict_tmp,
                                                      metric_list=[metric],
                                                      ref=ref)

        for m in measure_names:
            dict_res_tc[m].append(
                list(tc_validate_tmp['data_sources'][m].values())[0]
                                 )

    dict_ci = {m: {} for m in measure_names}
    for m in measure_names:
        dict_ci[m]['mean'] = np.mean(dict_res_tc[m])
        percentiles = np.percentile(dict_res_tc[m],
                                    [(100-conf)/2, conf + (100-conf)/2])
        dict_ci[m]['ci_l'] = percentiles[0]
        dict_ci[m]['ci_u'] = percentiles[1]

    return dict_ci
