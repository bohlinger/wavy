from wavy.satellite_module import satellite_class as sc
from wavy.insitu_module import insitu_class as ic
from wavy.model_module import model_class as mc
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


def triple_collocation(data,
                       metrics=['var', 'rmse', 'si',
                                'rho', 'mean', 'std'],
                       r2=0, 
                       ref=None,
                       dec=3):
    '''
    Runs the triple collocation given a dictionary
    containing three measurements, returns results
    in a dictionary.

    data: {'name of measurement':list of values}
    metrics: Str "all" or List of the metrics to return, among 'var',
    'rmse', 'si', 'rho', 'sensitivity', 'snr', 'snr_db', 'fmse', 'mean',
    'std'
    r2: representativeness error or cross correlation error between the 
    first two measurements in data. Default 0. 
    ref: Name of one of the measurements, must correspond
    to one key of data. Default first key from data. 
    dec: Number of decimals to round the results to. Default 3.

    returns: dict of dict of the metrics for each measurement
    {'name of measurement': {'metric name':metric}}
    '''
    measure_names = list(data.keys())

    if ref is None:
        ref = measure_names[0]
    
    if isinstance(data[measure_names[0]], (np.ndarray, list)): 
        measures = list(data.values())
        mean_ref = np.mean(data[ref])
    elif isinstance(data[measure_names[0]], (sc, ic, mc)):
        measures = [data[k].vars.Hs.values for k in measure_names]
        mean_ref = np.mean(data[ref].vars.Hs.values)
    
    results = {key: {} for key in measure_names}

    A = measures[0]
    B = measures[1]
    C = measures[2]

    cov_ab = np.cov(A, B)
    cov_bc = np.cov(B, C)
    cov_ac = np.cov(A, C)

    # Sensitivity
    sens = [(cov_ac[0][1]*(cov_ab[0][1] - r2))/cov_bc[0][1],
            ((cov_ab[0][1] - r2)*cov_bc[0][1])/cov_ac[0][1],
            (cov_bc[0][1]*cov_ac[0][1])/(cov_ab[0][1] - r2)]

    # Estimate of the variance of random error
    if any(m in metrics for m in ['var', 'rmse', 
                                  'si', 'snr',
                                  'snr_db']) or metrics == 'all':
        var = [cov_ab[0][0] - sens[0],
               cov_ab[1][1] - sens[1],
               cov_bc[1][1] - sens[2]]

    # Root Mean Square Error
    if any(m in metrics for m in ['rmse', 'si']) or metrics == 'all':
        rmse = [np.sqrt(var[0]),
                np.sqrt(var[1]),
                np.sqrt(var[2])]

    # Scatter Index
    if 'si' in metrics or metrics == 'all':
        si = [rmse[0]/mean_ref,
              rmse[1]/mean_ref,
              rmse[2]/mean_ref]

    # Signal to Noise Ratio
    if any(m in metrics for m in ['snr', 'snr_db']) or metrics == 'all':
        snr = [sens[0]/var[0],
               sens[1]/var[1],
               sens[2]/var[2]]

    # Fractional Mean Squared Error
    if 'fmse' in metrics or metrics == 'all':
        fmse = [1/(1+snr[0]),
                1/(1+snr[1]),
                1/(1+snr[2])]

    # Signal to Noise Ratio (dB)
    if 'snr_db' in metrics or metrics == 'all':
        snr_db = [10*np.log10(snr[0]),
                  10*np.log10(snr[1]),
                  10*np.log10(snr[2])]

    # Data truth correlation
    if 'rho' in metrics or metrics == 'all':
        rho = [sens[0]/cov_ab[0][0],
               sens[1]/cov_ab[1][1],
               sens[2]/cov_bc[1][1]]

    for i, k in enumerate(measure_names):
        if 'var' in metrics or metrics == 'all':
            results[k]['var'] = var[i]
        if 'rmse' in metrics or metrics == 'all':
            results[k]['rmse'] = rmse[i]
        if 'si' in metrics or metrics == 'all':
            results[k]['si'] = si[i]
        if 'sensitivity' in metrics or metrics == 'all':
            results[k]['sensitivity'] = sens[i]
        if 'rho' in metrics or metrics == 'all':
            results[k]['rho'] = rho[i]
        if 'snr' in metrics or metrics == 'all':
            results[k]['snr'] = snr[i]
        if 'fmse' in metrics or metrics == 'all':
            results[k]['fmse'] = fmse[i]
        if 'snr_db' in metrics or metrics == 'all':
            results[k]['snr_db'] = snr_db[i]
        if 'mean' in metrics or metrics == 'all':
            results[k]['mean'] = np.mean(measures[i])
        if 'std' in metrics or metrics == 'all':
            results[k]['std'] = np.std(measures[i])

    results = pd.DataFrame(results).transpose().round(dec)
    results.attrs['ref'] = ref
    
    return results


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
