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
import xarray as xr


def filter_collocation_distance(data, dist_max, name):
    """ 
    Filters the datasets according to a maximum collocation
    distance between satellite and in-situ.

    data (dict of wavy objects): wavy objects to filter given
                                 the maximum collocation distance.
                                 One of the objects must contain
                                 the collocation distance.
    dist_max (float): Maximum collocation distance in km.
    name (string): key from the dictionary that refers to the
                   wavy object containing the distance

    returns:
    data_filtered (dict of wavy objects): dictionary of the wavy objects 
                                 filtered using the maximum collocation 
                                 distance given.
    """
    data_filtered = {}

    dist_data = data[name].vars.colloc_dist.values
    
    idx_dist = (dist_data <= dist_max)

    for k in data.keys():

        wavy_obj_tmp = copy.deepcopy(data[k])
        wavy_obj_tmp.vars = wavy_obj_tmp.vars.where(idx_dist).\
                                              dropna(dim='time')

        data_filtered[k] = wavy_obj_tmp
        
    return data_filtered


def filter_values(data, ref_data, min=0.0, max=25.0, return_ref_data=False):
    '''
    Filters the values for each data serie given as input. 

    data (dict of arrays): data to filter
    ref_data (string or array): Either a string corresponding
                    to a key in data or an array. Values for all
                    data are filtered with respect to the ref_data.
    min (float): minimum value that ref_data should take. 
    max (float): maximum value that ref_data should take.
    '''

    if isinstance(ref_data, str):
        ref_data = data[ref_data]
    
    idx = (ref_data >= min) & (ref_data < max)

    data_filtered = {}

    for k in data.keys():
        data_filtered[k] = data[k][idx]
    
    if return_ref_data == False:
        return data_filtered
    else:
        ref_data_filtered=ref_data[idx]
        return data_filtered, ref_data_filtered


def filter_dynamic_collocation(data, mod_1, mod_2, max_rel_diff=0.05): 
    '''
    Filter data when the two given model data differ by more than a given
    percentage. Dynamical collocation filtering method for collocation
    refers to Dodet et al., 2025.

    data (dict of lists): data to filter 
    mod_1 (string or list): Either key from data for the first model data
                            or the list of values of the model directly
    mod_2 (string or list): Either key from data for the first model data
                            or the list of values of the model directly
    max_rel_diff (float): Maximum relative difference (abs(mod_1-mod_2)/mod_1)
                          allowed between values from mod_1 and mod_2.

    returns:
    data_filtered (dict of lists): filtered data
    '''

    if isinstance(mod_1, str):
        mod_1 = data[mod_1]

    if isinstance(mod_2, str):
        mod_2 = data[mod_2]

    mod_1 = np.array(mod_1)
    mod_2 = np.array(mod_2)
    
    idx = np.abs(mod_1 - mod_2)/mod_1 < 0.05

    data_filtered = {}

    for k in data.keys():
        data_filtered[k] = data[k][idx]

    return data_filtered


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


def get_CDF(data, step, 
            llim=None, ulim=None, 
            data_min=None, data_max=None, 
            dec=3, no_empty_bins=True):

      N = len(data)
      if data_max == None:
          data_max = np.ceil(np.max(data)) + 2 
      if data_min == None:
          data_min = np.floor(np.min(data)) - 2

      if llim==None and ulim==None:
          bins = np.arange(data_min,data_max+step,step)
      elif llim==None and ulim!=None:
          bins = np.concatenate([np.arange(data_min,ulim,step),
                                np.array([data_max])])
      elif llim!=None and ulim==None:
          bins = np.concatenate([np.array([data_min]), 
                                np.arange(llim,data_max+step,step)])
      else:
          bins = np.concatenate([np.array([data_min]), 
                                 np.arange(llim,ulim,step), 
                                 np.array([data_max])])
   
      count_bins = [np.sum((data > bins[i]) &\
                           (data <= bins[i+1])) for i in range(len(bins)-1)]  

      idx_null = [i for i, v in enumerate(count_bins) if v == 0]

      if no_empty_bins==True:
          if len(idx_null) > 0:
              list_bins_null = ['({},{}]'.format(round(bins[i],dec), 
                                                 round(bins[i+1],dec))\
                                                 for i in idx_null]
              print('Warning: bins ' + ', '.join(list_bins_null) +\
                    ' do not have any data.') 
    
              if idx_null[0]!=0:
                  print("Invalid CDF, contains empty bins!")
                  return None
      
      CDF = [np.sum(count_bins[0:i])/N for i in range(1,len(count_bins)+1)]

      df_CDF = pd.DataFrame({'lower bound':bins[:-1], 
                             'upper bound':bins[1:], 
                             'CDF':CDF, 
                             'count':count_bins})     

      return df_CDF

def CDF_matching_cal(old, CDF_old, CDF_new):

    new = []
    min_bound = np.min(CDF_new['lower bound'])
    
    for x_tmp in old: 

        row_old_tmp = CDF_old[(x_tmp > CDF_old['lower bound']) &\
                              (x_tmp <= CDF_old['upper bound'])].iloc[0,:]
        b_i = row_old_tmp['upper bound']
        b_i_1 = row_old_tmp['lower bound']
        a_x = (x_tmp - b_i_1)/(b_i - b_i_1)
        C_b_i = row_old_tmp['CDF']        
        C_b_i_1 = CDF_old[CDF_old['upper bound']==b_i_1]['CDF'].values[0]
        C_x = a_x * C_b_i_1 + (1 - a_x) * C_b_i
                
        row_new_tmp = CDF_new[CDF_new['CDF'] >= C_x].iloc[0,:]
        b_k = row_new_tmp['upper bound']
        b_k_1 = row_new_tmp['lower bound']
        C_tild_b_k = row_new_tmp['CDF']

        if b_k_1 <= min_bound:
            x_tild = min_bound
        else:
            C_tild_b_k_1 =\
                     CDF_new[CDF_new['upper bound'] == b_k_1].iloc[0,:]['CDF']
            a_C_x =  (C_x - C_tild_b_k_1)/(C_tild_b_k - C_tild_b_k_1)
            x_tild = a_C_x * b_k_1 + (1 - a_C_x) * b_k  

        new.append(x_tild)

    return np.array(new)
    

def calibration_triplets_cdf_matching(data, ref, step, seed=5):
    
    measure_names = list(data.keys())

    size = len(data[measure_names[0]]) 
    
    tc_res = triple_collocation(data, metrics=['var'], ref=ref)
    
    np.random.seed(seed)
    
    max_var = np.max(tc_res['var'])
    
    data_err = {}
    
    for k in data.keys():
    
        if tc_res.loc[k,'var'] == max_var:
            data_err[k] = copy.copy(data[k])
        else: 
            e_diff = np.random.normal(0,np.sqrt(max_var - tc_res.loc[k,'var']),
                                      size)
            data_err[k] =  copy.copy(data[k]) + e_diff
        
    data_max = np.ceil(np.max([np.max(d) for d in data_err.values()])) + 1
    data_min = np.ceil(np.min([np.min(d) for d in data_err.values()])) - 1
    
    CDF_dict = {}
    
    max_low_idx = 0
    min_up_idx = np.inf
    
    for k in data.keys():
    
        CDF_tmp = get_CDF(data_err[k], 
                          data_min=data_min,
                          data_max=data_max,
                          step=step, 
                          dec=3,
                          no_empty_bins=False)
        CDF_dict[k] = CDF_tmp
    
        list_null_count = [idx for idx in range(len(CDF_tmp)) if\
                           CDF_tmp['count'][idx]==0]
        
        list_low_idx = [j for j in list_null_count if j < len(CDF_tmp)/2]
        if len(list_low_idx) > 0:
            max_low_idx_tmp = np.max(list_low_idx)
        else: 
            max_low_idx_tmp = 0
    
        list_up_idx = [j for j in list_null_count if j > len(CDF_tmp)/2]
        if len(list_up_idx) > 0:
            min_up_idx_tmp = np.min(list_up_idx)
        else:
            min_up_idx_tmp = len(CDF_tmp) - 1
        
        if max_low_idx_tmp > max_low_idx:
             max_low_idx = max_low_idx_tmp
        
        if min_up_idx_tmp < min_up_idx:
             min_up_idx = min_up_idx_tmp 
    
    for k in CDF_dict.keys():
        
        CDF_tmp = CDF_dict[k] 
          
        nb_values = CDF_tmp['count'].sum()
        
        core_cdf = CDF_tmp.iloc[max_low_idx+1:min_up_idx, :]
        
        lower_bin_cdf = CDF_tmp.iloc[:max_low_idx+1, :]
        
        lower_bin_cdf_df = pd.DataFrame({
                           'lower bound':np.min(lower_bin_cdf['lower bound']), 
                           'upper bound':np.max(lower_bin_cdf['upper bound']), 
                           'CDF': np.sum(lower_bin_cdf['count'])/nb_values, 
                           'count':np.sum(lower_bin_cdf['count'])
                           },
                           index=[0])
        
        upper_bin_cdf = CDF_tmp.iloc[min_up_idx:, :]
        
        upper_bin_cdf_df = pd.DataFrame({
                           'lower bound':np.min(upper_bin_cdf['lower bound']), 
                           'upper bound':np.max(upper_bin_cdf['upper bound']), 
                           'CDF': CDF_tmp['CDF'][min_up_idx-1]+\
                                  np.sum(upper_bin_cdf['count'])/nb_values, 
                           'count':np.sum(upper_bin_cdf['count'])
                           },
                           index=[min_up_idx])
        
        final_cdf = pd.concat([lower_bin_cdf_df, core_cdf, upper_bin_cdf_df], 
                              ignore_index=True)
    
        CDF_dict[k] = final_cdf
    
    max_matching = final_cdf.iloc[-1]['lower bound'] 
    min_matching = final_cdf.iloc[0]['upper bound']

    print("Upper limit for CDF matching: ", round(max_matching,3))
    print("Lower limit for CDF matching: ", round(min_matching,3))
    
    idx_to_cal = (data[measure_names[0]] < max_matching) &\
                 (data[measure_names[1]] < max_matching) &\
                 (data[measure_names[2]] < max_matching) &\
                 (data[measure_names[0]] > min_matching) &\
                 (data[measure_names[1]] > min_matching) &\
                 (data[measure_names[2]] > min_matching)
    
    data_to_cal = {measure_names[i]:data[measure_names[i]][idx_to_cal] for\
                   i in range(3)}
    
    data_cal = {}
    
    for name in measure_names: 
    
         if name == ref: 
             data_cal[name] = data_to_cal[ref]
    
         else: 
             data_cal[name] = CDF_matching_cal(data_to_cal[name], 
                                               CDF_dict[name], 
                                               CDF_dict[ref])
    
    # Remove values that would have exceeded the upper and lower 
    # bounds after calibration
    idx_cal = (data_cal[measure_names[0]] < max_matching) &\
              (data_cal[measure_names[1]] < max_matching) &\
              (data_cal[measure_names[2]] <= max_matching) &\
              (data_cal[measure_names[0]] > min_matching) &\
              (data_cal[measure_names[1]] > min_matching) &\
              (data_cal[measure_names[2]] > min_matching)
    
    data_cal_final = {name:data_cal[name][idx_cal] for name in data_cal.keys()}

    return data_cal_final
    

def calibration_triplets_tc(data, ref, r2=0, return_cal_cst=False):
    '''
    Calibrate A and B relatively to R using triple collocation calibration
    constant estimates, following Gruber et al., 2016 method.
    
    data (dict of lists of floats): Dictionary of the data to calibrate.
    ref (string): Name of the reference data to use for calibration.
    r2 (float): Representativeness error
    cal_cst (bool): If True, returns a dictionary for the calibration
                    constantes in addition to the calibrated data. 

    returns:
    data_cal (dict of lists of floats): Dictionary of the calibrated 
                                        data series
    '''
    R = data[ref]

    measure_names = list(data.keys())

    if ref == measure_names[0]:
        A = data[measure_names[1]]
        B = data[measure_names[2]]
        idx_R = 0
        idx_A = 1
        idx_B = 2
    elif ref == measure_names[1]:
        A = data[measure_names[0]]
        B = data[measure_names[2]]
        idx_R = 1
        idx_A = 0
        idx_B = 2
    elif ref == measure_names[2]:
        A = data[measure_names[0]]
        B = data[measure_names[1]]
        idx_R = 2
        idx_A = 0
        idx_B = 1
    else:
        print("Invalid reference. {} does not appear\
               in the keys of the input data dictionary.")
   
    c_AB = np.cov(A, B)[0, 1] 
    c_RA = np.cov(R, A)[0, 1]
    c_RB = np.cov(R, B)[0, 1]

    a_A = c_AB/c_RB
    a_B = c_AB/(c_RA - a_A*r2)
    a_R = 1.0

    A_R = (a_R/a_A)*(A - np.mean(A)) + np.mean(R)
    B_R = (a_R/a_B)*(B - np.mean(B)) + np.mean(R)

    res = {idx_R:R, idx_A:A_R, idx_B:B_R}
    res_cst = {idx_R:a_R, idx_A:a_A, idx_B:a_B}
    
    data_cal = {measure_names[i]:res[i] for i in range(3)}
    cal_cst = {measure_names[i]:res_cst[i] for i in range(3)}

    if return_cal_cst==False:
        return data_cal
    else:
        return data_cal, cal_cst


def least_squares_merging(data, tc_results=None, return_var=False, **kwargs):
    '''
    Merges the three data series given as input following the least 
    squares merging method described in Yilmaz et al., 2012. 

    data (dict of lists of floats): Dictionary of the data to calibrate.
    tc_results (pandas DataFrame): table of the results of triple collocation
               for the given data. Must contain the variance. If None, the 
               triple collocation is performed using the data and kwargs given.
    return_var (bool): If True, returns the variance of the error of the merged 
               data in addition to the merged data.

    returns:
    least_squares_merge (numpy array): series of merged data
    least_squares_var (float): variance of error of the merged data    
    '''
    measure_names = list(data.keys())
    data_0 = data[measure_names[0]]
    data_1 = data[measure_names[1]]
    data_2 = data[measure_names[2]]

    if tc_results is None:
        r2 = kwargs.get('r2', 0)
        ref = kwargs.get('ref', None)
        dec = kwargs.get('dec', 3)
        tc_results = triple_collocation(data, 
                                        metrics=['var'],
                                        r2=r2,
                                        ref=ref,
                                        dec=dec)

    s_0_sq = tc_results['var'][measure_names[0]]
    s_1_sq = tc_results['var'][measure_names[1]]
    s_2_sq = tc_results['var'][measure_names[2]]
    
    s_12 = s_1_sq * s_2_sq
    s_02 = s_0_sq * s_2_sq
    s_01 = s_0_sq * s_1_sq
    w_0 = s_12/(s_01+s_02+s_12)
    w_1 = s_02/(s_01+s_02+s_12)
    w_2 = s_01/(s_01+s_02+s_12)
    least_squares_merge = np.array(w_0 * data_0 + w_1 * data_1 + w_2 * data_2)
    least_squares_var = w_0**2 * s_0_sq + w_1**2 * s_1_sq + w_2**2 * s_2_sq

    if return_var==False:
        return least_squares_merge
    else: 
        return least_squares_merge, least_squares_var


def spectra_cco(cco, fs, returns='average',nsample=64):
    """
    Calculate the power spectra for both model and observation series 
    of a collocation class object. The series are separated into 
    samples of a given length. Up to one missing value for each
    sample is filled with interpolation. A Hamming window is applied
    for each sample. 

    cco (collocation_class object or xarray dataset): collocation class object 
              for which the spectra are calculated. If xarray.dataset 
              is given, it must have a single time dimension and obs_values 
              and model_values variables.
    fs (float): frequency of the time series
    returns (str): If 'average' returns the average power spectra of the time series. 
                   If 'list' returns lists of power spectra of each samples of the given series.
    nsample (int): size of the samples to calculate the power spectra. 

    return:
    PS_mod (numpy array): average spectra of cco model values or individual spectra for each sample
    PS_obs (numpy array): average spectra of cco observation values or individual spectra for each sample
    f (numpy array): frequencies for the spectra
    """
    from scipy.signal import periodogram
    from wavy.collocation_module import collocation_class as cc
    
    if isinstance(cco, cc): 
        cco_vars = cco.vars.dropna(dim='time')
    elif isinstance(cco, xr.Dataset):
        cco_vars = cco.dropna(dim='time')
    
    step_cco = [(cco_vars.time.values[i+1] - cco_vars.time.values[i])/np.timedelta64(1,'s') for\
                i in range(len(cco_vars.time)-1)]
    step_cco.append(np.nan)
    median_step = np.median(step_cco)
        
    i_0 = 0
    i_1 = i_0 + nsample
    
    list_PS_obs = []
    list_PS_mod = []
    
    while i_1 < len(cco_vars.time.values): 
    
        sample_tmp = cco_vars.isel(time=range(i_0,i_1))
    
        step_tmp = [(sample_tmp.time.values[i+1] - sample_tmp.time.values[i])/np.timedelta64(1,'s') for\
                    i in range(len(sample_tmp.time)-1)]
        step_tmp.append(np.nan)  
        
        idx_1_na = np.argwhere((step_tmp > 1.5*median_step) & (step_tmp < 2.5*median_step))
        idx_2_na = np.argwhere(step_tmp > 2.5*median_step)
    
        if len(idx_2_na) > 0:
            i_0 = idx_2_na[-1][0]
            i_1 = i_0 + nsample
            continue
    
        for idx in idx_1_na:
            idx = idx[0]
        
            time_idx = sample_tmp.time.values[idx]
            obs_idx = sample_tmp.obs_values.values[idx]
            mod_idx = sample_tmp.model_values.values[idx]
            time_idx_1 = sample_tmp.time.values[idx+1]
            obs_idx_1 = sample_tmp.obs_values.values[idx+1]
            mod_idx_1 = sample_tmp.model_values.values[idx+1]
            
            time_between = time_idx + np.timedelta64(int(1000*median_step),'ms')
            obs_between= (obs_idx + obs_idx_1)/2
            mod_between = (mod_idx + mod_idx_1)/2 
        
            ds_between = xr.Dataset(
                                data_vars={'obs_values': (('time'), [obs_between]),
                                           'model_values': (('time'), [mod_between])},
                                coords={'time': [time_between]}
                                   ) 
    
            sample_tmp = xr.concat([sample_tmp, ds_between], dim='time').sortby('time')
    
        sample_tmp = sample_tmp.isel(time=range(0,nsample))
        
        obs_val_tmp = sample_tmp.obs_values.values
        mod_val_tmp = sample_tmp.model_values.values
        
        f, PS_obs_tmp =  periodogram(obs_val_tmp, fs=fs, window='hamming')
        f, PS_mod_tmp =  periodogram(mod_val_tmp, fs=fs, window='hamming')
    
        list_PS_obs.append(PS_obs_tmp)
        list_PS_mod.append(PS_mod_tmp)
    
        i_0 = i_1 - len(idx_1_na)
        i_1 = i_0 + nsample
    
    mean_PS_mod = np.mean(np.array(list_PS_mod), axis=0)
    mean_PS_obs = np.mean(np.array(list_PS_obs), axis=0)
    
    if returns == 'list':
        PS_mod = np.array(list_PS_mod)
        PS_obs = np.array(list_PS_obs)
    elif returns == 'average':
        PS_mod = mean_PS_mod
        PS_obs = mean_PS_obs
        
    return PS_mod, PS_obs, f
    
    
def integrate_r2(PS_mod, PS_obs, f, threshold=np.inf, threshold_type='inv_freq'):
    """
    Estimates the representativeness error r2 by integrating the difference
    between the average power spectra of the model and the observations. 
    Integrates up to a given threshold resolution (inverse of the frequency), 
    or from a given threshold frequency. 
    s
    PS_mod (numpy array): power spectra for the model
    PS_obs (numpy array): power spectra for the observations
    f (numpy array): frequencies of the power spectra
    threshold (float): either upper resolution to integrate to, or minimum 
                       frequency to integrate from
    threshold_type (str): indicate whether the threshold corresponds to a minimum 
                          frequency ('freq') or to a maximum resolution ('inv_freq') 

    return:
    r2 (float): representativeness error estimate
    """
    freq_step = f[1] - f[0]
    diff_PS = PS_obs - PS_mod
    weighted_diff_PS = freq_step*diff_PS

    if threshold_type == 'freq':
        if threshold==0:
            threshold=np.inf
        else:
            threshold = 1/threshold    
    

    f_1 = [1/f[i] if f[i] != 0 else np.inf for i in range(len(f)) ]
    idx_threshold = np.argwhere(np.array(f_1) <= threshold)[0][0]
    r2 = np.sum(weighted_diff_PS[idx_threshold:])
    
    return r2
