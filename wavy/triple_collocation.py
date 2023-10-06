from wavy.satellite_module import satellite_class as sc
from wavy.insitu_module import insitu_class as ic
from wavy.utils import parse_date, haversineA
from wavy.utils import find_included_times
from wavy.validationmod import validate, disp_validation
import numpy as np
import pandas as pd
import copy

def collocate_sat_and_insitu(sco, ico, twin=5):
    '''
     twin in minutes (as in find_included_times)
     create blocks as lists
     compute distances for all instances
     choose min dist for all instances
     return time series with one distance/value/lat/lon/time
      for each buoy time stamp
    '''
    import datetime
    # create blocks

    dts_src = [datetime.datetime.utcfromtimestamp(t.tolist()/1e9)
               for t in sco.vars['time'].values]
    lats_src = sco.vars['lats'].values
    lons_src = sco.vars['lons'].values
    dts_trgt = [datetime.datetime.utcfromtimestamp(t.tolist()/1e9)
                for t in ico.vars['time'].values]  # ico.vars
    dts_blck = []
    lons_blck = []
    lats_blck = []
    idx_blck = []
    blcks = []
    
    # i is in range of length of in situ data
    for i in range(len(dts_trgt)):
        """
        idx = list of indexes from sat data which falls in the
              time window given by twin (in minutes) in in-situ
              datetimes at index i
        """
        idx = find_included_times(dts_src, dts_trgt[i], twin=twin) 
        # if some indexes are returned
        if len(idx) > 0:
            # list of arrays of matching sat datetime
            dts_blck.append(np.array(dts_src)[idx])
            # list of arrays matching sat latitude
            lats_blck.append(np.array(lats_src)[idx])
            # list of arrays matching sat longitude
            lons_blck.append(np.array(lons_src)[idx])
            # list of list of indices from sat data which were
            # matched mby in situ 
            idx_blck.append(idx)
            # list of indices from in-situ data which matched some sat data
            blcks.append(i)
   
    blck_dict = {
                 # save the obtained matching list of arrays from
                 # sat data into the dictionnary
                 'src_datetime': dts_blck,
                 'src_latitude': lats_blck,
                 'src_longitude': lons_blck,
                 'src_idx': idx_blck,
                 # get datetime, latitude and longitude of in-situ data
                 # which matched and save it
                 'trgt_datetime': [ico.vars['time'].values[b] for b in blcks],
                 'trgt_latitude': [ico.vars['lats'].values[b] for b in blcks],
                 'trgt_longitude': [ico.vars['lons'].values[b] for b in blcks],
                 'trgt_idx': blcks
                }
   
    # compute all distances for all blocks and add to block dict
    # (list of dists added to blck_dict)
    blck_dict['dst'] = []
    blck_dict['min_dst'] = []
    blck_dict['min_dst_idx'] = []
    # For each datetime of in-situ data which matched sat data
    for i in range(len(blck_dict['trgt_datetime'])):
        # lat and lon lists from sat data which matched in-situ data i
        lon1 = blck_dict['src_longitude'][i]
        lat1 = blck_dict['src_latitude'][i]
        # create lists of lat and lon of the point i of in-situ data which
        # matched sat data, repeated at same length as sat lat, lon list.
        lon2 = [blck_dict['trgt_longitude'][i]]*len(lon1)
        lat2 = [blck_dict['trgt_latitude'][i]]*len(lat1)
        # compute distance of every point of sat data matched with the#
        # corresponding in-situ data point i
        dist = list(haversineA(lon1, lat1, lon2, lat2)[0])

        # save the distances lists, the minimum distance each time and the
        # corresponding index inside the sat list matched at point i
        blck_dict['dst'].append(dist)
        blck_dict['min_dst'].append(np.min(dist))
        blck_dict['min_dst_idx'].append(dist.index(np.min(dist)))

    sco_filter = copy.deepcopy(sco)
    ico_filter = copy.deepcopy(ico)

    sco_filter.vars = sco_filter.vars.isel(
                      time=[blck_dict['src_idx'][i]\
                               [blck_dict['min_dst_idx'][i]]\
                               for i in range(len(blck_dict['min_dst_idx']))])

    ico_filter.vars = ico_filter.vars.isel(
                      time=[idx for idx in blck_dict['trgt_idx']])

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
    return covariance(A, A) - (covariance(A, B)*covariance(A, C))/covariance(B, C)

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
    
def triple_collocation_validate(results_dict, ref=None):
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
    validate_dict = {}

    data_sources = results_dict.keys()

    if len(list(data_sources)) != 3:
        raise Exception("Exactly three data sources must be provided.")

    # Check if ref argument is, given and correct
    # if not, takes the first data source as reference
    if ref is None: 
        ref = list(data_sources)[0]
        print("Since no reference was given,", ref, "was taken as default.")

    if ref not in data_sources:
        print("Incorrect argument `ref` value given, should be one of the keys of `tc_validate`.")
        ref = list(data_sources)[0]
        print(ref, "was taken as reference instead.\n")

    # Save data sources into a list
    ds_list = list(results_dict.values())

    # Iterates on key:value from data_sources dict
    for i,k in enumerate(data_sources):

        # Temporary dictionary to save
        # results of current data source
        dict_tmp = {}

        # Rotates the data sources at each step i,
        # so A is the ith element of `results_dict.values()`
        ds_list.insert(0, ds_list.pop(i))
        A = ds_list[0]
        B = ds_list[1]
        C = ds_list[2]

        # Calculates all metrics for A which is, at step i,
        # the ith element of `results_dict.values()`
        dict_tmp["var_est"] = variance_estimates(A, B, C)
        dict_tmp["RMSE"] = RMSE(A, B, C)
        dict_tmp["SI"] = SI(A, B, C, results_dict[ref])
        dict_tmp["sensitivity"] = sensitivity_estimates(A, B, C)
        #dict_tmp["SNR"] = SNR(A,B,C, dict_tmp["var_est"])
        dict_tmp["fMSE"] = fMSE(A, B, C, dict_tmp["var_est"])
        dict_tmp["SNR_dB"] = SNR_dB(A, B, C, dict_tmp["var_est"])

        # Save the metric dictionary into validate_dict
        validate_dict[k] = dict_tmp

    validate_dict = {"data_sources": validate_dict}
    validate_dict["reference_dataset"] = ref
    
    return validate_dict      

def disp_tc_validation(tc_validate, dec=3):
    '''
    Displays results of triple collocation into a table
    '''

    ref = tc_validate["reference_dataset"]

    data_sources_res = tc_validate["data_sources"]
    data_sources_names = data_sources_res.keys()

    metrics_names = data_sources_res[ref].keys()
    
    data = np.array([np.round(list(data_sources_res[k].values()), dec) for k in data_sources_names])
    data = np.transpose(data)
    
    format_row = "{:>12}" * (len(data_sources_names) + 1)
    print(format_row.format("", *data_sources_names))
    for ds, row in zip(metrics_names, data):
        print(format_row.format(ds, *row))
    
    print("\n The reference for the SI is:", ref)

