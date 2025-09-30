import pytest
from wavy.insitu_module import insitu_class as ic
from wavy.satellite_module import satellite_class as sc
from wavy.model_module import model_class as mc
import xarray as xr
import wavy.triple_collocation as tc
import numpy as np
import pandas as pd

def test_triple_collocation(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd='2014-01-01',
             ed='2018-12-31',
             nID='history_cmems_NRT',
             name='Norne')
    ico.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_ico.nc")
                              )
    # Import satellite data
    sco = sc(sd='2014-01-01',
             ed='2018-12-31',
             nID='CCIv1_L3',
             name='multi')
    sco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_sco.nc")
                              )
    # Import model data
    mco = mc(sd='2014-01-01',
             ed='2018-12-31',
             nID='NORA3_hc_waves')
    mco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_mco.nc")
                              )
    # Create dictionary for triple collocation function
    dict_data = {'insitu': ico, 'satellite': sco, 'model': mco}
    # Apply triple collocation
    ref = 'insitu'
    tc_res = tc.triple_collocation(dict_data, ref=ref)

    assert isinstance(tc_res, pd.DataFrame)
    assert tc_res.attrs['ref'] in dict_data.keys()
    assert len(tc_res) == 3
    assert len(list(tc_res)) == 6
    
    # Simulated data
    n = 1000
    T = [np.sin(0.2 * i) + 1.2 for i in range(n)]
    b_x = 1.1
    b_y = 0.5
    b_z = 1.8

    np.random.seed(1)
    s_x = 0.1
    e_x = np.random.normal(0, s_x, n)

    np.random.seed(5)
    s_y = 0.2
    e_y = np.random.normal(0, s_y, n)

    np.random.seed(11)
    s_z = 0.5
    e_z = np.random.normal(0, s_z, n)

    X = [b_x * T[i] + e_x[i] for i in range(n)]
    Y = [b_y * T[i] + e_y[i] for i in range(n)]
    Z = [b_z * T[i] + e_z[i] for i in range(n)]

    dict_data = {'X': X, 'Y': Y, 'Z': Z}
    ref = 'X'

    tc_res = tc.triple_collocation(dict_data, ref=ref)

    assert isinstance(tc_res, pd.DataFrame)
    assert tc_res.attrs['ref'] in dict_data.keys()
    assert len(tc_res) == 3
    assert len(list(tc_res)) == 6
    

def test_calibration_triplets_cdf_matching(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd='2014-01-01',
             ed='2018-12-31',
             nID='history_cmems_NRT',
             name='Norne')
    ico.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_ico.nc")
                              )
    # Import satellite data
    sco = sc(sd='2014-01-01',
             ed='2018-12-31',
             nID='CCIv1_L3',
             name='multi')
    sco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_sco.nc")
                              )
    # Import model data
    mco = mc(sd='2014-01-01',
             ed='2018-12-31',
             nID='NORA3_hc_waves')
    mco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_mco.nc")
                              )
    # Create dictionary for triple collocation function
    dict_data = {'insitu': ico.vars.Hs.values, 
                 'satellite': sco.vars.Hs.values,
                 'model': mco.vars.Hs.values}
    # Apply triple collocation
    ref = 'insitu'
    
    data_cal = tc.calibration_triplets_cdf_matching(dict_data, 
                                                    ref=ref, 
                                                    step=0.04)
    
    assert list(dict_data.keys()) == list(data_cal.keys())


def test_calibration_triplets_tc(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd='2014-01-01',
             ed='2018-12-31',
             nID='history_cmems_NRT',
             name='Norne')
    ico.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_ico.nc")
                              )
    # Import satellite data
    sco = sc(sd='2014-01-01',
             ed='2018-12-31',
             nID='CCIv1_L3',
             name='multi')
    sco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_sco.nc")
                              )
    # Import model data
    mco = mc(sd='2014-01-01',
             ed='2018-12-31',
             nID='NORA3_hc_waves')
    mco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_mco.nc")
                              )
    # Create dictionary for triple collocation function
    dict_data = {'insitu': ico.vars.Hs.values, 
                 'satellite': sco.vars.Hs.values,
                 'model': mco.vars.Hs.values}
    # Apply triple collocation
    ref = 'insitu'
    
    data_cal = tc.calibration_triplets_tc(dict_data, ref=ref)
    
    assert list(dict_data.keys()) == list(data_cal.keys())

    
def test_least_squares_merging(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd='2014-01-01',
             ed='2018-12-31',
             nID='history_cmems_NRT',
             name='Norne')
    ico.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_ico.nc")
                              )
    # Import satellite data
    sco = sc(sd='2014-01-01',
             ed='2018-12-31',
             nID='CCIv1_L3',
             name='multi')
    sco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_sco.nc")
                              )
    # Import model data
    mco = mc(sd='2014-01-01',
             ed='2018-12-31',
             nID='NORA3_hc_waves')
    mco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_mco.nc")
                              )
    # Create dictionary for triple collocation function
    dict_data = {'insitu': ico.vars.Hs.values, 
                 'satellite': sco.vars.Hs.values,
                 'model': mco.vars.Hs.values}
    # Apply triple collocation
    ref = 'insitu'
    
    least_squares_merge = tc.least_squares_merging(dict_data)

    assert len(least_squares_merge) == len(dict_data['insitu'])

    
def test_spectra_cco(test_data):

    assert True
    
def test_integrate_r2(test_data):

    assert True
    
def test_filter_collocation_distance(test_data):

    assert True

def test_filter_values(test_data):

    # Wavy objects
    # Import in-situ data
    ico = ic(sd='2014-01-01',
             ed='2018-12-31',
             nID='history_cmems_NRT',
             name='Norne')
    ico.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_ico.nc")
                              )
    # Import satellite data
    sco = sc(sd='2014-01-01',
             ed='2018-12-31',
             nID='CCIv1_L3',
             name='multi')
    sco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_sco.nc")
                              )
    # Import model data
    mco = mc(sd='2014-01-01',
             ed='2018-12-31',
             nID='NORA3_hc_waves')
    mco.vars = xr.open_dataset(
                              str(test_data/"triple_collocation/Norne_mco.nc")
                              )
    # Create dictionary for triple collocation function
    dict_data = {'insitu': ico.vars.Hs.values, 
                 'satellite': sco.vars.Hs.values,
                 'model': mco.vars.Hs.values}
    # Apply triple collocation
    ref = 'insitu'

    data_filtered = tc.filter_values(dict_data, ref_data=ref)

    assert len(data_filtered['insitu']) <= len(dict_data['insitu'])

def test_filter_dynamic_collocation(test_data):

    assert True
