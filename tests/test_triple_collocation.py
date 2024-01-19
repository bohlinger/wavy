import pytest
from wavy.insitu_module import insitu_class as ic
from wavy.satellite_module import satellite_class as sc
from wavy.model_module import model_class as mc
import xarray as xr
import wavy.triple_collocation as tc
import numpy as np


def test_triple_collocation_simulated_data():

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

    tc_result = tc.triple_collocation_validate(dict_data, ref=ref)

    assert isinstance(tc_result, dict)

    keys_lvl_1 = tc_result.keys()
    assert 'data_sources' in keys_lvl_1
    assert 'reference_dataset' in keys_lvl_1
    assert tc_result['reference_dataset'] in dict_data.keys()

    keys_lvl_2 = tc_result['data_sources'].keys()
    assert set(keys_lvl_2) == set(dict_data.keys())

    assert ref in keys_lvl_2
    assert len(tc_result['data_sources'][ref].keys()) == 6


def test_ci():

    n = 1500
    T = [np.sin(0.2 * i) + 2 for i in range(n)]
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

    ci_95 = tc.bootstrap_ci(dict_data)

    assert isinstance(ci_95, dict)
    assert set(dict_data.keys()) == set(ci_95.keys())
    assert set(['mean', 'ci_l', 'ci_u']) == set(ci_95['X'].keys())


def test_triple_collocation_wavy_objects(test_data):

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
    tc_res = tc.tc_analysis_wavy_objects(dict_data, ref=ref)
    assert isinstance(tc_res, dict)

    keys_lvl_1 = tc_res.keys()
    assert 'data_sources' in keys_lvl_1
    assert 'reference_dataset' in keys_lvl_1
    assert tc_res['reference_dataset'] in dict_data.keys()

    keys_lvl_2 = tc_res['data_sources'].keys()
    assert set(keys_lvl_2) == set(dict_data.keys())

    assert ref in keys_lvl_2
    assert len(tc_res['data_sources'][ref].keys()) == 6
