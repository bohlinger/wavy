import pytest
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
    assert len(tc_result['data_sources'][ref].keys()) == 8


# def test_triple_collocation_real_data():

#    # use "real" data
#    # continue as above
