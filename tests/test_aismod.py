import wavy.ais_module as ais
import numpy as np
import pytest


@pytest.mark.need_credentials
def test_get_ais_data():

    bbound = ['5.89', '62.3', '6.5', '62.7']
    sd = '201701030800'
    ed = '201701030900'

    ais_ds = ais.get_AIS_data(bbound, sd, ed)
    assert len(ais_ds['time']) > 0
    assert not all(np.isnan(v) for v in ais_ds['time'])
    assert not all(np.isnan(v) for v in ais_ds['lons'])
    assert not all(np.isnan(v) for v in ais_ds['lats'])
