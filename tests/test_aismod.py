import wavy.ais_module as ais
import numpy as np
import pytest


@pytest.mark.need_credentials
def test_get_ais_data():

    bbox = ['5.89', '62.3', '6.5', '62.7']
    sd = '2017-01-03 08'
    ed = '2017-01-03 09'

    ais_ds = ais.get_AIS_data(bbox, sd, ed)
    assert len(ais_ds['time']) > 0
    assert not all(np.isnan(v) for v in ais_ds['time'])
    assert not all(np.isnan(v) for v in ais_ds['lons'])
    assert not all(np.isnan(v) for v in ais_ds['lats'])
