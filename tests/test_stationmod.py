from datetime import datetime, timedelta
import numpy as np

from wavy.stationmod import station_class as sc

sdate = datetime(2017, 11, 1)
edate = datetime(2017, 12, 1)
station = 'ekofiskL'
leadtime = 6


def test_sc_obj():
    sc_obj = sc(station, '1', sdate, edate)
    assert sc_obj.__class__.__name__ == 'station_class'
