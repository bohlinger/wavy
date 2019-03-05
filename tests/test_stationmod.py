import sys
import unittest
from datetime import datetime, timedelta
import numpy as np
sys.path.append(r'/home/patrikb/wavy/wavy')
from stationmod import station_class as sc
from stationmod import matchtime
from stationmod import get_loc_idx

class TestStationmod(unittest.TestCase):

    def setUp(self):
        self.sdate = datetime(2017,11,1)
        self.edate = datetime(2017,12,1)
        self.station = 'ekofiskL'
        self.leadtime = 6

    def test_sc_obj(self):
        sc_obj = sc(self.station,self.sdate,self.edate)
        self.assertEqual(sc_obj.__class__.__name__,'station_class')


if __name__ == '__main__':
    unittest.main()
