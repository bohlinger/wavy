import sys
import unittest
from datetime import datetime, timedelta
import numpy as np
sys.path.append(r'/home/patrikb/wavy/wavy')
from stationmod import station_class as sc
from stationmod import get_model, matchtime
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
        fc_date = self.sdate + timedelta(hours=self.leadtime)
        #mwam4
        mtime,mhs,mlons,mlats = \
                    get_model('mwam4',fc_date,leadtime=self.leadtime)
        ctime, cidx = \
            matchtime(fc_date,fc_date,sc_obj.time,sc_obj.basedate)
        idx,idy,distM,mlats_new,mlons_new = \
            get_loc_idx(mlats,mlons,sc_obj.lat,sc_obj.lon,mhs)
        self.assertEqual(mhs.__class__.__name__,'MaskedArray')
        self.assertEqual(np.isnan(np.nanmean(mhs)),False)
        self.assertGreater(len(idx),0)
        self.assertEqual(np.isnan(idx[0]),False)

    def test_get_model(self):
        # ARCMFC
        #20171104_MyWaveWam8r625_b20171101.nc
        fc_date = datetime(2017,11,4)
        init_date = datetime(2017,11,1)
        mtime,mhs,mlons,mlats = \
                    get_model('ARCMFC',fc_date,init_date=init_date)
        self.assertEqual(mhs.__class__.__name__,'MaskedArray')
        self.assertEqual(np.isnan(np.nanmean(mhs)),False)
        

if __name__ == '__main__':
    unittest.main()
