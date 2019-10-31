import sys
import unittest
from datetime import datetime
sys.path.append(r'../wavy')
import satmod

class TestSatmod(unittest.TestCase):

    def setUp(self):
        self.date = datetime(2017,8,1,12)
        self.timewin = 30
        self.distlim = 6
        self.region = 'ARCMFC'

    def test_get_remotefiles(self):
        import os
        import time
        from pathfinder import satpath_ftp_014_001, satpath_lustre
        satmod.get_remotefiles(satpath_ftp_014_001,satpath_lustre,
                    self.date,self.date,timewin=self.timewin,
                    corenum=8,download=True)
        creation_date = os.path.getctime(satpath_lustre 
                        + '2017/08/'
                        + 'global_vavh_l3_rt_s3a_C0020_P0583_'
                        + '20170801T072304_'
                        + '20170801T080956_20170801T111812.nc')
        now = time.time()
        limit = now - 60 # not older than 60 seconds
        self.assertLessEqual(limit,creation_date)

    def test_sa_obj(self):
        from datetime import datetime
        from satmod import sentinel_altimeter as sa
        sa_obj = sa(self.date,timewin=self.timewin,region=self.region)
        self.assertEqual(sa_obj.__class__.__name__,'sentinel_altimeter')
        sa_obj = sa(self.date,timewin=self.timewin,region=self.region,
                    mode="ARCMFC")
        self.assertEqual(sa_obj.__class__.__name__,'sentinel_altimeter')



if __name__ == '__main__':
    unittest.main()
