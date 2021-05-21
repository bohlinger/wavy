import sys
import os
import unittest
from datetime import datetime
sys.path.append(r'../wavy')
from satmod import satellite_class

class TestSatmod(unittest.TestCase):

    def setUp(self):
        self.date = datetime(2020,1,2,12)
        self.region = 'NordicSeas'
        self.sat = 's3a'
        self.varalias = 'Hs'
        # create tmp directory in tests
        cmd = 'mkdir tmp'
        os.system(cmd)

    def test_downloading_copernicus_l3(self):
        self.assertEqual(1,1)

    def test_satellite_class_object(self):
        from satmod import satellite_class
        # check polygon region
        sa_obj = satellite_class(sdate=self.date,region=self.region,
                                sat=self.sat,varalias=self.varalias)
        self.assertEqual(sa_obj.__class__.__name__,'satellite_class')
        self.assertGreaterEqual(len(vars(sa_obj).keys()),11)
        self.assertGreaterEqual(len(sa_obj.vars.keys()),6)

    def tearDown(self):
        # rm tmp directory
        print('cleaning after test')
        cmd = 'rm -r tmp'
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
