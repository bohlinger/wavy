import sys
import os
import unittest
from datetime import datetime
sys.path.append(r'../wavy')
import satmod
from satmod import satellite_class
import yaml

class TestSatmod(unittest.TestCase):

    def setUp(self):
        now = datetime.now()
        self.date = datetime(now.year,now.month,now.day)
        self.region = 'NordicSeas'
        self.sat = 's3a'
        self.varalias = 'Hs'
        self.twin = 60
        self.nproc = 2
        self.instr = 'altimeter'
        self.provider = 'cmems'
        self.path_local = 'tmp_unittest/'

        # create tmp directory in tests
        cmd = 'mkdir ' + self.path_local
        os.system(cmd)
        # backup of existing config files if exist
        cmd = ( 'cp ../config/*.yaml ' + self.path_local + '.' )
        os.system(cmd)
        # set relevant config files to default
        cmd = ( 'cp ../config/region_specs.yaml.default ' +
                '../config/region_specs.yaml')
        os.system(cmd)
        cmd = ( 'cp ../config/model_specs.yaml.default ' +
                '../config/model_specs.yaml')
        os.system(cmd)
        cmd = ( 'cp ../config/satellite_specs.yaml.default ' +
                '../config/satellite_specs.yaml')
        os.system(cmd)
        cmd = ( 'cp ../config/variable_info.yaml.default ' +
                '../config/variable_info.yaml')
        os.system(cmd)
        # parsing yaml files and make available as dict
        moddir = os.path.abspath(os.path.join(
                            os.path.dirname( __file__ ),
                            '..', 'config/region_specs.yaml'))
        with open(moddir,'r') as stream:
            region_dict=yaml.safe_load(stream)

        moddir = os.path.abspath(os.path.join(
                            os.path.dirname( __file__ ),
                            '..', 'config/model_specs.yaml'))
        with open(moddir,'r') as stream:
            model_dict=yaml.safe_load(stream)

        moddir = os.path.abspath(os.path.join(
                            os.path.dirname( __file__ ),
                            '..', 'config/satellite_specs.yaml'))
        with open(moddir,'r') as stream:
            satellite_dict=yaml.safe_load(stream)

        moddir = os.path.abspath(os.path.join(
                            os.path.dirname( __file__ ),
                            '..', 'config/variable_info.yaml'))
        with open(moddir,'r') as stream:
            variable_info=yaml.safe_load(stream)

        self.region_dict = region_dict
        self.satellite_dict = satellite_dict
        self.model_dict = model_dict
        self.variable_info = variable_info

    def test_ftp_files_and_init_satellite_class(self):
        # make satpath
        path_remote = self.satellite_dict[self.instr]\
                        [self.provider]['remote']['path']\
                        + self.sat + '/'
        # evoke fct get_remote_files
        satmod.get_remote_files(path_remote,self.path_local,
                                self.date,self.date,
                                self.twin,self.nproc,
                                self.instr,self.provider)
        # check if file were download to tmp directory
        filelist = os.listdir(self.path_local)
        nclist = [i for i in range(len(filelist))\
                  if '.nc' in filelist[i]]
        self.assertGreaterEqual(len(nclist),1)
        # init satellite_object and check for polygon region
        from satmod import satellite_class
        sa_obj = satellite_class(sdate=self.date,region=self.region,
                                sat=self.sat,varalias=self.varalias,
                                download_path=self.path_local)
        self.assertEqual(sa_obj.__class__.__name__,'satellite_class')
        self.assertGreaterEqual(len(vars(sa_obj).keys()),11)
        self.assertGreaterEqual(len(sa_obj.vars.keys()),6)
        self.assertEqual('error' in vars(sa_obj).keys(),False)
        # another one checking for error

    def tearDown(self):
        # restore config files
        cmd = ( 'cp ' + self.path_local + '*.yaml ../config/.' )
        os.system(cmd)
        # rm tmp directory
        print('cleaning after test')
        cmd = 'rm -r ' + self.path_local
        os.system(cmd)

if __name__ == '__main__':
    unittest.main()
