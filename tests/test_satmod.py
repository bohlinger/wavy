import sys
import os
from datetime import datetime
import pytest

from wavy.wconfig import load_or_default
import wavy.sat_collectors
from wavy.satmod import satellite_class as sc

satellite_dict = load_or_default('satellite_specs.yaml')

#def test_ftp_files_and_satellite_class_features(tmpdir):
@pytest.mark.need_credentials
def test_collectors_cmems_L3(tmpdir):
    sdate = "2023-2-1 12"
    edate = "2023-2-1 12"
    sdate_dt = datetime(2023, 2, 1, 12)
    edate_dt = datetime(2023, 2, 1, 12)
    twin = 30
    nproc = 1
    mission = 's3a'
    product = 'cmems_L3_NRT'
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission': mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt, edate=edate_dt,
                            twin=twin, nproc=nproc,
                            product=product, api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
              if '.nc' in filelist[i]]
    assert len(nclist) >= 1

@pytest.mark.need_credentials
def test_collectors_cmems_L3_s6a(tmpdir):
    sdate = "2023-2-1 12"
    edate = "2023-2-1 12"
    sdate_dt = datetime(2023,2,1,12)
    edate_dt = datetime(2023,2,1,12)
    twin = 30
    nproc = 1
    mission = 's6a'
    product = 'cmems_L3_s6a'
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt,edate=edate_dt,
                            twin=twin,nproc=nproc,
                            product=product,api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1

@pytest.mark.need_credentials
def test_collectors_cmems_MY(tmpdir):
    sdate = "2018-1-2 12"
    edate = "2018-1-2 12"
    sdate_dt = datetime(2018,1,2,12)
    edate_dt = datetime(2018,1,2,12)
    twin = 30
    nproc = 1
    mission = 'j3'
    product = 'cmems_L3_MY'
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt,edate=edate_dt,
                            twin=twin,nproc=nproc,
                            product=product,api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1

@pytest.mark.need_credentials
def test_collectors_cmems_cci_L2P(tmpdir):
    sdate = "2018-1-2 12"
    edate = "2018-1-2 12"
    sdate_dt = datetime(2018,1,2,12)
    edate_dt = datetime(2018,1,2,12)
    twin = 30
    nproc = 1
    mission = 'j3'
    product = 'cci_L2P'
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt,edate=edate_dt,
                            twin=twin,nproc=nproc,
                            product=product,api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1

@pytest.mark.need_credentials
def test_collectors_reader_cmems_cci_L3(tmpdir):
    sdate = "2018-12-1 01"
    edate = "2018-12-1 23"
    sdate_dt = datetime(2018,12,1,1)
    edate_dt = datetime(2018,12,1,23)
    twin = 30
    nproc = 1
    mission = 'multi'
    product = 'cci_L3'
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt,edate=edate_dt,
                            twin=twin,nproc=nproc,
                            product=product,api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1
    region = 'NorwegianSea'
    varalias = 'Hs'
    sco = sc(sdate=sdate,edate=edate,region=region,
             mission=mission,twin=twin,varalias=varalias,
             product=product,path_local=tmpdir)
    assert sco.__class__.__name__ == 'satellite_class'
    assert len(vars(sco).keys()) >= 11
    assert len(sco.vars.keys()) >= 6
    assert not 'error' in vars(sco).keys()

@pytest.mark.need_credentials
def test_collectors_aviso_cfo(tmpdir):
    sdate_dt = datetime(2022,3,10,10)
    edate_dt = datetime(2022,3,10,10)
    twin = 30
    nproc = 1
    mission = 'cfo'
    product = 'cfo_swim_L2P'
    # evoke fct get_remote_files
    api_url = None
    dict_for_sub = {'mission':mission}
    wavy.sat_collectors.get_remote_files(
                            path_local=tmpdir,
                            sdate=sdate_dt,edate=edate_dt,
                            twin=twin,nproc=nproc,
                            product=product,api_url=api_url,
                            mission=mission,
                            dict_for_sub=dict_for_sub)
    # check if file were download to tmp directory
    filelist = os.listdir(tmpdir)
    nclist = [i for i in range(len(filelist))\
                if '.nc' in filelist[i]]
    assert len(nclist) >= 1

def test_readers_aviso_cfo(test_data):
    # convertion to Tp
    sco = sc(sdate="2022-2-26",edate="2022-2-27",
            mission="cfo",product="cfo_swim_L2P",varalias='Lp',return_var='Tp',
            path_local=str(test_data/"SWIM_L2P"))
    assert len(vars(sco).keys()) >= 11

def test_readers_and_write_to_nc_cmems_NRT(tmpdir,test_data):
    # evoke fct get_remote_files
    sdate = "2020-11-1 12"
    edate = "2020-11-1 12"
    region = 'NordicSeas'
    mission = 's3a'
    varalias = 'Hs'
    twin = 30
    nproc = 1
    product = 'cmems_L3_NRT'
    # init satellite_object and check for polygon region
    sco = sc(sdate=sdate,edate=edate,region=region,
             mission=mission,twin=twin,varalias=varalias,
             product=product,path_local=str(test_data/"L3"))
    assert sco.__class__.__name__ == 'satellite_class'
    assert len(vars(sco).keys()) >= 11
    assert len(sco.vars.keys()) >= 6
    assert not 'error' in vars(sco).keys()
    # write to nc
    sco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
    # check if created -> assert
    # read nc
    # check if varalias assert
