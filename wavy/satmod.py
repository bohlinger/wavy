#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to acquire, read, and prepare
geophysical variables related to wave height from satellite
altimetry files for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import sys
import numpy as np
from datetime import datetime, timedelta
import os
from copy import deepcopy
import time
from urllib.request import urlretrieve, urlcleanup # python3
from ftplib import FTP
import netCDF4 as netCDF4
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed
import pandas as pd
import pyproj
import zipfile
import tempfile
from tqdm import tqdm
import xarray as xr

# own imports
from wavy.ncmod import ncdumpMeta
from wavy.ncmod import read_netcdfs, get_filevarname
from wavy.ncmod import find_attr_in_nc, dumptonc_ts_sat
from wavy.utils import find_included_times, NoStdStreams
from wavy.utils import sort_files, parse_date
from wavy.utils import make_pathtofile, make_subdict
from wavy.utils import finditem, haversineA
from wavy.credentials import get_credentials
from wavy.modelmod import make_model_filename_wrapper
from wavy.modelmod import read_model_nc_output_lru
from wavy.wconfig import load_or_default
from wavy.filtermod import filter_main,vardict_unique
from wavy.filtermod import rm_nan_from_vardict
# ---------------------------------------------------------------------#

# read yaml config files:
region_dict = load_or_default('region_specs.yaml')
model_dict = load_or_default('model_specs.yaml')
satellite_dict = load_or_default('satellite_specs.yaml')
variable_info = load_or_default('variable_info.yaml')

# --- global functions ------------------------------------------------#

def tmploop_get_remote_files(i,matching,user,pw,
                            server,remote_path,
                            path_local):
    """
    Function to download files using ftp. Tries 10 times before failing.
    """
    print("File: ",matching[i])
    dlstr=('ftp://' + user + ':' + pw + '@'
                + server + remote_path + matching[i])
    for attempt in range(10):
        print ("Attempt to download data: ")
        try:
            print ("Downloading file")
            urlretrieve(dlstr, os.path.join(path_local, matching[i]))
            urlcleanup()
        except Exception as e:
            print (e.__doc__)
            print (e.message)
            print ("Waiting for 10 sec and retry")
            time.sleep(10)
        else:
            break
    else:
        print ('An error was raised and I ' +
              'failed to fix problem myself :(')
        print ('Exit program')
        sys.exit()

def get_remote_files_cmems(\
sdate,edate,twin,nproc,sat,product,path_local,dict_for_sub):
    '''
    Download swath files from CMEMS and store them at defined
    location. Time stamps in file name stand for:

    from, to, creation
    '''
    # credentials
    server = satellite_dict[product]['src']['server']
    user, pw = get_credentials(remoteHostName = server)
    tmpdate = deepcopy(sdate)
    filesort = False
    while (tmpdate <= edate):
        # create remote path
        path_template = satellite_dict[product]['src']\
                                      ['path_template']
        strsublst = satellite_dict[product]['src']\
                                  ['strsub']
        subdict = make_subdict(strsublst,class_object_dict=dict_for_sub)
        path_remote = make_pathtofile(path_template,\
                                      strsublst,subdict,\
                                      date=tmpdate)
        if path_local is None:
            # create local path
            path_template = satellite_dict[product]['dst']\
                                          ['path_template']
            strsublst = satellite_dict[product]['dst']\
                                      ['strsub']
            path_local = make_pathtofile(path_template,\
                                     strsublst,subdict,\
                                     date=sdate)
            filesort = True
        print ('# ----- ')
        print ('Chosen source: ')
        print (sat + ' values from ' + product + ': ' + server)
        print ('# ----- ')
        # get list of accessable files
        ftp = FTP(server)
        ftp.login(user, pw)
        ftp.cwd(path_remote)
        content=FTP.nlst(ftp)
        #choose files according to sdate/edate
        tmplst=[]
        tmpdate_new = tmpdate-timedelta(minutes=twin)
        tmpdate_end = edate+timedelta(minutes=twin)
        while (tmpdate_new <= tmpdate_end):
            matchingtmp = [s for s in content
                            if tmpdate_new.strftime('%Y%m%dT%H')
                            in s ]
            tmplst = tmplst + matchingtmp
            tmpdate_new = tmpdate_new + timedelta(minutes=twin)
        matching = np.unique(tmplst)
        # check if download path exists if not create
        if not os.path.exists(path_local):
            os.makedirs(path_local,exist_ok=True)
        # Download matching files
        print ('Downloading ' + str(len(matching))
                + ' files: .... \n')
        print ("Used number of possible simultaneous downloads "
                + str(nproc) + "!")
        Parallel(n_jobs=nproc)(
                        delayed(tmploop_get_remote_files)(
                        i,matching,user,pw,server,
                        path_remote,path_local
                        ) for i in range(len(matching))
                        )
        # update time
        tmpdate = datetime((tmpdate + relativedelta(months=+1)).year,
                            (tmpdate + relativedelta(months=+1)).month,1)
    if filesort is True:
        # sort files
        print("Data is being sorted into subdirectories " \
            + "year and month ...")
        filelst = [f for f in os.listdir(path_local)
                    if os.path.isfile(os.path.join(path_local,f))]
        sort_files(path_local,filelst,product,sat)
    print ('Files downloaded to: \n', path_local)

def get_remote_files_aviso(\
sdate,edate,twin,nproc,sat,product,path_local,dict_for_sub):
    '''
    Download swath files from AVISO+ and store them at defined
    location.
    '''
    # credentials
    server = satellite_dict[product]['src']['server']
    user, pw = get_credentials(remoteHostName = server)
    tmpdate = deepcopy(sdate)
    filesort = False
    while (tmpdate <= edate):
        # create remote path
        path_template = satellite_dict[product]['src']\
                                      ['path_template']
        strsublst = satellite_dict[product]['src']\
                                  ['strsub']
        subdict = make_subdict(strsublst,class_object_dict=dict_for_sub)
        path_remote = make_pathtofile(path_template,\
                                      strsublst,subdict,\
                                      date=tmpdate)
        if path_local is None:
            # create local path
            path_template = satellite_dict[product]['dst']\
                                          ['path_template']
            strsublst = satellite_dict[product]['dst']\
                                      ['strsub']
            path_local = make_pathtofile(path_template,\
                                     strsublst,subdict,\
                                     date=sdate)
            filesort = True
        print ('# ----- ')
        print ('Chosen source: ')
        print (sat + ' values from ' + product + ': ' + server)
        print ('# ----- ')
        # get list of accessable files
        ftp = FTP(server)
        ftp.login(user, pw)
        ftp.cwd(path_remote)
        content=FTP.nlst(ftp)
        #choose files according to sdate/edate
        tmplst=[]
        tmpdate_new = tmpdate-timedelta(minutes=twin)
        tmpdate_end = edate+timedelta(minutes=twin)
        while (tmpdate_new <= tmpdate_end):
            matchingtmp = [s for s in content
                            if tmpdate_new.strftime('%Y%m%dT%H')
                            in s ]
            tmplst = tmplst + matchingtmp
            tmpdate_new = tmpdate_new + timedelta(minutes=twin)
        matching = np.unique(tmplst)
        # check if download path exists if not create
        if not os.path.exists(path_local):
            os.makedirs(path_local,exist_ok=True)
        # Download matching files
        print ('Downloading ' + str(len(matching))
                + ' files: .... \n')
        print ("Used number of possible simultaneous downloads "
                + str(nproc) + "!")
        Parallel(n_jobs=nproc)(
                        delayed(tmploop_get_remote_files)(
                        i,matching,user,pw,server,
                        path_remote,path_local
                        ) for i in range(len(matching))
                        )
        # update time
        tmpdate = datetime((tmpdate + relativedelta(years=+1)).year,
                            (tmpdate + relativedelta(years=+1)).month,1)
    if filesort is True:
        # sort files
        print("Data is being sorted into subdirectories " \
            + "year and month ...")
        filelst = [f for f in os.listdir(path_local)
                    if os.path.isfile(os.path.join(path_local,f))]
        sort_files(path_local,filelst,product,sat)
    print ('Files downloaded to: \n', path_local)

def get_remote_files_cci(\
sdate,edate,twin,nproc,sat,product,path_local,dict_for_sub):
    '''
    Download swath files from CCI and store them at defined
    location.
    '''
    # credentials
    server = satellite_dict[product]['src']['server']
    level = satellite_dict[product]['processing_level']
    user, pw = get_credentials(remoteHostName = server)
    tmpdate = deepcopy(sdate)
    filesort = False
    while (tmpdate <= edate):
        print(tmpdate)
        # create remote path
        path_template = satellite_dict[product]['src']\
                                      ['path_template']
        strsublst = satellite_dict[product]['src']\
                                  ['strsub']
        dict_for_sub['mission'] =\
                        satellite_dict[product]['mission'][sat]
        subdict = make_subdict(strsublst,class_object_dict=dict_for_sub)
        path_remote = make_pathtofile(path_template,\
                                      strsublst,subdict,\
                                      date=tmpdate)
        if path_local is None:
            # create local path
            subdict['mission'] = sat
            path_template = satellite_dict[product]['dst']\
                                          ['path_template']
            strsublst = satellite_dict[product]['dst']\
                                      ['strsub']
            path_local = make_pathtofile(path_template,\
                                     strsublst,subdict,\
                                     date=sdate)
            filesort = True
        print ('# ----- ')
        print ('Chosen source: ')
        print (sat + ' values from ' + product + ': ' + server)
        print ('# ----- ')
        # get list of accessable files
        ftp = FTP(server)
        ftp.login(user, pw)
        ftp.cwd(path_remote)
        content = FTP.nlst(ftp)
        #choose files according to sdate/edate
        tmplst = []
        tmpdate_new = tmpdate-timedelta(minutes=twin)
        tmpdate_end = edate+timedelta(minutes=twin)
        while (tmpdate_new <= tmpdate_end):
            if level == 'L2P':
                matchingtmp = [s for s in content
                            if tmpdate_new.strftime('-%Y%m%dT%H')
                            in s ]
            elif level == 'L3':
                matchingtmp = [s for s in content
                            if tmpdate_new.strftime('-%Y%m%d')
                            in s ]
            tmplst = tmplst + matchingtmp
            tmpdate_new = tmpdate_new + timedelta(minutes=twin)
        matching = np.unique(tmplst)
        # check if download path exists if not create
        if not os.path.exists(path_local):
            os.makedirs(path_local,exist_ok=True)
        # Download matching files
        print ('Downloading ' + str(len(matching))
                + ' files: .... \n')
        print ("Used number of simultaneous downloads "
                + str(nproc) + "!")
        Parallel(n_jobs=nproc)(
                        delayed(tmploop_get_remote_files)(
                        i,matching,user,pw,server,
                        path_remote,path_local
                        ) for i in range(len(matching))
                        )
        # update time
        tmpdate += timedelta(days=1)
    if filesort is True:
        # sort files
        print("Data is being sorted into subdirectories " \
            + "year and month ...")
        print(path_local)
        filelst = [f for f in os.listdir(path_local)
                    if os.path.isfile(os.path.join(path_local,f))]
        sort_files(path_local,filelst,product,sat)
    print ('Files downloaded to: \n', path_local)

def get_remote_files_eumetsat(\
product,sdate,edate,api_url,sat,path_local,dict_for_sub):
    '''
    Download swath files from EUMETSAT and store them at defined
    location. This fct uses the SentinelAPI for queries.
    '''
    import sentinelsat as ss
    products = None
    dates = (sdate.strftime('%Y-%m-%dT%H:%M:%SZ'),\
             edate.strftime('%Y-%m-%dT%H:%M:%SZ'))
    filesort = False
    if path_local is None:
        # create local path
        path_template = satellite_dict[product]['dst']\
                                      ['path_template']
        strsublst = satellite_dict[product]['dst']\
                                  ['strsub']
        subdict = make_subdict(strsublst,
                               class_object_dict=dict_for_sub)
        path_local = make_pathtofile(path_template,\
                                     strsublst,
                                     subdict,\
                                     date=sdate)
        filesort = True
    kwargs = make_query_dict(product,sat)
    if api_url is None:
        api_url_lst = \
            satellite_dict[product]['src']['api_url']
        for url in api_url_lst:
            print('Source:',url)
            try:
                user, pw = get_credentials(remoteHostName=url)
                api = ss.SentinelAPI(user, pw, url)
                products = api.query(area=None, date=dates,**kwargs)
                break
            except Exception as e:
                if isinstance(e,ss.exceptions.ServerError):
                    print(e)
    else:
        user, pw = get_credentials(remoteHostName = api_url)
        api = ss.SentinelAPI(user, pw, api_url)
        products = api.query(area=None, date=dates,**kwargs)
    if products is not None:
        # check if download path exists if not create
        if not os.path.exists(path_local):
            os.makedirs(path_local,exist_ok=True)
        api.download_all(products,directory_path=path_local)
        #api.download(product_id)
    else: print('No products found!')
    if filesort is True:
        # sort files
        print("Data is being sorted into subdirectories " \
            + "year and month ...")
        filelst = [f for f in os.listdir(path_local)
                    if os.path.isfile(os.path.join(path_local,f))]
        sort_files(path_local,filelst,product,sat)
    print ('Files downloaded to: \n', path_local)

def get_remote_files(path_local,sdate,edate,twin,
                    nproc,product,api_url,sat,dict_for_sub):
    '''
    Download swath files and store them at defined location.
    It is currently possible to download L3 altimeter data from
    CMEMS, L3 and L2P from CEDA CCI, and L2 from EUMETSAT,
    as well as L2P from aviso+ for cfosat swim data.
    '''
    if product=='cmems_L3':
        get_remote_files_cmems(sdate,edate,twin,nproc,\
                               sat,product,path_local,\
                               dict_for_sub)
    elif product=='cfo_swim_L2P':
        get_remote_files_aviso(sdate,edate,twin,nproc,\
                               sat,product,path_local,\
                               dict_for_sub)
    elif product=='eumetsat_L2':
        get_remote_files_eumetsat(product,sdate,edate,\
                                  api_url,sat,path_local,\
                                  dict_for_sub)
    elif product=='cci_L2P' or product=='cci_L3':
        get_remote_files_cci(sdate,edate,twin,nproc,\
                               sat,product,path_local,\
                               dict_for_sub)

def make_query_dict(product,sat):
    '''
    fct to setup queries of L2 data using SentinelAPI
    '''
    level = satellite_dict[product]['mission'].get('processing')
    SAT = satellite_dict[product]['mission'].get(sat)
    kwargs =  {'platformname': 'Sentinel-3',
               'instrumentshortname': 'SRAL',
               'productlevel': level,
               'filename': SAT + '*WAT*'}
    return kwargs

def get_local_files(sdate,edate,twin,product,dict_for_sub=None,
path_local=None):
    """
    Function to retrieve list of files/paths for available
    locally stored satellite data. This list is used for
    other functions to query and parsing.

    param:
        sdate - start date (datetime object)
        edate - end date (datetime object)
        twin - time window (temporal constraint) in minutes
        product - roduct as of satellite_specs.yaml
        dict_for_sub - dictionary for substitution in templates
        local_path - a path if defined

    return:
        pathlst - list of paths
        filelst - list of files
    """
    filelst = []
    pathlst = []
    tmpdate = sdate-timedelta(minutes=twin)
    if path_local is None:
        print('path_local is None -> checking config file')
        while (tmpdate <= edate + relativedelta(months=+1)):
            try:
                # create local path for each time
                path_template = \
                        satellite_dict[product]['dst'].get(
                                              'path_template')
                strsublst = \
                        satellite_dict[product]['dst'].get('strsub')
                subdict = \
                        make_subdict(strsublst,
                                     class_object_dict=dict_for_sub)
                path_local = make_pathtofile(path_template,\
                                             strsublst,subdict)
                path_local = (
                            os.path.join(
                            path_local,
                            tmpdate.strftime('%Y'),
                            tmpdate.strftime('%m'))
                            )
                print(path_local)
                if os.path.isdir(path_local):
                    tmplst = np.sort(os.listdir(path_local))
                    filelst.append(tmplst)
                    pathlst.append([os.path.join(path_local,e)
                                    for e in tmplst])
                tmpdate = tmpdate + relativedelta(months=+1)
                path_local = None
            except Exception as e:
                print(e)
                tmpdate = tmpdate + relativedelta(months=+1)
        filelst = np.sort(flatten(filelst))
        pathlst = np.sort(flatten(pathlst))
    else:
        filelst = np.sort(os.listdir(path_local))
        pathlst = [os.path.join(path_local,e) for e in filelst]
    idx_start,tmp = check_date(filelst, sdate - timedelta(minutes=twin))
    tmp,idx_end = check_date(filelst, edate + timedelta(minutes=twin))
    if idx_end == 0:
        idx_end = len(pathlst)-1
    del tmp
    pathlst = np.unique(pathlst[idx_start:idx_end+1])
    filelst = np.unique(filelst[idx_start:idx_end+1])
    print (str(int(len(pathlst))) + " valid files found")
    return pathlst, filelst

def read_local_ncfiles(pathlst,product,varalias,
sd,ed,twin,variable_info):
    """
    Wrapping function to read satellite netcdf files.

    param:
        pathlst - list of paths to be parsed
        product - product as specified in satellite_specs.yaml
        varalias
        sd - start date (datetime object)
        ed - start date (datetime object)
        twin - time window (temporal constraint) in minutes
        variable_info - from variable_info.yaml

    return:
        dictionary of variables for the satellite_class object
    """
    # adjust start and end
    sd = sd - timedelta(minutes=twin)
    ed = ed + timedelta(minutes=twin)
    # get meta data
    ncmeta = ncdumpMeta(pathlst[0])
    ncvar = get_filevarname(varalias,variable_info,
                            satellite_dict[product],ncmeta)
    # retrieve sliced data
    ds = read_netcdfs(pathlst)
    ds_sort = ds.sortby('time')
    ds_sliced = ds_sort.sel(time=slice(sd, ed))
    # make dict and start with stdvarname for varalias
    stdvarname = variable_info[varalias]['standard_name']
    var_sliced = ds_sliced[ncvar]
    vardict = {}
    vardict[stdvarname] = list(var_sliced.values)
    # add coords to vardict
    # 1. retrieve list of coordinates
    coords_lst = list(var_sliced.coords.keys())
    # 2. iterate over coords_lst
    for varname in coords_lst:
        stdcoordname = ds_sliced[varname].attrs['standard_name']
        if stdcoordname == 'longitude':
            vardict[stdcoordname] = \
                list(((ds_sliced[varname].values - 180) % 360) - 180)
        elif stdcoordname == 'time':
            # convert to unixtime
            df_time = ds_sliced[varname].to_dataframe()
            unxt = (pd.to_datetime(df_time[varname]).view(int) / 10**9)
            vardict[stdcoordname] = unxt.values
            vardict['time_unit'] = variable_info[stdcoordname]['units']
        else:
            vardict[stdcoordname] = list(ds_sliced[varname].values)
    return vardict

def unzip_eumetsat(pathlst,tmpdir):
    """
    Function to unzip eumetsat files prior to reading

    param:
        pathlst - list of paths to zipped files
        tmpdir - temporary folder to unzipped files

    return:
        pathlst_new - new list of paths to unzipped files
    """
    for count, f in enumerate(pathlst):
        zipped = zipfile.ZipFile(f)
        enhanced_measurement = zipped.namelist()[-1]
        extracted = zipped.extract(enhanced_measurement,
                                   path=tmpdir.name)
        fname = extracted.split('/')[-1]
        dignumstr = '_{num:0' + str(len(str(len(pathlst)))) + 'd}.nc'
        # cp extracted file to parent tmp
        cmdstr = ('cp ' + extracted + ' ' + tmpdir.name
                + '/' + fname[0:-3]
                + dignumstr.format(num=count))
                #+ '_{num:04d}.nc'.format(num=count))
        os.system(cmdstr)
        # delete subfolder
        cmdstr = ('rm -r '
                 + os.path.dirname(os.path.realpath(extracted)))
        os.system(cmdstr)
    flst = os.listdir(tmpdir.name)
    pathlst_new = []
    for f in flst:
        pathlst_new.append(os.path.join(tmpdir.name,f))
    return pathlst_new

def read_local_files_eumetsat(pathlst,product,varalias,satellite_dict):
    '''
    Read and concatenate all data to one timeseries for each variable.
    Fct is tailored to EUMETSAT files.
    '''
    # --- find variable cf names --- #
    print ("Processing " + str(int(len(pathlst))) + " files")
    print (pathlst[0])
    print (pathlst[-1])
    # --- find ncvar cf names --- #
    tmpdir = tempfile.TemporaryDirectory()
    zipped = zipfile.ZipFile(pathlst[0])
    enhanced_measurement = zipped.namelist()[-1]
    extracted = zipped.extract(enhanced_measurement, path=tmpdir.name)
    stdname = variable_info[varalias]['standard_name']
    ncmeta = ncdumpMeta(extracted)
    ncvar = get_filevarname(varalias,variable_info,
                            satellite_dict[product],ncmeta)
    tmpdir.cleanup()
    # --- create vardict --- #
    vardict = {}
    tmpdir = tempfile.TemporaryDirectory()
    print('tmp directory is established:',tmpdir.name)
    for f in tqdm(pathlst):
        #path = f[0:-len(f.split('/')[-1])]
        zipped = zipfile.ZipFile(f)
        enhanced_measurement = zipped.namelist()[-1]
        extracted = zipped.extract(enhanced_measurement,
                                   path=tmpdir.name)
        ds = xr.open_dataset(extracted)
        ds_var = ds[ncvar]
        if stdname in vardict.keys():
            vardict[stdname] += list(ds[ncvar])
        else:
            vardict[stdname] = list(ds[ncvar])
        coords_lst = list(ds_var.coords.keys())
        for varname in coords_lst:
            stdcoordname = ds[varname].attrs['standard_name']
            if stdcoordname == 'longitude':
                if stdcoordname in vardict.keys():
                    vardict[stdcoordname] += list(ds[varname].values)
                else:
                    vardict[stdcoordname] = list(ds[varname].values)
            elif stdcoordname == 'time':
                # convert to unixtime
                df_time = ds[varname].to_dataframe()
                unxt = (pd.to_datetime(df_time[varname])\
                        .view(int) / 10**9)
                if stdcoordname in vardict.keys():
                    vardict[stdcoordname] += list(unxt.values)
                else:
                    vardict[stdcoordname] = list(unxt.values)
            else:
                if stdcoordname in vardict.keys():
                    vardict[stdcoordname] += list(ds[varname].values)
                else:
                    vardict[stdcoordname] = list(ds[varname].values)
    # transform to -180 to 180 degrees
    tmp = np.array(vardict['longitude'])
    vardict['longitude'] = list(((tmp - 180) % 360) - 180)
    vardict['time_unit'] = variable_info['time']['units']
    tmpdir.cleanup()
    return vardict

def read_local_files(pathlst,product,varalias,
    sd,ed,twin,variable_info,satellite_dict):
    '''
    wrapping function to read altimetry files

    return:
        vardict - dictionary of variables for altimeter data
    '''
    # read local files depending on product
    if product == 'cmems_L3':
        vardict = read_local_ncfiles(pathlst,product,varalias,
                                     sd,ed,twin,variable_info)
    elif (product == 'cci_L2P' or product == 'cci_L3'):
        vardict = read_local_ncfiles(pathlst,product,varalias,
                                     sd,ed,twin,variable_info)
    elif (product == 'eumetsat_L2'):
        sys.exit('!!! eumetsat L2 temporarily not provided !!!')
        vardict = read_local_files_eumetsat(pathlst,product,\
                                            varalias,\
                                            satellite_dict)
    return vardict

def get_sat_ts(sdate,edate,twin,region,product,pathlst,
varalias,poi,distlim,variable_info,satellite_dict):
    """
    Main function to obtain data from satellite missions.
    reads files, apply region and temporal filter

    return: adjusted dictionary according to spatial and
            temporal contarinst
    """
    cvardict = read_local_files(pathlst,product,varalias,
                               sdate,edate,twin,variable_info,
                               satellite_dict)
    print('Total: ', len(cvardict['time']), ' footprints found')
    print('Apply region mask')
    ridx = match_region(cvardict['latitude'],
                        cvardict['longitude'],
                        region=region,
                        grid_date=sdate)
    print('Region mask applied')
    rvardict = {}
    for element in cvardict:
        if element != 'time_unit':
            rvardict[element] = list(np.array(
                                    cvardict[element]
                                    )[ridx])
        else:
            rvardict[element] = cvardict[element]
    del cvardict,ridx
    if len(rvardict['time'])>0:
        rvardict['datetime'] = netCDF4.num2date(
                                    rvardict['time'],
                                    rvardict['time_unit'])
        print('For chosen region and time: ',
                len(rvardict['time']),'footprints found')
        # convert to datetime object
        timedt = rvardict['datetime']
        rvardict['datetime'] = [datetime(t.year,t.month,t.day,
                                         t.hour,t.minute,t.second,\
                                         t.microsecond)\
                                for t in timedt]
    else:
        print('For chosen region and time: 0 footprints found!')
    if poi is not None:
        pvardict = {}
        pidx = match_poi(rvardict,twin,distlim,poi)
        for element in rvardict:
            if element != 'time_unit':
                pvardict[element] = list(np.array(
                                        rvardict[element]
                                        )[pidx])
            else:
                pvardict[element] = rvardict[element]
        rvardict = pvardict
        print('For chosen poi: ',
                len(rvardict['time']),'footprints found')
    # find variable name as defined in file
    if product == 'cmems_L3':
        ncdict = ncdumpMeta(pathlst[0])
    elif (product == 'cci_L2P' or product == 'cci_L3'):
        ncdict = ncdumpMeta(pathlst[0])
    elif product == 'eumetsat_L2':
        tmpdir = tempfile.TemporaryDirectory()
        zipped = zipfile.ZipFile(pathlst[0])
        enhanced_measurement = zipped.namelist()[-1]
        extracted = zipped.extract(enhanced_measurement,
                                   path=tmpdir.name)
        ncdict = ncdumpMeta(extracted)
        tmpdir.cleanup()
    rvardict['meta'] = ncdict
    return rvardict

def crop_vardict_to_period(vardict,sdate,edate):
    """
    Function to crop the variable dictionary to a given period
    """
    for key in vardict:
        if (key != 'time_unit' and key != 'meta' and key != 'datetime'):
            vardict[key] =  list(np.array(vardict[key])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ])
        else:
            vardict[key] = vardict[key]
    vardict['datetime'] =   list(np.array(vardict['datetime'])[ \
                                ( (np.array(vardict['datetime'])>=sdate)
                                & (np.array(vardict['datetime'])<=edate)
                                ) ])
    return vardict


# flatten all lists before returning them
# define flatten function for lists
''' fct does the following:
flat_list = [item for sublist in TIME for item in sublist]
or:
for sublist in TIME:
for item in sublist:
flat_list.append(item)
'''
flatten = lambda l: [item for sublist in l for item in sublist]

def check_date(filelst,date):
    '''
    Checks if str in lst according to desired date (sdate,edate)

    return: idx for file
    '''
    idx = []
    for i in range(len(filelst)):
        element = filelst[i]
        tmp = element.find(date.strftime('%Y%m%d'))
        if tmp>=0:
            idx.append(i)
    if len(idx)<=0:
        idx=[0]
    return idx[0],idx[-1]

# ---------------------------------------------------------------------#


class satellite_class():
    '''
    Class to handle netcdf files containing satellite data e.g.
    Hs[time], lat[time], lon[time]

    This class offers the following added functionality:
     - get swaths of desired days and read
     - get the closest time stamp(s)
     - get the location (lon, lat) for this time stamp
     - get Hs or 10m wind value for this time
     - region mask
    '''

    def __init__(
        self,sdate=None,mission='s3a',product='cmems_L3',
        edate=None,twin=None,download=False,path_local=None,
        region='mwam4',nproc=1,varalias='Hs',api_url=None,
        filterData=False,poi=None,distlim=None,
        **kwargs):
        print('# ----- ')
        print(" ### Initializing satellite_class object ###")
        print(" ")
        # parse and translate date input
        sdate = parse_date(sdate)
        edate = parse_date(edate)
        # check settings
        if (sdate is None and edate is None and poi is not None):
            sdate = poi['datetime'][0]
            edate = poi['datetime'][-1]
        if edate is None:
            print ("Requested time: ", str(sdate))
            edate = sdate
        else:
            print ("Requested time frame: " +
                str(sdate) + " - " + str(edate))
        if twin is None:
            twin = int(30)
        stdname = variable_info[varalias].get('standard_name')
        units = variable_info[varalias].get('units')
        # define some class variables
        self.sdate = sdate
        self.edate = edate
        self.varalias = varalias
        self.units = units
        self.stdvarname = stdname
        self.twin = twin
        self.region = region
        self.mission = mission
        self.obstype = 'satellite_altimeter'
        self.product = product
        self.provider = satellite_dict[product].get('provider')
        self.processing_level = \
                satellite_dict[product].get('processing_level')
        print('Chosen time window is:', twin, 'min')
        # make satpaths
        if path_local is None:
            path_template = satellite_dict[product]\
                                          ['dst']\
                                          ['path_template']
            self.path_local = path_template
        else:
            self.path_local = path_local
        # retrieve files
        if download is False:
            print ("No download initialized, checking local files")
        else:
            print ("Downloading necessary files ...")
            get_remote_files(path_local,sdate,edate,twin,
                            nproc,product,api_url,mission,vars(self))
        print(" ")
        print(" ## Find files ...")
        t0=time.time()
        pathlst, _ = get_local_files(sdate,edate,twin,
                                     product,vars(self),
                                     path_local=path_local)
        print(" ")
        print(" ## Read files ...")
        if len(pathlst) > 0:
#            for i in range(1):
            try:
                if filterData == True:
                    # extend time period due to filter
                    if 'stwin' not in kwargs.keys():
                        kwargs['stwin'] = 1 # needs to be changed
                    if 'etwin' not in kwargs.keys():
                        kwargs['etwin'] = 1
                    twin_tmp = twin + kwargs['stwin'] + kwargs['etwin']
                    # retrieve data
                    rvardict = get_sat_ts( sdate,edate,
                                           twin_tmp,region,
                                           product,pathlst,
                                           varalias,poi,distlim,
                                           variable_info,
                                           satellite_dict)
                    # filter data
                    rvardict = filter_main( rvardict,
                                            varalias = varalias,
                                            **kwargs )
                    # crop to original time period
                    sdate_tmp = sdate - timedelta(minutes=twin)
                    edate_tmp = sdate + timedelta(minutes=twin)
                    rvardict = crop_vardict_to_period(rvardict,
                                                      sdate_tmp,
                                                      edate_tmp)
                    self.filter = True
                    self.filterSpecs = kwargs
                else:
                    rvardict = get_sat_ts( sdate,edate,twin,region,
                                           product,pathlst,
                                           varalias,poi,distlim,
                                           variable_info,
                                           satellite_dict )
                    # make ts in vardict unique
                    rvardict = vardict_unique(rvardict)
                    # rm NaNs
                    rvardict = rm_nan_from_vardict(varalias,rvardict)
                # find variable name as defined in file
                if product == 'cmems_L3':
                    ncdict = ncdumpMeta(pathlst[0])
                elif (product == 'cci_L2P' or product == 'cci_L3'):
                    ncdict = ncdumpMeta(pathlst[0])
                elif product == 'eumetsat_L2':
                    tmpdir = tempfile.TemporaryDirectory()
                    zipped = zipfile.ZipFile(pathlst[0])
                    enhanced_measurement = zipped.namelist()[-1]
                    extracted = zipped.extract(enhanced_measurement,
                                               path=tmpdir.name)
                    ncdict = ncdumpMeta(extracted)
                    tmpdir.cleanup()
                with NoStdStreams():
                    filevarname = get_filevarname(varalias,
                                              variable_info,
                                              satellite_dict[product],
                                              ncdict)
                rvardict['meta'] = ncdict
                # define more class variables
                self.vars = rvardict
                self.varname = filevarname
                t1=time.time()
                print(" ")
                print( '## Summary:')
                print(str(len(self.vars['time'])) + " footprints retrieved.")
                print("Time used for retrieving satellite data:",\
                        round(t1-t0,2),"seconds")
                print(" ")
                print ("### Satellite object initialized ###")
                print ('# ----- ')
            except Exception as e:
                print(e)
                print('Error encountered')
                print('No satellite_class object initialized')
        else:
            print('No satellite data found')
            print('No satellite_class object initialized')
            print ('# ----- ')

    def get_item_parent(self,item,attr):
        """
        Offers possibility to explore netcdf meta info.
        by specifying what you are looking for (item),
        e.g. part of a string, and in which attribute (attr),
        e.g. standard_name, this function returns the
        parent parameter name of the query string.

        param:
            item - (partial) string e.g. [m]
            attr - attribute e.g. units

        return: list of matching parameter strings

        e.g. for satellite_class object sco:

        .. code ::

            sco.get_item_parent('m','units')
        """

        ncdict = self.vars['meta']
        lst = [i for i in ncdict.keys() \
                if (attr in ncdict[i].keys() \
                and item in ncdict[i][attr]) \
                ]
        if len(lst) >= 1:
            return lst
        else: return None

    def get_item_child(self,item):
        """
        Gets all attributes connected to given parameter name.

        param:
            item - (partial) string e.g. [m]

        return: matching parameter string

        e.g. for satellite_class object sco:

        .. code ::

            sco.get_item_child('time')
        """

        ncdict = self.vars['meta']
        parent = finditem(ncdict,item)
        return parent

    def quicklook(self,m=False,ts=False,projection=None):
        """
        Enables to explore the class object (and retrieved results)
        by plotting time series and map.

        param:
            m - map figure (True/False)
            ts - time series (True/False)
            projection - specified projection for cartopy

        return:
            figures
        """
        if m:
            import cartopy.crs as ccrs
            import cmocean
            import matplotlib.pyplot as plt
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            lons = self.vars['longitude']
            lats = self.vars['latitude']
            var = self.vars[self.stdvarname]
            if projection is None:
                projection = ccrs.PlateCarree()
            lonmax,lonmin = np.max(lons),np.min(lons)
            latmax,latmin = np.max(lats),np.min(lats)
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_extent(  [lonmin, lonmax,latmin, latmax],
                            crs = projection )
            sc = ax.scatter(lons,lats,s=10,
                            c = var,
                            marker='o', edgecolor = 'face',
                            cmap=cmocean.cm.amp,
                            transform=ccrs.PlateCarree())
            axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
            fig.colorbar(sc, cax=axins, label=self.varalias
                                        + ' [' + self.units + ']')
            ax.coastlines()
            gl = ax.gridlines(draw_labels=True,crs=projection,
                              linewidth=1, color='grey', alpha=0.4,
                              linestyle='-')
            gl.top_labels = False
            gl.right_labels = False
            plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
            ax.set_title(self.mission + ' (' + self.provider + ')\n'
                      + 'from ' + str(self.vars['datetime'][0])
                      + ' to ' + str(self.vars['datetime'][-1]))
            #fig.suptitle('', fontsize=16) # unused
            plt.show()
        if ts:
            import matplotlib.pyplot as plt
            import matplotlib.dates as mdates
            fig = plt.figure(figsize=(9,3.5))
            ax = fig.add_subplot(111)
            colors = ['k']
            ax.plot(self.vars['datetime'],
                    self.vars[self.stdvarname],
                    linestyle='None',color=colors[0],
                    label=self.mission,
                    marker='o',alpha=.5,ms=2)
            plt.ylabel(self.varalias + '[' + self.units + ']')
            plt.legend(loc='best')
            plt.tight_layout()
            #ax.set_title()
            plt.show()

    def write_to_nc(self,pathtofile=None,file_date_incr=None):
        """
        Write class variables to netcdf files.

        param:
            pathtofile
            file_date_incr - what date increment to use for files
        """
        if 'error' in vars(self):
            print('Erroneous satellite_class file detected')
            print('--> dump to netCDF not possible !')
        else:
            tmpdate = self.sdate
            edate = self.edate
            while tmpdate <= edate:
                if pathtofile is None:
                    path_template = satellite_dict[self.product]\
                                                  ['dst']\
                                                  ['path_template']
                    file_template = satellite_dict[self.product]\
                                                  ['dst']\
                                                  ['file_template']
                    strsublst = satellite_dict[self.product]\
                                              ['dst']['strsub']
                    if 'filterData' in vars(self).keys():
                        file_template = 'filtered_' + file_template
                    tmppath = os.path.join(path_template,file_template)
                    subdict = make_subdict(strsublst,
                                           class_object_dict=vars(self))
                    pathtofile = make_pathtofile(tmppath,strsublst,
                                                 subdict,
                                                 date=tmpdate)
                title = (self.obstype
                       + ' observations from '
                       + self.mission)
                dumptonc_ts_sat(self,pathtofile,title)
                # determine date increment
                if file_date_incr is None:
                    file_date_incr = satellite_dict[self.product]\
                                    ['dst'].get('file_date_incr','m')
                if file_date_incr == 'm':
                    tmpdate += relativedelta(months = +1)
                elif file_date_incr == 'Y':
                    tmpdate += relativedelta(years = +1)
                elif file_date_incr == 'd':
                    tmpdate += timedelta(days = +1)
        return

def poi_sat(indict,twin,distlim,poi,ridx,i):
    """
    return: indices for values matching the spatial and
            temporal constraints
    """
    tidx = find_included_times(
                list(np.array(indict['datetime'])[ridx]),
                target_t=poi['datetime'][i],
                twin=twin)
    slons = list(np.array(indict['longitude'])[ridx][tidx])
    slats = list(np.array(indict['latitude'])[ridx][tidx])
    plons = [poi['longitude'][i]]*len(slons)
    plats = [poi['latitude'][i]]*len(slats)
    dists_tmp = haversineA( slons,slats,plons,plats )
    sidx = np.argwhere(np.array(dists_tmp)<=distlim).flatten()
    dists = list(np.array(dists_tmp)[sidx])
    return list(np.array(tidx)[sidx])

def match_poi(indict, twin, distlim, poi):
    """
    return: idx that match to region
    """
    region={'llcrnrlat':np.min(poi['latitude']),
            'urcrnrlat':np.max(poi['latitude']),
            'llcrnrlon':np.min(poi['longitude']),
            'urcrnrlon':np.max(poi['longitude'])}
    ridx = match_region_rect(indict['latitude'],
                             indict['longitude'],
                             region=region)
    sat_dict = deepcopy(indict)
    idx = [poi_sat(sat_dict,twin,distlim,poi,ridx,i) \
                  for i in range(len(poi['datetime']))]
    idx = list(np.array(ridx)[flatten(idx)])
    return idx

def match_region(LATS,LONS,region,grid_date):
    """
    Function to filter satellite data according to region

    return:
        indices that match the region
    """
    # region in region_dict[poly]:
    # find values for given region
    if (region not in region_dict['poly'] and \
    region not in model_dict and \
    region not in region_dict['geojson']):
        if region is None:
            region = 'global'
        else:
            if (region not in region_dict['rect'] \
            and region not in region_dict['geojson'] \
            and isinstance(region,dict)==False):
                sys.exit("Region is not defined")
            elif isinstance(region,dict):
                print("Region bound by min/max of poi coordinates")
            else:
                print("Specified region: " + region + "\n"
                + " --> Bounds: " + str(region_dict['rect'][region]))
                region = region_dict['rect'][region]
        ridx = match_region_rect(LATS,LONS,region=region)
    elif region in region_dict['geojson']:
        print("Region is defined as geojson")
        ridx = match_region_geojson(LATS,LONS,region=region)
    else:
        ridx = match_region_poly(LATS,LONS,region=region,
                                grid_date=grid_date)
    return ridx

def match_region_rect(LATS,LONS,region):
    """
    Takes care of regular grid regions
    """
    if (region is None or region == "global"):
        region = "global"
        ridx = range(len(LATS))
    else:
        llcrnrlat = region["llcrnrlat"]
        urcrnrlat = region["urcrnrlat"]
        llcrnrlon = region["llcrnrlon"]
        urcrnrlon = region["urcrnrlon"]
        ridx = np.where(\
                    (np.array(LATS)>llcrnrlat)\
                    &(np.array(LATS)<urcrnrlat)\
                    &(np.array(LONS)>llcrnrlon)\
                    &(np.array(LONS)<urcrnrlon)\
                    )[0]
    print (len(ridx), " values found for chosen region and time frame.")
    return ridx

def match_region_geojson(LATS,LONS,region):
    """
    Takes care of regions defines as geojson
    """
    import geojson
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    LATS = list(LATS)
    LONS = list(LONS)
    lats = np.array(LATS).ravel()
    lons = np.array(LONS).ravel()
    points = np.c_[lons,lats]
    ridx = []
    fstr = region_dict['geojson'][region]['fstr']
    with open(fstr) as f:
        gj = geojson.load(f)
    fidx = region_dict['geojson'][region].get('fidx')
    if fidx is not None:
        geo = {'type': 'Polygon',
           'coordinates':\
                gj['features'][fidx]['geometry']['coordinates'][0]}
        poly = Polygon([tuple(l)
                        for l in geo['coordinates'][0]],
                        closed=True)
        hits = Path(poly.xy).contains_points(points,radius=1e-9)
        ridx_tmp = list(np.array(range(len(LONS)))[hits])
        ridx += ridx_tmp
    else:
        for i in range(len(gj['features'])):
            geo = {'type': 'Polygon',
               'coordinates':\
                    gj['features'][i]['geometry']['coordinates'][0]}
            poly = Polygon([tuple(l)
                            for l in geo['coordinates'][0]],
                            closed=True)
            hits = Path(poly.xy).contains_points(points,radius=1e-9)
            ridx_tmp = list(np.array(range(len(LONS)))[hits])
            ridx += ridx_tmp
    return ridx

def match_region_poly(LATS,LONS,region,grid_date):
    """
    Takes care of region defined as polygon
    """
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    import numpy as np
    if (region not in region_dict['poly'] \
        and region not in model_dict):
        sys.exit("Region polygone is not defined")
    elif isinstance(region,dict)==True:
        print ("Manuall specified region: \n"
            + " --> Bounds: " + str(region))
        poly = Polygon(list(zip(region['lons'],
            region['lats'])), closed=True)
    elif (isinstance(region,str)==True and region in model_dict):
        try:
            print('Use date for retrieving grid: ', grid_date)
            filestr = make_model_filename_wrapper(\
                                    region,grid_date,'best')
            meta = ncdumpMeta(filestr)
            flon = get_filevarname('lons',variable_info,\
                                model_dict[region],meta)
            flat = get_filevarname('lats',variable_info,\
                                model_dict[region],meta)
            time = get_filevarname('time',variable_info,\
                                model_dict[region],meta)
            #M = xa.open_dataset(filestr, decode_cf=True)
            #model_lons = M[flon].data
            #model_lats = M[flat].data
            model_lons, model_lats, _ = \
                read_model_nc_output_lru(filestr,flon,flat,time)
        except (KeyError,IOError,ValueError) as e:
            print(e)
            if 'grid_date' in model_dict[region]:
                grid_date = model_dict[region]['grid_date']
                print('Trying default date ', grid_date)
            else:
                grid_date = datetime(
                                    datetime.now().year,
                                    datetime.now().month,
                                    datetime.now().day
                                    )
            filestr = make_model_filename_wrapper(\
                                    region,grid_date,'best')
            meta = ncdumpMeta(filestr)
            flon = get_filevarname('lons',variable_info,\
                                model_dict[region],meta)
            flat = get_filevarname('lats',variable_info,\
                                model_dict[region],meta)
            time = get_filevarname('time',variable_info,\
                                model_dict[region],meta)
            #M = xa.open_dataset(filestr, decode_cf=True)
            #model_lons = M[flon].data
            #model_lats = M[flat].data
            model_lons, model_lats, _ = \
                read_model_nc_output_lru(filestr,flon,flat,time)
        if (len(model_lons.shape)==1):
            model_lons, model_lats = np.meshgrid(
                                    model_lons,
                                    model_lats
                                    )
        print('Check if footprints fall within the chosen domain')
        ncdict = ncdumpMeta(filestr)
        try:
            proj4 = find_attr_in_nc('proj',ncdict=ncdict,
                                    subattrstr='proj4')
        except IndexError:
            print('proj4 not defined in netcdf-file')
            print('Using proj4 from model config file')
            proj4 = model_dict[region]['proj4']
        proj_model = pyproj.Proj(proj4)
        Mx, My = proj_model(model_lons,model_lats,inverse=False)
        Vx, Vy = proj_model(LONS,LATS,inverse=False)
        xmax, xmin = np.max(Mx), np.min(Mx)
        ymax, ymin = np.max(My), np.min(My)
        ridx = list( np.where((Vx>xmin) & (Vx<xmax) &
                              (Vy>ymin) & (Vy<ymax))[0] )
    elif isinstance(region,str)==True:
        print ("Specified region: " + region + "\n"
        + " --> Bounded by polygon: \n"
        + "lons: " + str(region_dict['poly'][region]['lons']) + "\n"
        + "lats: " + str(region_dict['poly'][region]['lats']))
        poly = Polygon(list(zip(region_dict['poly'][region]['lons'],
            region_dict['poly'][region]['lats'])), closed=True)
        # check if coords in region
        LATS = list(LATS)
        LONS = list(LONS)
        lats = np.array(LATS).ravel()
        lons = np.array(LONS).ravel()
        points = np.c_[lons,lats]
        # radius seems to be important to correctly define polygone
        # see discussion here:
        # https://github.com/matplotlib/matplotlib/issues/9704
        hits = Path(poly.xy).contains_points(points,radius=1e-9)
        ridx = list(np.array(range(len(LONS)))[hits])
    if not ridx:
        print ("No values for chosen region and time frame!!!")
    else:
        print ("Values found for chosen region and time frame.")
    return ridx
