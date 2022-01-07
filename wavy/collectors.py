#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module comprises functions to collect remote satellite data.
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
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed

# own imports
from wavy.utils import sort_files
from wavy.utils import make_pathtofile, make_subdict
from wavy.credentials import get_credentials
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_specs.yaml')

# --- def functions ---------------------------------------------------#

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
        print(path_remote)
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

        print(path_local)
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

def get_remote_files(path_local,sdate,edate,twin,
                    nproc,product,api_url,sat,dict_for_sub):
    '''
    Download swath files and store them at defined location.
    It is currently possible to download L3 altimeter data from
    CMEMS, L3 and L2P from CEDA CCI, and L2 from EUMETSAT,
    as well as L2P from aviso+ for cfosat swim data.
    '''
    if (product=='cmems_L3_NRT' or product=='cmems_L3_MY'):
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
