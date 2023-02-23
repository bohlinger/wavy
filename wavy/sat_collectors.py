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
from urllib.parse import quote
from ftplib import FTP
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed

# own imports
from wavy.utils import sort_files
from wavy.utils import make_pathtofile, make_subdict
from wavy.utils import date_dispatcher
from wavy.credentials import get_credentials
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default('satellite_specs.yaml')

# --- def functions ---------------------------------------------------#

def tmploop_get_remote_files(i: int, matching: str,
                             user: str, pw: str,
                             server: str, remote_path: str,
                             path_local: str):
    """
    Function to download files using ftp. Tries 10 times before failing.
    """
    print("File: ",matching[i])
    print("src path: ", remote_path)
    pw = quote(pw) # to escape special characters
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

def get_remote_files_cmems(**kwargs):
    '''
    Download swath files from CMEMS and store them at defined
    location. Time stamps in file name stand for:

    from, to, creation
    '''
    product = kwargs.get('product')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    twin = kwargs.get('twin',30)
    nproc = kwargs.get('nproc',1)
    mission = kwargs.get('mission','s3a')
    path_local = kwargs.get('path_local')
    dict_for_sub = kwargs.get('dict_for_sub')
    # check if search str template
    file_search_template = \
        satellite_dict[product]['src'].get('file_search_template',\
                                            '%Y%m%dT%H')
    # credentials
    server = satellite_dict[product]['src']['server']
    user, pw = get_credentials(remoteHostName = server)
    tmpdate = deepcopy(sdate)
    filesort = False
    path_template_src = satellite_dict[product]['src']\
                                  ['path_template']
    strsublst_src = satellite_dict[product]['src']\
                              ['strsub']
    subdict_src = make_subdict(strsublst_src,
                               class_object_dict=dict_for_sub)
    while (tmpdate <= edate):
        # create remote path
        path_remote = make_pathtofile(path_template_src,\
                                      strsublst_src,subdict_src,\
                                      date=tmpdate)
        if path_local is None:
            # create local path
            path_template_dst = satellite_dict[product]['dst']\
                                          ['path_template']
            strsublst_dst = satellite_dict[product]['dst']\
                                          ['strsub']
            subdict_dst = make_subdict(strsublst_dst,
                                           class_object_dict=dict_for_sub)
            path_local = make_pathtofile(path_template_dst,\
                                         strsublst_dst,subdict_dst,\
                                         date=tmpdate)
            filesort = True

        print ('# ----- ')
        print ('Chosen source: ')
        print (mission + ' values from ' + product + ': ' + server)
        print(path_remote)
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
                            if tmpdate_new.strftime(file_search_template)
                            in s ]
            tmplst = tmplst + matchingtmp
            tmpdate_new = tmpdate_new + timedelta(minutes=twin)
        matching = np.unique(tmplst)
        print(matching)
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
        #tmpdate = datetime((tmpdate + relativedelta(months=+1)).year,
        #                    (tmpdate + relativedelta(months=+1)).month,1)
        date_incr = satellite_dict[product]['src'].get('date_incr', 'm')
        tmpdate = date_dispatcher(tmpdate, date_incr=date_incr)
    if filesort is True:
        # sort files
        print("Data is being sorted into subdirectories " \
            + "year and month ...")
        filelst = [f for f in os.listdir(path_local)
                    if os.path.isfile(os.path.join(path_local,f))]
        sort_files(path_local,filelst,product,mission)
    print ('Files downloaded to: \n', path_local)

def get_remote_files_aviso(**kwargs):
    '''
    Download swath files from AVISO+ and store them at defined
    location.
    '''
    product = kwargs.get('product')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    twin = kwargs.get('twin',30)
    nproc = kwargs.get('nproc',1)
    mission = kwargs.get('mission','cfo')
    path_local = kwargs.get('path_local')
    dict_for_sub = kwargs.get('dict_for_sub')
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
        print (mission + ' values from ' + product + ': ' + server)
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
        sort_files(path_local,filelst,product,mission)
    print ('Files downloaded to: \n', path_local)

def get_remote_files_cci(**kwargs):
    '''
    Download swath files from CCI and store them at defined
    location.
    '''
    product = kwargs.get('product')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    twin = kwargs.get('twin',30)
    nproc = kwargs.get('nproc',1)
    mission = kwargs.get('mission','multi')
    path_local = kwargs.get('path_local')
    dict_for_sub = kwargs.get('dict_for_sub')
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
                        satellite_dict[product]['mission'][mission]
        subdict = make_subdict(strsublst,class_object_dict=dict_for_sub)
        path_remote = make_pathtofile(path_template,\
                                      strsublst,subdict,\
                                      date=tmpdate)
        if path_local is None:
            # create local path
            subdict['mission'] = mission
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
        print (mission + ' values from ' + product + ': ' + server)
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
        sort_files(path_local,filelst,product,mission)
    print ('Files downloaded to: \n', path_local)

def get_remote_files_eumetsat(**kwargs):
    '''
    Download swath files from EUMETSAT and store them at defined
    location. This fct uses the SentinelAPI for queries.
    '''
    product = kwargs.get('product')
    sdate = kwargs.get('sdate')
    edate = kwargs.get('edate')
    mission = kwargs.get('mission','s3a')
    path_local = kwargs.get('path_local')
    dict_for_sub = kwargs.get('dict_for_sub')
    api_url = kwargs.get('api_url')
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
    query_dict = make_query_dict(product,mission)
    print(query_dict)
    if api_url is None:
        api_url_lst = \
            satellite_dict[product]['src']['api_url']
        for url in api_url_lst:
            print('Source:',url)
            try:
                user, pw = get_credentials(remoteHostName=url)
                api = ss.SentinelAPI(user, pw, url)
                products = api.query(area=None, date=dates,**query_dict)
                break
            except Exception as e:
                print(e)
    else:
        user, pw = get_credentials(remoteHostName = api_url)
        api = ss.SentinelAPI(user, pw, api_url)
        products = api.query(area=None, date=dates,**query_dict)
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
        sort_files(path_local,filelst,product,mission)
    print ('Files downloaded to: \n', path_local)

def make_query_dict(product: str,mission: str) -> dict:
    '''
    fct to setup queries of L2 data using SentinelAPI
    '''
    level = satellite_dict[product]['mission'].get('processing_level')
    SAT = satellite_dict[product]['mission'].get(mission)
    kwargs =  {'platformname': 'Sentinel-3',
               'instrumentshortname': 'SRAL',
               'productlevel': level,
               'filename': SAT + '*WAT*'}
    return kwargs

def get_remote_files(**kwargs):
    '''
    Download swath files and store them at defined location.
    It is currently possible to download L3 altimeter data from
    CMEMS, L3 and L2P from CEDA CCI, and L2 from EUMETSAT,
    as well as L2P from aviso+ for cfosat swim data.
    '''
    dispatch_collector = {
                'cmems_L3_NRT': get_remote_files_cmems,
                'cmems_L3_s6a': get_remote_files_cmems,
                'cmems_L3_MY': get_remote_files_cmems,
                'cfo_swim_L2P': get_remote_files_aviso,
                'eumetsat_L2': get_remote_files_eumetsat,
                'cci_L2P': get_remote_files_cci,
                'cci_L3': get_remote_files_cci,
                }
    product = kwargs.get('product')
    # check if product available in dispatcher
    if product in dispatch_collector.keys():
        pass
    else:
        product = 'cmems_L3_NRT'

    dispatch_collector[product](**kwargs)
