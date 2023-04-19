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
obs_dict = load_or_default('insitu_specs.yaml')

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
    nID = kwargs.get('nID','gullfaks_local')
    path_local = kwargs.get('path_local')
    dict_for_sub = kwargs.get('dict_for_sub')
    # check if search str template
    file_search_template = \
        obs_dict[nID]['src'].get('file_search_template',\
                                            '%Y%m%dT%H')
    # credentials
    server = obs_dict[nID]['src']['server']
    user, pw = get_credentials(remoteHostName = server)
    tmpdate = deepcopy(sdate)
    filesort = False
    path_template_src = obs_dict[nID]['src']\
                                  ['path_template_remote']
    strsublst_src = obs_dict[nID]['src']['strsub']
    subdict_src = make_subdict(strsublst_src,
                               class_object_dict=dict_for_sub)
    while (tmpdate <= edate):
        # create remote path
        path_remote = make_pathtofile(path_template_src,\
                                      strsublst_src,subdict_src,\
                                      date=tmpdate)
        if path_local is None:
            # create local path
            path_template_dst = obs_dict[nID]['dst']\
                                          ['path_template']
            strsublst_dst = obs_dict[nID]['dst']\
                                          ['strsub']
            subdict_dst = make_subdict(strsublst_dst,
                                           class_object_dict=dict_for_sub)
            path_local = make_pathtofile(path_template_dst,\
                                         strsublst_dst,subdict_dst,\
                                         date=tmpdate)
            filesort = True

        print ('# ----- ')
        print ('Chosen source: ')
        print (nID + ' values from ' + obs_dict[nID]['product'] + ': ' + server)
        print(path_remote)
        print ('# ----- ')
        # get list of accessable files
        print('HERE1')
        ftp = FTP(server)
        print('HERE2')
        ftp.login(user, pw)
        print('HERE3')
        print(ftp)
        print(ftp.login)
        print(path_remote)
        ftp.cwd(path_remote)
        print('HERE4')
        content=FTP.nlst(ftp)
        print('HERE5')
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
        date_incr = obs_dict[product]['src'].get('date_incr', 'm')
        tmpdate = date_dispatcher(tmpdate, date_incr=date_incr)
    if filesort is True:
        # sort files
        print("Data is being sorted into subdirectories " \
            + "year and month ...")
        filelst = [f for f in os.listdir(path_local)
                    if os.path.isfile(os.path.join(path_local,f))]
        sort_files(path_local,filelst,product,mission)
    print ('Files downloaded to: \n', path_local)


def get_remote_files(**kwargs):
    '''
    Download insitu files and store them at defined location.
    '''
    dispatch_collector = {
                'cmems': get_remote_files_cmems,
                }
    product = kwargs.get('product')
    # check if product available in dispatcher
    if product in dispatch_collector.keys():
        pass
    else:
        product = 'cmems_L3_NRT'

    dispatch_collector[product](**kwargs)
