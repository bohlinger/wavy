#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module comprises functions to collect remote insitu data.
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import sys
import numpy as np
from datetime import datetime, timedelta
import os
from copy import deepcopy
import time
from urllib.request import urlretrieve, urlcleanup
from urllib.parse import quote
from ftplib import FTP
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed
import logging
import copernicusmarine as cmc
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=30)
logger = logging.getLogger(__name__)

# own imports
from wavy.utils import sort_files
from wavy.utils import make_pathtofile, make_subdict
from wavy.utils import date_dispatcher
from wavy.credentials import get_credentials
from wavy.wconfig import load_or_default
# ---------------------------------------------------------------------#

# read yaml config files:
insitu_dict = load_or_default('insitu_cfg.yaml')

# --- def functions ---------------------------------------------------#

def tmploop_get_remote_files(i: int, matching: str,
                             user: str, pw: str,
                             server: str, remote_path: str,
                             path_local: str):
    """
    Function to download files using ftp. Tries 10 times before failing.
    """
    print("File: ", matching[i])
    print("src path: ", remote_path)
    pw = quote(pw)  # to escape special characters
    dlstr = ('ftp://' + user + ':' + pw + '@'
             + server + remote_path + matching[i])
    for attempt in range(10):
        print(attempt, "attempt to download data: ")
        try:
            print("Downloading file")
            urlretrieve(dlstr, os.path.join(path_local, matching[i]))
            urlcleanup()
        except Exception as e:
            print(e.__doc__)
            print(e.message)
            print("Waiting for 10 sec and retry")
            time.sleep(10)
        else:
            break
    else:
        print('An error was raised and I ' +
              'failed to fix problem myself :(')
        print('Exit program')
        sys.exit()

def get_remote_files_cmems(**kwargs):
    '''
    Insitu files from CMEMS and store them at defined
    location. Time stamps in file name stand for:

    from, to, creation
    '''
    cfg = kwargs.get('cfg')
    product = kwargs.get('nID')
    sdate = kwargs.get('sd')
    edate = kwargs.get('ed')
    twin = int(np.max([kwargs.get('twin', 30), 30]))
    nproc = kwargs.get('nproc', 1)
    name = kwargs.get('name')
    dict_for_sub = kwargs

    # define path
    path = kwargs.get('path', None)

    # check if search str template
    file_search_template = cfg.download['ftp']\
        .get('search_str', '%Y%m%dT').replace('name',
                                              name)

    # credentials
    server = insitu_dict[product]['download']['ftp']['server']
    user, pw = get_credentials(remoteHostName=server)

    # create paths
    tmpdate = deepcopy(sdate)
    path_template_src = cfg.download['ftp']['src_tmplt']
    strsublst_src = cfg.download['ftp']['strsub']
    subdict_src = make_subdict(strsublst_src,
                               class_object_dict=dict_for_sub)
    while (tmpdate <= edate):
        try:
            # create remote path
            path_remote = make_pathtofile(path_template_src,
                                          strsublst_src, subdict_src,
                                          date=tmpdate)

            if path is None:
                # create local path
                path_template_dst = cfg.download['ftp']['trgt_tmplt']
                strsublst_dst = cfg.download['ftp']['strsub']
                subdict_dst = make_subdict(strsublst_dst,
                                           class_object_dict=dict_for_sub)
                path_local = make_pathtofile(path_template_dst,
                                             strsublst_dst, subdict_dst,
                                             date=tmpdate)
            else:
                path_local = path

            print('# ----- ')
            print('Chosen source: ')
            print(name + ' values from ' + product + ': ' + server)
            print(path_remote)
            print('# ----- ')
            # get list of accessable files
            ftp = FTP(server)
            ftp.login(user, pw)
            ftp.cwd(path_remote)
            content = FTP.nlst(ftp)

            # choose files according to sdate/edate
            tmplst = []
            tmpdate_new = tmpdate-timedelta(minutes=twin)
            tmpdate_end = edate+timedelta(minutes=twin)
            while (tmpdate_new <= tmpdate_end):
                matchingtmp = [s for s in content
                               if tmpdate_new.strftime(file_search_template)
                               in s]

                tmplst = tmplst + matchingtmp
                tmpdate_new = tmpdate_new + timedelta(minutes=twin)
            matching = np.unique(tmplst)
            print(matching)

            # check if download path_local exists if not create
            if not os.path.exists(path_local):
                os.makedirs(path_local, exist_ok=True)

            # Download matching files
            print('Downloading ' + str(len(matching))
                  + ' files: .... \n')
            print("Used number of possible simultaneous downloads "
                  + str(nproc) + "!")
            Parallel(n_jobs=nproc)(
                            delayed(tmploop_get_remote_files)(
                                i, matching, user, pw, server,
                                path_remote, path_local
                                ) for i in range(len(matching))
                            )
        except Exception as e:
            logger.exception(e)
        # update time
        path_date_incr_unit = cfg.download['ftp']\
            .get('path_date_incr_unit', 'm')
        path_date_incr = cfg.download['ftp']\
            .get('path_date_incr', 1)
        tmpdate = date_dispatcher(tmpdate,
                                  path_date_incr_unit, path_date_incr)
        print('####################################')
        print(path_local)
        print('####################################')

    print('Files downloaded to: \n', path_local)
    
    
def get_remote_files_copernicusmarine(**kwargs):
    '''
    Download swath files from CMEMS using copernicusmarine parckage
    and store them at defined location. Time stamps in file name stand for:

    from, to, creation
    '''
    product = kwargs.get('nID')
    sdate = kwargs.get('sd')
    edate = kwargs.get('ed')
    nproc = kwargs.get('nproc', 1)
    name = kwargs.get('name', 'Draugen')
    dict_for_sub = kwargs
    # define path
    path = kwargs.get('path', None)
    # Get time increment
    time_incr = insitu_dict[product]['download']['copernicus']\
                .get('time_incr','h')
    
    # Chose search template for time given time_incr
    if time_incr=='m':
        file_search_template = '%Y%m'
        dataset_part = 'monthly'
    print('Date search format:', file_search_template)

    # Get dataset_id
    dataset_id = insitu_dict\
                            [product]['download']['copernicus']\
                            ['dataset_id']
    strsublst_src = insitu_dict[product]['download']\
                            ['copernicus']['strsub']
    subdict_src = make_subdict(strsublst_src,
                               class_object_dict=dict_for_sub)

    # replace name of the mission in dataset_id
    dataset_id = make_pathtofile(dataset_id,
                                  strsublst_src, 
                                  subdict_src)
    
    # Initialize start date to match original files time increment
    tmpdate = deepcopy(sdate)
    tmpdate_end = deepcopy(edate)


    # tmpdate = tmpdate - timedelta(hours=1)

    # tmpdate_end = tmpdate_end + timedelta(hours=1)
 
    try:
        
        print('# ----- ')
        print('Chosen source: ')
        print(name + ' values from ' + product + ': ' + 'copernicusmarine')
        print('# ----- ')
        
        while (tmpdate <= tmpdate_end):

            if path is None:
	        # create local path
                path_template_dst = insitu_dict[product]['download']\
	                                ['copernicus']['trgt_tmplt']
                strsublst_dst = insitu_dict[product]['download']\
	                                ['copernicus']['strsub']
                subdict_dst = make_subdict(strsublst_dst,
	                                   class_object_dict=dict_for_sub)
                path_local = make_pathtofile(path_template_dst,
	                                     strsublst_dst, subdict_dst,
	                                     date=tmpdate)
            else:
                path_local = path
	
            print('* --------------')
            print('Downloading for date:', tmpdate)
            print('* --------------')

            # check if download path_local exists if not create
            if not os.path.exists(path_local):
                os.makedirs(path_local, exist_ok=True)
    
            # Create regexp filter
            tmpdate_str = tmpdate.strftime(file_search_template)
            regexp_tmp = "*MO_{}_{}*.nc".format(name, tmpdate_str)
            # Fetch data corresponding to tmp date 
            print(regexp_tmp)
            try:
                cmc.get(
                        dataset_id = dataset_id,
                        dataset_part = 'monthly',
                        filter = regexp_tmp,
                        no_directories = True,
                        output_directory = path_local,
                        force_download=True,
                        overwrite_output_data=True)
            except Exception as e:
                print(e)
                pass
            
            if time_incr=='m':
                tmpdate = tmpdate + relativedelta(months=+1)    
     
    except Exception as e:
        logger.exception(e)
   
    print('# -----------------------------------')
    print('Files downloaded to: \n', path_local)        

