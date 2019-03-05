#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from satellites. I try to mostly follow the PEP 
convention for python code style. Constructive comments on style and 
effecient programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# progress bar and other stuff
import sys

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os

# get_altim
if sys.version_info <= (3, 0):
    from urllib import urlretrieve, urlcleanup # python2
else:
    from urllib.request import urlretrieve, urlcleanup # python3
import gzip
import ftplib
from ftplib import FTP

# read_altim
import netCDF4 as netCDF4

# libraries for quip and quim imported in fct for speedup

# create_file
import calendar

# libraries for parallel computing
from joblib import Parallel, delayed
import multiprocessing as mp

# bintime
import math

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

import time

# get necessary paths for module
import pathfinder

# region spec
from region_specs import poly_dict, region_dict
# --- global functions ------------------------------------------------#

def progress(count, total, status=''):
    "from: https://gist.github.com/vladignatyev/06860ec2040cb497f0f3"
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()

def credentials_from_netrc():
    import netrc
    import os.path
    # get user home path
    usrhome = os.getenv("HOME")
    netrc = netrc.netrc()
    remoteHostName = "nrt.cmems-du.eu"
    user = netrc.authenticators(remoteHostName)[0]
    pw = netrc.authenticators(remoteHostName)[2]
    return user, pw

def credentials_from_txt():
    import os.path
    print("try local file credentials.txt")
    # get user home path
    usrhome = os.getenv("HOME")
    my_dict = {}
    with open(usrhome + "/credentials.txt", 'r') as f:
        for line in f:
            items = line.split('=')
            key, values = items[0], items[1]
            my_dict[key] = values
    # 1:-2 to remove the linebreak and quotation marks \n
    user = my_dict['user'][1:-2]
    pw = my_dict['pw'][1:-2]
    return user, pw

def get_credentials():
    # get credentials from .netrc
    import os.path
    usrhome = os.getenv("HOME")
    if os.path.isfile(usrhome + "/.netrc"):
        try:
            print ('Attempt to obtain credentials from .netrc')
            user, pw = credentials_from_netrc()
            return user, pw
        except AttributeError:
            print ("std copernicus user in netrc file not registered")
            print ("try local file credentials.txt")
            print ("file must contain:")
            print ("user='username'")
            print ("pw='yourpassword'")
            except_key = 1
    elif (os.path.isfile(usrhome + "/.netrc") == False
    or except_key == 1):
        try:
            user, pw = credentials_from_txt()
            return user, pw
        except:
            print("Credentials could not be obtained")

def tmploop_get_remotefiles(i,matching,user,pw,
                            server,path,sdate,
                            satpath_lustre):
    #for element in matching:
    print("File: ",matching[i])
    dlstr=('ftp://' + user + ':' + pw + '@' 
                + server + path + matching[i])
    for attempt in range(2):
        print ("Attempt to download data: ")
        try:
            print ("Downloading file")
            urlretrieve(dlstr, satpath_lustre + matching[i])
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

def get_remotefiles(satpath_ftp_014_001,destination,sdate,edate,timewin,
                    corenum,download):
    '''
    Download swath files and store them at defined location
    Example file:
    fname='global_vavh_l3_rt_s3a' 
          + '_C0028_P0624_' 
          + '20180306T174941_20180306T183759_20180306T203441.nc.gz'
    time stamps in file name stand for:
        from,to,creation
    '''
    if download is None:
        print ("No download initialized, checking local files")
        return
    # credentials
    user, pw = get_credentials()
    tmpdate = deepcopy(sdate)
    pathlst = []
    filelst = []
    while (tmpdate <= edate):
        # server and path
        if sdate >= datetime(2017,7,9):
            server='nrt.cmems-du.eu'
            path=(satpath_ftp_014_001
                 + '/'
                 + str(tmpdate.year)
                 + '/'
                 + tmpdate.strftime('%m')
                 + '/')
            print ('# ----- ')
            print ('Chosen product: ')
            print ('WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001')
            print ('# ----- ')
        else:
            sys.exit("Product not available for chosen date!")
        # get list of accessable files
        ftp = FTP(server)
        ftp.login(user, pw)
        ftp.cwd(path)
        content=FTP.nlst(ftp)
        #choose files according to verification date
        tmplst=[]
        tmpdate=sdate-timedelta(minutes=timewin)
        tmpend=edate+timedelta(minutes=timewin)
        while (tmpdate.strftime("%Y%m%dT%H")
            <= tmpend.strftime("%Y%m%dT%H")):
            matchingtmp = [s for s in content
                            if ('_'
                            + str(tmpdate.year)
                            + str(tmpdate)[5:7]
                            + str(tmpdate)[8:10]
                            + 'T'
                            + str(tmpdate)[11:13])
                        in s]
            tmplst = tmplst + matchingtmp
            tmpdate = tmpdate + timedelta(minutes=timewin)
        matching = tmplst
        print ("Download initiated")
        # Download and gunzip matching files
        print ('Downloading ' + str(len(matching)) + ' files: .... \n')
        print ("Used number of cores " + str(corenum) + "!")
        Parallel(n_jobs=corenum)(
                        delayed(tmploop_get_remotefiles)(
                        i,matching,user,pw,server,
                        path,sdate,destination
                        ) for i in range(len(matching))
                        )
        # update time
        tmpdate = tmpdate + relativedelta(months=+1)
    print ('Organizing downloaded files')
    os.system('cd /lustre/storeA/project/fou/om/altimeter '
               + '&& pwd '
               + '&& ./organize.sh')
    print ('Files downloaded to: \n' + destination)

def tmploop_get_model(
    j,sat_time_dt,model_time_dt_valid,timewin,distlim,model,
    sat_rlats,sat_rlons,sat_rHs,
    model_rlats,model_rlons,model_rHs,
    moving_win
    ):
    from utils import haversine
    lat_win = 0.1
    if (sat_time_dt[j]>=model_time_dt_valid[0]
    -timedelta(minutes=timewin) 
    and sat_time_dt[j]<=model_time_dt_valid[0]
    +timedelta(minutes=timewin)):
        if (model=='mwam8' or model=='mwam4' or model=='ARCMFC' 
        or model=='ARCMFCnew'):
            sat_rlat=sat_rlats[j]
            sat_rlon=sat_rlons[j]
            # constraints to reduce workload
            model_rlats_new=model_rlats[
                        (model_rlats>=sat_rlat-lat_win)
                        & 
                        (model_rlats<=sat_rlat+lat_win)
                        &
                        (model_rlons>=sat_rlon-moving_win)
                        &
                        (model_rlons<=sat_rlon+moving_win)
                        ]
            model_rlons_new=model_rlons[
                        (model_rlats>=sat_rlat-lat_win)
                        & 
                        (model_rlats<=sat_rlat+lat_win)
                        &
                        (model_rlons>=sat_rlon-moving_win)
                        &
                        (model_rlons<=sat_rlon+moving_win)
                        ]
            tmp=range(len(model_rlats))
            tmp_idx=np.array(tmp)[
                        (model_rlats>=sat_rlat-lat_win)
                        & 
                        (model_rlats<=sat_rlat+lat_win)
                        &
                        (model_rlons>=sat_rlon-moving_win)
                        &
                        (model_rlons<=sat_rlon+moving_win)
                        ]
            # compute distances
            distlst=map(
                    haversine,
                    [sat_rlon]*len(model_rlons_new),
                    [sat_rlat]*len(model_rlons_new),
                    model_rlons_new,model_rlats_new
                    )
            tmp_idx2 = distlst.index(np.min(distlst))
            idx_valid = tmp_idx[tmp_idx2]
            if distlst[tmp_idx2]<=distlim:
                nearest_all_dist_matches=distlst[tmp_idx2]
                nearest_all_date_matches=sat_time_dt[j]
                nearest_all_model_Hs_matches=\
                               model_rHs[0][idx_valid]
                nearest_all_sat_Hs_matches=sat_rHs[j]
                nearest_all_sat_lons_matches=sat_rlon
                nearest_all_sat_lats_matches=sat_rlat
                nearest_all_model_lons_matches=\
                                model_rlons[idx_valid]
                nearest_all_model_lats_matches=\
                                model_rlats[idx_valid]
                return nearest_all_date_matches,nearest_all_dist_matches,\
                    nearest_all_model_Hs_matches,nearest_all_sat_Hs_matches,\
                    nearest_all_sat_lons_matches, nearest_all_sat_lats_matches,\
                    nearest_all_model_lons_matches, nearest_all_model_lats_matches
            else:
                return

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
    returns idx for file
    '''
    # check if str in lst according to wished date (sdate,edate)
    idx = []
    for i in range(len(filelst)):
        element = filelst[i]
        tmp = element.find(date.strftime('%Y%m%d'))
        if tmp>=0:
            #return first index available
            idx.append(i)
    return idx[0],idx[-1]

# ---------------------------------------------------------------------#


class sentinel_altimeter():
    '''
    class to handle netcdf files containing satellite altimeter derived
    level 3 data i.e. Hs[time], lat[time], lon[time] 
    This class offers the following added functionality:
     - get swaths of desired days and read
     - get the closest time stamp(s)
     - get the location (lon, lat) for this time stamp
     - get Hs value for this time
    '''
    satpath_lustre = pathfinder.satpath_lustre
    satpath_copernicus = pathfinder.satpath_copernicus
    satpath_ftp_014_001 = pathfinder.satpath_ftp_014_001 
    from region_specs import region_dict

    def __init__(self,sdate,edate=None,timewin=None,download=None,region=None,
                corenum=None,mode=None,polyreg=None,tiling=None):
        print ('# ----- ')
        print (" ### Initializing sentinel_altimeter instance ###")
        print ('# ----- ')
        if corenum is None:
            corenum=1
        if edate is None:
            print ("Requested time: ", str(sdate))
            edate=sdate
            if timewin is None:
                timewin = int(30)
        else:
            print ("Requested time frame: " + 
                str(sdate) + " - " + str(edate))
        get_remotefiles(self.satpath_ftp_014_001,self.satpath_lustre,
                        sdate,edate,timewin,corenum,download)
        pathlst, filelst = self.get_localfilelst(sdate,edate,timewin,mode,region)
        fLATS,fLONS,fTIME,fVAVHS,fMAXS,fVAVHS_smooth = \
                                        self.read_localfiles(pathlst,mode)
        idx = np.array(range(len(fVAVHS)))[~np.isnan(fVAVHS)]
        fLATS = list(np.array(fLATS)[idx])
        fLONS = list(np.array(fLONS)[idx])
        fTIME = list(np.array(fTIME)[idx])
        fVAVHS = list(np.array(fVAVHS)[idx])
        fVAVHS_smooth = list(np.array(fVAVHS_smooth)[idx])
        gloc = [fLATS,fLONS]
        gHs = fVAVHS
        gHs_smooth = fVAVHS_smooth
        # find values for give time constraint
        ctime,cidx,timelst = self.matchtime(sdate,edate,fTIME,timewin)
        cloc = [np.array(fLATS)[cidx],np.array(fLONS)[cidx]]
        cHs = list(np.array(gHs)[cidx])
        cHs_smooth = list(np.array(gHs_smooth)[cidx])
        cTIME = list(np.array(fTIME)[cidx])
        # find values for given region
        latlst,lonlst,rlatlst,rlonlst,ridx,region = \
            self.matchregion(cloc[0],cloc[1],region=region,\
            polyreg=polyreg, tiling=tiling)
        rloc = [cloc[0][ridx],cloc[1][ridx]]
        rHs = list(np.array(cHs)[ridx])
        rHs_smooth = list(np.array(cHs_smooth)[ridx])
        rtime = list(np.array(ctime)[ridx])
        rTIME = list(np.array(cTIME)[ridx])
        self.edate = edate
        self.sdate = sdate
        self.basetime = datetime(2000,1,1)
        self.loc = rloc # regional coords [lats,lons]
        self.Hs = rHs # regional HS
        self.Hs_smooth = rHs_smooth # region smoothed ts
        self.idx = ridx # region indices
        self.dtime = rtime # region time steps as datetime obj
        self.time = rTIME # region time steps in seconds from basedate
        self.timewin = timewin
        self.region = region
        print ("Sentinel object initialized including " 
                + str(len(self.Hs)) + " footprints.")

    def get_localfilelst(self,sdate,edate,timewin,mode,region):
        if mode == 'ARCMFC':
            tmpdatestr=(self.satpath_lustre + '/monthly/')
            tmplst = np.sort(os.listdir(tmpdatestr))
            pathlst = []
            filelst = []
            tmpdate = deepcopy(sdate)
            while (tmpdate <= edate + timedelta(minutes=timewin)):
                for element in tmplst:
                    if (element.find(tmpdate.strftime("%Y%m"))>0 and \
                    element.find(region))>0:
                        pathlst.append(tmpdatestr + element)
                        filelst.append(element)
                tmpdate = tmpdate + relativedelta(months=+1)
        else:
            print ("Time window: ", timewin)
            tmpdate=deepcopy(sdate-timedelta(minutes=timewin))
            pathlst = []
            filelst = []
            while (tmpdate <= edate + timedelta(minutes=timewin)):
                tmpdatestr=(self.satpath_lustre
                            + str(tmpdate.year)
                            + '/' + tmpdate.strftime('%m')
                            + '/')
                tmplst = np.sort(os.listdir(tmpdatestr))
                filelst.append(tmplst)
                for element in tmplst:
                    pathlst.append(tmpdatestr + element)
                if (edate is not None and edate!=sdate):
                    tmpdate = tmpdate + timedelta(hours=1)
                else:
                    tmpdate = tmpdate + relativedelta(months=+1)
            filelst=np.sort(flatten(filelst))
            pathlst=np.sort(pathlst)
            idx_start,tmp = check_date(pathlst,sdate-timedelta(minutes=timewin))
            tmp,idx_end = check_date(pathlst,edate+timedelta(minutes=timewin))
            del tmp
            pathlst = np.unique(pathlst[idx_start:idx_end+1])
            filelst = np.unique(filelst[idx_start:idx_end+1])
            print (str(int(len(pathlst))) + " valid files found")
        return pathlst,filelst

    def read_localfiles(self,pathlst,mode):
        '''
        read and concatenate all data to one timeseries for each variable
        '''
        from utils import runmean
        # --- open file and read variables --- #
        LATS = []
        LONS = []
        VAVHS = []
        TIME = []
        MAXS = []
        count = 0
        print ("Processing " + str(int(len(pathlst))) + " files")
        print (pathlst[0])
        print (pathlst[-1])
        for element in pathlst:
            progress(count,str(int(len(pathlst))),'')
            try:
                # file includes a 1-D dataset with dimension time
                f = netCDF4.Dataset(element,'r')
                if mode=='ARCMFC':
                    TIME.append(f.variables['rtime'][:])
                    LATS.append(f.variables['rlats'][:])
                    LONS.append(f.variables['rlons'][:])
                    VAVH=f.variables['rHs'][:]
                    VAVHS.append(f.variables['rHs'][:])
                    MAXS.append(np.max(VAVH.astype('float')))
                    f.close()
                else:
                    tmp = f.variables['longitude'][:]
                    # transform
                    lons = ((tmp - 180) % 360) - 180
                    lats = f.variables['latitude'][:]
                    time = f.variables['time'][:]
                    VAVH = f.variables['VAVH'][:]
                    f.close()
                    LATS.append(lats)
                    LONS.append(lons)
                    TIME.append(time)
                    VAVHS.append(VAVH)
                    MAXS.append(np.max(VAVH.astype('float')))
            except (IOError):
                print ("No such file or directory")
            count = count + 1
        print ('\n')
        # apply flatten to all lists
        fTIME,fLATS,fLONS,fVAVHS,fMAXS=\
            flatten(TIME),flatten(LATS),flatten(LONS),\
            flatten(VAVHS),MAXS
        fTIME_unique,indices=np.unique(fTIME,return_index=True)
        fLATS=list(np.array(fLATS)[indices])
        fLONS=list(np.array(fLONS)[indices])
        fTIME=list(np.array(fTIME)[indices])
        fVAVHS=list(np.array(fVAVHS)[indices])
        # smooth Hs time series
        fVAVHS_smooth,fVAVHS_std = runmean(fVAVHS,5,'centered')
        return fLATS, fLONS, fTIME, fVAVHS, fMAXS, fVAVHS_smooth

    def bintime(self,binframe=None):
        '''
        fct to return frequency of occurrence per chose time interval
        e.g. days, weeks, months, years ...
        currently only daily bins
        '''
        listofdatetimes=self.rtime
        listofdatetimes_sec=self.rTIME
        sdate=listofdatetimes[0]
        edate=listofdatetimes[-1]
        trange=int(math.ceil((edate-sdate).total_seconds()/60./60./24.))
        freqlst=[]
        datelst=[]
        dincr=0
        timewin = 0
        for i in range(trange):
            lodt=sdate+timedelta(days=dincr)
            ctime,cidx,timelst = self.matchtime(
                                    datetime(lodt.year,lodt.month,lodt.day),
                                    datetime(lodt.year,lodt.month,lodt.day)
                                    +timedelta(days=1),
                                    listofdatetimes_sec,
                                    timewin = timewin)
            count=len(np.array(self.rHs)[
                                cidx][
                                ~np.isnan(np.array(self.rHs)[cidx])
                                     ]
                    )
            freqlst.append(count)
            dincr=dincr+1
            datelst.append(lodt)
        return freqlst,datelst

    def matchtime(self,sdate,edate,fTIME,timewin=None):
        '''
        fct to obtain the index of the time step closest to the 
        requested time step from buoy and forecast including the 
        respective time stamp(s). Similarily, indices are chosen
        for the time and defined region.
        '''
        if timewin is None:
            timewin = 0
        basetime=datetime(2000,1,1)
        # create list of datetime instances
        timelst=[]
        ctime=[]
        cidx=[]
        idx=0
        print ('Time window is: ', timewin)
        if (edate is None or sdate==edate):
            for element in fTIME:
                tmp=basetime + timedelta(seconds=element)
                timelst.append(tmp)
                # choose closest match within window of win[minutes]
                if (tmp >= sdate-timedelta(minutes=timewin) 
                and tmp < sdate+timedelta(minutes=timewin)):
                    ctime.append(tmp)
                    cidx.append(idx)
                del tmp
                idx=idx+1
        if (edate is not None and edate!=sdate):
            for element in fTIME:
                tmp=basetime + timedelta(seconds=element)
                timelst.append(tmp)
                if (tmp >= sdate-timedelta(minutes=timewin)
                and tmp < edate+timedelta(minutes=timewin)):
                    ctime.append(tmp)
                    cidx.append(idx)
                del tmp
                idx=idx+1
        return ctime, cidx, timelst

    def matchregion(self,LATS,LONS,region,polyreg,tiling):
        # find values for given region
        if polyreg is None:
            if region is None:
                region = 'Global'
            if ~isinstance(region,str)==True:
                print ("Manuall specified region "
                    + [llcrnrlat,urcrnrlat,llcrnrlon,urcrnrlon] + ": \n"
                    + " --> Bounds: " + str(region))
            else:
                if region not in region_dict.keys():
                    sys.exit("Region is not defined")
                else:
                    print ("Specified region: " + region + "\n"
                      + " --> Bounds: " + str(self.region_dict[region]))
            latlst,lonlst,rlatlst,rlonlst,ridx = \
                self.matchregion_prim(LATS,LONS,region=region)
        else:
            region = polyreg
            latlst,lonlst,rlatlst,rlonlst,ridx = \
                self.matchregion_poly(LATS,LONS,region=region,tiling=tiling)
        return latlst, lonlst, rlatlst, rlonlst, ridx, region

    def matchregion_prim(self,LATS,LONS,region):
        if (region is None or region == "Global"):
            region = "Global"
            latlst = LATS
            lonlst = LONS
            rlatlst= LATS
            rlonlst= LONS
            ridx = range(len(LATS))
        else:
            if isinstance(region,str)==True:
                if (region=='ARCMFC' or region=='Arctic'):
                    boundinglat = self.region_dict[region]["boundinglat"]
                else:
                    llcrnrlat = self.region_dict[region]["llcrnrlat"]
                    urcrnrlat = self.region_dict[region]["urcrnrlat"]
                    llcrnrlon = self.region_dict[region]["llcrnrlon"]
                    urcrnrlon = self.region_dict[region]["urcrnrlon"]
            else:
                llcrnrlat = region[0]
                urcrnrlat = region[1]
                llcrnrlon = region[2]
                urcrnrlon = region[3]
            # check if coords in region
            latlst=[]
            lonlst=[]
            rlatlst=[]
            rlonlst=[]
            ridx=[]
            for i in range(len(LATS)):
                latlst.append(LATS[i])
                lonlst.append(LONS[i])
                if (
                (region=='ARCMFC' or region=='Arctic')
                and LATS[i] >= boundinglat
                    ):
                    rlatlst.append(LATS[i])
                    rlonlst.append(LONS[i])
                    ridx.append(i)
                elif(
                region != 'ARCMFC' and region != 'Arctic'
                and LATS[i] >= llcrnrlat
                and LATS[i] <= urcrnrlat
                and LONS[i] >= llcrnrlon
                and LONS[i] <= urcrnrlon
                    ):
                    rlatlst.append(LATS[i])
                    rlonlst.append(LONS[i])
                    ridx.append(i)
        if not ridx:
            print ("No values for chosen region and time frame!!!")
        else:
            print ("Values found for chosen region and time frame.")
        return latlst, lonlst, rlatlst, rlonlst, ridx

    def matchregion_poly(self,LATS,LONS,region,tiling):
        from matplotlib.patches import Polygon
        from matplotlib.path import Path
        import numpy as np
        from model_specs import model_dict
        if (region not in poly_dict.keys() \
            and region not in model_dict.keys()):
            sys.exit("Region is not defined")
        elif isinstance(region,dict)==True:
            print ("Manuall specified region: \n"
                + " --> Bounds: " + str(region))
            poly = Polygon(list(zip(region['lons'],
                region['lats'])), closed=True)
        elif (isinstance(region,str)==True and region in model_dict.keys()):
            from modelmod import get_model
            if region == 'mwam8':
                grid_date = datetime(2019,2,1,6)
            elif (region == 'MoskNC' or region == 'MoskWC'):
                grid_date = datetime(2018,3,1)
            else:
                grid_date = datetime(2019,2,1)
            if region == 'ARCMFC':
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc", model=region, fc_date=grid_date,
                    init_date=grid_date)
            else:
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(simmode="fc", model=region, fc_date=grid_date,
                    leadtime=0)
            def ggt(a, b):
                if a < b: a,b = b,a
                while a%b != 0:
                    a,b = b,a%b
                return b
            # tiling of model domain
            in1,in2 = model_lats.shape[0],model_lats.shape[1]
            if isinstance(tiling,int):
                g = tiling
            else:
                if ggt(in1,in2)==1:
                    g = ggt(in1-1,in2-1)
                else:
                    g = ggt(in1,in2)
            xtilesize = model_lats.shape[0]/g
            ytilesize = model_lats.shape[1]/g
            xidx = []
            for i in range(g):
                xidx.append(xtilesize*i)
            xidx.append(model_lats.shape[0])
            yidx = []
            for i in range(g):
                yidx.append(ytilesize*i)
            yidx.append(model_lats.shape[1])
            print(xidx)
            print(yidx)
            tiles = []
            for i in range(g):
                for j in range(g):
                    tiles.append([model_lons[xidx[i]:xidx[i+1],
                                yidx[j]:yidx[j+1]],
                                model_lats[xidx[i]:xidx[i+1],
                                yidx[j]:yidx[j+1]]])
            # create polygon for each tile
            rlatlst = []
            rlonlst = []
            ridx = []
            for i in range(len(tiles)):
                model_lons = tiles[i][0]
                model_lats = tiles[i][1]
                # create polygon of model domain
                idx = 1
                lontmp1=list(model_lons[:,0])
                lontmp2=list(model_lons[:,-1])
                lontmp3=list(model_lons[0,:])
                lontmp4=list(model_lons[-1,:])
                mlons = (lontmp3[::-1][::idx] + lontmp1[::idx] \
                    + lontmp4[::idx] + lontmp2[::-1][::idx] \
                    + [lontmp3[::-1][0]] + [lontmp3[::-1][0]])
                lattmp1=list(model_lats[:,0])
                lattmp2=list(model_lats[:,-1])
                lattmp3=list(model_lats[0,:])
                lattmp4=list(model_lats[-1,:])
                mlats = (lattmp3[::-1][::idx] + lattmp1[::idx] \
                    + lattmp4[::idx] + lattmp2[::-1][::idx] \
                    + [lattmp3[::-1][0]] + [lattmp3[::-1][0]])
                poly = Polygon(list(zip(mlons,mlats)), closed=True)
                # check if coords in region
                LATS = list(LATS)
                LONS = list(LONS)
                latlst = LATS
                lonlst = LONS
                lats = np.array(LATS).ravel()
                lons = np.array(LONS).ravel()
                points = np.c_[lons,lats]
                hits = Path(poly.xy).contains_points(points)
                rlatlst.append(list(np.array(LATS)[hits]))
                rlonlst.append(list(np.array(LONS)[hits]))
                ridx.append(list(np.array(range(len(LONS)))[hits]))
            rlatlst = flatten(rlatlst)
            rlonlst = flatten(rlonlst)
            ridx = flatten(ridx)
        elif isinstance(region,str)==True:
            print ("Specified region: " + region + "\n"
              + " --> Bounded by polygon: \n"
              + "lons: " + str(poly_dict[region]['lons']) + "\n"
              + "lats: " + str(poly_dict[region]['lats']))
            poly = Polygon(list(zip(poly_dict[region]['lons'],
                poly_dict[region]['lats'])), closed=True)
            # check if coords in region
            LATS = list(LATS)
            LONS = list(LONS)
            latlst = LATS
            lonlst = LONS
            lats = np.array(LATS).ravel()
            lons = np.array(LONS).ravel()
            points = np.c_[lons,lats]
            hits = Path(poly.xy).contains_points(points)
            rlatlst = list(np.array(LATS)[hits])
            rlonlst = list(np.array(LONS)[hits])
            ridx = list(np.array(range(len(LONS)))[hits])
        if not ridx:
            print ("No values for chosen region and time frame!!!")
        else:
            print ("Values found for chosen region and time frame.")
        return latlst, lonlst, rlatlst, rlonlst, ridx


def get_pointsat(sa_obj,station=None,lat=None,lon=None,distlim=None):
    from utils import haversine
    from stationlist import locations
    if ((lat is None or lon is None) and (station is None)):
        print ("location is missing")
    if (lat is None or lon is None):
        lat=locations[station][0]
        lon=locations[station][1]
    lats = sa_obj.loc[0]
    lons = sa_obj.loc[1]
    Hs = sa_obj.Hs
    time = sa_obj.time
    sample = []
    dists = []
    for i in range(len(lats)):
        if ((lats[i] < lat+1) and (lats[i] > lat-1)):
            if ((lons[i] < lon+1) and (lons[i] > lon-1)):
                dist=haversine(lons[i],lats[i],lon,lat)
                if (dist<=distlim):
                    sample.append(Hs[i])
                    dists.append(dist)
    return sample, dists
