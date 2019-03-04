#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and process wave
field related data from satellites. I try to mostly follow the PEP 
convention for python code style. Constructive comments on style and 
effecient programming are most welcome!

Future plans involve expanding the help function and adding 
functionalities such that the program can be executed in the shell
directly by only writing e.g. 
./satmod.py [-option] [sdate] [edate] [model]
where the options can be e.g.   -h for help     (existing)
                                -p for plot     (wish list)
                                -s for save     (wish list)
                                -d for download (wish list)
                                - ...
'''
__version__ = "0.5.0"
__author__="Patrik Bohlinger, Norwegian Meteorological Institute"
__maintainer__ = "Patrik Bohlinger"
__email__ = "patrikb@met.no"
__status__ = "under development with operation ARCMFC branch"

# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# progress bar and other stuff
import sys

# ignore irrelevant warnings from matplotlib for stdout
import warnings
#warnings.filterwarnings("ignore")

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

def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    '''
    Add a vertical color bar to an image plot.
    '''
    from mpl_toolkits import axes_grid1
    import matplotlib.pyplot as plt
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, 
                                    aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)

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

def prod_orbcyc_nc(outpath=None,sdate=None,edate=None,to_nc=None):
    """
    Read and store one orbital repeat cycle and dump to file
    Repeat cycle is 27 days
    """
    if (sdate is None or edate is None):
        sdate=datetime(2017,9,1)
        edate=sdate+timedelta(days=27)
    sa_obj = sa(sdate,edate=edate,timewin=30,region="Global",\
                mode="ARCMFC")
    if to_nc==True:
        sa_obj.dumptonc(outpath,ncmode='orbcyc')
    return sa_obj

def learn_orbpos(sa_obj):
    """
    - Read in the orbital cycle netcdf file
    - Model the cycle using a Gaussian Process model comparable to Mauna 
        Loa example on scikit page
    """
    # learn lats
    # learn lons
    # retrun scikit object
    return

def pred_relative_date(sa_obj,datelst):
    from stationmod import matchtime
    #sa_obj.rtime[0]
    return

def pred_setpos(sa_obj,sdate,edate,discr=1):
    """
    predict lats and lon for a given period (sdate to edate) and 
    discretization in seconds (disrc)
    """
    tmpdate = sdate
    t_discr = []
    while tmpdate<=edate:
        t_discr.append(tmpdate)
        tmpdate = tmpdate + timedelta(seconds=1)
    pred_relative_date(sa_obj,t)
    return

def get_sentpos(time):
    """
    - Fct to obtain position of Sentinel 3a footprint based on
    the orbital cycle, information from the manual under:
    https://sentinel.esa.int/documents/247904/685236/ + 
        Sentinel-3_User_Handbook
    --> Orbital repeat cycle: 27 days 
        (14+7/27 orbits per day, 385 orbits per cycle)
    - Predicted position will be in netcdf file together with
      the associated i,j indices for model product as function of time
    """
    return lat,lon
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
    # from 20170709, depricated
    satpath_ftp_008_052 = pathfinder.satpath_ftp_008_052
    # from 20180320, now valid from start!
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
        self.rloc = rloc # regional coords [lats,lons]
        self.cloc = cloc # close time coords [lats,lons]
        self.gloc = gloc # global coords [lats,lons]
        self.fTIME = fTIME # custom time steps in seconds
        self.cTIME = cTIME # custom time steps in seconds
        self.ctime = ctime # custom time steps as datetime obj
        self.cidx = cidx # adjacent indices
        self.gHs = gHs
        self.cHs = cHs
        self.rHs = rHs
        self.rHs_smooth = rHs_smooth
        self.ridx = ridx # region indices
        self.rtime = rtime # region time steps as datetime obj
        self.rTIME = rTIME # region time steps in seconds
        self.timewin = timewin
        self.gHsMax = np.array(fMAXS)
        self.region = region
        print ("Sentinel object initialized including " 
                + str(len(self.rHs)) + " footprints.")

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

    def quim(self,region=None):
        # ignore irrelevant warnings from matplotlib for stdout
        import warnings
        warnings.filterwarnings("ignore")
        from mpl_toolkits.basemap import Basemap, cm
        if (region is None or region == 'Global'):
            # Mollweide
            m = Basemap(
                    projection='moll',lon_0=0,resolution='c',
                    area_thresh=10000
                    )
            # labels = [left,right,top,bottom]
            m.drawparallels(
                        np.arange(-80.,81.,20.),
                        labels=[True,False,False,False]
                        )
        elif (region=='Arctic' or region=='ARCMFC'):
            # Polar Stereographic Projection
            m = Basemap(
                projection='npstere',
                boundinglat=self.region_dict[region]["boundinglat"],
                lon_0=0,
                resolution='l',area_thresh=100
                )
            m.drawparallels(
                    np.arange(40.,81.,10.),
                    labels=[True,False,False,False]
                    )
            m.drawmeridians(
                    np.arange(-180.,181.,10.),
                    labels=[False,False,False,True]
                    )
        elif (region=='Moskenes' or region=='Sulafj' or 
            region=='mwam4' or region=='mwam8' or region=='man' or
            region=='Mosk_dom'):   
            # Lambert Conformal Projection over Moskenes
            m = Basemap(
                llcrnrlon=self.region_dict[region]["llcrnrlon"],
                llcrnrlat=self.region_dict[region]["llcrnrlat"],
                urcrnrlon=self.region_dict[region]["urcrnrlon"],
                urcrnrlat=self.region_dict[region]["urcrnrlat"],
                projection='lcc', resolution='f',area_thresh=1,
                lat_1=67.2,lat_2=68,lat_0=67.6,lon_0=12.5
                )
            if (region=='Moskenes' or region=='Sulafj'):
                m.drawparallels(
                    np.arange(40.,81.,1.),
                    labels=[True,False,False,False]
                    )
                m.drawmeridians(
                    np.arange(-180.,181.,1.),
                    labels=[False,False,False,True]
                    )
            elif (region=='mwam4' or region=='mwam8' or region=='man'):
                m.drawparallels(
                    np.arange(40.,81.,5.),
                    labels=[True,False,False,False]
                    )
                m.drawmeridians(
                    np.arange(-180.,181.,5.),
                    labels=[False,False,False,True]
                    )
        m.drawcoastlines()
        #m.fillcontinents(color='navajowhite')
        return m

    def quip(self,region=None,show=None,save=None, outpath=None):
        # ignore irrelevant warnings from matplotlib for stdout
        import warnings
        warnings.filterwarnings("ignore")
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap, cm
        if region is None:
            region = 'Global'
            loc = self.cloc
            Hs = self.cHs
        else:
            loc = self.rloc
            Hs = self.rHs
        if outpath is None:
            outpath = 'outpath/'
        cmap=cm.GMT_haxby
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        m = self.quim(region)
        x, y = m(loc[1],loc[0])
        sizedict = {'Global':10,'ARCMFC':10,'Moskenes':20,'Sulafj':20,
                    'mwam4':20,'mwam8':20,'man':20,'Mosk_dom':20}
        sc=m.scatter(
                    x,y,sizedict[region],c=Hs,marker='o',linewidth='.0',
                    cmap=cmap, vmin=0, vmax=np.nanmax(Hs)
                    )
        cbar=add_colorbar(sc)
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label('m',size=10)
        ax.set_title('Altimeter Hs for: \n' 
                    + str(self.sdate) 
                    + ' to '
                    + str(self.edate),
                    fontsize=10
                    )
        if save == True:
            os.system('mkdir -p ' + outpath)
            plt.savefig(
                        outpath
                        + 'altimeter_' 
                        + region 
                        + "_" 
                        + self.sdate.strftime("%Y%m%d%H%M%S") 
                        + "_" 
                        + self.edate.strftime("%Y%m%d%H%M%S")
                        + '.pdf', 
                        format='pdf'
                        )
        if show == True:
            plt.show()

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

    def discont_on_swath():
        '''
        check for differences of waves during a swath to check for 
        anomalies connected to Moskstraumen
        '''

    def plotavail(self,btr_dates,btr_freq,show=None,save=None):
        '''
        btr_ = bintime result
        '''
        # ignore irrelevant warnings from matplotlib for stdout
        import warnings
        warnings.filterwarnings("ignore")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(20,5))
        fig.suptitle("Number of footprints in "
                    + self.region 
                    + " region from " 
                    + str(btr_dates[0])  
                    + " to " 
                    + str(btr_dates[-1]))
        ax = fig.add_subplot(111)
        ax.plot(btr_dates,btr_freq,'ro')
        ax.plot(btr_dates,btr_freq,'r',linewidth=0.5)
        # format your data to desired format. Here I chose YYYY-MM-DD but 
        # you can set it to whatever you want.
        import matplotlib.dates as mdates
        ax.xaxis_date()
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
        ax.set_ylabel('occurrence frequency')
        # rotate and align the tick labels so they look better
        fig.autofmt_xdate()
        if save == True:
            plt.savefig(
                    'altimeter_occfreq_'
                    + self.region
                    + "_" 
                    + btr_dates[0].strftime("%Y%m%d%H%M%S") 
                    + "_" 
                    + btr_dates[-1].strftime("%Y%m%d%H%M%S")
                    + '.pdf', 
                    format='pdf'
                    )
        if show == True:
            plt.show()

    def quipdiag(self):
        '''
        Create diagnostic figures tailored to chosen region.
        To come in future!
        '''
        # ignore irrelevant warnings from matplotlib for stdout
        import warnings
        warnings.filterwarnings("ignore")
        import matplotlib.pyplot as plt

        basetime=datetime(2000,1,1)
        timelst=[]
        for i in range(len(self.rtime)):
            timelst.append(basetime+timedelta(seconds=self.rtime[i]))
        #bin dates
        #plt.plot()

    def grid(self):
        '''
        Grid data for a given time frame. Result is something like a 
        mean wave field based on satellite coverage.
        To come in future!
        '''

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
                    tiles.append([model_lons[xidx[i]:xidx[i+1],yidx[j]:yidx[j+1]],model_lats[xidx[i]:xidx[i+1],yidx[j]:yidx[j+1]]])
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
#            # create polygon of model domain
#            idx = 1 # sparseness
#            lontmp1=list(model_lons[:,0])
#            lontmp2=list(model_lons[:,-1])
#            lontmp3=list(model_lons[0,:])
#            lontmp4=list(model_lons[-1,:])
#            mlons = (lontmp3[::-1][::idx] + lontmp1[::idx] \
#                    + lontmp4[::idx] + lontmp2[::-1][::idx] \
#                    + [lontmp3[::-1][0]] + [lontmp3[::-1][0]])
#            lattmp1=list(model_lats[:,0])
#            lattmp2=list(model_lats[:,-1])
#            lattmp3=list(model_lats[0,:])
#            lattmp4=list(model_lats[-1,:])
#            mlats = (lattmp3[::-1][::idx] + lattmp1[::idx] \
#                    + lattmp4[::idx] + lattmp2[::-1][::idx] \
#                    + [lattmp3[::-1][0]] + [lattmp3[::-1][0]])
#            poly = Polygon(list(zip(mlons,mlats)), closed=True)
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

    def dumptonc(self,outpath,ncmode=None,timeframe=None):
        """
        1. check if nc file already exists
        2. - if so use append mode
           - if not create file
        3. make sure files are only for one single month
        """
        # 1. check if nc file already exists
        win = int(30)
        if self.edate is None:
            edate=self.sdate + timedelta(minutes=win)
            sdate=edate-timedelta(minutes=win)
        else:
            sdate=self.sdate
            edate=self.edate
        if (ncmode is None or ncmode == 'man'):
            filename = ("altimeter_" 
                    + self.region
                    + "_" 
                    + sdate.strftime("%Y%m%d%H%M%S") 
                    + "_" 
                    + edate.strftime("%Y%m%d%H%M%S")
                    + ".nc")
        elif (ncmode == 'auto'):
            if (timeframe == 'monthly' or timeframe is None):
                filename = ("global_vavh_l3_rt_s3a_"
                    + self.region
                    + "_" 
                    + sdate.strftime("%Y%m")
                    + ".nc")   
            if timeframe == 'daily':
                filename = ("global_vavh_l3_rt_s3a_"
                    + self.region
                    + "_" 
                    + sdate.strftime("%Y%m%d")
                    + ".nc")
        elif (ncmode == 'orbcyc'):
            filename = ("orbcyc_"
                    + self.region
                    + "_"
                    + sdate.strftime("%Y%m%d%H%M%S")
                    + "_"
                    + edate.strftime("%Y%m%d%H%M%S")
                    + ".nc")
        fullpath = outpath + filename
        if os.path.isfile(fullpath):
            print ('Dump altimeter wave data to file: ' + fullpath)
            nc = netCDF4.Dataset(
                            fullpath,mode='a',
                            clobber=False
                            )
            # variables
            #ncrtime=nc.variables['rtime'][:]
            startidx = len(nc['rtime'])
            endidx = len(nc['rtime'])+len(self.rTIME)
            nc.variables['rtime'][startidx:endidx] = self.rTIME[:]
            nc.variables['rHs'][startidx:endidx] = self.rHs[:]
            nc.variables['rlons'][startidx:endidx] = self.rloc[1][:]
            nc.variables['rlats'][startidx:endidx] = self.rloc[0][:]
        else:
            os.system('mkdir -p ' + outpath)
            print ('Dump altimeter wave data to file: ' + fullpath)
            nc = netCDF4.Dataset(
                            fullpath,mode='w',
                            format='NETCDF4'
                            )
            nc.title = 'altimeter significant wave height'
            rtimerange=len(self.ridx)
            dimsize = None
            # dimensions
            dimrtime = nc.createDimension(
                                    'rtime', 
                                    size=dimsize
                                    )
            # variables
            ncrtime = nc.createVariable(
                                   'rtime',
                                   np.float64, 
                                   dimensions=('rtime')
                                   )   
            ncrlats = nc.createVariable(
                                   'rlats',
                                   np.float64, 
                                   dimensions=('rtime')
                                   )        
            ncrlons = nc.createVariable(
                                   'rlons',
                                   np.float64, 
                                   dimensions=('rtime')
                                   )        
            ncrHs = nc.createVariable(
                                   'rHs',
                                   np.float64, 
                                   dimensions=('rtime')
                                   )        
 
            # generate time for netcdf file
            basetime=datetime(2000,1,1)
            ncrtime.units = 'seconds since 2000-01-01 00:00:00'
            ncrtime[:] = self.rTIME
            ncrHs.units = 'm'
            ncrHs[:] = self.rHs
            ncrlons.units = 'degrees east'
            ncrlons[:] = self.rloc[1]
            ncrlats.units = 'degrees north'
            ncrlats[:] = self.rloc[0]
        nc.close()

    def get_model(self,model,init_date,fc_date,
        timewin=None,distlim=None,corenum=None,expname=None,simmode=None):
        """ 
        Get model data.
        Hindcast AnaC:
        /lustre/storeB/users/anac/HINDCAST2017
        /lustre/storeB/users/anac/HINDCAST2017/BETAMAX1.20
        """
        from utils import haversine
        from model_specs import model_dict
        print ("Get model data according to date ....")
        if timewin is None:
            timewin = int(30)
        if distlim is None:
            distlim = int(6)
        if model == 'ARCMFC':
            filestr = (model_dict[model]['path']
                  + fc_date.strftime('%Y%m%d')
                  + init_date.strftime(model_dict[model]['file_template']))
        elif model == 'ARCMFCnew':
            filestr = (model_dict[model]['path']
                  + expname
                  + init_date.strftime(model_dict[model]['file_template']))
        elif model == 'mwam4':
            filestr = (init_date.strftime(model_dict[model]['path_template'])
                  + init_date.strftime(model_dict[model]['file_template']))
        print (filestr)
        f = netCDF4.Dataset(filestr,'r')
        model_lons = f.variables[model_dict[model]['lons']][:]
        model_lats = f.variables[model_dict[model]['lats']][:]
        model_time = f.variables[model_dict[model]['time']][:]
        # Hs [time,lat,lon]
        model_Hs = f.variables[model_dict[model]['Hs']][:].squeeze()
        f.close()
        model_basetime = model_dict[model]['basetime']
        model_time_dt=[]
        for element in model_time:
            model_time_dt.append(model_basetime 
                            + timedelta(seconds=element))
        sat_time_dt=self.rtime
        if simmode is None:
            model_time_dt_valid=[model_time_dt[model_time_dt.index(fc_date)]]
            print ("date matches found:")
            print (model_time_dt_valid)
        elif simmode == 'cont':
            model_time_dat_valid = [model_time_dt]
        """
        get stellite time steps close to model time step in a given 
        time frame. 
        """
        # Constrain to region
        if ((model=='mwam8' or model=='mwam4' or model=='ARCMFC')
        and (self.region!='Arctic' and self.region!='ARCMFC')):
            model_rlats = model_lats[
                        (model_lats 
                        >= self.region_dict[self.region]["llcrnrlat"])
                      & (model_lats
                        <= self.region_dict[self.region]["urcrnrlat"])
                      & (model_lons
                        >= self.region_dict[self.region]["llcrnrlon"])
                      & (model_lons
                        <= self.region_dict[self.region]["urcrnrlon"])
                                ]
            model_rlons = model_lons[
                        (model_lats 
                        >= self.region_dict[self.region]["llcrnrlat"])
                      & (model_lats
                        <= self.region_dict[self.region]["urcrnrlat"])
                      & (model_lons
                        >= self.region_dict[self.region]["llcrnrlon"])
                      & (model_lons
                        <= self.region_dict[self.region]["urcrnrlon"])
                                ]
            model_rHs=[]
            for i in range(len(model_Hs)):
                tmpA = model_Hs[i,:]
                tmpB = tmpA[
                        (model_lats 
                        >= self.region_dict[self.region]["llcrnrlat"])
                      & (model_lats
                        <= self.region_dict[self.region]["urcrnrlat"])
                      & (model_lons
                        >= self.region_dict[self.region]["llcrnrlon"])
                      & (model_lons
                        <= self.region_dict[self.region]["urcrnrlon"])
                                ]
                model_rHs.append(tmpB)
                del tmpA, tmpB
        elif ((model=='mwam8' or model=='mwam4' 
        or model=='ARCMFC' or model=='ARCMFCnew')
        and (self.region=='Arctic' or self.region=='ARCMFC')):
            model_rlats = model_lats[
                        (model_lats 
                        >= self.region_dict[self.region]["boundinglat"])
                                ]
            model_rlons = model_lons[
                        (model_lats 
                        >= self.region_dict[self.region]["boundinglat"])
                                ]
            model_rHs=[]
            for i in range(len(model_Hs)):
                tmpA = model_Hs[i,:]
                tmpB = tmpA[
                        (model_lats 
                        >= self.region_dict[self.region]["boundinglat"])
                                ]
                model_rHs.append(tmpB)
                del tmpA, tmpB
        # Compare wave heights of satellite with model with 
        # constraint on distance and time frame
        nearest_all_date_matches=[]
        nearest_all_dist_matches=[]
        nearest_all_model_Hs_matches=[]
        nearest_all_sat_Hs_matches=[]
        nearest_all_sat_lons_matches=[]
        nearest_all_sat_lats_matches=[]
        nearest_all_model_lons_matches=[]
        nearest_all_model_lats_matches=[]
        # create local variables before loop
        sat_rlats=self.rloc[0]
        sat_rlons=self.rloc[1]
        sat_rHs=self.rHs
        # moving window compensating for increasing latitudes
        try:
            moving_win = round(
                    (distlim / 
                     haversine(0,
                        np.max(sat_rlats),
                        1,
                        np.max(sat_rlats))
                    ),
                    2)
        except (ValueError):
            moving_win = .6
        print ("Searching for matches with moving window of degree:",\
                moving_win)
        for j in range(len(sat_time_dt)):
            progress(j,str(int(len(sat_time_dt))),'')
            try:
                resultlst=tmploop_get_model(
                    j,sat_time_dt,model_time_dt_valid,timewin,distlim,model,
                    sat_rlats,sat_rlons,sat_rHs,
                    model_rlats,model_rlons,model_rHs,moving_win)
                nearest_all_date_matches.append(resultlst[0])
                nearest_all_dist_matches.append(resultlst[1])
                nearest_all_model_Hs_matches.append(resultlst[2])
                nearest_all_sat_Hs_matches.append(resultlst[3])
                nearest_all_sat_lons_matches.append(resultlst[4])
                nearest_all_sat_lats_matches.append(resultlst[5])
                nearest_all_model_lons_matches.append(resultlst[6])
                nearest_all_model_lats_matches.append(resultlst[7])
            except:
                pass
        results_dict = {
            'valid_date':np.array(model_time_dt_valid),
            'date_matches':np.array(nearest_all_date_matches),
            'dist_matches':np.array(nearest_all_dist_matches),
            'model_Hs_matches':np.array(nearest_all_model_Hs_matches),
            'sat_Hs_matches':np.array(nearest_all_sat_Hs_matches),
            'sat_lons_matches':np.array(nearest_all_sat_lons_matches),
            'sat_lats_matches':np.array(nearest_all_sat_lats_matches),
            'model_lons_matches':np.array(nearest_all_model_lons_matches),
            'model_lats_matches':np.array(nearest_all_model_lats_matches)
            }
        return results_dict

def validate(results_dict,boot=None):
    import numpy as np
    from utils import rmsd, corr, mad, bias, scatter_index
    """
    produced metrics:
    mean of product --> mop
    mean of reference --> mor
    mean square difference --> msd
    number of data values --> nov
    scatter index --> SI
    """
    flatten = lambda l: [item for sublist in l for item in sublist]
    date_matches = np.array(results_dict['date_matches'])
    model_Hs_matches = np.array(results_dict['model_Hs_matches'])
    sat_Hs_matches = np.array(results_dict['sat_Hs_matches'])
    if (boot is None or boot ==  False):
        mop = np.nanmean(model_Hs_matches)
        mor = np.nanmean(sat_Hs_matches)
        msd, rmsd = rmsd(model_Hs_matches,sat_Hs_matches)
        nov = len(sat_Hs_matches)
        mad = mad(model_Hs_matches,sat_Hs_matches)
        corr = corr(model_Hs_matches,sat_Hs_matches)
        bias = bias(model_Hs_matches,sat_Hs_matches)
        SI = scatter_index(model_Hs_matches,sat_Hs_matches)
        arcmfc_validation_dict = {
            'mop':mop,
            'mor':mor,
            'msd':msd,
            'nov':nov,
            'rmsd':rmsd,
            'corr':corr,
            'mad':mad,
            'bias':bias,
            'SI':SI}
    elif boot is True:
        from utils import bootstr, marginalize
        reps=1000
        newmodel,newsat,newidx = marginalize(model_Hs_matches,sat_Hs_matches)
        sat_boot,boot_idx=bootstr(newsat,reps)
        print (len(sat_boot[np.isnan(sat_boot)]))
        RMSD=np.zeros(reps)*np.nan
        MSD=np.zeros(reps)*np.nan
        BIAS=np.zeros(reps)*np.nan
        CORR=np.zeros(reps)*np.nan
        for i in range(reps):
            results_dict = {'date_matches':date_matches[newidx[boot_idx[:,i]]],
                        'model_Hs_matches':newmodel[boot_idx[:,i]],
                        'sat_Hs_matches':newsat[boot_idx[:,i]]}
            try:
                RMSD[i]=validate(results_dict)['rmsd']
                MSD[i]=validate(results_dict)['mad']
                BIAS[i]=validate(results_dict)['bias']
                CORR[i]=validate(results_dict)['corr']
            except IndexError as e:
                print (e)
        arcmfc_validation_dict = {'rmsd':RMSD,'mad':MSD,'bias':BIAS,'corr':CORR}
    return arcmfc_validation_dict

def get_pointsat(sa_obj,station=None,lat=None,lon=None,distlim=None):
    from utils import haversine
    from stationlist import locations
    if ((lat is None or lon is None) and (station is None)):
        print ("location is missing")
    if (lat is None or lon is None):
        lat=locations[station][0]
        lon=locations[station][1]
    proxim = 1
    lats = sa_obj.rloc[0]
    lons = sa_obj.rloc[1]
    Hs = sa_obj.rHs
    time = sa_obj.rtime
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

# --- help ------------------------------------------------------------#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
Module encompassing classes and methods to read and 
process wave field related data from satellites.\n
Usage example in python: 
# load modules
from satmod import sentinel_altimeter as sa
from datetime import datetime, timedelta
from satmod import validate\n
# assume 12h leadtime
init_date = datetime(2018,5,1,0,0,0) - timedelta(hours=12)
fc_date = datetime(2018,5,1,0,0,0) 
# get satellite waves for fc_date
timewin = 30 # units = minutes
sa_obj = sa(fc_date,timewin=timewin,
            region="ARCMFC"[,download=True])\n
# instead you can also use a time period
sdate = init_date
edate = fc_date
# if region is within ARCMFC the argument mode="ARCMFC"
# is recommended for quicker execution:
sa_obj = sa(sdate,edate=edate,timewin=timewin,
            region="ARCMFC",mode="ARCMFC")\n
# possible to save data for region in netcdf
sa_obj.dumptonc("outpath/")\n
# possible to have a quick look at the swath
sa_obj.quip("ARCMFC",show=True[,save=True])\n
# model/sentinel collocation:
# get according model data with time and space constraints
results = sa_obj.get_model('ARCMFC',init_date,fc_date,\
                            timewin=30,distlim=6)\n
# validate the model
valid_dict = validate(results)
# or with bootstrap
valid_dict = validate(results, boot=True)
        """,
        formatter_class = RawTextHelpFormatter
        )
    args = parser.parse_args()
