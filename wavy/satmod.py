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
import sys

# progress bar
from utils import progress, sort_files

# all class
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os
import yaml
import time

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

# create_file
import calendar

# libraries for parallel computing
from joblib import Parallel, delayed

# get_remote
from dateutil.relativedelta import relativedelta
from copy import deepcopy

# credentials
from credentials import get_credentials
# --- global functions ------------------------------------------------#

# read yaml config files:
with open("../config/region_specs.yaml", 'r') as stream:
    region_dict=yaml.safe_load(stream)
with open("../config/model_specs.yaml", 'r') as stream:
    model_dict=yaml.safe_load(stream)
with open("../config/satellite_specs.yaml",'r') as stream:
    satellite_dict=yaml.safe_load(stream)
with open("../config/variable_shortcuts.yaml",'r') as stream:
    shortcuts_dict=yaml.safe_load(stream)

def tmploop_get_remote_files(i,matching,user,pw,
                            server,path,
                            path_local):
    #for element in matching:
    print("File: ",matching[i])
    dlstr=('ftp://' + user + ':' + pw + '@' 
                + server + path + matching[i])
    for attempt in range(10):
        print ("Attempt to download data: ")
        try:
            print ("Downloading file")
            urlretrieve(dlstr, path_local + '/' + matching[i])
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

def get_remote_files(path_remote,path_local,sdate,edate,timewin,
                    corenum,download,instr,provider):
    '''
    Download swath files and store them at defined location
    time stamps in file name stand for: from, to, creation
    '''
    if download is False:
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
            server = satellite_dict[instr][provider]['remote']['server']
            path = (path_remote
                 + '/'
                 + str(tmpdate.year)
                 + '/'
                 + tmpdate.strftime('%m')
                 + '/')
            print ('# ----- ')
            print ('Chosen source: ')
            print (instr + ' from ' + provider + ': ' + server)
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
        tmpdate_new=tmpdate-timedelta(minutes=timewin)
        tmpend=edate+timedelta(minutes=timewin)
        while (tmpdate_new.strftime("%Y%m%dT%H")
            <= tmpend.strftime("%Y%m%dT%H")):
            matchingtmp = [s for s in content
                            if ('_'
                            + str(tmpdate_new.year)
                            + str(tmpdate_new)[5:7]
                            + str(tmpdate_new)[8:10]
                            + 'T'
                            + str(tmpdate_new)[11:13])
                            in s
                            ]
            tmplst = tmplst + matchingtmp
            tmpdate_new = tmpdate_new + timedelta(minutes=timewin)
        matching = tmplst
        print ("Download initiated")
        # Download and gunzip matching files
        print ('Downloading ' + str(len(matching)) + ' files: .... \n')
        print ("Used number of cores " + str(corenum) + "!")
        Parallel(n_jobs=corenum)(
                        delayed(tmploop_get_remote_files)(
                        i,matching,user,pw,server,
                        path,path_local
                        ) for i in range(len(matching))
                        )
        # update time
        tmpdate = datetime((tmpdate + relativedelta(months=+1)).year,
                            (tmpdate + relativedelta(months=+1)).month,1)
    print ('Files downloaded to: \n' + path_local)

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


class satellite_class():
    '''
    class to handle netcdf files containing satellite data i.e. 
    Hs[time], lat[time], lon[time] 
    This class offers the following added functionality:
     - get swaths of desired days and read
     - get the closest time stamp(s)
     - get the location (lon, lat) for this time stamp
     - get Hs value for this time
    '''

    def __init__(
        self,sdate,sat='s3a',instr='altimeter',provider='cmems',
        edate=None,timewin=None,download=False,region=None,
        corenum=1,varlst=None
        ):
        print ('# ----- ')
        print (" ### Initializing satellite_class object ###")
        print ('# ----- ')
        # check settings
        if edate is None:
            print ("Requested time: ", str(sdate))
            edate = sdate
            if timewin is None:
                timewin = int(30)
        else:
            print ("Requested time frame: " + 
                str(sdate) + " - " + str(edate))
        # make satpaths
        path_local = satellite_dict[instr][provider]['local']['path']\
                        + '/' + sat + '/'
        path_remote = satellite_dict[instr][provider]['remote']['path']\
                        + '/' + sat + '/'
        self.path_local = path_local
        self.path_remote = path_remote
        # retrieve files
        get_remote_files(path_remote,path_local,
                        sdate,edate,timewin,corenum,download,
                        instr,provider)
        pathlst, filelst = self.get_local_filelst(
                                sdate,edate,timewin,region
                                )
        if pathlst > 0:
            vardict = self.read_local_files(pathlst,varlst)
            print('Total: ', len(vardict['time']), ' footprints found')
            # find values for give time constraint
            cidx,dtimelst = self.matchtime(
                                sdate,edate,vardict['time'],
                                vardict['time_unit'],timewin
                                )
            cvardict = {}
            for element in vardict:
                if element != 'time_unit':
                    cvardict[element] = list(np.array(vardict[element])[cidx])
            del vardict
            print('In chosen time period: ', len(cvardict['time']), 
                ' footprints found')
            # find values for given region
            ridx = self.matchregion(cvardict['latitude'],
                                cvardict['longitude'],
                                region=region,grid_date=sdate)
            rvardict = {}
            for element in cvardict:
                if element != 'time_unit':
                    rvardict[element] = list(np.array(cvardict[element])[ridx])
            del cvardict
            print('For chosen region and time: ', len(rvardict['time']), 
                ' footprints found')
            # define class variables
            self.edate = edate
            self.sdate = sdate
            self.vars = rvardict
            self.timewin = timewin
            self.region = region
            self.sat = sat
            print ("Satellite object initialized including " 
                + str(len(self.vars['time'])) + " footprints.")
        else:
            print('No satellite_class object initialized')

    def get_local_filelst(self,sdate,edate,timewin,region):
        print ("Time window: ", timewin)
        tmpdate=deepcopy(sdate-timedelta(minutes=timewin))
        pathlst = []
        filelst = []
        try:
            while (tmpdate <= edate + timedelta(minutes=timewin)):
                tmpdatestr=(self.path_local
                        + str(tmpdate.year)
                        + '/' + tmpdate.strftime('%m')
                        + '/')
                tmplst = np.sort(os.listdir(tmpdatestr))
                filelst.append(tmplst)
                pathlst.append([(tmpdatestr + e) for e in tmplst])
                if (edate is not None and edate!=sdate):
                    tmpdate = tmpdate + timedelta(hours=1)
                else:
                    tmpdate = tmpdate + relativedelta(months=+1)
            filelst=np.sort(flatten(filelst))
            pathlst=np.sort(flatten(pathlst))
            idx_start,tmp = check_date(pathlst,sdate-timedelta(minutes=timewin))
            tmp,idx_end = check_date(pathlst,edate+timedelta(minutes=timewin))
            del tmp
            pathlst = np.unique(pathlst[idx_start:idx_end+1])
            filelst = np.unique(filelst[idx_start:idx_end+1])
            print (str(int(len(pathlst))) + " valid files found")
        except FileNotFoundError as e:
            print(e)
        return pathlst,filelst

    def read_local_files(self,pathlst,varlst=['Hs']):
        '''
        read and concatenate all data to one timeseries for each variable
        '''
        # --- open file and read variables --- #
        varlst = varlst + ['lons','lats','time']
        varlst_cf = []
        for var in varlst:
            varlst_cf.append(shortcuts_dict[var])
        count = 0
        print ("Processing " + str(int(len(pathlst))) + " files")
        print (pathlst[0])
        print (pathlst[-1])
        vardict = {}
        for element in pathlst:
            progress(count,str(int(len(pathlst))),'')
            stdname_lst = []
            try:
                # file includes a 1-D dataset with dimension time
                f = netCDF4.Dataset(element,'r')
                ncvars = [v for v in f.variables]
                for ncvar in ncvars:
                    stdname = f.variables[ncvar].getncattr('standard_name')
                    if stdname in varlst_cf:
                        if stdname in vardict:
                            if stdname in stdname_lst:
                                print("Caution: variable " +
                                        "standard_name is not unique !!!")
                                print("only 1. appearance is used.")
                                print("variable " + ncvar + " is neglected")
                            else:
                                tmp = list(f.variables[ncvar][:])
                                vardict[stdname] += tmp
                        else:
                            tmp = list(f.variables[ncvar][:])
                            vardict[stdname] = tmp
                        stdname_lst.append(stdname)
            except (IOError):
                print ("No such file or directory")
            count = count + 1
        print ('\n')
        # remove redundant entries
        time_unique,indices=np.unique(vardict['time'],return_index=True)
        for key in vardict:
            vardict[key]=list(np.array(vardict[key])[indices])
        if (f.variables[shortcuts_dict['lons']].getncattr('valid_min') == 0):
            # transform to -180 to 180 degrees
            tmp = np.array(vardict[shortcuts_dict['lons']])
            vardict[shortcuts_dict['lons']] = list(((tmp - 180) % 360) - 180)
        # add reference time from netcdf
        vardict['time_unit'] = f.variables['time'].units
        # close nc-file
        f.close()
        return vardict

    def bintime(self,binframe=None):
        '''
        fct to return frequency of occurrence per chose time interval
        e.g. days, weeks, months, years ...
        currently only daily bins
        '''
        return freqlst,datelst

    def matchtime(self,sdate,edate,time,reftime,timewin=None):
        '''
        fct to obtain the index of the time step closest to the 
        requested time step from buoy and forecast including the 
        respective time stamp(s).
        '''
        if timewin is None:
            timewin = 0
        # create list of datetime instances
        dtimelst=[]
        cidx=[]
        idx=0
        print ('Time window is: ', timewin)
        if (edate is None or sdate==edate):
            for element in time:
                tmp = netCDF4.num2date(element,reftime)
                dtimelst.append(tmp)
                # choose closest match within window of win[minutes]
                if (tmp >= sdate-timedelta(minutes=timewin) 
                and tmp < sdate+timedelta(minutes=timewin)):
                    cidx.append(idx)
                del tmp
                idx=idx+1
        if (edate is not None and edate!=sdate):
            for element in time:
                tmp = netCDF4.num2date(element,reftime)
                dtimelst.append(tmp)
                if (tmp >= sdate-timedelta(minutes=timewin)
                and tmp < edate+timedelta(minutes=timewin)):
                    cidx.append(idx)
                del tmp
                idx=idx+1
        return cidx, list(np.array(dtimelst)[cidx])

    def matchregion(self,LATS,LONS,region,grid_date):
        # region in region_dict[poly]:
        # find values for given region
        if (region not in region_dict['poly'] and region not in model_dict):
            if region is None:
                region = 'global'
            if ~isinstance(region,str)==True:
                print ("Manually specified region "
                    + [llcrnrlat,urcrnrlat,llcrnrlon,urcrnrlon] + ": \n"
                    + " --> Bounds: " + str(region))
            else:
                if region not in region_dict['rect']:
                    sys.exit("Region is not defined")
                else:
                    print ("Specified region: " + region + "\n"
                      + " --> Bounds: " + str(region_dict['rect'][region]))
            ridx = self.matchregion_rect(LATS,LONS,region=region)
        else:
            ridx = self.matchregion_poly(LATS,LONS,region=region,
                                    grid_date=grid_date)
        return ridx

    def matchregion_rect(self,LATS,LONS,region):
        if (region is None or region == "global"):
            region = "global"
            ridx = range(len(LATS))
        else:
            if isinstance(region,str)==True:
                if "boundinglat" in region_dict['rect'][region]:
                    boundinglat = region_dict['rect'][region]["boundinglat"]
                else:
                    llcrnrlat = region_dict['rect'][region]["llcrnrlat"]
                    urcrnrlat = region_dict['rect'][region]["urcrnrlat"]
                    llcrnrlon = region_dict['rect'][region]["llcrnrlon"]
                    urcrnrlon = region_dict['rect'][region]["urcrnrlon"]
            else:
                llcrnrlat = region[0]
                urcrnrlat = region[1]
                llcrnrlon = region[2]
                urcrnrlon = region[3]
            # check if coords in region
            ridx=[]
            for i in range(len(LATS)):
                if (
                "boundinglat" in region_dict['rect'][region]
                and LATS[i] >= boundinglat
                    ):
                    ridx.append(i)
                elif(
                not "boundinglat" in region_dict['rect'][region]
                and LATS[i] >= llcrnrlat
                and LATS[i] <= urcrnrlat
                and LONS[i] >= llcrnrlon
                and LONS[i] <= urcrnrlon
                    ):
                    ridx.append(i)
        if not ridx:
            print ("No values for chosen region and time frame!!!")
        else:
            print ("Values found for chosen region and time frame.")
        return ridx

    def matchregion_poly(self,LATS,LONS,region,grid_date):
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
            from modelmod import get_model
            import pyproj
            try:
                grid_date = grid_date
                print('Use date for retrieving grid: ', grid_date)
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(model=region, fc_date=grid_date)
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
                model_Hs,model_lats,model_lons,model_time,model_time_dt = \
                    get_model(model=region, fc_date=grid_date)
            if (len(model_lats.shape)==1):
                model_lons, model_lats = np.meshgrid(model_lons, model_lats)
            if (region=='global'):
                rlatlst, rlonlst = LATS, LONS
            else:
                proj4 = model_dict[region]['proj4']
                proj_model = pyproj.Proj(proj4)
                Mx, My = proj_model(model_lons,model_lats,inverse=False)
                Vx, Vy = proj_model(LONS,LATS,inverse=False)
                xmax, xmin = np.max(Mx), np.min(Mx)
                ymax, ymin = np.max(My), np.min(My)
                ridx = list(np.where((Vx>xmin) & (Vx<xmax) & 
                                (Vy>ymin) & (Vy<ymax))[0]
                            )
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
            latlst = LATS
            lonlst = LONS
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


def get_pointsat(sa_obj,station=None,lat=None,lon=None,distlim=None):
    from utils import haversine
    # read yaml config files:
    with open("../config/stationlist.yaml", 'r') as stream:
        locations=yaml.safe_load(stream)
    if ((lat is None or lon is None) and (station is None)):
        print ("location is missing")
    if (lat is None or lon is None):
        lat=locations[station][0]
        lon=locations[station][1]
    print('Get footprints near lat:', lat, ' lon:', lon)
    lats = sa_obj.loc[0]
    lons = sa_obj.loc[1]
    Hs = sa_obj.Hs
    time = sa_obj.time
    sample = []
    distsp = []
    lonsp = []
    latsp = []
    timep = []
    idx = []
    for i in range(len(lats)):
        dist=haversine(lons[i],lats[i],lon,lat)
        if (dist<=distlim):
            sample.append(Hs[i])
            distsp.append(dist)
            lonsp.append(lons[i])
            latsp.append(lats[i])
            timep.append(time[i])
            idx.append(i)
    return sample, distsp, lonsp, latsp, timep, idx
