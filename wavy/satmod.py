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
# standard library imports
import sys
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os
from copy import deepcopy
import yaml
import time
if sys.version_info <= (3, 0):
    from urllib import urlretrieve, urlcleanup # python2
else:
    from urllib.request import urlretrieve, urlcleanup # python3
import gzip
import ftplib
from ftplib import FTP
import netCDF4 as netCDF4
import calendar
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed
import xarray as xa
import pyresample

# own imports
from .ncmod import ncdumpMeta, get_varname_for_cf_stdname_in_ncfile
from .ncmod import find_attr_in_nc
from utils import progress, sort_files, collocate_times
from credentials import get_credentials
from modelmod import get_filevarname
from modelmod import model_class as mc
from modelmod import make_model_filename_wrapper
import pyproj

# ---------------------------------------------------------------------#

# read yaml config files:
moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/region_specs.yaml'))
with open(moddir,'r') as stream:
    region_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/model_specs.yaml'))
with open(moddir,'r') as stream:
    model_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/satellite_specs.yaml'))
with open(moddir,'r') as stream:
    satellite_dict=yaml.safe_load(stream)

moddir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'config/variable_info.yaml'))
with open(moddir,'r') as stream:
    variable_info=yaml.safe_load(stream)

# --- global functions ------------------------------------------------#

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

def get_remote_files(path_remote,path_local,sdate,edate,twin,
                    nproc,instr,provider):
    '''
    Download swath files and store them at defined location
    time stamps in file name stand for: from, to, creation
    '''
    # credentials
    user, pw = get_credentials()
    tmpdate = deepcopy(sdate)
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
        tmpdate_new=tmpdate-timedelta(minutes=twin)
        tmpend=edate+timedelta(minutes=twin)
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
            tmpdate_new = tmpdate_new + timedelta(minutes=twin)
        matching = tmplst
        # Download and gunzip matching files
        print ('Downloading ' + str(len(matching))
                + ' files: .... \n')
        print ("Used number of processes " + str(nproc) + "!")
        Parallel(n_jobs=nproc)(
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
    # check if str in lst according to desired date (sdate,edate)
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
    Class to handle netcdf files containing satellite data i.e.
    Hs[time], lat[time], lon[time]
    This class offers the following added functionality:
     - get swaths of desired days and read
     - get the closest time stamp(s)
     - get the location (lon, lat) for this time stamp
     - get Hs value for this time
    '''

    def __init__(
        self,sdate,sat='s3a',instr='altimeter',provider='cmems',
        edate=None,twin=None,download=False,download_path=None,
        remote_ftp_path=None,region=None,nproc=1,varalias='Hs'
        ):
        print ('# ----- ')
        print (" ### Initializing satellite_class object ###")
        # check settings
        if edate is None:
            print ("Requested time: ", str(sdate))
            edate = sdate
        else:
            print ("Requested time frame: " +
                str(sdate) + " - " + str(edate))
        if twin is None:
            twin = int(30)
        print('Please wait ...')
        print('Chosen time window is:', twin, 'min')
        # make satpaths
        if download_path is None:
            path_local = satellite_dict[instr][provider]\
                        ['local']['path'] + '/' + sat + '/'
        else:
            path_local = download_path
        if remote_ftp_path is None:
            path_remote = satellite_dict[instr][provider]\
                        ['remote']['path'] + sat + '/'
        else:
            path_remote = remote_ftp_path
        self.path_local = path_local
        self.path_remote = path_remote
        # retrieve files
        if download is False:
            print ("No download initialized, checking local files")
        else:
            print ("Downloading necessary files ...")
            get_remote_files(path_remote,path_local,
                        sdate,edate,twin,nproc,
                        instr,provider)

        t0=time.time()
        pathlst, filelst = self.get_local_filelst(
                            sdate,edate,twin,region)
        if len(pathlst) > 0:
            vardict = self.read_local_files(pathlst,provider,varalias)
            print('Total: ', len(vardict['time']), ' footprints found')
            # find values for give time constraint
            dtime = list( netCDF4.num2date(vardict['time'],
                                           vardict['time_unit']) )
            cidx = collocate_times(dtime,sdate=sdate,edate=edate,twin=twin)
            cvardict = {}
            for element in vardict:
                if element != 'time_unit':
                    cvardict[element] = list(np.array(
                                            vardict[element]
                                            )[cidx])
                else:
                    cvardict[element] = vardict[element]
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
                    rvardict[element] = list(np.array(
                                            cvardict[element]
                                            )[ridx])
                else:
                    rvardict[element] = cvardict[element]
            del cvardict
            if len(rvardict['time'])>0:
                rvardict['datetime'] = netCDF4.num2date(
                                            rvardict['time'],
                                            rvardict['time_unit'])
                print('For chosen region and time: ',
                        len(rvardict['time']),'footprints found')
                # convert to datetime object
                timedt = rvardict['datetime']
                rvardict['datetime'] = [datetime(t.year,t.month,t.day,
                          t.hour,t.minute,t.second) for t in timedt]
            else:
                print('For chosen region and time: 0 footprints found!')
            # find variable name as defined in file
            stdname = variable_info[varalias]['standard_name']
            ncdict = ncdumpMeta(pathlst[0])
            filevarname = get_varname_for_cf_stdname_in_ncfile(
                                                ncdict,stdname)
            if (len(filevarname) or filename is None) > 1:
                filevarname = satellite_dict[instr][provider]\
                            ['misc']['vardef'][stdname]
            else:
                filevarname = get_varname_for_cf_stdname_in_ncfile(
                                                    ncdict,stdname)[0]
            rvardict['model_meta'] = ncdict
            # define class variables
            self.edate = edate
            self.sdate = sdate
            self.vars = rvardict
            self.varalias = varalias
            self.varname = filevarname
            self.stdvarname = stdname
            self.twin = twin
            self.region = region
            self.sat = sat
            t1=time.time()
            print("Time used for retrieving satellite data:",\
                    round(t1-t0,2),"seconds")
            print ("Satellite object initialized including "
                + str(len(self.vars['time'])) + " footprints.")
            #print (" ### satellite_class object initialized ###")
            print ('# ----- ')
        else:
            print('No satellite_class object initialized')
            print ('# ----- ')

    def get_local_filelst(self,sdate,edate,twin,region):
        print ("Time window: ", twin)
        tmpdate = deepcopy(sdate-timedelta(minutes=twin))
        pathlst = []
        filelst = []
        try:
            while (tmpdate <= edate + timedelta(minutes=twin)):
                if self.path_local == 'tmp_unittest/':
                    tmpdatestr = self.path_local
                else:
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
            idx_start,tmp = check_date(pathlst,sdate-timedelta(minutes=twin))
            tmp,idx_end = check_date(pathlst,edate+timedelta(minutes=twin))
            del tmp
            pathlst = np.unique(pathlst[idx_start:idx_end+1])
            filelst = np.unique(filelst[idx_start:idx_end+1])
            print (str(int(len(pathlst))) + " valid files found")
        except FileNotFoundError as e:
            print(e)
        return pathlst,filelst

    def read_local_files(self,pathlst,provider,varalias):
        '''
        read and concatenate all data to one timeseries for each variable
        '''
        # --- open file and read variables --- #
        varlst = [varalias] + ['lons','lats','time']
        varlst_cf = []
        for var in varlst:
            varlst_cf.append(variable_info[var]['standard_name'])
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
                    if satellite_dict['altimeter'][provider]\
                    ['misc']['vardef'] is not None:
                        if (stdname in satellite_dict['altimeter'][provider]\
                        ['misc']['vardef']
                        and ncvar != satellite_dict['altimeter'][provider]\
                        ['misc']['vardef'][stdname]):
                            ncvar = satellite_dict['altimeter'][provider]\
                                ['misc']['vardef'][stdname]
                    if stdname in varlst_cf:
                        if stdname in vardict:
                            if stdname in stdname_lst:
                                print("Caution: variable " +
                                        "standard_name is not unique !!!")
                                if satellite_dict['altimeter'][provider]\
                                ['misc']['vardef'] is not None:
                                    if stdname in satellite_dict['altimeter']\
                                    [provider]['misc']['vardef']:
                                        print( "As defined in "
                                        + "satellite_specs.yaml, "
                                        + "the following "
                                        + "nc-variable name "
                                        + "is chosen:\n", ncvar, "for "
                                        + "stdname: ", stdname )
                                else:
                                    print("Only 1. appearance is used.")
                                    print("variable " + ncvar + " is neglected")
                            else:
                                tmp = list(f.variables[ncvar][:])
                                vardict[stdname] += tmp
                        else:
                            tmp = list(f.variables[ncvar][:])
                            vardict[stdname] = tmp
                        stdname_lst.append(stdname)
            except (IOError) as e:
                print ("No such file or directory")
                print (e)
            count = count + 1
        print ('\n')
        # remove redundant entries
        time_unique,indices=np.unique(vardict['time'],return_index=True)
        for key in vardict:
            vardict[key]=list(np.array(vardict[key])[indices])
        if (f.variables[variable_info['lons']\
        ['standard_name']].getncattr('valid_min') == 0):
            # transform to -180 to 180 degrees
            tmp = np.array(vardict[variable_info['lons']['standard_name']])
            vardict[variable_info['lons']['standard_name']] \
                = list(((tmp - 180) % 360) - 180)
        # add reference time from netcdf
        vardict['time_unit'] = f.variables['time'].units
        # close nc-file
        f.close()
        return vardict

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
            try:
                print('Use date for retrieving grid: ', grid_date)
                filestr = make_model_filename_wrapper(\
                                        region,grid_date,'best')
                model_meta = ncdumpMeta(filestr)
                flon = get_filevarname(region,\
                                    'lons',variable_info,\
                                    model_dict,model_meta)
                flat = get_filevarname(region,\
                                    'lats',variable_info,\
                                    model_dict,model_meta)
                M = xa.open_dataset(filestr, decode_cf=True)
                model_lons = M[flon].data
                model_lats = M[flat].data
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
                model_meta = ncdumpMeta(filestr)
                flon = get_filevarname(region,\
                                    'lons',variable_info,\
                                    model_dict,model_meta)
                flat = get_filevarname(region,\
                                    'lats',variable_info,\
                                    model_dict,model_meta)
                M = xa.open_dataset(filestr, decode_cf=True)
                model_lons = M[flon].data
                model_lats = M[flat].data
            if (len(model_lons.shape)==1):
                model_lons, model_lats = np.meshgrid(
                                        model_lons,
                                        model_lats
                                        )
            print('Check if footprints fall within the chosen domain')
            if (region=='global'):
                rlatlst, rlonlst = LATS, LONS
            else:
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
