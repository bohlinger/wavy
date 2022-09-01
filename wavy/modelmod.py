#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
The main task of this module is to acquire, read, and prepare
geophysical variables from model output files for further use.
'''
# --- import libraries ------------------------------------------------#
# standard library imports
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import time
from functools import lru_cache
from tqdm import tqdm

# own imports
from wavy.utils import hour_rounder, make_fc_dates
from wavy.utils import finditem, parse_date, NoStdStreams
from wavy.utils import convert_meteorologic_oceanographic
from wavy.utils import find_direction_convention
from wavy.ncmod import check_if_ncfile_accessible
from wavy.ncmod import ncdumpMeta, get_filevarname
from wavy.wconfig import load_or_default

# ---------------------------------------------------------------------#

def get_model_filedate(model, fc_date, leadtime):
    '''
    get init_date for latest model output file and checks if available

    param:
        model - modelname type(str)
        fc_date - datetime object
        leadtime - integer in hours

    return:
        suitable datetime to create model filename
    '''
    if ('init_times' in model_dict[model].keys()
            and model_dict[model]['init_times'] is not None):
        init_times = \
            np.array(model_dict[model]['init_times']).astype('float')
    else:
        print('init_times for chosen model not specified in config file')
        print('Assuming continuous simulation with hourly values')
        init_times = np.array(range(25)).astype('float')
    date = fc_date - timedelta(hours=leadtime)
    date_hour = hour_rounder(date).hour
    if date_hour in init_times:
        print('Leadtime', leadtime , \
              'available for date' ,fc_date)
        init_diffs = date_hour - init_times
        init_diffs[init_diffs < 0] = np.nan
        h_idx = np.where(init_diffs==\
                        np.min(init_diffs[~np.isnan(init_diffs)]))
        h = int(init_times[h_idx[0][0]])
        return datetime(date.year, date.month, date.day, h)
    else:
        #print('Leadtime', leadtime , \
        #      'not available for date' ,fc_date, '!')
        return False


def make_model_filename(model, fc_date, leadtime):
    """
    creates/returns filename based on fc_date,leadtime

        param:
        model - modelname type(str)
        fc_date - datetime object
        leadtime - integer in hours

    return:
        filename (consists of path + filename)

    comment:
            - special characters are escaped by adding "\\"
            - the escapes need to be removed for certain libraries
              like xarray and netCDF4
    """
    if model in model_dict:
        if 'xtra_h' in model_dict[model]:
            filedate = get_model_filedate(model, fc_date, leadtime)
            pathdate = filedate + timedelta(hours=leadtime) \
                                * model_dict[model]['lt_switch_p']
            tmpstr = model_dict[model]['file_template']
            for i in range(model_dict[model]['nr_filedates']):
                filedatestr = model_dict[model]['filedate_formats'][i]
                replacestr = (filedate \
                            + timedelta(hours = leadtime \
                                    - (leadtime % \
                                        model_dict[model]['init_step']))
                                    * model_dict[model]['lt_switch_f'][i]
                            + timedelta(hours = \
                                    model_dict[model]['xtra_h'][i])).\
                            strftime(filedatestr)
                tmpstr = tmpstr.replace('filedate', replacestr, 1)
            filename = (
                        pathdate.strftime(\
                            model_dict[model]['path_template'])
                        + tmpstr)
        else:
            filedate = get_model_filedate(model, fc_date, leadtime)
            filename =\
                (filedate.strftime(model_dict[model]['path_template'])\
                + filedate.strftime(model_dict[model]['file_template']))
    else:
        raise ValueError("Chosen model is not specified in model_specs.yaml")
    # replace/escape special characters
    filename = filename.replace(" ", "\\ ")\
                       .replace("?", "\\?")\
                       .replace("&", "\\&")\
                       .replace("(", "\\(")\
                       .replace(")", "\\)")\
                       .replace("*", "\\*")\
                       .replace("<", "\\<")\
                       .replace(">", "\\>")
    return filename


def make_model_filename_wrapper(model, fc_date, leadtime, max_lt = None):
    """
    Wrapper function of make_model_filename. Organizes various cases.

    param:
        model - modelname type(str)
        fc_date - datetime object
        leadtime - integer in hours
        max_lt - maximum lead time allowed

    return:
        filename
    """
    if leadtime is None:
        leadtime = 'best'
    if (isinstance(fc_date, datetime) and leadtime != 'best'):
        filename = make_model_filename(model, fc_date, leadtime)
    elif (isinstance(fc_date, datetime) and leadtime == 'best'):
        switch = False
        leadtime = generate_bestguess_leadtime(model, fc_date)
        while switch is False:
            filename = make_model_filename(model, fc_date, leadtime)
            # check if file is accessible
            switch = check_if_ncfile_accessible(filename)
            if (switch is False):
                print("Desired file:", filename , " not accessible")
                print("Continue to look for date with extended leadtime")
                leadtime = leadtime + model_dict[model]['init_step']
            if max_lt is not None and leadtime > max_lt:
                print("Leadtime:",leadtime,"is greater as maximum allowed leadtime:",max_lt)
                return None
    elif (isinstance(fc_date, list) and isinstance(leadtime, int)):
        filename = [make_model_filename(model,date,leadtime) \
                    for date in fc_date]
    elif (isinstance(fc_date, list) and leadtime == 'best'):
        leadtime = generate_bestguess_leadtime(model, fc_date)
        filename = [make_model_filename(model,fc_date[i],leadtime[i]) \
                    for i in range(len(fc_date))]
    return filename

def make_list_of_model_filenames(model,fc_dates,lt):
    """
    return: flst - list of model files to be opened
            dlst - list of dates to be chosen within each file
    """
    #fn = make_model_filename_wrapper('mwam4',datetime(2021,1,1,1),1)
    flst = []
    for d in fc_dates:
        fn = make_model_filename_wrapper(model,d,lt)
        flst.append(fn)
    return flst

@lru_cache(maxsize=32)
def read_model_nc_output_lru(filestr,lonsname,latsname,timename):
    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    filestr=filestr.replace('\\','')
    f = netCDF4.Dataset(filestr, 'r')
    # get coordinates and time
    model_lons = f.variables[lonsname][:]
    model_lats = f.variables[latsname][:]
    model_time = f.variables[timename]
    model_time_dt = list(
        netCDF4.num2date(model_time[:], units=model_time.units))
    f.close()
    return model_lons,model_lats,model_time_dt

def read_model_nc_output(filestr,lonsname,latsname,timename):
    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    filestr=filestr.replace('\\','')
    f = netCDF4.Dataset(filestr, 'r')
    # get coordinates and time
    model_lons = f.variables[lonsname][:]
    model_lats = f.variables[latsname][:]
    model_time = f.variables[timename]
    model_time_dt = list(
        netCDF4.num2date(model_time[:], units=model_time.units))
    f.close()
    return model_lons,model_lats,model_time_dt

def get_model_fc_mode(filestr, model, fc_date, varalias=None, **kwargs):
    """
    fct to retrieve model data for correct time
    """
    vardict = {}
    print("Get model data according to selected date(s) ....")
    print(filestr)
    meta = ncdumpMeta(filestr)
    stdvarname = variable_info[varalias]['standard_name']
    lonsname = get_filevarname('lons',variable_info,
                               model_dict[model],meta)
    latsname = get_filevarname('lats',variable_info,
                               model_dict[model],meta)
    timename = get_filevarname('time',variable_info,
                               model_dict[model],meta)
    # get other variables e.g. Hs [time,lat,lon]
    filevarname = get_filevarname(varalias,variable_info,
                                  model_dict[model],meta)
    try:
        model_lons,model_lats,model_time_dt = \
        read_model_nc_output_lru(filestr,lonsname,latsname,timename)
    except Exception as e:
        print(e)
        print('continue with uncached retrieval')
        model_lons,model_lats,model_time_dt = \
        read_model_nc_output(filestr,lonsname,latsname,timename)

    vardict[variable_info['lons']['standard_name']] = model_lons
    vardict[variable_info['lats']['standard_name']] = model_lats

    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    filestr = filestr.replace('\\','')
    f = netCDF4.Dataset(filestr, 'r')
    model_time = f.variables[timename]
    l = kwargs.get('vertical_level',0)
    if (type(filevarname) is dict):
        print(
            'Target variable can be computed from vector \n'
            'components with the following aliases: ', filevarname)
        model_time_dt_valid = model_time_dt[model_time_dt.index(fc_date)]
        model_time_valid = float(model_time[model_time_dt.index(fc_date)])
        model_time_unit = model_time.units
        vardict[variable_info['time']['standard_name']] = model_time_valid
        vardict['datetime'] = model_time_dt_valid
        vardict['time_unit'] = model_time_unit
        for key in filevarname.keys():
            filevarname_dummy = get_filevarname(\
                                    filevarname[key][0],
                                    variable_info,
                                    model_dict[model],meta)
            if filevarname_dummy is not None:
                print(filevarname[key][0], 'exists')
                break
        print('Use aliases:', filevarname[key])
        model_var_dummy = f.variables[filevarname_dummy]
        if len(model_var_dummy.dimensions) == 4:
            model_var_dummy = model_var_dummy[:,l,:,:].squeeze()
        if len(model_var_dummy.shape) == 3: # for multiple time steps
            model_var_valid_tmp = \
                model_var_dummy[model_time_dt.index(fc_date),:,:].squeeze()**2
            for i in range(1, len(filevarname[key])):
                filevarname_dummy = get_filevarname(\
                                        filevarname[key][i],
                                        variable_info,
                                        model_dict[model],meta)
                if len(f.variables[filevarname_dummy].dimensions) == 4:
                    model_var_valid_tmp += \
                        f.variables[filevarname_dummy][
                            model_time_dt.index(fc_date),l,:,:
                            ].squeeze()**2
                elif len(f.variables[filevarname_dummy].dimensions) == 3:
                    model_var_valid_tmp += \
                        f.variables[filevarname_dummy][
                            model_time_dt.index(fc_date),:,:
                            ].squeeze()**2
            model_var_valid = np.sqrt(model_var_valid_tmp)
        elif len(model_var_dummy.dimensions) == 2:
            model_var_valid_tmp = model_var_dummy[:, :]**2
            for i in range(1, len(filevarname[key])):
                filevarname_dummy = get_filevarname(\
                                            filevarname[key][i],
                                            variable_info,
                                            model_dict[model],
                                            meta)
                model_var_valid_tmp += \
                    f.variables[filevarname_dummy][:,:]**2
            model_var_valid = np.sqrt(model_var_valid_tmp)
        else:
            print('Dimension mismatch!')
        vardict[stdvarname] = model_var_valid
    else:
        model_time_dt_valid = model_time_dt[model_time_dt.index(fc_date)]
        model_time_valid = float(model_time[model_time_dt.index(fc_date)])
        model_time_unit = model_time.units
        vardict[variable_info['time']['standard_name']] = model_time_valid
        vardict['datetime'] = model_time_dt_valid
        vardict['time_unit'] = model_time_unit
        model_var_link = f.variables[filevarname]
        if len(model_var_link.dimensions) == 4:
            model_var_link = model_var_link[:,l,:,:].squeeze()
        if len(model_var_link.shape) == 3:  # for multiple time steps
            model_var_valid = \
                model_var_link[model_time_dt.index(fc_date),:,:].squeeze()
        elif len(model_var_link.dimensions) == 2:
            model_var_valid = model_var_link[:, :].squeeze()
        else:
            print('Dimension mismatch!')
        vardict[variable_info[varalias]['standard_name']] = \
                                                    model_var_valid
    # transform masked array to numpy array with NaNs
    f.close()
    vardict['longitude'] = vardict['longitude'].filled(np.nan)
    vardict['latitude'] = vardict['latitude'].filled(np.nan)
    vardict[variable_info[varalias]['standard_name']] = \
        vardict[variable_info[varalias]['standard_name']].filled(np.nan)
    vardict['meta'] = meta
    # make lats,lons 2D if only 1D (regular grid)
    if len(vardict['longitude'].shape)==1:
        LATS,LONS = np.meshgrid(vardict['latitude'],
                                vardict['longitude'])
        vardict['longitude']=np.transpose(LONS)
        vardict['latitude']=np.transpose(LATS)
    return vardict, filevarname


def generate_bestguess_leadtime(model, fc_date, lidx=None):
    """
    fct to return leadtimes for bestguess
    """
    if isinstance(fc_date, list):
        leadtime = \
            [generate_bestguess_leadtime(model,date) for date in fc_date]
    else:
        init_times = \
            np.array(model_dict[model]['init_times']).astype('float')
        diffs = fc_date.hour - np.array(init_times)
        gtz = diffs[diffs >= 0]
        if len(gtz) == 0:
            leadtime = int(np.abs( np.min(np.abs(diffs)) \
                                 - model_dict[model]['init_step'] ))
        elif (len(gtz)>0 and lidx is not None):
            leadtime = int(np.sort(diffs[diffs >= 0])[lidx])
        else:
            leadtime = int(np.min(diffs[diffs >= 0]))
    return leadtime

def get_model(model=None,
              sdate=None,
              edate=None,
              date_incr=None,
              fc_date=None,
              leadtime=None,
              varalias=None,
              st_obj=None,
              transform_lons=None):
    """
    toplevel function to get model data
    """
    if st_obj is not None:
        sdate = st_obj.sdate
        edate = st_obj.edate
        varalias = st_obj.varalias
    if (sdate is not None
        and edate is not None
        and date_incr is not None):
        fc_date = make_fc_dates(sdate, edate, date_incr)
    filestr = make_model_filename_wrapper(model=model,
                                          fc_date=fc_date,
                                          leadtime=leadtime)
    if (isinstance(filestr, list) and st_obj is None):
        vardict, \
        filevarname = get_model_fc_mode(filestr=filestr[0],model=model,
                                    fc_date=fc_date[0],varalias=varalias)
        vardict[variable_info[varalias]['standard_name']]=\
                    [vardict[variable_info[varalias]['standard_name']]]
        vardict['time'] = [vardict['time']]
        vardict['datetime'] = [vardict['datetime']]
        vardict['leadtime'] = leadtime
        for i in tqdm(range(1, len(filestr))):
            tmpdict, \
            filevarname = get_model_fc_mode(\
                            filestr=filestr[i],model=model,
                            fc_date=fc_date[i],varalias=varalias)
            vardict[variable_info[varalias]['standard_name']].append(
                tmpdict[variable_info[varalias]['standard_name']])
            vardict['time'].append(tmpdict['time'])
            vardict['datetime'].append(tmpdict['datetime'])
        vardict[variable_info[varalias]['standard_name']]=\
            np.array(vardict[variable_info[varalias]['standard_name']])
    else:
        vardict, \
        filevarname = get_model_fc_mode(filestr=filestr,model=model,
                                    fc_date=fc_date,varalias=varalias)
        vardict['time'] = [vardict['time']]
        vardict['datetime'] = [vardict['datetime']]
        vardict['leadtime'] = leadtime
    # transform to datetime
    with NoStdStreams():
        tmpd = [parse_date(str(d)) for d in vardict['datetime']]
    vardict['datetime'] = tmpd
    del tmpd
    if transform_lons==180:
        # transform longitudes from -180 to 180
        vardict['longitude'] = ((vardict['longitude'] - 180) % 360) - 180
    elif transform_lons==360:
        print('not yet implemented !!')
    # adjust conventions
    # check if variable is one with conventions
    if 'convention' in variable_info[varalias].keys():
        convention_set = False
        print('Chosen variable is defined with conventions')
        print('... checking if correct convention is used ...')
        # 1. check if clear from standard_name
        file_stdvarname = find_direction_convention(filevarname,vardict['meta'])
        if "to_direction" in file_stdvarname:
            print('Convert from oceanographic to meteorologic convention')
            vardict[variable_info[varalias]['standard_name']] = \
                    convert_meteorologic_oceanographic(\
                        vardict[variable_info[varalias]['standard_name']])
            convention_set = True
        elif "from_direction" in file_stdvarname:
            print('standard_name indicates meteorologic convention')
            convention_set = True
            pass
        # 2. overwrite convention from config file
        if ('convention' in model_dict[model].keys() and
        model_dict[model]['convention'] == 'oceanographic' and
        convention_set is False):
            print('Convention is set in config file')
            print('This will overwrite conventions from standard_name in file!')
            print('\n')
            print('Convert from oceanographic to meteorologic convention')
            vardict[variable_info[varalias]['standard_name']] = \
                    convert_meteorologic_oceanographic(\
                        vardict[variable_info[varalias]['standard_name']])
            convention_set = True
    return vardict, fc_date, leadtime, filestr, filevarname


# ---------------------------------------------------------------------#

# read yaml config files:
model_dict = load_or_default('model_specs.yaml')
variable_info = load_or_default('variable_info.yaml')


class model_class():
    '''
    class to read and process model data
    model: e.g. Hs[time,lat,lon], lat[rlat,rlon], lon[rlat,rlon]
    This class should communicate with the satellite, model, and
    station classes.
    '''
    def __init__(self,
                 model='mwam4',
                 sdate=None,
                 edate=None,
                 date_incr=1,
                 fc_date=None,
                 leadtime=None,
                 varalias='Hs',
                 st_obj=None,
                 transform_lons=None):
        print('# ----- ')
        print(" ### Initializing model_class object ###")
        print(' ')
        # parse and translate date input
        sdate = parse_date(sdate)
        edate = parse_date(edate)
        fc_date = parse_date(fc_date)
        if st_obj is not None:
            sdate = st_obj.sdate
            edate = st_obj.edate
            varalias = st_obj.varalias
        if fc_date is not None:
            print("Requested time: ", str(fc_date))
        elif (edate is None and fc_date is None and sdate is not None):
            fc_date = sdate
            print("Requested time: ", str(fc_date))
        elif (sdate is None and edate is None and fc_date is None):
            now = datetime.now()
            fc_date = datetime(now.year, now.month, now.day, now.hour)
            print("Requested time: ", str(fc_date))
        elif (sdate is not None and edate is not None
              and date_incr is not None):
            fc_date = sdate
            print("Requested time frame: " + str(sdate) + " - " + str(edate))

        print(" ")
        print(" ## Read files ...")
        t0 = time.time()
        vardict, \
        fc_date, leadtime, \
        filestr, \
        filevarname = get_model(model=model,sdate=sdate,edate=edate,
                            date_incr=date_incr,fc_date=fc_date,
                            leadtime=leadtime,varalias=varalias,
                            st_obj=st_obj,transform_lons=transform_lons)
        stdname = variable_info[varalias]['standard_name']
        units = variable_info[varalias]['units']
        varname = filevarname
        # define class variables
        self.fc_date = fc_date
        self.sdate = sdate
        self.edate = edate
        self.leadtime = leadtime
        self.model = model
        self.varalias = varalias
        self.varname = varname
        self.stdvarname = stdname
        self.units = units
        self.vars = vardict
        self.filestr = filestr
        # create label for plotting
        self.label = self.model
        t1 = time.time()
        print(" ")
        print( '## Summary:')
        print("Time used for retrieving model data:", round(t1 - t0, 2),
              "seconds")
        print(" ")
        print(" ### model_class object initialized ###")
        print('# ----- ')

    def get_item_parent(self,item,attr):
        ncdict = self.vars['meta']
        lst = [i for i in ncdict.keys() \
                if (attr in ncdict[i].keys() \
                and item in ncdict[i][attr]) \
                ]
        if len(lst) >= 1:
            return lst
        else: return None

    def get_item_child(self,item):
        ncdict = self.vars['meta']
        parent = finditem(ncdict,item)
        return parent

    def quicklook(self,a=True,projection=None,date=None,**kwargs):
        m = kwargs.get('m',a)
        if m is True:
            import cartopy.crs as ccrs
            import cmocean
            import matplotlib.pyplot as plt
            import matplotlib.cm as mplcm
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            dt = parse_date(date)
            lons = self.vars['longitude']
            lats = self.vars['latitude']
            var = self.vars[self.stdvarname]
            fc_date = self.fc_date
            if date is not None:
                # find date in self.vars and give a snapshot
                idx = self.vars['datetime'].index(dt)
                var = self.vars[self.stdvarname][idx,:]
                fc_date = self.fc_date[idx]
            elif (date is None and len(var.shape)==3):
                # find date in self.vars and give a snapshot
                var = self.vars[self.stdvarname][0,:]
                fc_date = self.fc_date[0]
            if projection is None:
                projection = ccrs.PlateCarree()
            # parse kwargs
            cflevels = kwargs.get('cflevels',10)
            clevels = kwargs.get('clevels',10)
            vartype = variable_info[self.varalias].get('type','default')
            if kwargs.get('cmap') is None:
                if vartype == 'cyclic':
                    cmap = mplcm.twilight
                    cflevels = kwargs.get('cflevels',range(0,365,5))
                    clevels = kwargs.get('clevels',range(0,365,5))
                else:
                    cmap = cmocean.cm.amp
            else:
                cmap = kwargs.get('cmap')
            lonmax,lonmin = np.max(lons),np.min(lons)
            latmax,latmin = np.max(lats),np.min(lats)
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection=projection)
            ax.set_extent(  [lonmin, lonmax,latmin, latmax],
                            crs = projection )
            cf = ax.contourf(lons,lats, var,
                            cmap=cmap,levels=cflevels,
                            transform=ccrs.PlateCarree())
            c = ax.contour(lons,lats, var,
                           cmap=cmap,levels=clevels,
                           transform=ccrs.PlateCarree())
            axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="100%",  # height : 50%
                       loc='lower left',
                       bbox_to_anchor=(1.01, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
            fig.colorbar(cf, cax=axins, label=self.varalias
                                        + ' [' + self.units + ']')
            ax.coastlines()
            ax.gridlines()
            plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
            ax.set_title(self.model + ' for ' + str(fc_date))
            plt.show()

    def write_to_pickle(self,pathtofile=None):
        import pickle
        # writing
        pickle.dump( self, open( pathtofile, "wb" ) )
        print('model_class object written to:',pathtofile)
        # for reading
        # mco = pickle.load( open( pathtofile, "rb" ) )
