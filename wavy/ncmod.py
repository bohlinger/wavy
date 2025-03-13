#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses classes and methods to read and write to netcdf
files from model, station, or satellite output. I try to mostly follow
the PEP convention for python code style. Constructive comments on style
and effecient programming are most welcome!
'''
# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this class. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# standard library imports
import netCDF4
import xarray as xr
import numpy as np
from datetime import datetime
import os
import sys
from copy import deepcopy
from functools import lru_cache
from tqdm import tqdm
import zipfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
#from hidefix import xarray

# own imports
from wavy.wconfig import load_or_default
from wavy.utils import find_included_times, finditem
from wavy.utils import flatten

# read yaml config files:
insitu_dict = load_or_default('insitu_cfg.yaml')
satellite_dict = load_or_default('satellite_cfg.yaml')
variable_info = load_or_default('variable_def.yaml')

# --- global functions ------------------------------------------------#
"""
definition of some global functions
"""
# currently None
# ---------------------------------------------------------------------#

def check_if_ncfile_accessible(fstr):
    # remove escape character because xarray handles white spaces
    # but cannot handle escape characters (apparently)
    fstr_repl = fstr.replace('\\', '')
    try:
        xs = xr.open_dataset(fstr_repl)
        return True
    except (OSError, FileNotFoundError) as e:
        print("Desired file not accessible")
        print(e)
        return False

def read_netcdfs(paths, dim='time', decode_times=None, use_cftime=None):
    @lru_cache(maxsize=128)
    def process_one_path(path):
        # use a context manager, to ensure the file gets closed after use
        with xr.open_dataset(path,
                decode_times=decode_times,
                use_cftime=use_cftime) as ds:
            # use it after closing each original file
            ds.load()
            return ds
    datasets = [process_one_path(p) for p in tqdm(paths)]
    print("Concatenate ...")
    combined = xr.concat(datasets, dim,
                         coords='minimal',
                         data_vars='minimal',
                         compat='override',
                         combine_attrs='override',
                         join='override')
    print("... done concatenating")
    return combined

def build_usr_pw_path(path, remoteHostName, usr, pw):
    searchstr = '/thredds/dodsC'
    idx = path.find(searchstr)
    newpath = f'https://{usr}:{pw}@my.cmems-du.eu/' + path[idx::]
    return newpath

def read_netcdfs_with_credentials_aggregated(
path, remoteHostName, usr, pw, decode_times=None, use_cftime=None):
    tmp_path = build_usr_pw_path(path, remoteHostName, usr, pw)
    ds = xr.open_dataset(tmp_path,
                         decode_times=decode_times,
                         use_cftime=use_cftime)
    return ds

def read_netcdfs_naive(paths, varnames, varalias):
    # read all paths
    print('read all data files')
    darrays = [process_one_path_netCDF4(p, varnames, varalias)\
               for p in tqdm(paths)]
    # consolidate all data
    print('consolidate data')
    darray = consolidate_darrays(darrays)
    # create xr dataset
    print('build xr dataset')
    ds = build_xr_ds(darray, varnames, varalias)
    return ds

def process_one_path_netCDF4(path: str, varnames: dict, varalias: str):
    @lru_cache(maxsize=128)
    def read_nc(path):
        nc = netCDF4.Dataset(path)
        time_var = nc.variables[varnames['time']]
        dtime = netCDF4.num2date(time_var[:], time_var.units)
        lons = nc.variables[varnames['lons']][:]
        lats = nc.variables[varnames['lats']][:]
        var = nc.variables[varnames[varalias]][:]
        nc.close()
        return var, lons, lats, dtime
    var, lons, lats, dtime = read_nc(path)
    return np.array([var, lons, lats, dtime])

def consolidate_darrays(darrays: list) -> np.ndarray:
    A = darrays[0]
    for i in range(1, len(darrays)):
        A = np.append(A, darrays[i], axis=1)
    return A

def build_xr_ds(darray: np.ndarray, varnames: dict, varalias: str):
    import xarray as xr
    ds = xr.Dataset({
            varnames[varalias]: xr.DataArray(
                    data=darray[0],
                    dims=[varnames['time']],
                    coords={varnames['time']: darray[3]}
                    ),
            varnames['lons']: xr.DataArray(
                    data=darray[1],
                    dims=[varnames['time']],
                    coords={varnames['time']: darray[3]}
                    ),
            varnames['lats']: xr.DataArray(
                    data=darray[2],
                    dims=[varnames['time']],
                    coords={varnames['time']: darray[3]}
                    ),
            varnames['time']: xr.DataArray(
                    data=darray[3],
                    dims=[varnames['time']],
                    coords={varnames['time']: darray[3]}
                    )
                },
            attrs={'title': 'wavy dataset'}
        )
    return ds

def build_xr_ds_from_dict(dict_var, var_name_ref):
    ds = xr.Dataset({
        var_name: xr.DataArray(
            data=dict_var[var_name],
            dims=[var_name_ref],
            coords={var_name_ref: dict_var[var_name_ref]}
            ) for var_name in dict_var.keys()},
                attrs={'title': 'wavy dataset'})
    return ds

#def read_netcdfs_hidefix(paths):
#    ds = xr.open_mfdataset(paths, engine='hidefix')
#    return ds

#def read_netcdf_hidefix(path):
#    return xr.open_mfdataset(path, engine='hidefix')

#def tpe_hidefix(paths):
#    with ThreadPoolExecutor() as tpe:
#        results = list(tpe.map(read_netcdf_hidefix, paths))
#    return xr.concat(results, dim='time',
#                     coords='minimal',
#                     data_vars='minimal',
#                     compat='override',
#                     combine_attrs='override')

def read_netcdfs_KF(paths, dim='time'):
    # https://github.com/knutfrode/concepts/blob/main/Open_MFDataset_overlap.ipynb
    datasets = [xr.open_dataset(p, chunks='auto') for p in paths]
    print('Concatenating...')
    ds = xr.concat(datasets, dim=dim,
                   compat='override',
                   combine_attrs='override',
                   join='override',
                   coords='all')
    print("... done concatenating")
    return ds

def read_mf_netcdfs(paths):
    ds = xr.open_mfdataset(paths, parallel=True)
    return ds

@lru_cache(maxsize=32)
def process_one_path_lru(path, t, varname):
    with xr.open_dataset(path) as ds:
        da = ds.sel(time=t)[varname]
        da.load()
        return da

def process_one_path(path,t,varname):
    with xr.open_dataset(path) as ds:
        da = ds.sel(time=t)[varname]
        da.load()
        return da

def read_netcdfs_sel_lru(paths,dlst,varname,dim='time'):
    dataarr = [process_one_path_lru(paths[i],dlst[i],varname)\
                for i in tqdm(range(len(paths)))]
    print("Concatenate ...")
    combined = xr.concat(dataarr,dim,
                         coords='minimal',
                         compat='override',
                         combine_attrs='override')
    combined = combined.to_dataset()
    print("... done concatenating")
    return combined

def read_netcdfs_sel(paths,dlst,varname,dim='time'):
    dataarr = [process_one_path(paths[i],dlst[i],varname)\
                for i in tqdm(range(len(paths)))]
    print("Concatenate ...")
    combined = xr.concat(dataarr, dim,
                         coords='minimal',
                         compat='override',
                         combine_attrs='override')
    combined = combined.to_dataset()
    print("... done concatenating")
    return combined

def get_swim_var_coords(varalias):
    varname = satellite_dict['cfo_swim_L2P']['vardef'][varalias]
    varidx = satellite_dict['cfo_swim_L2P']['vardef'][varalias + '_idx']
    lonname = satellite_dict['cfo_swim_L2P']['vardef']['lons']
    lonidx = satellite_dict['cfo_swim_L2P']['vardef']['lons_idx']
    latname = satellite_dict['cfo_swim_L2P']['vardef']['lats']
    latidx = satellite_dict['cfo_swim_L2P']['vardef']['lats_idx']
    timename = satellite_dict['cfo_swim_L2P']['vardef']['time']
    timeidx = satellite_dict['cfo_swim_L2P']['vardef']['time_idx']

    varnamedict = {
            'varname':varname,
            'varidx':varidx,
            'lonname':lonname,
            'lonidx':lonidx,
            'latname':latname,
            'latidx':latidx,
            'timename':timename,
            'timeidx':timeidx
            }

    return varnamedict

def read_swim_nc(path,varnamedict):
    ds = xr.open_dataset(path)
    var = eval("ds[varnamedict['varname']].values"+varnamedict['varidx'])
    time = eval("ds[varnamedict['timename']].values"+varnamedict['timeidx'])
    lons = eval("ds[varnamedict['lonname']].values"+varnamedict['lonidx'])
    lats = eval("ds[varnamedict['latname']].values"+varnamedict['latidx'])
    return var, time, lons, lats

def read_swim_netcdfs(pathlst,varalias):
    varnamedict = get_swim_var_coords(varalias)

    varlst = []
    timelst = []
    lonslst = []
    latslst = []

    for f in pathlst:
        var,time,lon,lat = \
            read_swim_nc(f,varnamedict) 
        varlst.append(var)
        timelst.append(time)
        lonslst.append(lon)
        latslst.append(lat)

    var = flatten(varlst)
    time = flatten(timelst)
    lons = flatten(lonslst)
    lats = flatten(latslst)

    vardict = {
            varalias: var,
            'time': time,
            'longitude': lons,
            'latitude': lats
            }
    return vardict

def get_arcmfc_ts(pathtofile):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        sys.exit('File does not exist')
    else:
        nc = netCDF4.Dataset(
            pathtofile, mode='r',
            )
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:], time_var.units)
        sHs = nc.variables['obs_values'][:]
        mHs = nc.variables['model_values'][:]
        nc.close()
    return dtime, sHs, mHs

def get_arcmfc_stats(pathtofile):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        print('File does not exist')
        return
    else:
        nc = netCDF4.Dataset(
            pathtofile, mode='r',
            )
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:], time_var.units)
        mop = nc.variables['mop'][:]
        mor = nc.variables['mor'][:]
        rmsd = nc.variables['rmsd'][:]
        msd = nc.variables['msd'][:]
        corr = nc.variables['corr'][:]
        mad = nc.variables['mad'][:]
        bias = nc.variables['bias'][:]
        SI = nc.variables['SI'][:]
        nov = nc.variables['nov'][:]
        nc.close()
        valid_dict = {
            'mop': mop,
            'mor': mor,
            'msd': msd,
            'nov': nov,
            'rmsd': rmsd,
            'msd': msd,
            'corr': corr,
            'mad': mad,
            'bias': bias,
            'SI': SI}
        return valid_dict, dtime

def get_filevarname(varalias, variable_info, srcdict, ncdict):
    stdname = variable_info[varalias]['standard_name']
    print(' Get filevarname for \n' + 'stdvarname:', stdname,
          '\n' + 'varalias:', varalias)
    filevarname = get_varname_for_cf_stdname_in_ncfile(ncdict, stdname)
    if (filevarname is None and 'alias' in variable_info[varalias]):
        filevarname = get_varname_for_cf_stdname_in_ncfile(
            ncdict, variable_info[varalias]['alias'])
    if (filevarname is not None and len(filevarname) > 1):
        print(' !!! standard_name: ', stdname, ' is not unique !!!',
              '\nThe following variables have the same standard_name:\n',
              filevarname)
        print(' Searching *_cfg.yaml config file for definition')
        filevarname = None
    if filevarname is not None:
        return filevarname[0]
    tmpdict = finditem(srcdict, 'vardef')
    if len(tmpdict)==0:
        tmpdict = [{'None':None}]
    vardefdict = tmpdict[0]
    if (filevarname is None and varalias in vardefdict.keys()):
        filevarname = vardefdict[varalias]
        print(' Variable defined in *_cfg.yaml is:')
        print(varalias, '=', filevarname)
        return filevarname
    elif (filevarname is None
          and varalias not in vardefdict.keys()
          and 'aliases_of_vector_components' in variable_info[varalias]):
        print(
          'Checking variable_info if variable can be ' +
          'computed from vector components')
        filevarname = variable_info[varalias]['aliases_of_vector_components']
        return filevarname
    else:
        print(' !!! variable not defined nor ' +
              'available in nc-file !!!')

def get_nc_ts(pathtofile, varlst):
    import os.path
    indicator = os.path.isfile(pathtofile)
    if indicator is False:
        dtime = False
        sys.exit('File does not exist')
    else:
        vardict = {}
        for name in varlst:
            nc = netCDF4.Dataset(
                pathtofile, mode='r',
                )
            var = nc.variables[name][:]
            vardict[name] = var
        time_var = nc.variables['time']
        dtime = netCDF4.num2date(time_var[:], time_var.units)
        vardict['dtime'] = dtime
        vardict['time'] = time_var[:]
        vardict['time_unit'] = time_var.units
        nc.close()
    return vardict

def get_varlst_from_nc_1D(pathtofile, varlst, sdate, edate):
    # retrieve
    vardict = {}
    nc = netCDF4.Dataset(pathtofile)
    time_var = nc.variables['time']
    vardict['time'] = time_var[:]
    vardict['time_unit'] = time_var.units
    dtvar = netCDF4.num2date(time_var[:], time_var.units)
    dtime = [datetime(dt.year, dt.month, dt.day,
             dt.hour, dt.minute, dt.second)
             for dt in dtvar]
    idx = find_included_times(dtime, sdate=sdate, edate=edate)
    for name in varlst:
        var = nc.variables[name][idx]
        vardict[name] = var
    vardict['dtime'] = np.array(dtime)[idx]
    nc.close()
    return vardict

def dumptonc_stats(pathtofile, title, time_dt, time_unit, valid_dict):
    """
    1. check if nc file already exists
    2. - if so use append mode
       - if not create file
    """
    # create time vector in seconds since first date
    time = netCDF4.date2num(time_dt, time_unit)
    mop = np.array(valid_dict['mop'])
    mor = np.array(valid_dict['mor'])
    rmsd = np.array(valid_dict['rmsd'])
    msd = np.array(valid_dict['msd'])
    corr = np.array(valid_dict['corr'])
    mad = np.array(valid_dict['mad'])
    bias = np.array(valid_dict['bias'])
    SI = np.array(valid_dict['SI'][1])
    nov = np.array(valid_dict['nov'])
    print('Dump data to file: ' + pathtofile)
    if os.path.isfile(pathtofile):
        nc = netCDF4.Dataset(
                        pathtofile, mode='a',
                        clobber=False
                        )
        # variables
        startidx = len(nc['time'])
        endidx = len(nc['time'])+1
        nc.variables['time'][startidx:endidx] = time
        nc.variables['mop'][startidx:endidx] = mop
        nc.variables['mor'][startidx:endidx] = mor
        nc.variables['rmsd'][startidx:endidx] = rmsd
        nc.variables['msd'][startidx:endidx] = msd
        nc.variables['corr'][startidx:endidx] = corr
        nc.variables['mad'][startidx:endidx] = mad
        nc.variables['bias'][startidx:endidx] = bias
        nc.variables['SI'][startidx:endidx] = SI
        nc.variables['nov'][startidx:endidx] = nov
    else:
        outpath = pathtofile[0:-len(pathtofile.split('/')[-1])]
        os.makedirs(outpath, exist_ok=True)
        nc = netCDF4.Dataset(
                        pathtofile, mode='w',
                        # format='NETCDF4'
                        )
        nc.title = title
        dimsize = None
        # dimensions
        dimtime = nc.createDimension(
                                'time',
                                size=dimsize
                                )
        # variables
        nctime = nc.createVariable(
                               'time',
                               np.float64,
                               dimensions=('time')
                               )
        ncmop = nc.createVariable(
                               'mop',
                               np.float64,
                               dimensions=('time')
                               )
        ncmor = nc.createVariable(
                               'mor',
                               np.float64,
                               dimensions=('time')
                               )
        ncrmsd = nc.createVariable(
                               'rmsd',
                               np.float64,
                               dimensions=('time')
                               )
        ncmsd = nc.createVariable(
                               'msd',
                               np.float64,
                               dimensions=('time')
                               )
        nccorr = nc.createVariable(
                               'corr',
                               np.float64,
                               dimensions=('time')
                               )
        ncmad = nc.createVariable(
                               'mad',
                               np.float64,
                               dimensions=('time')
                               )
        ncbias = nc.createVariable(
                               'bias',
                               np.float64,
                               dimensions=('time')
                               )
        ncSI = nc.createVariable(
                               'SI',
                               np.float64,
                               dimensions=('time')
                               )
        ncnov = nc.createVariable(
                               'nov',
                               np.float64,
                               dimensions=('time')
                               )
        # generate time for netcdf file
        # time
        nctime.standard_name = 'time matches'
        nctime.long_name = \
            'associated time steps between model and observation'
        nctime.units = time_unit
        nctime[:] = time
        # mop
        ncmop.standard_name = 'mop'
        ncmop.long_name = 'mean of product (wave model)'
        ncmop.units = 'm'
        ncmop[:] = mop
        # mor
        ncmor.standard_name = 'mor'
        ncmor.long_name = 'mean of reference (observations)'
        ncmor.units = 'm'
        ncmor[:] = mor
        # rmsd
        ncrmsd.standard_name = 'rmsd'
        ncrmsd.long_name = 'root mean square deviation'
        ncrmsd.units = 'm'
        ncrmsd[:] = rmsd
        # msd
        ncmsd.standard_name = 'msd'
        ncmsd.long_name = 'mean square deviation'
        ncmsd.units = 'm^2'
        ncmsd[:] = msd
        # corr
        nccorr.standard_name = 'corr'
        nccorr.long_name = 'correlation coefficient'
        nccorr.units = 'none'
        nccorr[:] = corr
        # mad
        ncmad.standard_name = 'mad'
        ncmad.long_name = 'mean absolute deviation'
        ncmad.units = 'm'
        ncmad[:] = mad
        # bias
        ncbias.standard_name = 'bias'
        ncbias.long_name = 'Bias (mean error)'
        ncbias.units = 'm'
        ncbias[:] = bias
        # SI
        ncSI.standard_name = 'SI'
        ncSI.long_name = 'scatter index'
        ncSI.units = 'none'
        ncSI[:] = SI
        # nov
        ncnov.standard_name = 'nov'
        ncnov.long_name = 'number of values'
        ncnov.units = 'none'
        ncnov[:] = nov
    nc.close()

@lru_cache(maxsize=32)
def ncdump(nc_fid, verb=True):
    '''
    Function from:
    http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
    #
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.
    #
    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed
    #
    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("WARNING: %s does not contain variable attributes" % key)
    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

@lru_cache(maxsize=32)
def ncdumpMeta(pathtofile):
    '''
    Returns dict of netcdf-file content
    Input: str pointing to netcdf-file
    Output: dict of attributes
    '''
    # remove escape character because netCDF4 handles white spaces
    # but cannot handle escape characters (apparently)
    pathtofile = pathtofile.replace('\\', '')
    # init netCDF4 instance
    nc = netCDF4.Dataset(pathtofile, mode='r')
    # init empty dict
    ncdict = {}
    # retrieve variable attributes
    for var in nc.variables:
        ncattrs = nc.variables[var].ncattrs()
        ncdict[var] = {}
        for ncattr in ncattrs:
            ncdict[var][ncattr] = nc.variables[var].getncattr(ncattr)
    # retrieve global attributes
    ncdict['global'] = {}
    nc_attrs = nc.ncattrs()
    for nc_attr in nc_attrs:
        ncdict['global'][nc_attr] = nc.getncattr(nc_attr)
    return ncdict


def find_attr_in_nc(attrstr, pathtofile=None, ncdict=None, subattrstr=None):
    """
    fct to find a specific attribute with its value in a netcdf-file
    when only a fraction of attribute name is given, can also search
    in a deeper layer of attribute hierarchy.
    input:  path to the nc-file or ncdict
            string of desired attribute
            string of desired sub-attribute (optional)
    output: dict
    """
    if pathtofile is not None and ncdict is None:
        # get content of netcdf-file
        ncdict = ncdumpMeta(pathtofile)
    # browse for str using list comprehension
    res1 = [i for i in ncdict.keys() if attrstr in i]
    if subattrstr is None:
        return ncdict[res1[0]]
    else:
        res2 = [i for i in ncdict[res1[0]] if subattrstr in i]
        return ncdict[res1[0]][res2[0]]


def get_varname_for_cf_stdname_in_ncfile(ncdict, stdname):
    lst = [i for i in ncdict.keys()
           if ('standard_name' in ncdict[i].keys()
           and stdname in ncdict[i]['standard_name'])
           ]
    if len(lst) >= 1:
        return lst
    else: return None
