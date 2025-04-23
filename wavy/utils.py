"""
utility fcts for the verification
"""
import numpy as np
import netCDF4
import xarray as xr
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from math import radians, cos, sin, asin, sqrt, floor
import sys
import subprocess
import os
from dateutil.parser import parse
import math
from wavy.wconfig import load_or_default
import pandas as pd

# ---------------------------------------------------------------------#
variable_def = load_or_default('variable_def.yaml')
# ---------------------------------------------------------------------#

def grab_PID():
    """
    Retrieves PID and prints it
    """
    import os
    # retrieve PID
    PID = os.getpid()
    print("\n")
    print("PID - with the license to kill :) ")
    print(str(PID))
    print("\n")
    return

def haversineP(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    km = 6367 * c
    return km

def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (np.sin(dlat/2.0)**2
       + np.cos(lat1)
       * np.cos(lat2)
       * np.sin(dlon/2.0)**2)
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km

def haversineA(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    Note: lon1,lat1,lon2, and lat2 can be lists
    """
    # convert decimal degrees to radians
    rads = np.deg2rad(np.array([lon1, lat1, lon2, lat2]))
    # haversine formula
    if isinstance(lon1, list):
        dlon = rads[2, :] - rads[0, :]
        dlat = rads[3, :] - rads[1, :]
        a = np.sin(dlat/2)**2 \
            + np.cos(rads[1, :]) * np.cos(rads[3, :]) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        km = 6367 * c
        return list(km)
    else:
        dlon = rads[2] - rads[0]
        dlat = rads[3] - rads[1]
        a = np.sin(dlat/2)**2 \
            + np.cos(rads[1]) * np.cos(rads[3]) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        km = 6367 * c
        return [km]

def runmean_old(vec, win, mode=None, weights=None) -> tuple:
    """
    Computes the running mean with various configurations.

    Args:
        vec (numpy.ndarray | list): array of values to me smoothed
        win (int): window length
        mode (str): string: left, centered, right
        weights (numpy.ndarray | list): weights (same size as win)

    Returns:
        tuple (out (numpy.ndarray), std (numpy.ndarray)):
                array of smoothed values and std deviation
    """
    win = int(win)
    if mode is None:
        mode = 'centered'
    out = np.zeros(len(vec))*np.nan
    std = np.zeros(len(vec))*np.nan
    length = len(vec)-win+1
    if mode == 'left':
        count = win-1
        start = win-1
        for i in range(length):
            out[count] = np.mean(vec[count-start:count+1])
            std[count] = np.std(vec[count-start:count+1])
            count = count+1
    elif mode == 'centered':
        start = int(floor(win/2))
        for i in range(start, length):
            if win % 2 == 0:
                sys.exit("window length needs to be odd!")
            else:
                sidx = int(i-start)
                eidx = int(i+start+1)
                if weights is not None:
                    out[i] = np.sum(vec[sidx:eidx]*weights)
                else:
                    out[i] = np.mean(vec[sidx:eidx])
                std[i] = np.std(vec[sidx:eidx])
    elif mode == 'right':
        count = int(0)
        for i in range(length):
            out[count] = np.mean(vec[i:i+win])
            std[count] = np.std(vec[i:i+win])
            count = count+1
    return out, std

def runmean(vec, win, mode=None, weights=None) -> tuple:
    """
    Computes the running mean with various configurations.

    Args:
        vec (numpy.ndarray | list): array of values to me smoothed
        win (int): window length
        mode (str): string: left, centered, right
        weights (numpy.ndarray | list): weights (same size as win)

    Returns:
        tuple (out (numpy.ndarray), std (numpy.ndarray)):
                array of smoothed values and std deviation
    """
    win = int(win)
    if mode is None:
        mode = 'centered'
    out = np.zeros(len(vec))*np.nan
    std = np.zeros(len(vec))*np.nan
    if mode == 'left':
        length = len(vec)
        start = win-1
        for i in range(start, length):
            out[i] = np.mean(vec[i-win+1:i+1])
            std[i] = np.std(vec[i-win+1:i+1])
    elif mode == 'centered':
        length = len(vec)-floor(win/2)
        start = int(floor(win/2))
        for i in range(start, length):
            if win % 2 == 0:
                sys.exit("window length needs to be odd!")
            else:
                sidx = int(i-start)
                eidx = int(i+start+1)
                if weights is not None:
                    out[i] = np.sum(vec[sidx:eidx]*weights)
                else:
                    out[i] = np.mean(vec[sidx:eidx])
                std[i] = np.std(vec[sidx:eidx])
    elif mode == 'right':
        length = len(vec)
        for i in range(length-win+1):
            out[i] = np.mean(vec[i:i+win])
            std[i] = np.std(vec[i:i+win])
    return out, std

def runmean_conv(x: np.ndarray, win: int, mode='flat') -> np.ndarray:
    """
    running mean using convolution

    Args:
        x (numpy.ndarray): array of values to me smoothed
        win (int): window length
        mode (str): which type of smoothing window to pic

    Notes:
        https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

    Returns:
        out (numpy array): array of smoothed values

    Raises:
        ValueError: for wrong dimension of x and wrong windowsize
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < win:
        raise ValueError("Input vector needs to be bigger than window size.")
    if win < 3:
        print("window length too small -> returning original signal")
        return x
    s = np.r_[x[win-1:0:-1], x, x[-2:-win-1:-1]]
    if mode == 'flat':  # moving average
        w = np.ones(win, 'd')
    else:
        # ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
        w = eval('numpy.' + mode + '(win)')
    out = np.convolve(w/w.sum(), s, mode='valid')
    return out

def bootstr(a, reps):
    """
    Conducts a simple naive bootstrap:

    input:    - a is time series of length n
              - reps (number of repetitions)
    output:   - an array of dim n x m where
                m is the number of repetitions
              - indices of draws
    """
    n = len(a)
    b = np.random.choice(a, (n, reps))
    bidx = np.zeros(b.shape) * np.nan
    for i in range(len(a)):
        tmp = np.where(b == a[i])
        bidx[tmp[0], tmp[1]] = i
        del tmp
    return b, bidx.astype('int')

def marginalize(a, b=None):
    """
    Removes entries in both time series that are NaN.

    input: np.arrays with np.nan for invalids
    """
    if b is None:
        return a[~np.isnan(a)]
    else:
        comb = a + b
        idx = np.array(range(len(a)))[~np.isnan(comb)]
        a1 = a[idx]
        b1 = b[idx]
        return a1, b1, idx

def hour_rounder(t, method='nearest'):
    '''
    Rounds to nearest hour adding a timedelta hour if minute >= 30 (default), 
    or to the integer hour before (floor) or the integer hour after (ceil) the 
    given time.
    '''
    if method=='nearest':
        add_hour=t.minute//30
    elif method=='floor':
        add_hour=0
    elif method=='ceil':
        add_hour=1

    t = (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
              + timedelta(hours=add_hour))

    return t 

def hour_rounder_pd(times):
    '''
    Rounds to nearest hour by adding a timedelta hour if minute >= 30
    '''
    df = pd.DataFrame(columns=['time'], data=times)
    rounded = df.time.round('H').values
    return rounded

def sort_files(dirpath, filelst, product, sat):
    """
    mv files to sub-folders of year and month
    """
    if product == 'cmems_L3_NRT':
        sort_cmems_l3_nrt(dirpath, filelst, sat)
    elif product == 'cmems_L3_s6a':
        sort_cmems_l3_s6a(dirpath, filelst, sat)
    elif product == 'cmems_L3_MY':
        sort_cmems_l3_my(dirpath, filelst, sat)
    elif (product == 'cci_L2P' or product == 'cci_L3'):
        sort_cci(dirpath, filelst)
    elif product == 'eumetsat_L2':
        sort_eumetsat_l2(dirpath, filelst)
    elif product == 'cfo_swim_L2P':
        sort_aviso_l2p(dirpath, filelst)

def sort_aviso_l2p(dirpath: str, filelst: list):
    '''
    Sort AVISO files according to year and month.
    '''
    for e in filelst:
        if os.path.isfile(os.path.join(dirpath, e)):
            tmp = e.split('_')
            d1 = parse(tmp[-2])
            #d2 = parse(tmp[-1].split('.')[0])
            year, month = d1.strftime("%Y"), d1.strftime("%m")
            folder = os.path.join(dirpath, year, month)
            os.makedirs(folder, exist_ok=True)
            cmd = 'mv ' + dirpath + '/' + e + ' ' + folder
            os.system(cmd)

def sort_cmems_l3_nrt(dirpath: str, filelst: list, sat: str):
    '''
    Sort L3 files according to year and month.
    '''
    for e in filelst:
        if os.path.isfile(os.path.join(dirpath, e)):
            tmp = 'global_vavh_l3_rt_' + sat + '_'
            year, month = e[len(tmp):len(tmp)+4], e[len(tmp)+4:len(tmp)+6]
            folder = os.path.join(dirpath, year, month)
            os.makedirs(folder, exist_ok=True)
            cmd = 'mv ' + dirpath + '/' + e + ' ' + folder
            os.system(cmd)

def sort_cmems_l3_s6a(dirpath: str, filelst: list, sat: str):
    '''
    Sort L3 s6a files according to year and month.
    '''
    for e in filelst:
        if os.path.isfile(os.path.join(dirpath, e)):
            tmp = 'global_vavh_l3_rt_' + sat + '_lr_'
            year, month = e[len(tmp):len(tmp)+4], e[len(tmp)+4:len(tmp)+6]
            folder = os.path.join(dirpath, year, month)
            os.makedirs(folder, exist_ok=True)
            cmd = 'mv ' + dirpath + '/' + e + ' ' + folder
            os.system(cmd)

def sort_cmems_l3_my(dirpath: str, filelst: list, sat: str):
    '''
    Sort L3 files according to year and month.
    '''
    for e in filelst:
        if os.path.isfile(os.path.join(dirpath, e)):
            tmp = 'global_vavh_l3_rep_' + sat + '_'
            year, month = e[len(tmp):len(tmp)+4], e[len(tmp)+4:len(tmp)+6]
            folder = os.path.join(dirpath, year, month)
            os.makedirs(folder, exist_ok=True)
            cmd = 'mv ' + dirpath + '/' + e + ' ' + folder
            os.system(cmd)

def sort_cci(dirpath: str, filelst: list):
    '''
    Sort L2P and L3 files according to year and month.
    '''
    for e in filelst:
        if os.path.isfile(os.path.join(dirpath, e)):
            tmp = e.split('-')[-2]
            year, month = tmp[0:4], tmp[4:6]
            folder = os.path.join(dirpath, year, month)
            os.makedirs(folder, exist_ok=True)
            cmd = 'mv ' + dirpath + '/' + e + ' ' + folder
            os.system(cmd)

def sort_eumetsat_l2(dirpath: str, filelst: list):
    '''
    Sort L2 files according to year and month.
    '''
    for e in filelst:
        splits = e.split('____')
        if os.path.isfile(os.path.join(dirpath, e)):
            year, month = splits[1][0:4], splits[1][4:6]
            print(year, month)
            folder = os.path.join(dirpath, year, month)
            print(folder)
            os.makedirs(folder, exist_ok=True)
            cmd = 'mv ' + dirpath + '/' + e + ' ' + folder
            os.system(cmd)

def get_size(obj, seen=None):
    """
    Recursively finds size of objects

    From:
    https://goshippo.com/blog/measure-real-size-any-python-object/
    """
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size

def find_included_times_pd(unfiltered_t: list,
                           sdate: datetime, edate: datetime) -> list:
    idx = np.array(range(len(unfiltered_t)))
    df = pd.to_datetime(unfiltered_t)
    mask = ((df >= sdate.isoformat()) & (df <= edate.isoformat()))
    return list(idx[mask])

def find_included_times(
    unfiltered_t: list, target_t=None,
    sdate=None, edate=None, twin=0) -> list:
    """
    Find index/indices of unfiltered time series that fall
    within a tolerance time window around the target time
    or within a time window specified by sdate and edate
    """
    if (sdate is None and edate is None):
        idx = [i for i in range(len(unfiltered_t))
               if (unfiltered_t[i] >= target_t-timedelta(minutes=twin)
               and unfiltered_t[i] <= target_t+timedelta(minutes=twin))]
    else:
        idx = [i for i in range(len(unfiltered_t))
               if (unfiltered_t[i] >= sdate-timedelta(minutes=twin)
               and unfiltered_t[i] <= edate+timedelta(minutes=twin))]
    return idx

def collocate_times(
    unfiltered_t: list, target_t=None,
    sdate=None, edate=None, twin=None) -> list:
    """
    fct for collocating times within a given twin as tolerance
    target_t and unfiltered_t need to be lists of datetime objects
    twin is in minutes.

    returns idx
    """
    if twin is None:
        twin = 0
    if ((twin is None or twin == 0) and (sdate is None and edate is None)):
        idx = [unfiltered_t.index(t) for t in target_t if (t in unfiltered_t)]
    else:
        if (sdate is None and edate is None):
            idx = [find_included_times(unfiltered_t, target_t=t,
                                       sdate=sdate, edate=edate, twin=twin)
                   for t in target_t]
            idx = flatten(idx)
        else:
            idx = find_included_times(unfiltered_t, sdate=sdate,
                                      edate=edate, twin=twin)
    return idx

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

def make_fc_dates(sdate: datetime, edate: datetime,
date_incr_unit: str, date_incr: int) -> list:
    '''
    fct to create forecast date vector
    '''
    sdate = parse_date(str(sdate))
    edate = parse_date(str(edate))
    fc_dates = []
    while sdate <= edate:
        fc_dates.append(sdate)
        tmp_date = parse_date(str(sdate))
        sdate = date_dispatcher(tmp_date, date_incr=date_incr_unit,
                                incr=date_incr)
    return fc_dates

def system_call(command: str):
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read()

def make_subdict(strsublst, class_object=None, class_object_dict=None):
    if class_object_dict is None:
        class_object_dict = vars(class_object)
    subdict = {}
    if strsublst is None:
        pass
    else:
        for strsub in strsublst:
            if strsub in class_object_dict:
                subdict[strsub] = class_object_dict[strsub]
            else:
                print(strsub, 'is not available and not substituted')

    return subdict

def get_pathtofile(pathlst, strsublst, subdict, date):
    '''
    Finds and returns path of file given templates and keywords and date.
    '''
    i = 0
    switch = False
    if isinstance(pathlst, list):
        pass
    else:
        pathlst = [pathlst]
    while switch is False:
        try:
            pathtofile = date.strftime(pathlst[i])
        except IndexError as e:
            print(e)
            print('Index to large for pathlst')
            print('-> returning None')
            return None
        for strsub in strsublst:
            pathtofile = pathtofile.replace(strsub, subdict[strsub])
            #pathtofile = pathtofile.replace(strsub,kwargs[strsub])
        # check if thredds and if accessible using netCDF4a
        if ('thredds' in pathtofile and pathtofile[-3::] == '.nc'):
            # check if available
            try:
                nc = netCDF4.Dataset(pathtofile)
                nc.close()
                switch = True
            except Exception as e:
                print(e)
                print(pathtofile, 'not accessible')
        else:
            if os.path.isfile(pathtofile) is not False:
                switch = True
        if switch is False:
            print(pathtofile, 'does not exist!')
            i += 1
    return pathtofile

def finditem(search_dict, field):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []
    for key, value in search_dict.items():
        if key == field:
            fields_found.append(value)
        elif isinstance(value, dict):
            results = finditem(value, field)
            for result in results:
                fields_found.append(result)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = finditem(item, field)
                    for another_result in more_results:
                        fields_found.append(another_result)
    return fields_found

def make_pathtofile(tmppath, strsublst, subdict, date=None):
    '''
    Creates a path given templates and keywords and date.
    '''
    if date is not None:
        pathtofile = date.strftime(tmppath)
    else: pathtofile = tmppath
    if strsublst is None:
        pass
    else:
        for strsub in strsublst:
            if strsub in subdict:
                pathtofile = pathtofile.replace(strsub, subdict[strsub])
            else:
                print(strsub,
                      'in substitutables not needed for destination path')
    return pathtofile

def find_direction_convention(filevarname, ncdict):
    file_stdvarname = get_item_child(ncdict, filevarname)[0]['standard_name']
    return file_stdvarname

def convert_meteorologic_oceanographic(alpha):
    """
    fct to convert angles from meteorological convention to
    oceanographic and vice versa.
    """
    return (alpha+180)%360

class NoStdStreams(object):
    '''
    Suppress stdout.
    if argument is verbose stdout is shown,
    e.g.: with NoStdStreams(verbose)

    https://codereview.stackexchange.com/questions/25417/
    is-there-a-better-way-to-make-a-function-silent-on-need
    '''
    def __init__(self,stdout = None, stderr = None):
        self.devnull = open(os.devnull,'w')
        self._stdout = stdout or self.devnull or sys.stdout
        self._stderr = stderr or self.devnull or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        self.devnull.close()

def get_item_parent(ncdict,item,attr):
    lst = [i for i in ncdict.keys() \
            if (attr in ncdict[i].keys() \
            and item in ncdict[i][attr]) \
            ]
    if len(lst) >= 1:
        return lst
    else: return None

def get_item_child(ncdict,item):
    parent = finditem(ncdict,item)
    return parent

def parse_date(indate):
    #print("Parsing date")
    if isinstance(indate,datetime):
        return indate
    elif isinstance(indate,str):
        #print('Translate to datetime')
        return parse(indate)
    else:
        print('Not able to parse input return as is')
        return indate

def dispersion_deep_water(T=None,k=None,l=None,cp=None,cg=None):
    """
    computes requested variable from dispersion relation in deep water
    """
    return

def dispersion_shallow_water(l=None,h=None,T=None):
    return

def dispersion_intermediate_water(l=None,h=None,T=None):
    return

def calc_deep_water_T(l=None):
    g = 9.81
    return np.sqrt(l*2*math.pi/g)

def calc_shallow_water_T(l=None,h=None):
    g = 9.81
    return l/np.sqrt(g*h)

def wave_length_mask_swim(ds,llim=50,ulim=2000):
    """
    remove all results for wavelengths below given llim
    """
    res = ds.where((calc_deep_water_T(ds)>llim)&(calc_deep_water_T(ds)<ulim))
    mask = ~np.isnan(res)
    return mask

def compute_quantiles(ts,lq):
    """
    fct to compute quantiles for given ts

    param: 
        ts - iterable of ts
        lq - iterable of quantiles

    return:
        qA - numpy array of quantiles
    """
    ts = marginalize(ts)
    return np.array([np.quantile(ts, q) for q in lq])

def get_obsdict(obstype):
    if obstype == 'insitu':
        obsdict = load_or_default('insitu_specs.yaml')
    elif obstype == 'satellite_altimeter':
        obsdict = load_or_default('satellite_specs.yaml')
    else:
        print("obstype", obstype, "is not applicable")
        obsdict = None
    return obsdict

def find_tagged_obs(tags, obstype):
    d = get_obsdict(obstype)
    l = []
    for t in tags:
        l += [k for k in d if t in d[k].get('tags', [''])]
    return list(np.unique(l))

def expand_nID_for_sensors(nID, obstype):
    obsdict = get_obsdict(obstype)
    sensors = list(obsdict[nID]['sensor'])
    return sensors

def date_dispatcher(date, date_incr='d', incr=1):
    dispatch_date = {
                'h': date_next_hour,
                'd': date_next_day,
                'm': date_next_month,
                'y': date_next_year
                }
    return dispatch_date[date_incr](date, incr)

def date_next_hour(date, incr):
    date += timedelta(hours=incr)
    return date

def date_next_day(date, incr):
    date += timedelta(days=incr)
    return date

def date_next_month(date, incr):
    return datetime((date + relativedelta(months=+incr)).year,
                    (date + relativedelta(months=+incr)).month,
                    1)

def date_next_year(date, incr):
    return datetime((date + relativedelta(years=+incr)).year,
                    (date + relativedelta(years=+incr)).month,
                    1)

def footprint_pulse_limited_radius(Hs: float, h: float, tau: float) -> float:
    """
    Pulse limited footprint radius according to Chelton et al. 2001
    as referenced in coastal altimtery book p. 458, EQ 17.1

    Footprint size across track depends on:
    - significant wave height
    - satellite specs:
        - height over ground (h)
        - pulse duration (tau)
    """

    # constants
    c = 299792458# m/s, speed of light
    R = 6371*10**3 # m, radius Earth

    # equation
    r = np.sqrt(((c * tau + 2*Hs)*h) / (1+(h / R)))

    return r


def build_xr_ds(var: tuple, varnames: tuple):
    ds = xr.Dataset({
            varnames[0]: xr.DataArray(
                    data=var[0],
                    dims=[varnames[3]],
                    coords={'time': var[3]},
                    attrs=variable_def[varnames[0]],
                    ),
            varnames[1]: xr.DataArray(
                    data=var[1],
                    dims=[varnames[3]],
                    coords={'time': var[3]},
                    attrs=variable_def[varnames[1]],
                    ),
            varnames[2]: xr.DataArray(
                    data=var[2],
                    dims=[varnames[3]],
                    coords={'time': var[3]},
                    attrs=variable_def[varnames[2]],
                    ),
            varnames[3]: xr.DataArray(
                    data=var[3],
                    dims=[varnames[3]],
                    coords={'time': var[3]},
                    attrs=variable_def[varnames[3]],
                    )
                },
            attrs={'title': 'wavy dataset'}
        )
    return ds
