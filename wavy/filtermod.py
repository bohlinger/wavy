import numpy as np
from copy import deepcopy
import yaml
import os
from datetime import datetime, timedelta
import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import roaring_landmask
import shapely.geometry as sgeom
from shapely.ops import unary_union
import shapely.vectorized

# own imports
from wavy.utils import find_included_times, collocate_times
from wavy.wconfig import load_or_default
from wavy.utils import runmean_conv

ROAR = None

variable_info = load_or_default('variable_info.yaml')

flatten = lambda l: [item for sublist in l for item in sublist]

def filter_main(vardict_in,varalias='Hs',**kwargs):
    """ Governing function of filtermod

    Tasks:
        - check if prior/post transforms are needed
        - check if cleaning is needed
        - check if filter is needed
        - check if land mask is needed
            - if so apply cleaning/filters to subsets
              i.e. each chunk will be fed into filter_data
              and consolidated when finished with all chunks

    Args:
        vardict

    Returns:
        vardict
    """
    stdvarname = variable_info[varalias]['standard_name']
    # clone vardict_in
    vardict = deepcopy(vardict_in)
    # rm NaNs
    vardict = rm_nan_from_vardict(varalias,vardict)
    # apply physical limits
    if kwargs.get('limits') is not None:
        vardict = apply_limits(varalias,vardict)

    # start main filter section
    if ( kwargs.get('land_mask') is not None \
    and kwargs['land_mask'] == True ):
        del kwargs['land_mask']
        vardict,sea_mask = apply_land_mask(vardict,**kwargs)
        indices = start_stop(sea_mask, True)
        # initialize newvardict with same keys
        newvardict = deepcopy(vardict)
        for key in vardict:
            if (key != 'time_unit' and key != 'meta'):
                newvardict[key] = []
        for start_idx, stop_idx in indices:
            lenofchunk = len(list(range(start_idx, stop_idx)))
            print('Length of chunk:',lenofchunk)
            tmpdict = {}
            for key in vardict:
                if (key != 'time_unit' and key != 'meta'):
                    tmpdict[key] = vardict[key][start_idx:stop_idx]
                else:
                    tmpdict[key] = vardict[key]
            newtmpdict = filter_main(tmpdict,
                                     varalias=varalias,
                                     **kwargs)
            # append to newvardict
            for key in tmpdict:
                if (key != 'time_unit' and key != 'meta'):
                    newvardict[key].append(newtmpdict[key])
        # flatten newvardict lists
        for key in newvardict:
            if (key != 'time_unit' and key != 'meta'):
                newvardict[key] = flatten(newvardict[key])
        vardict = newvardict

    else:
        if kwargs.get('slider') is not None:
            # create chunks with size of slider
            vardict = filter_slider(vardict,varalias,**kwargs)
        else:
            if kwargs.get('priorOp') is not None:
                method = kwargs.get('priorOp')
                vardict = apply_priorOp(varalias,vardict,
                                        method = method)
            if kwargs.get('cleaner') is not None:
                #output_dates = kwargs.get('output_dates')
                method = kwargs.get('cleaner')
                #date_incr = kwargs.get('date_incr')
                vardict = apply_cleaner(varalias,vardict,
                                        method = method,
                                        **kwargs)
            if kwargs.get('smoother') is not None:
                output_dates = kwargs.get('output_dates')
                method = kwargs.get('smoother')
                vardict = apply_smoother(varalias,vardict,
                                         output_dates = output_dates,
                                         method = method,
                                         **kwargs)
            if kwargs.get('postOp') is not None:
                method = kwargs.get('postOp')
                vardict = apply_postOp(varalias,vardict,method = method)

    return vardict

def filter_slider(vardict,varalias,**kwargs):
        slider = kwargs['slider']
        del kwargs['slider']
        newvardict = deepcopy(vardict)
        for key in vardict:
            if (key != 'time_unit' and key != 'meta'):
                newvardict[key] = []
        for i in range(0,len(vardict['time']),slider):
            start_idx = i
            stop_idx = i + slider
            if i == range(0,len(vardict['time']),slider)[-1]:
                stop_idx = len(vardict['time']) -1
            tmpdict = {}
            for key in vardict:
                if (key != 'time_unit' and key != 'meta'):
                    tmpdict[key] = vardict[key][start_idx:stop_idx]
                else:
                    tmpdict[key] = vardict[key]
            newtmpdict = filter_main(tmpdict,
                                     varalias=varalias,
                                     **kwargs)
            # append to newvardict
            for key in tmpdict:
                if (key != 'time_unit' and key != 'meta'):
                    newvardict[key].append(newtmpdict[key])
        # flatten newvardict lists
        for key in newvardict:
            if (key != 'time_unit' and key != 'meta'):
                newvardict[key] = flatten(newvardict[key])
        return newvardict

def rm_nan_from_vardict(varalias,vardict):
    stdvarname = variable_info[varalias]['standard_name']
    nanmask = ~np.isnan(vardict[stdvarname])
    for key in vardict.keys():
        if (key != 'time_unit' and key != 'meta'):
            vardict[key] = list(np.array(vardict[key])[nanmask])
    return vardict

def start_stop(a, trigger_val):
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False,np.equal(a, trigger_val),False]
    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])
    # Get the start and end indices with slicing along the shifting ones
    return zip(idx[::2], idx[1::2]-1)

def apply_land_mask(vardict,**kwargs):
    """ Mask out parts covering land

    Args:
        land_mask_resolution (str): resolution of the GSHHG land mask in
                                    cartopy
        plot (bool): Plot the extent, satellite swath, and masked points
                        if True.
        conservative_mask (int): Number of points to extend the land mask
                                    in the ocean

    Returns:
        land_mask

    """
    global ROAR

    print('Apply land mask')
    if ROAR is None:
        ROAR = roaring_landmask.RoaringLandmask.new()

    longitudes = np.array(vardict['longitude'])
    latitudes = np.array(vardict['latitude'])
    land_mask = ROAR.contains_many(longitudes, latitudes)
    sea_mask = np.invert(land_mask)

    for key in vardict.keys():
        if (key != 'time_unit' and key != 'meta'):
            vardict[key] = list(np.array(vardict[key])[sea_mask])
    print('Number of disregarded values:', len(sea_mask[sea_mask==False]))
    return vardict, sea_mask

def apply_limits(varalias,vardict):
    print('Apply limits')
    print('Crude cleaning using valid range defined in variable_info.yaml')
    stdvarname = variable_info[varalias]['standard_name']
    clean_dict = deepcopy(vardict)
    y = vardict[stdvarname]
    llim = variable_info[varalias]['valid_range'][0]
    ulim = variable_info[varalias]['valid_range'][1]
    tmpdict = {'y':y}
    df = pd.DataFrame(data = tmpdict)
    dfmask = df['y'].between(llim, ulim, inclusive=True)
    for key in vardict:
        if (key != 'time_unit' and key != 'meta'):
            clean_dict[key] = list(np.array(vardict[key])[dfmask.values])
    return clean_dict

def square_data(varalias,vardict):
    stdvarname = variable_info[varalias]['standard_name']
    newdict = deepcopy(vardict)
    var_squared = list(np.array(vardict[stdvarname])**2)
    newdict[stdvarname] = var_squared
    return newdict

def root_data(varalias,vardict):
    stdvarname = variable_info[varalias]['standard_name']
    newdict = deepcopy(vardict)
    var_root = list(np.sqrt(np.array(vardict[stdvarname])))
    newdict[stdvarname] = var_root
    return newdict

def apply_priorOp(varalias,vardict,method=None):
    print("Prepare data prior to further treatment")
    print('Apply ',method)
    if method is None:
        newdict = vardict
    if method == 'square':
        newdict = square_data(varalias,vardict)
    return newdict

def apply_postOp(varalias,vardict,method=None):
    print("Transform data back after treatment")
    print('Apply ',method)
    if method is None:
        newdict = vardict
    if method == 'root':
        newdict = root_data(varalias,vardict)
    return newdict

def apply_cleaner(varalias,vardict,method='linearGAM',**kwargs):
    # rigorous data cleaning use techniques like:
    # blockVariance, GP, GAM, (quantile regression) random forest, ...
    print('Apply cleaner')
    print('Cleaning data using method:', method)
    stdvarname = variable_info[varalias]['standard_name']
    if kwargs.get('itr') is not None:
        itr = kwargs['itr']
    else: itr = 1
    for i in range(itr):
        clean_dict = rm_nan_from_vardict(varalias,vardict)
        dt = vardict['datetime']
        x = vardict['time']
        y = vardict[stdvarname]
        if method=='linearGAM':
            idx = cleaner_linearGAM(x,y,varalias,**kwargs)
            ts_clean = np.array(y)
            ts_clean[idx] = np.nan
        if method=='GP':
            idx = cleaner_GP(x,y,varalias,**kwargs)
            ts_clean = np.array(y)
            ts_clean[idx] = np.nan
        if method=='expectileGAM':
            idx = cleaner_expectileGAM(x,y,varalias,**kwargs)
            ts_clean = np.array(y)
            ts_clean[idx] = np.nan
        #clean_dict['indices'] = idx
        clean_dict[stdvarname] = ts_clean
        vardict = clean_dict
    return vardict

def apply_smoother(varalias,vardict,output_dates=None,method=None,date_incr=None,**kwargs):
    """
    Applies a smoother to the data
    **kwargs includes method specific input for chosen method
    Methods are:
            block-average
            running mean using convolution
            GP
            GAM
            Lanczos
            ...
    Caution:    for some smoothers much more of time series has
                to be included.
    """
    print('Apply smoother')
    print('Smooth data using method:',method)
    stdvarname = variable_info[varalias]['standard_name']
    newdict = deepcopy(vardict)
    # determine the output grid
    if (isinstance(date_incr,int) and output_dates is None):
    # increments are in #hours
    # create output grid --> list of time stamps depending on choice
        sd = vardict['datetime'][0]
        ed = vardict['datetime'][-1]
        steps = int((ed-sd).total_seconds()/(date_incr*60*60))+1
        tmpd = sd
        output_dates = [tmpd + timedelta(hours=i) \
                        for i in range(0,steps,date_incr) \
                        if (tmpd + timedelta(hours=i) <= ed)]
        del tmpd
    elif output_dates is None: # original datetimes are used
        output_dates = vardict['datetime']
    output_grid = netCDF4.date2num(output_dates,\
                                   units=vardict['time_unit'])
    smoothed_ts, idx = smoothing(varalias,newdict,output_grid,\
                                    output_dates,method=method,\
                                    date_incr=date_incr,\
                                    **kwargs)
    newdict[stdvarname] = list(smoothed_ts)
    for key in newdict:
        if (key != stdvarname and key != 'time_unit' and key != 'meta'):
            newdict[key] = list(np.array(newdict[key])[idx])
    return newdict

def smoothing(varalias,vardict,output_grid,\
output_dates, method='linearGAM', date_incr=None,**kwargs):
    stdvarname = variable_info[varalias]['standard_name']
    dt = vardict['datetime']
    x = vardict['time']
    y = vardict[stdvarname]
    X = output_grid
    idx = collocate_times(list(dt),output_dates)
    if method=='linearGAM':
        # NaNs need to be removed before gam
        tmpvar = np.array(y)
        tmptime = np.array(x)
        tmpdtime = np.array(dt)
        tmptime = tmptime[~np.isnan(tmpvar)]
        tmpdtime = tmpdtime[~np.isnan(tmpvar)]
        tmpvar = tmpvar[~np.isnan(tmpvar)]
        y = tmpvar
        x = tmptime
        dt = tmpdtime
        smoothed_ts = smoother_linearGAM(x,y,X,varalias,**kwargs)
    elif method=='expectileGAM':
        # NaNs need to be removed before gam
        tmpvar = np.array(y)
        tmptime = np.array(x)
        tmpdtime = np.array(dt)
        tmptime = tmptime[~np.isnan(tmpvar)]
        tmpdtime = tmpdtime[~np.isnan(tmpvar)]
        tmpvar = tmpvar[~np.isnan(tmpvar)]
        y = tmpvar
        x = tmptime
        dt = tmpdtime
        smoothed_ts = smoother_expectileGAM(x,y,X,varalias,**kwargs)
    elif method=='GP':
        # NaNs need to be removed before gp
        tmpvar = np.array(y)
        tmptime = np.array(x)
        tmpdtime = np.array(dt)
        tmptime = tmptime[~np.isnan(tmpvar)]
        tmpdtime = tmpdtime[~np.isnan(tmpvar)]
        tmpvar = tmpvar[~np.isnan(tmpvar)]
        y = tmpvar
        x = tmptime
        dt = tmpdtime
        smoothed_ts = smoother_GP(x,y,X,varalias,**kwargs)
    elif method=='NIGP':
        # NaNs need to be removed before gp
        tmpvar = np.array(y)
        tmptime = np.array(x)
        tmpdtime = np.array(dt)
        tmptime = tmptime[~np.isnan(tmpvar)]
        tmpdtime = tmpdtime[~np.isnan(tmpvar)]
        tmpvar = tmpvar[~np.isnan(tmpvar)]
        y = tmpvar
        x = tmptime
        dt = tmpdtime
        smoothed_ts = smoother_NIGP(x,y,X,varalias,**kwargs)
    elif method=='blockMean':
        # blocks are means from date_incr in hours
        # For each grid_input time_stamp compute mean of hour
        # if at least half of values are valid
        # else attribute NaN
        smoothed_ts = smoother_blockMean(dt,x,y,output_dates,date_incr)
    elif method=='lanczos':
        y = vardict[stdvarname]
        window = kwargs['window']
        cutoff = kwargs['cutoff']
        smoothed_ts = smoother_lanczos(y,window,cutoff)
        smoothed_ts = smoothed_ts[idx]
    else: print('Method not defined, please enter valid method')
    return smoothed_ts, idx

def lanczos_weights(window,cutoff):
    """ Calculate weights for a low pass Lanczos filter

    args:
        window: (integer) the length of the filter window
        cutoff: (float) the cutoff frequency in inverse time steps

    returns: weights

    example: https://scitools.org.uk/iris/docs/v1.2/examples/
             graphics/SOI_filtering.html
    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]

def smoother_lanczos(y,window,cutoff):
    from utils import runmean
    weights = lanczos_weights(window,cutoff)
    ts, std = runmean(y,window,mode='centered',weights=weights)
    return ts

def smoother_blockMean(dt,x,y,output_dates,date_incr):
    means = []
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    for i in range(len(output_dates)):
        # check if more than 50% of values are valid
        # if so compute mean
        idx = find_included_times(dt,\
                                  sdate = output_dates[i] - \
                                          timedelta(hours=date_incr),\
                                  edate=output_dates[i], twin=0)
        block = y[idx]
        nominator = len(block[np.isnan(block)])
        denominator = len(block)
        if denominator == 0:
            ratio = 1
        else:
            ratio = nominator/float(denominator)
        if ratio < 0.5:
            means.append(np.nanmean(block))
        else:
            means.append(np.nan)
    means = np.array(means)
    return means

def smoother_linearGAM(x,y,X,varalias,**kwargs):
    from pygam import LinearGAM, l, s
    if isinstance(x,list):
        x = np.array(x)
    x = x.reshape(len(x),1)
    if isinstance(y,list):
        y = np.array(y)
    if isinstance(X,list):
        X = np.array(X)
    if X is None:
        X = x.reshape(len(x),1)
    else:
        X = X.reshape(len(X),1)
    #if 'n_splines' in kwargs.keys():
    #    n_splines = kwargs['n_splines']
    #else:
    #    # This is because the automatic approach is too smooth
    #    n_splines = int(len(y)/5)
    #gam = LinearGAM(n_splines=n_splines,\
    #                terms=s(0,basis='ps')\
    #                ).gridsearch(x, y)
    gam = LinearGAM( terms=s(0,basis='ps')\
                    ).gridsearch(x, y )
    # sample on the input grid
    means = gam.predict(X)
    return means

def smoother_expectileGAM(x,y,X,varalias,**kwargs):
    from pygam import s, ExpectileGAM
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    if X is None:
        X = deepcopy(x)
    x = x.reshape(len(x),1)
    X = X.reshape(len(X),1)
    #if 'n_splines' in kwargs.keys():
    #    n_splines = kwargs['n_splines']
    #else:
    #    # This is because the automatic approach is too smooth
    #    n_splines = int(len(y)/5)
    if 'expectile' in kwargs.keys():
        expectile = kwargs['expectile']
    else:
        expectile = .5
    #gam50 = ExpectileGAM(expectile=expectile,terms=s(0),\
    #                    n_splines=n_splines).gridsearch(x, y)
    gam50 = ExpectileGAM(expectile=expectile,terms=s(0),\
                        ).gridsearch(x, y)
    # This practice of copying makes the models
    # less likely to cross and much faster
    # https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html
    # and copy the smoothing to the other models
    pred = gam50.predict(X)
    return pred

def smoother_GP(x,y,X,varalias,**kwargs):
    from sklearn import gaussian_process
    from sklearn.gaussian_process.kernels import RBF
    from sklearn.gaussian_process.kernels import WhiteKernel
    from sklearn.gaussian_process.kernels import RationalQuadratic
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    if isinstance(X,list):
        X = np.array(X)
    if X is None:
        X = x.reshape(-1,1)
    else:
        X = X.reshape(-1,1)
    x = x.reshape(-1,1)
    # create a zero mean process
    Y = y.reshape(-1,1) - np.nanmean(y)
    # define the kernel based on kwargs
    if 'kernel' in kwargs.keys():
        print('kernel is defined by user')
        kernel = kwargs['kernel']
    elif 'kernel_lst' in kwargs.keys():
        print('kernel constituents given by user')
        kernel = WhiteKernel(noise_level=1)
        if 'RBF' in kwargs['kernel_lst']:
            kernel += 1 * RBF(length_scale=1)
        if 'RationalQuadratic' in kwargs['kernel_lst']:
            kernel += 1 * RationalQuadratic(alpha=1,\
                                        length_scale=1)
    else:
        print('default kernel')
        kernel =  WhiteKernel(noise_level=1) +  1 * RBF(length_scale=1)
    gp = gaussian_process.GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=10)
    gp.fit(x, Y)
    print(gp.kernel_)
    y_pred, sigma = gp.predict(X, return_std=True)
    y_pred = y_pred + np.nanmean(y)
    return y_pred

def smoother_NIGP(x,y,X,varalias,**kwargs):
    import numpy as np
    from GPfcts import nll_fn_nigp
    from GPfcts import posterior_predictive_nigp
    from scipy.optimize import minimize
    from scipy.optimize import Bounds
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    if isinstance(X,list):
        X = np.array(X)
    if X is None:
        X = x.reshape(-1,1)
    else:
        X = X.reshape(-1,1)
    x = x.reshape(-1,1)
    # create a zero mean process
    Y = y.reshape(-1,1) - np.nanmean(y)
    # initialize using using standard GP
    mu = smoother_GP(x,Y,x,varalias,**kwargs)
    # define inits
    inits = np.array([1,1,1,1])
    # define bounds
    ulim = 1000
    bounds = Bounds([1, .001, .001, .001],[ulim,ulim,ulim,ulim])
    # continue with NIGP (depends on number of interations)
    if 'iter' in kwargs.keys():
        for i in range(kwargs['iter']):
            # get gradient
            fgrad = np.gradient(mu.ravel())
            # interpolate to points of interest
            fgrad_opt = np.interp(x.ravel(), X.ravel(), fgrad.ravel())
            fgrad_opt = fgrad_opt.reshape(-1,1)
            print(nll_fn_nigp(x, Y,Grad_fmean=fgrad_opt))
            # optimization
            res = minimize(
                    nll_fn_nigp(x, Y,Grad_fmean=fgrad_opt),
                    inits,
                    bounds=bounds,
                    method='L-BFGS-B')
                    #method='SLSQP')
            l_opt, sigma_f_opt, sigma_y_opt, sigma_x_opt = res.x
            # compute statistics
            mu, cov = posterior_predictive_nigp(
                                    x,x,Y,
                                    l = l_opt,
                                    sigma_f = sigma_f_opt,
                                    sigma_y = sigma_y_opt,
                                    sigma_x = sigma_x_opt,
                                    Grad_fmean = fgrad_opt )
    print(  'l:',l_opt,'sigma_f:',sigma_f_opt,
            'sigma_y:',sigma_y_opt,'sigma_x:',sigma_x_opt )
    # last step is to compute predictive posterior statistics
    # for output grid and add previously substracted mean
    mu, cov = posterior_predictive_nigp(
                                    X,x,Y,
                                    l = l_opt,
                                    sigma_f = sigma_f_opt,
                                    sigma_y = sigma_y_opt,
                                    sigma_x = sigma_x_opt,
                                    Grad_fmean = fgrad_opt )
    mu += np.mean(y)
    return mu

def cleaner_expectileGAM(x,y,varalias,**kwargs):
    from pygam import s, ExpectileGAM
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    X = x.reshape(len(x),1)
    #if 'n_splines' in kwargs.keys():
    #    n_splines = kwargs['n_splines']
    #else:
    #    # This is because the automatic approach is too smooth
    #    n_splines = int(len(y)/5)
    #gam50 = ExpectileGAM(expectile=.5,terms=s(0),\
    #                    n_splines=n_splines).gridsearch(X, y)
    gam50 = ExpectileGAM(expectile=.5,terms=s(0),\
                        ).gridsearch(X, y)
    # This practice of copying makes the models
    # less likely to cross and much faster
    # https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html
    # and copy the smoothing to the other models
    lam = gam50.lam
    # now fit a few more models
    if 'expectile_ulim' in kwargs.keys():
        expectile_ulim = kwargs['expectile_ulim']
    else:
        expectile_ulim = .95
    if 'expectile_llim' in kwargs.keys():
        expectile_llim = kwargs['expectile_llim']
    else:
        expectile_llim = .05
    #gam_ulim = ExpectileGAM(expectile=expectile_ulim, lam=lam,
    #                    terms=s(0),n_splines=n_splines).fit(X, y)
    #gam_llim = ExpectileGAM(expectile=expectile_llim, lam=lam,
    #                    terms=s(0),n_splines=n_splines).fit(X, y)
    gam_ulim = ExpectileGAM(expectile=expectile_ulim, lam=lam,
                        terms=s(0)).fit(X, y)
    gam_llim = ExpectileGAM(expectile=expectile_llim, lam=lam,
                        terms=s(0)).fit(X, y)
    ulim = gam_ulim.predict(X)
    llim = gam_llim.predict(X)
    idx = [i for i in range(len(y)) \
            if (y[i]>ulim[i] or y[i]<llim[i])]
    return idx

def cleaner_linearGAM(x,y,varalias,**kwargs):
    from pygam import LinearGAM, l, s
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    X = x.reshape(len(x),1)
    #if 'n_splines' in kwargs.keys():
    #    n_splines = kwargs['n_splines']
    #else:
    #    # This is because the automatic approach is too smooth
    #    #n_splines = int(len(y)/5)
    #gam = LinearGAM(n_splines=n_splines,\
    #                terms=s(0,basis='ps')\
    #                ).gridsearch(X, y)
    gam = LinearGAM(terms=s(0,basis='ps')).gridsearch(X, y)
    #gam = LinearGAM(n_splines=n_splines,terms=s(0)).gridsearch(X, y)
    # sample on the input grid
    means = gam.predict(X)
    bounds = gam.prediction_intervals(X, width=.95)
    idx = [i for i in range(len(y)) \
            if (y[i]>bounds[i,1] or y[i]<bounds[i,0])]
    return idx

def cleaner_GP(x,y,varalias,**kwargs):
    from sklearn import gaussian_process
    from sklearn.gaussian_process.kernels import RBF
    from sklearn.gaussian_process.kernels import WhiteKernel
    from sklearn.gaussian_process.kernels import RationalQuadratic
    if isinstance(x,list):
        x = np.array(x)
    if isinstance(y,list):
        y = np.array(y)
    X = x.reshape(-1,1)
    # create a zero mean process
    Y = y.reshape(-1,1) - np.nanmean(y)
    # define the kernel based on kwargs
    if 'kernel' in kwargs.keys():
        print('kernel is defined by user')
        kernel = kwargs['kernel']
    elif 'kernel_lst' in kwargs.keys():
        print('kernel constituents given by user')
        kernel = WhiteKernel(noise_level=1)
        if 'RBF' in kwargs['kernel_lst']:
            kernel += RBF(length_scale=1)
        if 'RationalQuadratic' in kwargs['kernel_lst']:
            kernel += RationalQuadratic(alpha=1,\
                                        length_scale=1)
    else:
        print('default kernel')
        kernel =  WhiteKernel(noise_level=1) +  1 * RBF(length_scale=1)
    gp = gaussian_process.GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=10)
    gp.fit(X, Y)
    print(gp.kernel_)
    y_pred, sigma = gp.predict(X, return_std=True)
    uplim = y_pred + (2*sigma).reshape(-1,1)
    lowlim = y_pred - (2*sigma).reshape(-1,1)
    idx = [i for i in range(len(Y)) \
            if (Y[i]>uplim[i] or Y[i]<lowlim[i])]
    return idx
