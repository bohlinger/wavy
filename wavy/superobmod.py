from utils import runmean_conv
from pygam import LinearGAM, l, s, ExpectileGAM
import numpy as np
import yaml


def superobbing(obs_obj,superob='block',outlier_detection='gam',\
missing_data_treatment='marginalize',**kwargs):
    """
    Applies a smoothing filter to create a super-observed ts
    **kwargs includes method specific input for chosen smoother
    Smoother on wish list are:
            block-average
            running mean using convolution
            GP
            GAM
            ...
    Caution:    for some smoothers much more of time series has 
                to be included.
    """
    print('under construction')
    # !if for satellites first a landmask has to be created!
    # remove outliers
    ol_dict = detect_outliers(obs_obj,method=outlier_detection)
    # create list of time stamps depending on choice: e.g. hourly timelist
    sd = obs_obj.sdate
    ed = obs_obj.edate
    
    # super observations are computed from cleaned time series
    # for chosen list of time stamps
    if missing_data_treatment == 'marginalize':
        # fill with blanks
        pass
    elif missing_data_treatment == 'impute':
        # fill with interpolation
        pass
    return sobs_ts

def detect_outliers(obs_obj,method=None,**kwargs):
    ol_dict={}
    dt = obs_obj.vars['datetime']
    x = obs_obj.vars['time']
    y = obs_obj.vars[obs_obj.stdvarname]
    # coars use approximate limits (in future use range from yaml)
    if (obs_obj.varalias == 'Hs' or obs_obj.varalias == 'Tm02'):
        ulim,llim = 30,0
    if (obs_obj.varalias == 'Mdir'):
        ulim,llim = 360,0
    # rogorous removal use techniques like: 
    # blockVariance, GP, GAM, (quantile regression) random forest, ...
    if method=='gam':
        idx = ol_linearGAM(x,y,**kwargs)
        ts_clean = y
        ts_clean[idx] = np.nan
    if method=='expectile':
        idx = ol_expectileGAM(x,y,**kwargs)
        ts_clean = y
        ts_clean[idx] = np.nan
    ol_dict['indices']=idx
    ol_dict['ts_clean']=ts_clean
    return ol_dict

def ol_expectileGAM(x,y,**kwargs):
    X = np.array(x).reshape(len(obs_obj.vars['time']),1)
    y = np.array(y)
    if 'n_splines' in kwargs.keys():
        n_splines = kwargs['n_splines']
    else:
        # This is because the automatic approach is too smooth
        n_splines = int(len(y)/10)
    gam50 = ExpectileGAM(expectile=0.50,terms=s(0),\
                        n_splines=int(len(y)/10)).gridsearch(X, y)
    # This practice of copying makes the models 
    # less likely to cross and much faster
    # https://pygam.readthedocs.io/en/latest/notebooks/tour_of_pygam.html
    # and copy the smoothing to the other models
    lam = gam50.lam
    # now fit a few more models
    gam95 = ExpectileGAM(expectile=0.95, lam=lam, terms=s(0),\
                        n_splines=int(len(y)/10)).fit(X, y)
    gam05 = ExpectileGAM(expectile=0.05, lam=lam, terms=s(0),\
                        n_splines=int(len(y)/10)).fit(X, y)
    ulim = gam95.predict(XX)
    llim = gam05.predict(XX)
    idx = [i for i in range(len(y)) \
            if (y[i]>ulim[i,1] or y[i]<llim[i])]
    return idx

def ol_linearGAM(x,y,**kwargs):
    X = np.array(x).reshape(len(obs_obj.vars['time']),1)
    y = np.array(y)
    if 'n_splines' in kwargs.keys():
        n_splines = kwargs['n_splines']
    else:
        # This is because the automatic approach is too smooth
        n_splines = int(len(y)/10)
    gam = LinearGAM(n_splines=n_splines,terms=s(0)).gridsearch(X, y)
    # sample on the input grid
    XX = X
    means = gam.predict(XX)
    bounds = gam.prediction_intervals(XX, width=.95)
    idx = [i for i in range(len(y)) \
            if (y[i]>bounds[i,1] or y[i]<bounds[i,0])]
    return idx
