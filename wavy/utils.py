"""
utility fcts for the verification
"""
import numpy as np
from datetime import datetime, timedelta
from math import radians, cos, sin, asin, sqrt, floor
import sys
from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

def block_detection(time,deltalim=None):
    if deltalim is None:
        deltalim = 1
    # forward check
    idx_a = []
    for i in range(1,len(time)):
        delta_t = time[i]-time[i-1]
        if delta_t>deltalim:
            idx_a.append(i)
    # backward check
    idx_b = []
    for i in range(0,len(time)-1):
        delta_t = time[i+1]-time[i]
        if delta_t>deltalim:
            idx_b.append(i)
    blocklst = []
    for i in range(len(idx_a)):
        if i == 0:
            tmp = [0,idx_b[i]]
            blocklst.append(tmp)
        if i < len(idx_a)-1:
            tmp = [idx_a[i],idx_b[i+1]]
            blocklst.append(tmp)
        if i == len(idx_a)-1:
            tmp = [idx_a[i],len(time)-1]
            blocklst.append(tmp)
    return idx_a, idx_b, blocklst

def identify_outliers_GP(x,y,mag):
    # mag -> magnitude in units std
    # if len(y)<11 flag all
    idx = []
    if len(y)<3:
        print ("block too short with length " + str(int(len(y))) + ", values flagged")
        x_pred = []
        y_pred = []
        sigma = []
        for i in range(len(y)):
            idx.append(i)
            x_pred.append(np.nan)
            y_pred.append(np.nan)
            sigma.append(np.nan)
    else:
        kernel = 1* RBF(length_scale=1) + WhiteKernel(noise_level=1)
        gp = gaussian_process.GaussianProcessRegressor(kernel=kernel)
        gp.fit(x.reshape(-1,1), y.reshape(-1,1))
        print (gp.kernel_)
        x_pred = x.reshape(-1,1)
        y_pred, sigma = gp.predict(x_pred, return_std=True)
        for i in range(len(y)):
            if ((y[i]>(y_pred[i]+(mag*sigma[i]))[0]) or (y[i]<(y_pred[i]-(mag*sigma[i]))[0])):
                idx.append(i)
    return idx,x_pred,y_pred,sigma

def identify_outliers(time,ts,ts_ref=None,hs_ll=None,hs_ul=None,dt=None):
    """
    time -> time series to check neighbour values
    ts -> time series to be checked for outliers
    ts_ref -> time series to compare to (optional)
    hs_ll -> hs lower limit over which values are checked
    hs_ul -> values over limit are being rejected
    """
    if hs_ll is None:
        hs_ll = 1.
    if hs_ul is None:
        hs_ul = 30.
    std_ts = np.nanstd(ts)
    mean_ts = np.nanmean(ts)
    # forward check
    idx_a = []
    for i in range(1,len(ts)):
        # transform to z
        if len(ts)<25:
            z = (ts[i] - np.nanmean(ts[:]))/np.nanstd(ts[:])
        elif i<13:
            z = (ts[i] - np.nanmean(ts[0:25]))/np.nanstd(ts[0:25])
        elif (i>=13 and i<(len(ts)-12)):
            z = (ts[i] - np.nanmean(ts[i-12:i+12]))/np.nanstd(ts[i-12:i+12])
        elif i>len(ts)-12:
            z = ((ts[i] - np.nanmean(ts[(len(ts-1)-25):-1]))
                /np.nanstd(ts[(len(ts-1)-25):-1]))
        if dt == True:
            delta_t = (time[i]-time[i-1]).total_seconds()
        else:
            delta_t = time[i]-time[i-1]
        if delta_t<2:
            #reject if value triples compared to neighbor
            # reject if greater than twice std (>2z)
            if ( ts[i] > hs_ll and ((ts[i-1] >= 3. * ts[i]) or (z>2)) ):
                idx_a.append(i)
        elif (ts[i] > hs_ll and z>2):
            idx_a.append(i)
    print (len(idx_a))
    # backward check
    idx_b = []
    for i in range(0,len(ts)-1):
        # transform to z
        if len(ts)<25:
            z = (ts[i] - np.nanmean(ts[:]))/np.nanstd(ts[:])
        elif i<13:
            z = (ts[i] - np.nanmean(ts[0:25]))/np.nanstd(ts[0:25])
        elif (i>=13 and i<len(ts)-12):
            z = (ts[i] - np.nanmean(ts[i-12:i+12]))/np.nanstd(ts[i-12:i+12])
        elif i>len(ts)-12:
            z = ((ts[i] - np.nanmean(ts[(len(ts-1)-25):-1]))
                /np.nanstd(ts[(len(ts-1)-25):-1]))
        if dt == True:
            delta_t = (time[i+1]-time[i]).total_seconds()
        else:
            delta_t = time[i+1]-time[i]
        if delta_t<2:
            #reject if value triples compared to neighbor
            # reject if greater than twice std (>2z)
            if ( ts[i] > hs_ll and ((ts[i+1] <= 1/3. * ts[i]) or (z>2)) ):
                idx_b.append(i)
        elif (ts[i] > hs_ll and z>2):
            idx_b.append(i)
    print (len(idx_b))
    idx_c = []
    for i in range(len(ts)):
        # reject if hs>hs_ul
        if ts[i]>hs_ul:
            idx_c.append(i)
    idx = np.unique(np.array(idx_a + idx_b + idx_c))
    if len(idx)>0:
        print(str(len(idx)) 
                + ' outliers detected of ' 
                + str(len(time)) 
                + ' values')
        return idx
    else:
        print('no outliers detected')
        return []

def grab_PID():
    """
    function to retrieve PID and display it to be able to kill the 
    python program that was just started
    """
    import os
    # retrieve PID
    PID = os.getpid()
    print ("\n")
    print ("PID - with the license to kill :) ")
    print (str(PID))
    print ("\n")
    return

def haversine(lon1, lat1, lon2, lat2):
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

def rmsd(a,b):
    '''
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    '''
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    n = len(a1)
    diff2 = (a1-b1)**2
    msd = diff2.sum()/n
    rmsd = np.sqrt(msd)
    return msd, rmsd

#def scatter_index(obs,model):
#    msd,rmsd = rmsd(obs,model)
#    SI = rmsd/np.nanmean(obs)*100.
#    return SI

def scatter_index(obs,model):
    from utils import rmsd
    msd,rmsd = rmsd(obs,model)
    stddiff = np.std(obs-model)
    SIrmse = rmsd/np.nanmean(obs)*100.
    SIstd = stddiff/np.nanmean(obs)*100.
    return SIrmse,SIstd

def corr(a,b):
    '''
    root mean square deviation
    if nans exist the prinziple of marginalization is applied
    '''
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    corr = np.corrcoef(a1,b1)[1,0]
    return corr

def bias(a,b):
    """
    if nans exist the prinziple of marginalization is applied
    """
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    bias = np.sum(a1-b1)/N
    return bias

def mad(a,b):
    """
    mean absolute deviation
    if nans exist the prinziple of marginalization is applied
    """
    a,b = np.array(a),np.array(b)
    comb = a + b
    idx = np.array(range(len(a)))[~np.isnan(comb)]
    a1=a[idx]
    b1=b[idx]
    N = len(a1)
    mad = np.sum(np.abs(a1-b1))/N
    return mad

def runmean(vec,win,mode=None):
    """
    input:  vec = vector of values to me smoothed
            win = window length
            mode= string: left, centered, right
    """
    win = int(win)
    if mode is None:
        mode='centered'
    out = np.zeros(len(vec))*np.nan
    std = np.zeros(len(vec))*np.nan
    length = len(vec)-win+1
    if mode=='left':
        count = int(win-1)
        start = int(win-1)
        for i in range(length):
            out[count] = np.mean(vec[count-start:count+1])
            std[count] = np.std(vec[count-start:count+1])
            count = count+1
    elif mode=='centered':
        count = int(floor(win/2))
        start = int(floor(win/2))
        for i in range(length):
            if win%2==0:
                sys.exit("window length needs to be odd!")
            else:
                sidx = int(count-start)
                eidx = int(count+start+1)
                out[count] = np.mean(vec[sidx-start:eidx])
                std[count] = np.std(vec[sidx:eidx])
                count = count+1
    elif mode=='right':
        count = int(0)
        for i in range(length):
            out[count] = np.mean(vec[i:i+win])
            std[count] = np.std(vec[i:i+win])
            count = count+1
    return out, std

def bootstr(a,reps):
    """
    input:    - is a time series of length n
              - reps (number of repetitions)
    output:   - an array of dim n x m where
                m is the number of repetitions
              - indices of draws
    """
    n = len(a)
    b = np.random.choice(a, (n, reps))
    bidx = np.zeros(b.shape) * np.nan
    for i in range(len(a)):
        tmp = np.where(b==a[i])
        bidx[tmp[0],tmp[1]] = i
        del tmp
    return b, bidx.astype('int')

def marginalize(a,b=None):
    if b is None:
        a = np.array(a)
        return a[~np.isnan[a]]
    else:
        a,b = np.array(a),np.array(b)
        comb = a + b
        idx = np.array(range(len(a)))[~np.isnan(comb)]
        a1=a[idx]
        b1=b[idx]
        return a1,b1,idx

def disp_validation(valid_dict):
    print('\n')
    print('# ---')
    print('Validation stats')
    print('# ---')
    print('Correlation Coefficient: ' + '{:0.2f}'.format(valid_dict['corr']))
    print('Root Mean Squared Error: ' + '{:0.2f}'.format(valid_dict['rmsd']))
    print('Mean Absolute Error: ' + '{:0.2f}'.format(valid_dict['mad']))
    print('Bias: ' + '{:0.2f}'.format(valid_dict['bias']))
    print('Scatter Index: ' + '{:0.2f}'.format(valid_dict['SI'][1]))
    print('Mean of Model: ' + '{:0.2f}'.format(valid_dict['mop']))
    print('Mean of Observations: ' + '{:0.2f}'.format(valid_dict['mor']))
    print('Number of Collocated Values: ' + str(valid_dict['nov']))
    print('\n')
    pass
