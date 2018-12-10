#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
This module encompasses methods for statistical analysis and validation 
of wave models. I try to mostly follow the PEP convention for python 
code style. Constructive comments on style and effecient programming 
are most welcome!
Future plans involve expanding the help functione.g. 
program [-option] [sdate] [edate]
where the options can be e.g.   -h for help
                                -p for plot
                                -s for save
                                -d for download
                                - ...
'''
__version__ = "0.0.1"
__author__="Patrik Bohlinger, Norwegian Meteorological Institute"
__maintainer__ = "Patrik Bohlinger"
__email__ = "patrikb@met.no"
__contributors__ = "Patrik Bohlinger (MET), Johannes Roehrs (MET)"
__status__ = "Development in start phase. Currently only 
            duplicate of Johannes' dataanalysis.py. Is to be 
            extended in the future."

# --- import libraries ------------------------------------------------#
'''
List of libraries needed for this module. Sorted in categories to serve
effortless orientation. May be combined at some point.
'''
# ignore irrelevant warnings from matplotlib for stdout
import warnings
warnings.filterwarnings("ignore")

# all
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import os
import netCDF4 as netCDF4
import calendar

# libraries for parallel computing
from joblib import Parallel, delayed

# Johannes' functions
import scipy as sp
import scipy.ndimage as nd
from numpy import ma
import pylab as pl
from scipy import stats
# ---------------------------------------------------------------------#

parser = argparse.ArgumentParser(
        description="""
        Module encompassing classes and methods for statistical
        analysis and valdiation of wave models.\n
        Usage: \n
        ....
        """,
        formatter_class = RawTextHelpFormatter
        )
parser.add_argument('-g', choices=['a', 'b'], default='a',
        help="""Place holder options, where \n
        a = alpha\n
        b = beta\n
        """)
args = parser.parse_args()

# ---------------------------------------------------------------------#
#
# basic statistics
#
def rmsd(a,b):
    '''
    root mean square deviation
    '''
    a,b = sp.array(a),sp.array(b)
    n = len(a)
    diff2 = (a-b)**2
    return sp.sqrt(diff2.sum()/n)

def mea(a,b):
    '''
    absolute mean error
    '''
    a,b = sp.array(a),sp.array(b)
    n = len(a)
    diff = abs(a-b)
    return diff.sum()/n

def bias(a,b):
    '''
    bias
    '''
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a, b = a[mask], b[mask]
    return a.mean()-b.mean()

def pearsonr(a,b):
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a, b = a[mask], b[mask]
    return stats.pearsonr(a,b)[0]

def slope(a,b):
    return stats.linregress(returnclean(a,b))[0]

def returnclean(a,b):
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a, b = a[mask], b[mask]
    return a,b

def returnclean_withzeros(a,b):
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a[mask==False]=0.
    b[mask==False]=0.
    return a,b


def correlate_vectors(u1,v1, u2,v2):
    ''' 
    u1,v1 : series of vector 1 components
    u2,v2 : series of vector 2 components
    
    returns:
     1 if vectors are always the same
     0 if they are uncorrelated
    -1 if they are opposite'''
    u1 = sp.array(u1)
    u2 = sp.array(u2)
    v1 = sp.array(v1)
    v2 = sp.array(v2)
    upper = u1*u2 + v1*v2
    lower = u1**2 + v1**2 + u2**2 + v2**2
    r =  2*upper.sum() / lower.sum()
    return r

def bootstrap2(sampleA, sampleB, samplesize = None, nsamples = 500, statfunc = sp.mean):
    """
    Arguments:
       sampleA - input sample of values
       sampleB - input of a paired sample, same lenght as sampleA
       nsamples - number of samples to generate
       samplesize - sample size of each generated sample
       statfunc- statistical function to apply to each generated sample.
 
    Performs resampling from sample with replacement, gathers
    statistic in a list computed by statfunc on the each generated sample.
    """
    if samplesize is None:                                                                   
        samplesize=len(sampleA)
    n = len(sampleA)
    X = []
    for i in range(nsamples):
        #print "i = ",  i, 
        selection = stats.randint.rvs(0, n-1, size=samplesize)
        resampleA = [sampleA[j] for j in selection] 
        resampleB = [sampleB[j] for j in selection] 
        x = statfunc(resampleA, resampleB)
        X.append(x)
    return X

#
# vector operations
#

def uv_path(u,v,time):
    '''integrate u,v to get a path

    time: datetime list for u,v
    u and v is m/s
    time as list of datetime objects
    path is returned with units in meters
    '''
    xpos=[0.]
    ypos=[0.]
    for t in range(len(time)-1):
        if sp.isscalar(u[t]):
            dt = time[t+1]-time[t]
            dx,dy = u[t] * dt.seconds , v[t] * dt.seconds
            xpos.append(xpos[-1]+dx)
            ypos.append(ypos[-1]+dy)
    return sp.array(xpos), sp.array(ypos)


def DD_FF(u,v):
    ''' calculates wind/current speed and direction from u and v components
    #
    if u and v are easterly and northerly components,
    DD is heading direction of the wind. 
    to get meteorological-standard, call DD_FF(-u, -v)
    '''
    DD = ma.arctan2(u, v)*180/sp.pi
    DD[DD < 0] = 360 + DD[DD < 0]
    FF = ma.sqrt(u**2 + v**2)
    return DD, FF

def UV(DD,FF,met=False):
    '''
    #requires testing
    DD is the heading direction of the wind
    to get met. standart, call UV(DD,FF,met=True)
    '''
    u = FF * sp.sin(DD*sp.pi/180)
    v = FF * sp.cos(DD*sp.pi/180)
    if met==True:
        u,v=-u,-v
    return u,v

#
# filters
#

def GodinTypeFilter(data, n, axis=0):
    ''' perform 3 times moving average over the specified array axis.
    suitable for time averaging
    '''
    weights = sp.ones((n),sp.float32)/n
    weights2 = sp.ones((n+1),sp.float32)/(n+1)
    data=nd.convolve1d(data, weights, axis=axis, mode='constant')
    data=nd.convolve1d(data, weights, axis=axis, mode='constant')
    data=nd.convolve1d(data, weights2, axis=axis, mode='constant')
    return data

def smooth2d(I,N):
    """Box average filter 2d"""
    kernel=sp.ones((N,N),sp.float32)/(N**2)
    return nd.convolve(I, kernel,mode='reflect')

def smooth1d(I,N):
    """running mean filter 1d"""
    kernel=sp.ones((N),sp.float32)/N
    return nd.convolve(I, kernel,mode='reflect')

def distance(p1,p2):
    ''' calculate distance between two points
    
    p1 and p2 are 2d coordinate tuples 
    
    returns distance in coordinate units'''
    return sp.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)


#
# vector field tools
#

def rotate_vector(v,alpha,deg='False'):
    '''
    conterclockwise rotation of a vector v by the angle alpha
    '''
    alpha = -alpha # ensure clockwise rotation (JR, July 2013)
    if not (deg == 'False'): alpha=alpha*sp.pi/180
    
    R=sp.array([[sp.cos(alpha),-sp.sin(alpha)],[sp.sin(alpha),sp.cos(alpha)]])
    rotatedv=sp.dot(R,v)
    return rotatedv

def rotate_vectorfield(U,V,alpha):
    '''rotate wind vectors clockwise. alpha may be a scalar or an array
	alpha is in degrees
	returns u,v '''
    alpha = sp.array(alpha)*sp.pi/180
    alpha = alpha.flatten()
    R = sp.array([[sp.cos(alpha), -sp.sin(alpha)], [sp.sin(alpha), sp.cos(alpha)] ])
    shpe = U.shape
    origwind = sp.array((U.flatten(), V.flatten()))
    if len(R.shape)==2:
        rotwind = dot(R, origwind) # for constant rotation angle
    else:
        # for rotation angle given as array with same dimensions as U and V:
        # k-loop with rotwind(k) = dot(R(i,j,k), origwind(j,k)) (einstein summation indices)
        rotwind = sp.einsum("ijk,ik -> jk", R, origwind)  # einstein summation indices
    Urot ,Vrot = rotwind[0,:], rotwind[1,:]
    Urot = Urot.reshape(shpe)
    Vrot = Vrot.reshape(shpe)
    return Urot, Vrot

def north_direction(lat):
    '''get the north direction relative to image positive y coordinate'''
    dlatdx = nd.filters.sobel(lat,axis=1,mode='constant',cval=sp.nan) #gradient in x-direction
    dlatdy = nd.filters.sobel(lat,axis=0,mode='constant',cval=sp.nan)
    ydir = lat[-1,0] -lat[0,0] # check if latitude is ascending or descending in y axis
    # same step might have to be done with x direction.
    return sp.arctan2(dlatdx,dlatdy*sp.sign(ydir) )*180/sp.pi

def anorm(a):
    '''calculate the norm of each row vector
    result[i]=norm(a[i,:])'''
    n=a.shape[0]
    betrag=sp.zeros(n)
    for i in range(n):
         betrag[i]=norm(a[i,:])
    return betrag


def timefind(timearray,time):
    '''find time index in timearray that is closest to time'''
    timea=timearray.copy()
    timea[timea>time]=time-(timea[timea>time]-time) # create function with a maximum at the desired position (time)
    return list(timea).index(max(timea))

def find_pos(lon,lat,lon_0,lat_0,f=1):
    '''find the best fitting position of a coordinate pair in an image
    lon: longitude image
    lat: latitude image
    lon_0, lat_0: coordinates of desired position; f: scale factor to speed up the calculation, but reduces the accurance of the output position; returns a tuple x,y, which gives the best position in lon,lat of lon_0,lat_0'''
    best_lat_i=lat[::f,::f].copy()
    best_lat_i[best_lat_i>lat_0]=2*lat_0-best_lat_i[best_lat_i>lat_0]
    best_lon_i=lon[::f,::f].copy()
    best_lon_i[best_lon_i>lon_0]=2*lon_0-best_lon_i[best_lon_i>lon_0]
    best_position=best_lat_i/lat_0.__abs__()+best_lon_i/lon_0.__abs__()
    best_y,best_x=nd.measurements.maximum_position(best_position)
    best_y,best_x=f*best_y,f*best_x
    return best_x,best_y

def find_pos1d(x,y,x0,y0):
    '''find index in x,y that suits the point x0,y0 most
    x, y are one-dimensional arrays or sequences of x and y coordinates
    returns the best fitting index of x,y as integer'''
    x,y = sp.array(x), sp.array(y)
    xdist=x-x0
    ydist=y-y0
    f=xdist**2+ydist**2 #calculate distance**2
    flist=list(f[sp.isfinite(f)]) # make a list without NaNs
    flist.sort() # sort the list, so that the smallest distance**2 is at the top
    return list(f).index(flist[0]) # return index of closest distance**2

def find_pos1d1(x,x0):
    x=sp.array(x)
    xdist=(x-x0)**2
    flist=list(xdist[sp.isfinite(xdist)])
    flist.sort()
    return list(xdist).index(flist[0])

def findtime(time,time0):
    time0 = pl.date2num(time0)
    timelist = pl.date2num(time)
    return find_pos1d1(timelist, time0)


#
# image tools
#

def downsample(img, n):
    '''resample image to lower resolution
    img: 2d-array
    n: factor, by which the image resolution will be reduced
    outputs the image on lower resolution
    '''
    Nin = img.shape
    Nout = Nin[0]/n, Nin[1]/n # get new image dimenstions
    # resize image in a way that all chunks to be averaged get stacked into a new dimension
    r_img = sp.zeros((n**2,Nout[0], Nout[1]), dtype=img.dtype) # initialize resized image
    for i in range(n):
        for j in range(n):
            r_img[i*n+j] = img[i:n*Nout[0]:n,j:n*Nout[1]:n]
    return r_img.mean(axis=0)


#
# interpolations
#
