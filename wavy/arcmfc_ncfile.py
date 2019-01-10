#!/usr/bin/env python
import sys
import time
from datetime import datetime, timedelta
import calendar
from copy import deepcopy
import numpy as np
from dateutil.relativedelta import relativedelta
from calendar import monthrange
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from utils import grab_PID
from stationmod import matchtime
import os
import argparse
from argparse import RawTextHelpFormatter
from custom_nc import get_arcmfc_stats

parser = argparse.ArgumentParser(
    description="""
Main program to run the monthly ARCMFC validation
with Sentinel.\n
Usage example in unix command line: 
./validate_arcmfc.py -d 201808\n
The argument consists of the year and month to be validated
If no date is given the last month is validated.
    """,
    formatter_class = RawTextHelpFormatter
    )
parser.add_argument('-d',help="validation for given month",
                    type=int)
args = parser.parse_args()
#args, rest = parser.parse_known_args()

if len(sys.argv)<=2:
    now = datetime.now()-relativedelta(months=1)
else:
    nowstr = str(sys.argv[2])
    now = datetime(int(nowstr[0:4]),int(nowstr[4:6]),1)

# retrieve PID
grab_PID()

flatten = lambda l: [item for sublist in l for item in sublist]
# ---
# Create Array 
# dim[time = days of month x 4,
#            forecasts=10,
#            surface=1,
#            metrics=4,
#            areas=1]
# ---
pathstr=('/lustre/storeA/project/fou/om/waveverification/Arc-MFC/' 
        + 'monthly/' 
        + now.strftime('%Y') 
        + '_' 
        + now.strftime('%b') + '/')
filestr=('product_quality_stats_ARCTIC_ANALYSIS_FORECAST_WAV_002_006_' 
        + now.strftime('%Y') 
        + now.strftime('%m') 
        + '01-' 
        + now.strftime('%Y') + now.strftime('%m') 
        + str(monthrange(now.year, now.month)[1]) + '.nc')

# cp original validation file to new file that can be changed
filestr_new = pathstr + filestr + ".test"
os.system("cp " + pathstr + filestr + " " + filestr_new)

nc = Dataset(pathstr + filestr , 'r')
#time:units = "days since 2001-01-01 12:00:00 UTC" ;
basetime = datetime(2001,1,1,12)
nctime = nc.variables['time'][:]
nc.close()
nc_start_time = basetime + timedelta(days=nctime[0])
nc_end_time = basetime + timedelta(days=nctime[-1])

# Retrieve time stamp for performance check
loop_time_start = time.time()

start_date = nc_start_time
end_date = nc_end_time
tmp_date = deepcopy(start_date)
forecasts = [12, 36, 60, 84, 108, 132, 156, 180, 204, 228]
M=np.ones([len(nctime),10,1,4,1])*9999.
dictlst_all=[]
excepts_all=[]
count1 = 0

while (tmp_date <= end_date):
    print (count1)
    if sys.version_info <= (3, 0):
        init_dates = map(lambda x: tmp_date - timedelta(hours=x),forecasts)
    else:
        init_dates = list(map(lambda x: tmp_date 
                        - timedelta(hours=x),forecasts))
    fc_date = init_dates[0]
    dictlst = []
    excepts = []
    count2 = 0
    for element in init_dates:
        print ("Validation for init_date: " 
        + str(element) 
        + "\n" 
        + " and fc_date: " 
        + str(fc_date))
        try:
            init_date = element
            start_time = time.time()
            # ---
            # Get stats ts
            inpath='/lustre/storeB/project/fou/om/ARCMFC/S3a/ValidationFiles/'
            filename_stats = fc_date.strftime("ARCMFC_val_ts_lt"
                                + "{:0>3d}".format(element)
                                + "h_%Y%m.nc")
            valid_dict, dtime = get_arcmfc_stats(inpath + filename_stats)
            # ---
            dictlst.append(valid_dict)
            if np.isnan(valid_dict['msd']):
                print ("msd is np.nan --> no values in range")
                print ("For init_date: ", init_date)
                print ("and fc_date: ", fc_date)
                M[count1,count2,0,3,0]=0
                excepts.append([fc_date,init_date])
                dictlst.append(9999.)
            else:
                dictnames=['mop','mor','msd','nov']
                for i in range(len(dictnames)):
                    M[count1,count2,0,i,0]=valid_dict[dictnames[i]]
            count2=count2+1
        except IOError as error:
            #print "No pass for date: ", str(init_date)
            print (error, " for date: ", str(init_date))
            print ("!!! Model run missing !!!")
            excepts.append([fc_date,init_date])
            dictlst.append(9999.)
            M[count1,count2,0,3,0]=0
            count2=count2+1
        except IndexError as error:
            print (error, " for date: ", str(init_date))
            print ("!!! No Sentinel pass !!!")
            excepts.append([fc_date,init_date])
            dictlst.append(9999.)
            M[count1,count2,0,3,0]=0
            count2=count2+1
        except ValueError as error:
            print (error, " for date: ", str(init_date))
            print ("!!! Date not in model file !!!")
            excepts.append([fc_date,init_date])
            dictlst.append(9999.)
            M[count1,count2,0,3,0]=0
            count2=count2+1
    count1=count1+1
    dictlst_all.append(dictlst)
    excepts_all.append(excepts)
    tmp_date = tmp_date + timedelta(hours=6)
loop_time = time.time() - loop_time_start
print ("Seconds needed for entire loop: ", loop_time)

#print ("\nAppending results to existing netcdf validation file ...")
#
#nc = Dataset(filestr_new, 'r+')
#nc.renameVariable('stats_VHM0','stats_VHM0_platform')
#nc_stats_VHM0_altimeter = nc.createVariable(
#                        'stats_VHM0_altimeter',
#                        np.float64,
#                        dimensions=('time',
#                        'forecasts',
#                        'surface',
#                        'metrics',
#                        'areas',),
#                        fill_value=9999.)
#nc_stats_VHM0_altimeter[:] = M
#nc_stats_VHM0_altimeter.standard_name = \
#                        "sea_surface_wave_significant_height"
#nc_stats_VHM0_altimeter.parameter = "stats_VHM0_altimeter"
#nc_stats_VHM0_altimeter.units = "m"
#nc_stats_VHM0_altimeter.reference = \
#                        "wave data from Sentinel-3a altimeter"
#nc_stats_VHM0_altimeter.reference_source = \
#                        "WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001"
#nc.close()
#print ("Data appended")
#print ("### VALIDATION FINISHED ###")
