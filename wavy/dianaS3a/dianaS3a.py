#!/usr/bin/env python

import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/wavy/dianaS3a')

# import libraries
from datetime import datetime, timedelta
from satmod import sentinel_altimeter as sa
from custom_nc import dumptonc_S3a
import os

# setup
now = datetime.now()
sdate = datetime(now.year,now.month,now.day,now.hour)- timedelta(hours=24)
edate = datetime(now.year,now.month,now.day,now.hour)
timewin = 0
region = 'mwam4'
outpath = '/lustre/storeB/project/fou/om/Diana/'

# get data
sa_obj = sa(sdate,edate=edate,timewin=timewin,region=region)

# dump to .ncfile
dumptonc_S3a(sa_obj,outpath,mode = 'diana')
