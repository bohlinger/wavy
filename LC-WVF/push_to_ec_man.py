#!/usr/bin/env python3
"""
needs: pysftp==0.2.8
bug in version 0.2.9
login for test if files are being pushed:
    lftp enmi_lcwfv@acquisition.ecmwf.int
    cd lcwfv_enmi
"""

import pysftp
import sys
import os
import numpy as np
sys.path.append(r'/home/patrikb/wavy/wavy')
from satmod import credentials_from_netrc

# obtain credentials
host = 'acquisition.ecmwf.int'
user, pw = credentials_from_netrc(host)
destination = 'lcwfv_enmi'

# obtain filename
path = '/lustre/storeA/project/fou/om/LC-WVF/'
filenamelst = [
                'wave_enmi_2019121606_prod_fc.grib2',
                'wave_enmi_2019121718_prod_fc.grib2',
                ]

for filename in filenamelst:
    with pysftp.Connection(host=host, username=user, password=pw) as sftp:
        fullname = path + filename
        print('transfer: ' + fullname)
        sftp.cwd(destination)
        sftp.put(fullname)
