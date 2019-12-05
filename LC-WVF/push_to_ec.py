#!/usr/bin/env python3
import pysftp
import sys
import os
sys.path.append(r'/home/patrikb/wavy/wavy')
from satmod import credentials_from_netrc

# obtain credentials
host = 'acquisition.ecmwf.int'
user, pw = credentials_from_netrc(host)

# obtain filename
filelst = np.sort(os.listdir('/lustre/storeA/project/fou/om/LC-WVF'))
filename = filelst[-1]

with pysftp.Connection(host=host, username=user, password=pw) as sftp:
    sftp.cwd('lcwf_enmi')   # alternativt sftp.cd
    sftp.put(filename)
