#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
obtain credentials for retrieving remote file
'''
# --- import libraries ------------------------------------------------#
import os.path

def credentials_from_netrc(remoteHostName=None):
    import netrc
    if remoteHostName is None:
        remoteHostName = "nrt.cmems-du.eu"
    netrc = netrc.netrc()
    user = netrc.authenticators(remoteHostName)[0]
    pw = netrc.authenticators(remoteHostName)[2]
    return user, pw

def credentials_from_txt(remoteHostName=None):
    print("Try local file credentials.txt")
    if remoteHostName is None:
        remoteHostName = "nrt.cmems-du.eu"
    # get user home path
    usrhome = os.getenv("HOME")
    with open(usrhome + "/credentials.txt", 'r') as f:
        for line in f:
            items = line.split(',')
            if remoteHostName in items[0]:
                user, pw = items[1].split('=')[1], \
                        items[2].split('=')[1][0:-1]
    return user, pw

def get_credentials(remoteHostName=None):
    # get credentials from .netrc
    usrhome = os.getenv("HOME")
    if os.path.isfile(usrhome + "/.netrc"):
        try:
            print('Obtaining credentials from .netrc')
            user, pw = \
                credentials_from_netrc(remoteHostName=remoteHostName)
            return user, pw
        except AttributeError:
            print("Credentials in netrc file not registered")
            print("Try local file credentials.txt")
            print("File must contain:")
            print("user='username'")
            print("pw='yourpassword'")
            except_key = 1
    elif (os.path.isfile(usrhome + "/.netrc") is False
    or except_key == 1):
        try:
            user, pw = credentials_from_txt()
            return user, pw
        except:
            print("Credentials could not be obtained")

