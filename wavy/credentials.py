#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
'''
obtain credentials for retrieving remote file
'''
# --- import libraries ------------------------------------------------#
import os.path
import netrc

def credentials_from_netrc(remoteHostName=None):
    import netrc
    if remoteHostName is None:
        remoteHostName = "nrt.cmems-du.eu"
    # get user home path
    usrhome = os.getenv("HOME")
    netrc = netrc.netrc()
    user = netrc.authenticators(remoteHostName)[0]
    pw = netrc.authenticators(remoteHostName)[2]
    return user, pw

def credentials_from_txt():
    print("try local file credentials.txt")
    # get user home path
    usrhome = os.getenv("HOME")
    my_dict = {}
    with open(usrhome + "/credentials.txt", 'r') as f:
        for line in f:
            items = line.split('=')
            key, values = items[0], items[1]
            my_dict[key] = values
    # 1:-2 to remove the linebreak and quotation marks \n
    user = my_dict['user'][1:-2]
    pw = my_dict['pw'][1:-2]
    return user, pw

def get_credentials():
    # get credentials from .netrc
    usrhome = os.getenv("HOME")
    if os.path.isfile(usrhome + "/.netrc"):
        try:
            print ('Attempt to obtain credentials from .netrc')
            user, pw = credentials_from_netrc()
            return user, pw
        except AttributeError:
            print ("std copernicus user in netrc file not registered")
            print ("try local file credentials.txt")
            print ("file must contain:")
            print ("user='username'")
            print ("pw='yourpassword'")
            except_key = 1
    elif (os.path.isfile(usrhome + "/.netrc") == False
    or except_key == 1):
        try:
            user, pw = credentials_from_txt()
            return user, pw
        except:
            print("Credentials could not be obtained")

