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
    usrhome = os.path.expanduser("~")
    cred_path = os.path.join(usrhome, "credentials.txt")
    with open(cred_path, 'r') as f:
        for line in f:
            items = line.split(',')
            if remoteHostName in items[0]:
                user, pw = items[1].split('=')[1], \
                        items[2].split('=')[1][0:-1]
    return user, pw

def get_credentials(remoteHostName=None):
    usrhome = os.path.expanduser("~")

    # Windows sometimes uses _netrc instead of .netrc
    netrc_paths = [
        os.path.join(usrhome, ".netrc"),
        os.path.join(usrhome, "_netrc")
    ]

    for cred_path in netrc_paths:
        if os.path.isfile(cred_path):
            try:
                print("Obtaining credentials from netrc")
                return credentials_from_netrc(remoteHostName=remoteHostName)
            except AttributeError:
                print("Credentials in netrc file not registered")
                break
            except Exception as e:
                print(f"Failed reading netrc: {e}")
                break

    # Fallback to credentials.txt
    try:
        print("Trying local file credentials.txt")
        return credentials_from_txt(remoteHostName=remoteHostName)
    except Exception as e:
        raise RuntimeError("Credentials could not be obtained") from e
