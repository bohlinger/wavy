#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------#
"""
This module comprises functions to collect remote satellite data.
"""

# --- import libraries ------------------------------------------------#
# standard library imports
import sys
import numpy as np
from datetime import datetime, timedelta
import os
from copy import deepcopy
import time
from urllib.request import urlretrieve, urlcleanup
from urllib.parse import quote
from ftplib import FTP
from dateutil.relativedelta import relativedelta
from joblib import Parallel, delayed
import logging
import copernicusmarine as cmc

# own imports
from wavy.utils import sort_files
from wavy.utils import make_pathtofile, make_subdict
from wavy.utils import date_dispatcher
from wavy.credentials import get_credentials
from wavy.wconfig import load_or_default

# ---------------------------------------------------------------------#

# read yaml config files:
satellite_dict = load_or_default("satellite_cfg.yaml")


# --- def functions ---------------------------------------------------#
def reorganize_nc_files(path_local):
    """Move .nc file from nested directory to parent and rename to dataset name (without coords)."""
    try:
        for item in os.listdir(path_local):
            item_path = os.path.join(path_local, item)
            if os.path.isdir(item_path):
                for file in os.listdir(item_path):
                    if file.endswith(".nc"):
                        src = os.path.join(item_path, file)
                        # Extract base name by removing coordinates and depth suffix
                        # e.g., cmems_obs-wave_glo_phy-swh_nrt_s3a-l3_PT1S_VAVH-WIND_SPEED_28.83E-52.17E_32.00S-12.41S_0.00m
                        # becomes: cmems_obs-wave_glo_phy-swh_nrt_s3a-l3_PT1S_VAVH-WIND_SPEED
                        base_name = item.rsplit("_", 1)[
                            0
                        ]  # Remove last part after last underscore
                        dst = os.path.join(path_local, f"{base_name}.nc")
                        os.rename(src, dst)
                        logger.info(f"Moved to {base_name}.nc")
                        print(f"✓ Moved: {base_name}.nc")
                os.rmdir(item_path)
    except Exception as e:
        logger.error(f"Failed to reorganize: {e}")
        print(f"Reorganization failed: {e}")


def tmploop_get_remote_files(
    i: int,
    matching: str,
    user: str,
    pw: str,
    server: str,
    remote_path: str,
    path_local: str,
    **kwargs,
):
    """
    Function to download files using ftp. Tries 10 times before failing.
    """
    logger = logging.getLogger(__name__)
    log_level = str(kwargs.get("logging", "WARNING").upper())
    logger.setLevel(getattr(logging, log_level, logging.WARNING))

    logger.info("File: " + str(matching[i]))
    logger.info("src path: " + str(remote_path))
    pw = quote(pw)  # to escape special characters
    dlstr = "ftp://" + user + ":" + pw + "@" + server + remote_path + matching[i]
    for attempt in range(10):
        logger.info(str(attempt) + "Attempt to download data: ")
        try:
            logger.info("Downloading file")
            urlretrieve(dlstr, os.path.join(path_local, matching[i]))
            urlcleanup()
        except Exception as e:
            logger.warning("Exception in tmploop_get_remote_files:")
            logger.warning(e)
            logger.warning("Waiting for 10 sec and retry")
            time.sleep(10)
        else:
            break
    else:
        logger.critical(
            "An error was raised and I " + "failed to fix problem myself :("
        )
        logger.critical("Exit program")
        sys.exit()


def get_remote_files_ftp(**kwargs):
    """
    Download swath files from CMEMS and store them at defined
    location. Time stamps in file name stand for:

    from, to, creation
    """
    logger = logging.getLogger(__name__)
    log_level = str(kwargs.get("logging", "WARNING").upper())
    logger.setLevel(getattr(logging, log_level, logging.WARNING))

    product = kwargs.get("nID")
    sdate = kwargs.get("sd")
    edate = kwargs.get("ed")
    twin = int(np.max([kwargs.get("twin", 30), 30]))
    nproc = kwargs.get("nproc", 1)
    name = kwargs.get("name", "s3a")
    # dict_for_sub = kwargs.get('dict_for_sub')
    dict_for_sub = kwargs

    # define path
    path = kwargs.get("path", None)
    # check if search str template
    file_search_template = kwargs.get("search_str")
    if file_search_template is None:
        file_search_template = satellite_dict[product]["download"]["ftp"].get(
            "search_str", "%Y%m%dT%H"
        )

    # credentials
    server = satellite_dict[product]["download"]["ftp"]["server"]
    user, pw = get_credentials(remoteHostName=server)
    tmpdate = deepcopy(sdate)
    path_template_src = satellite_dict[product]["download"]["ftp"]["src_tmplt"]
    strsublst_src = satellite_dict[product]["download"]["ftp"]["strsub"]
    subdict_src = make_subdict(strsublst_src, class_object_dict=dict_for_sub)
    while tmpdate <= edate:
        try:
            # create remote path
            path_remote = make_pathtofile(
                path_template_src, strsublst_src, subdict_src, date=tmpdate
            )

            if path is None:
                # create local path
                path_template_dst = satellite_dict[product]["download"]["ftp"][
                    "trgt_tmplt"
                ]
                strsublst_dst = satellite_dict[product]["download"]["ftp"]["strsub"]
                subdict_dst = make_subdict(
                    strsublst_dst, class_object_dict=dict_for_sub
                )
                path_local = make_pathtofile(
                    path_template_dst, strsublst_dst, subdict_dst, date=tmpdate
                )
            else:
                path_local = path

            print("# ----- ")
            print("Chosen source: ")
            print(name + " values from " + product + ": " + server)
            print(path_remote)
            print("# ----- ")
            # get list of accessable files
            ftp = FTP(server)
            ftp.login(user, pw)
            ftp.cwd(path_remote)
            content = FTP.nlst(ftp)

            # choose files according to sdate/edate
            tmplst = []
            tmpdate_new = tmpdate - timedelta(minutes=twin)
            tmpdate_end = edate + timedelta(minutes=twin)
            while tmpdate_new <= tmpdate_end:
                matchingtmp = [
                    s
                    for s in content
                    if tmpdate_new.strftime(file_search_template) in s
                ]
                tmplst = tmplst + matchingtmp
                tmpdate_new = tmpdate_new + timedelta(minutes=twin)
            matching = np.unique(tmplst)
            print(matching)

            # check if download path_local exists if not create
            if not os.path.exists(path_local):
                os.makedirs(path_local, exist_ok=True)

            # Download matching files
            print("Downloading " + str(len(matching)) + " files: .... \n")
            print("Used number of possible simultaneous downloads " + str(nproc) + "!")
            Parallel(n_jobs=nproc)(
                delayed(tmploop_get_remote_files)(
                    i, matching, user, pw, server, path_remote, path_local, **kwargs
                )
                for i in range(len(matching))
            )
        except Exception as e:
            logger.warning("Exception was raised during downloading in")
            logger.warning("get_remote_files_ftp")
            logger.exception(e)

        # update time
        path_date_incr_unit = satellite_dict[product]["download"]["ftp"].get(
            "path_date_incr_unit", "m"
        )
        path_date_incr = satellite_dict[product]["download"]["ftp"].get(
            "path_date_incr", 1
        )
        tmpdate = date_dispatcher(tmpdate, path_date_incr_unit, path_date_incr)
        print("####################################")
        print(path_local)
        print("####################################")

    print("Files downloaded to: \n", path_local)


def get_remote_files_copernicusmarine(**kwargs):
    """
    Download swath files from CMEMS using copernicusmarine parckage
    and store them at defined location. Time stamps in file name stand for:

    from, to, creation
    """
    logger = logging.getLogger(__name__)
    log_level = str(kwargs.get("logging", "WARNING").upper())
    logger.setLevel(getattr(logging, log_level, logging.WARNING))

    product = kwargs.get("nID")
    sdate = kwargs.get("sd")
    edate = kwargs.get("ed")
    name = kwargs.get("name", "s3a")

    # if CMEMS credentials are defined in environment other options
    # are overwritten
    if "COPERNICUSMARINE_SERVICE_USERNAME" in os.environ:
        username = os.getenv("COPERNICUSMARINE_SERVICE_USERNAME")
    else:
        username = kwargs.get("COPERNICUSMARINE_SERVICE_USERNAME")
    if "COPERNICUSMARINE_SERVICE_PASSWORD" in os.environ:
        password = os.getenv("COPERNICUSMARINE_SERVICE_PASSWORD")
    else:
        password = kwargs.get("COPERNICUSMARINE_SERVICE_PASSWORD")

    if username is None or password is None:
        print("--> Please provide complete credentials!")

    dict_for_sub = kwargs

    # define path
    path = kwargs.get("path", None)

    # Get time increment
    time_incr = satellite_dict[product]["download"]["copernicus"].get("time_incr", "h")

    # Chose search template for time given time_incr
    if time_incr == "h":
        file_search_template = "%Y%m%dT%H"
    elif time_incr == "d":
        file_search_template = "%Y%m%dT"
    elif time_incr == "m":
        file_search_template = "%Y%m"
    print("Date search format:", file_search_template)

    # Get dataset_id
    dataset_id = satellite_dict[product]["download"]["copernicus"]["dataset_id"]
    strsublst_src = satellite_dict[product]["download"]["copernicus"]["strsub"]
    subdict_src = make_subdict(strsublst_src, class_object_dict=dict_for_sub)

    # replace name of the mission in dataset_id
    dataset_id = make_pathtofile(dataset_id, strsublst_src, subdict_src)

    # Initialize start date to match original files time increment
    tmpdate = deepcopy(sdate)
    tmpdate_end = deepcopy(edate)

    while tmpdate.hour % 3 != 0:
        tmpdate = tmpdate - timedelta(hours=1)

    try:

        print("# ----- ")
        print("Chosen source: ")
        print(name + " values from " + product + ": " + "copernicusmarine")
        print("# ----- ")

        while tmpdate <= tmpdate_end:

            if path is None:
                # create local path
                path_template_dst = satellite_dict[product]["download"]["copernicus"][
                    "trgt_tmplt"
                ]
                strsublst_dst = satellite_dict[product]["download"]["copernicus"][
                    "strsub"
                ]
                subdict_dst = make_subdict(
                    strsublst_dst, class_object_dict=dict_for_sub
                )
                path_local = make_pathtofile(
                    path_template_dst, strsublst_dst, subdict_dst, date=tmpdate
                )
            else:
                path_local = path

            print("* --------------")
            print("Downloading for date:", tmpdate)
            print("* --------------")

            # check if download path_local exists if not create
            if not os.path.exists(path_local):
                os.makedirs(path_local, exist_ok=True)

            # Create regexp filter
            tmpdate_str = tmpdate.strftime(file_search_template)
            regexp_tmp = "*{}*_*_*.nc".format(tmpdate_str)
            # Fetch data corresponding to tmp date
            try:
                cmc.get(
                    dataset_id=dataset_id,
                    filter=regexp_tmp,
                    no_directories=True,
                    output_directory=path_local,
                    username=username,
                    password=password,
                    overwrite=True,
                )
                # force_download=True,
                # overwrite_output_data=True,
                # no_metadata_cache=no_metadata_cache
            except Exception as e:
                logger.warning("Exception was raised during downloading in")
                logger.warning("get_remote_files_copernicusmarine")
                logger.exception(e)
                pass

            if time_incr == "h":
                tmpdate = tmpdate + timedelta(hours=3)
            elif time_incr == "d":
                tmpdate = tmpdate + timedelta(days=1)
            elif time_incr == "m":
                tmpdate = tmpdate + relativedelta(months=+1)

    except Exception as e:
        logger.exception(e)

    print("# -----------------------------------")
    print("Files downloaded to: \n", path_local)


def get_remote_files_copernicusmarine_subset(**kwargs):
    """
    Download swath files from CMEMS using copernicusmarine parckage
    and store them at defined location. Time stamps in file name stand for:

    from, to, creation

    If mx_lt, mn_lt, mx_ln, mn_ln are provided, cmc.subset() is used to
    spatially crop the download; otherwise cmc.get() fetches global data.
    """
    product = kwargs.get("nID")
    sdate = kwargs.get("sd")
    edate = kwargs.get("ed")
    name = kwargs.get("name", "s3a")

    mx_lt = kwargs.get("mx_lt", None)
    mn_lt = kwargs.get("mn_lt", None)
    mx_ln = kwargs.get("mx_ln", None)
    mn_ln = kwargs.get("mn_ln", None)
    use_subset = all(v is not None for v in (mx_lt, mn_lt, mx_ln, mn_ln))

    # if CMEMS credentials are defined in environment other options
    # are overwritten
    if "COPERNICUSMARINE_SERVICE_USERNAME" in os.environ:
        username = os.getenv("COPERNICUSMARINE_SERVICE_USERNAME")
    else:
        username = kwargs.get("COPERNICUSMARINE_SERVICE_USERNAME")
    if "COPERNICUSMARINE_SERVICE_PASSWORD" in os.environ:
        password = os.getenv("COPERNICUSMARINE_SERVICE_PASSWORD")
    else:
        password = kwargs.get("COPERNICUSMARINE_SERVICE_PASSWORD")

    if username is None or password is None:
        print("--> Please provide complete credentials!")

    dict_for_sub = kwargs

    # define path
    path = kwargs.get("path", None)

    # Get time increment
    time_incr = satellite_dict[product]["download"]["copernicus"].get("time_incr", "h")

    # Chose search template for time given time_incr
    if time_incr == "h":
        file_search_template = "%Y%m%dT%H"
    elif time_incr == "d":
        file_search_template = "%Y%m%dT"
    elif time_incr == "m":
        file_search_template = "%Y%m"
    print("Date search format:", file_search_template)

    # Get dataset_id
    dataset_id = satellite_dict[product]["download"]["copernicus"]["dataset_id"]
    strsublst_src = satellite_dict[product]["download"]["copernicus"]["strsub"]
    subdict_src = make_subdict(strsublst_src, class_object_dict=dict_for_sub)

    # replace name of the mission in dataset_id
    dataset_id = make_pathtofile(dataset_id, strsublst_src, subdict_src)

    # Initialize start date to match original files time increment
    tmpdate = deepcopy(sdate)
    tmpdate_end = deepcopy(edate)

    while tmpdate.hour % 3 != 0:
        tmpdate = tmpdate - timedelta(hours=1)

    try:

        print("# ----- ")
        print("Chosen source: ")
        print(name + " values from " + product + ": " + "copernicusmarine")
        print("# ----- ")

        if use_subset:
            # Download entire spatial block in one call (no time looping)
            if path is None:
                path_template_dst = satellite_dict[product]["download"]["copernicus"][
                    "trgt_tmplt"
                ]
                strsublst_dst = satellite_dict[product]["download"]["copernicus"][
                    "strsub"
                ]
                subdict_dst = make_subdict(
                    strsublst_dst, class_object_dict=dict_for_sub
                )
                path_local = make_pathtofile(
                    path_template_dst, strsublst_dst, subdict_dst, date=sdate
                )
            else:
                path_local = path

            if not os.path.exists(path_local):
                os.makedirs(path_local, exist_ok=True)

            print("Downloading spatial subset block:")
            print(f"  Lat: {mn_lt}° to {mx_lt}°, Lon: {mn_ln}° to {mx_ln}°")
            print(f"  Time: {sdate} to {edate}")

            try:
                cmc.subset(
                    dataset_id=dataset_id,
                    output_directory=path_local,
                    username=username,
                    password=password,
                    minimum_longitude=mn_ln,
                    maximum_longitude=mx_ln,
                    minimum_latitude=mn_lt,
                    maximum_latitude=mx_lt,
                    start_datetime=sdate.strftime("%Y-%m-%dT%H:%M:%S"),
                    end_datetime=edate.strftime("%Y-%m-%dT%H:%M:%S"),
                    file_format="netcdf",
                )
                reorganize_nc_files(path_local)
            except Exception as error:
                print(error)
                pass
        else:
            # Loop through time increments for global download
            while tmpdate <= tmpdate_end:

                if path is None:
                    # create local path
                    path_template_dst = satellite_dict[product]["download"][
                        "copernicus"
                    ]["trgt_tmplt"]
                    strsublst_dst = satellite_dict[product]["download"]["copernicus"][
                        "strsub"
                    ]
                    subdict_dst = make_subdict(
                        strsublst_dst, class_object_dict=dict_for_sub
                    )
                    path_local = make_pathtofile(
                        path_template_dst, strsublst_dst, subdict_dst, date=tmpdate
                    )
                else:
                    path_local = path

                print("* --------------")
                print("Downloading for date:", tmpdate)
                print("* --------------")

                # check if download path_local exists if not create
                if not os.path.exists(path_local):
                    os.makedirs(path_local, exist_ok=True)

                # Create regexp filter
                tmpdate_str = tmpdate.strftime(file_search_template)
                regexp_tmp = "*{}*_*_*.nc".format(tmpdate_str)
                # Fetch data corresponding to tmp date
                try:
                    cmc.get(
                        dataset_id=dataset_id,
                        filter=regexp_tmp,
                        no_directories=True,
                        output_directory=path_local,
                        username=username,
                        password=password,
                        overwrite=True,
                    )
                except Exception as error:
                    print(error)
                    pass

                if time_incr == "h":
                    tmpdate = tmpdate + timedelta(hours=3)
                elif time_incr == "d":
                    tmpdate = tmpdate + timedelta(days=1)
                elif time_incr == "m":
                    tmpdate = tmpdate + relativedelta(months=+1)

    except Exception as e:
        logger.exception(e)

    if use_subset:
        reorganize_nc_files(path_local)

    print("# -----------------------------------")
    print("Files downloaded to: \n", path_local)
