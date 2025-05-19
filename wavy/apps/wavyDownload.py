#!/usr/bin/env python3

"""
download satellite data from Copernicus
"""
# --- imports -------------------------------------------------------- #
# standard library imports
import click
from wavy.satellite_module import satellite_class as sc
from datetime import datetime, timedelta
import time
from wavy.utils import parse_date
from pathlib import Path
from wavy.wconfig import load_or_default
# -------------------------------------------------------------------- #

@click.command(context_settings={"ignore_unknown_options": True})
@click.option('--sd', type=str, default=None,
        help='starting date and time of your query e.g.: 2023-10-1T00')
@click.option('--ed', type=str, default=None,
        help='ending date and time of your query e.g.: 2023-10-10T00')
@click.option('--nID', type=str, default='cmems_L3_NRT',
        help='nID as specified in satellite_cfg.yaml')
@click.option('--name', type=str, default=None,
        help='name as specified in satellite_cfg.yaml, if name equals "all", all names from chosen nID are considered')
@click.option('--nproc', type=int, default=None,
        help='chosen number of simultaneous processes (only valid for FTP downloads)')
@click.option('--path', type=str, default=None, help='custom specified target path')

def main(sd, ed, nid, name, path, nproc):
    """
    Wrapper for command line use of the wavy downloading functions.\n

    The following most common near real time (NRT) missions are available from the CMEMS webpage.\n
            \ncmems_L3_NRT:\n
            \n s3a - Sentinel-3A\n
            \n s3b - Sentinel-3B\n
            \n j3 - Jason-3 (deprecated reference mission)\n
            \n c2 - Cryosat-2\n
            \n al - SARAL/AltiKa\n
            \n cfo - CFOSAT\n
            \n h2b - HaiYang-2B\n
            \n s6a - Sentinel-6A Michael Freilich (reference mission)\n
            \n swon - SWOT nadir altimeter\n
            \n

    Note that these examples are not exclusive,\n
    almost any other source could be added and exploited.
    """

    # read yaml config files:
    satellite_dict = load_or_default('satellite_cfg.yaml')

    # settings
    now = datetime.now()

    if sd is None:
        sdate = now-timedelta(hours=24)
    else:
        sdate = parse_date(sd)

    if ed is None:
        edate = now
    else:
        edate = parse_date(ed)

    if name is None:
        namelst = [list(satellite_dict[nid]['name'].keys())[0]]
    elif name == 'all':
        namelst = list(satellite_dict[nid]['name'].keys())
    else:
        namelst = [name]

    if nproc is None:
        nproc = 1

    print(sdate)
    print(edate)
    print(nid)

    for n in namelst:
        print(n)

    for name in namelst:

        print("Attempting to download data for:", name)
        print("Time period:", str(sdate), "to", str(edate))

        start_time = time.time()

        sco = sc(sd=sdate, ed=edate,
                 nID=nid, name=name)
        sco.download(path=path, nproc=nproc)

        time1 = time.time() - start_time
        print("Time used for collecting data: ", time1, " seconds")

if __name__ == "__main__":
    main()
