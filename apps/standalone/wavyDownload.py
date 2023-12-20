#!/usr/bin/env python3
"""
download satellite data from Copernicus
"""
# --- imports -------------------------------------------------------- #
# standard library imports
import click
from wavy.satellite_module import satellite_class as sc
from datetime import datetime, timedelta
from pathlib import Path
# -------------------------------------------------------------------- #

# make sure that if name is name="all" then download all names for given nID

@click.command(context_settings={"ignore_unknown_options": True})
@click.option('--sd', type=str, default=datetime.now()-timedelta(hours=24),
        help='starting date and time of your query e.g.: 2023-10-1 00')
@click.option('--ed', type=str, default=datetime.now(),
        help='ending date and time of your query e.g.: 2023-10-10 00')
@click.option('--nID', type=str,
        help='nID as specified in satellite_cfg.yaml')
@click.option('--name', type=str,
        help='name as specified in satellite_cfg.yaml,\
        if name equals "all", all names from chosen nID are considered')
@click.option('--nproc', type=int, default=4,
        help='chosen number of simultaneous processes')
@click.argument('--path',
        help='target path for your download',
        type=click.Path(
            exists=True,
            readable=True,
            path_type=Path,
        ),)
@click.option('--search_str', type=str,
        help='identifyer string to search for in remote directory')
def main():
    """
    Wrapper for command line use of the wavy downloading functions.\n

    Here are some examples of supported files...
    The following most common missions are available from the CMEMS webpage:
            \ncmems_L3_NRT:\
            \n s3a - Sentinel-3A\
            \n s3b - Sentinel-3B\
            \n j3 - Jason-3 (reference mission)\
            \n c2 - Cryosat-2\
            \n al - SARAL/AltiKa\
            \n cfo - CFOSAT\
            \n h2b - HaiYang-2B\
            \n s6a - Sentinel-6A Michael Freilich\
            \n
    The following most common missions are available from CCIv1:
            \n
            \ncci_L2P:\
            \n j1 - Jason-1\
            \n j2 - Jason-2\
            \n j3 - Jason-3\
            \n c2 - Cryosat-2\
            \n envisat - Envisat\
            \n ers1 - European Remote-Sensing Satellite-1\
            \n ers2 - European Remote-Sensing Satellite-2\
            \n topex - TOPEX/Poseidon\
            \n al - SARAL/AltiKa\
            \n gfo - GEOSAT Follow-On\
            \n\
            \ncci_L3:\
            \n multi - multimission product 1991-2018\
            \n\

    Note that these examples are not exclusive,\n
    almost any other source could be added and exploited.
    """

    # settings
    now = datetime.now()
    if args.sat == 'all':
        namelst = list(satellite_dict[]['name'].keys())
    else:
        namelst = [name]

    if args.sd is None:
        sdate = datetime(now.year, now.month, now.day,
                         now.hour) - timedelta(hours=24)
    else:
        sdate = datetime(int(args.sd[0:4]), int(args.sd[4:6]),
                         int(args.sd[6:8]), int(args.sd[8:10]))

    if args.ed is None:
        edate = datetime(now.year, now.month, now.day, now.hour, now.minute)
    else:
        edate = datetime(int(args.ed[0:4]), int(args.ed[4:6]),
                         int(args.ed[6:8]), int(args.ed[8:10]))

    if args.nproc is None:
        args.nproc = 1

    if args.product is None:
        args.product = 'cmems_L3_NRT'

    for n in namelst:


if __name__ == "__main__":
    main()
