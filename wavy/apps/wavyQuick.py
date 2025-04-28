#!/usr/bin/env python

# import standard libraries
import numpy as np
from datetime import datetime, timedelta
from wavy.satmod import satellite_class as sa
import argparse
from argparse import RawTextHelpFormatter
import os
import yaml

# own imports
from wavy.modelmod import model_class as mc
from wavy.validationmod import validate, disp_validation
from wavy.quicklookmod import comp_fig, plot_sat
from wavy.collocmod import collocation_class as coll
from wavy.wconfig import load_or_default


def main():
    variable_info = load_or_default('variable_info.yaml')
    model_dict = load_or_default('model_specs.yaml')

# parser
    parser = argparse.ArgumentParser(description="""
    Check availability of satellite SWH data. Example:
    ./wavyQuick.py -sat s3a -reg mwam4 -mod mwam4 -sd 2020100112 -lt 0 -twin 30 --col --show
    ./wavyQuick.py -sat s3a -reg mwam4 -sd 2020100100 -ed 2020101000 --show
    ./wavyQuick.py -sat s3a -reg mwam4 -sd 2020100100 -ed 2020101000 -dump ./test.nc
        """,
                                        formatter_class=RawTextHelpFormatter)
    parser.add_argument("-reg", metavar='region', help="region to check")
    parser.add_argument("-sat",
                        metavar='satellite',
                        help="satellite mission, currently available\n \
            \ncmems_L3_NRT:\
            \n s3a - Sentinel-3A\
            \n s3b - Sentinel-3B\
            \n j3 - Jason-3 (reference mission)\
            \n c2 - Cryosat-2\
            \n al - SARAL/AltiKa\
            \n cfo - CFOSAT\
            \n h2b - HaiYang-2B\
            \n\
            \ncmems_L3_s6a:\
            \n s6a - Sentinel-6A Michael Freilich\
            \n\
            \neumetsat_L2:\
            \n s3a - Sentinel-3A\
            \n s3b - Sentinel-3B\
            \n\
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
            \ncfo_swim_L2P:\
            \n cfo - CFOSAT\
            \n\
            \nother options are:\
            \n all - for all cmems satellites\
            \n list - a list of chosen satellites using -l\n \
            \n")
    parser.add_argument("-product",
                        metavar='product',
                        help="available products as specified in *_specs.yaml")
    parser.add_argument('-l',
                        metavar='satellite list',
                        help='delimited list input for sats',
                        type=str)
    parser.add_argument("-sd",
                        metavar='startdate',
                        help="start date of time period to check")
    parser.add_argument("-ed",
                        metavar='enddate',
                        help="end date of time period to check")
    parser.add_argument("-mod", metavar='model', help="chosen wave model")
    parser.add_argument("-var",
                        metavar='varalias',
                        help="alias for chosen variable")
    parser.add_argument("-lt",
                        metavar='lead time',
                        type=int,
                        help="lead time from initialization")
    parser.add_argument("-twin",
                        metavar='time window',
                        type=int,
                        help="time window for collocation")
    parser.add_argument("-dist",
                        metavar='distance limit',
                        type=int,
                        help="distance limit for collocation")
    parser.add_argument("--col",
                        metavar="collocation",
                        help="collocation",
                        action='store_const',
                        const=True)
    parser.add_argument("--show",
                        help="show figure",
                        action='store_const',
                        const=True)
    parser.add_argument("-savep",
                        metavar="savepath",
                        help="save figure to path")
    parser.add_argument("-dump",
                        metavar="outpath",
                        help="dump data to .nc-file")

    args = parser.parse_args()
    print("Parsed arguments: ", args)

    flatten = lambda l: [item for sublist in l for item in sublist]

# setup
    if args.var is None:
        args.var = 'Hs'

    sdate = datetime(int(args.sd[0:4]), int(args.sd[4:6]), int(args.sd[6:8]),
                        int(args.sd[8:10]))

    if args.product is None:
        args.product = 'cmems_L3_NRT'

    if args.twin is None:
        args.twin = 30

    if args.dist is None:
        args.dist = 10

    if args.ed is None:
        edate = sdate
    else:
        edate = datetime(int(args.ed[0:4]), int(args.ed[4:6]),
                            int(args.ed[6:8]), int(args.ed[8:10]))
        twin = 0

    if (args.reg is None and args.mod is not None):
        args.reg = args.mod

# get data
    if args.sat == 'all':
        satlist = ['s3a', 's3b', 'j3', 'c2', 'al', 'cfo', 'h2b']
        lats = []
        lons = []
        var = []
        time = []
        dtime = []
        sats = []
        satnamelst = []
        for sat in satlist:
            try:
                sa_obj_tmp = sa(sdate=sdate,
                                mission=sat,
                                edate=edate,
                                twin=args.twin,
                                region=args.reg,
                                varalias=args.var,
                                product=args.product)
                if ('vars' in vars(sa_obj_tmp).keys()
                        and len(sa_obj_tmp.vars['time']) > 0):
                    sa_obj = sa_obj_tmp
                    lats.append(sa_obj.vars['latitude'])
                    lons.append(sa_obj.vars['longitude'])
                    var.append(sa_obj.vars[variable_info[args.var]\
                                                ['standard_name']])
                    time.append(sa_obj.vars['time'])
                    dtime.append(sa_obj.vars['datetime'])
                    sats.append(sat)
                    satnamelst.append([sat] * len(sa_obj.vars['time']))
            except Exception as e:
                print(e)
                print(sat + ' not available')
                pass
        lats = flatten(lats)
        lons = flatten(lons)
        var = flatten(var)
        time = flatten(time)
        dtime = flatten(dtime)
        satnames = flatten(satnamelst)
        sa_obj.vars['latitude'] = np.array(lats)
        sa_obj.vars['longitude'] = np.array(lons)
        sa_obj.vars[variable_info[args.var]['standard_name']] = np.array(var)
        sa_obj.vars['time'] = time
        sa_obj.vars['datetime'] = dtime
        sa_obj.region = args.reg
        sa_obj.mission = str(sats)
        sa_obj.mission_ts = satnames
    elif args.sat == 'list':
        satlist = args.l.split(',')
        lats = []
        lons = []
        var = []
        time = []
        sats = []
        satnamelst = []
        for sat in satlist:
            try:
                sa_obj_tmp = sa(sdate=sdate,
                                mission=sat,
                                edate=edate,
                                twin=args.twin,
                                region=args.reg,
                                varalias=args.var,
                                product=args.product)
                if ('vars' in vars(sa_obj_tmp).keys()
                        and len(sa_obj_tmp.vars['time']) > 0):
                    sa_obj = sa_obj_tmp
                    lats.append(sa_obj.vars['latitude'])
                    lons.append(sa_obj.vars['longitude'])
                    var.append(sa_obj.vars[variable_info[args.var]\
                                                ['standard_name']])
                    time.append(sa_obj.vars['time'])
                    sats.append(sat)
                    satnamelst.append([sat] * len(sa_obj.vars['time']))
            except Exception as e:
                print(e)
                print(sat + ' not available')
                pass
        lats = flatten(lats)
        lons = flatten(lons)
        var = flatten(var)
        time = flatten(time)
        satnames = flatten(satnamelst)
        sa_obj.vars['time'] = time
        sa_obj.vars['latitude'] = np.array(lats)
        sa_obj.vars['longitude'] = np.array(lons)
        sa_obj.vars[variable_info[args.var]['standard_name']] = np.array(var)
        sa_obj.region = args.reg
        sa_obj.mission = str(sats)
        sa_obj.mission_ts = satnames
    else:
        sa_obj = sa(sdate=sdate,
                    mission=args.sat,
                    edate=edate,
                    twin=args.twin,
                    region=args.reg,
                    varalias=args.var,
                    product=args.product)

# plot
    if (args.mod is None and sa_obj.region not in model_dict):
        plot_sat(sa_obj, savepath=args.savep, showfig=args.show)
    elif (args.mod is None and sa_obj.region in model_dict):
        print('Chosen region is a specified model domain')
        mc_obj = mc(model=sa_obj.region,
                    fc_date=model_dict[sa_obj.region]['grid_date'],
                    varalias=args.var)
        plot_sat(sa_obj, mc_obj=mc_obj, savepath=args.savep,
                showfig=args.show)
    elif (args.mod is not None and args.col is True and sdate==edate):
        # get model collocated values
        mc_obj = mc(model=args.mod,
                    fc_date=edate,
                    leadtime=args.lt,
                    varalias=args.var)
        #collocation
        coll_obj = coll(mc_obj_in=mc_obj,obs_obj_in=sa_obj,
                        distlim=args.dist,leadtime=args.lt)
        valid_dict = validate(coll_obj.vars)
        disp_validation(valid_dict)
        comp_fig(sa_obj=sa_obj,
                    mc_obj=mc_obj,
                    coll_obj=coll_obj,
                    savepath=args.savep,
                    showfig=args.show)
    elif (args.mod is not None and args.col is False):
        # get model collocated values
        mc_obj = mc(model=sa_obj.region,
                    fc_date=edate,
                    leadtime=args.lt,
                    varalias=args.var)
        results_dict = {'valid_date':[edate],
                        'date_matches':[edate-timedelta(minutes=args.twin),
                                        edate+timedelta(minutes=args.twin)],
                        'model_lons_matches':sa_obj.vars['longitude'],
                        'model_lats_matches':sa_obj.vars['latitude'],
                        'sat_matches':sa_obj.vars[variable_info[args.var]\
                                                        ['standard_name']]}
        comp_fig(args.mod,
                    sa_obj,
                    mc_obj.vars[variable_info[args.var]['standard_name']],
                    mc_obj.vars['longitude'],
                    mc_obj.vars['latitude'],
                    results_dict,
                    variable_info[args.var]['standard_name'],
                    savepath=args.savep,
                    showfig=args.show)
    elif (args.mod is not None and args.col is True and args.sat is not None):
        #collocation
        date_incr = model_dict[args.mod].get('date_incr',1)
        coll_obj = coll(model=args.mod, obs_obj_in=sa_obj, distlim=args.dist, leadtime = args.lt, date_incr = date_incr)
        valid_dict = coll_obj.validate_collocated_values()
        if args.show is True:
            print("Figure not shown! Request is ambiguous for multiple model time steps.")

# dump to .ncfile
    if args.dump is not None:
        if args.col is not None:
            coll_obj.write_to_nc(pathtofile=args.dump)
        else:
            sa_obj.write_to_nc(pathtofile=args.dump)

if __name__ == "__main__":
    main()
