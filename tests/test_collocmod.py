import pytest

from wavy.satellite_module import satellite_class as sc
from wavy.collocation_module import collocation_class as cc
from wavy.insitu_module import insitu_class as ic
from wavy.insitu_module import poi_class as pc

# include possibility for collocating different variable
# varalias = 'Hs', 'U', aso...

def test_sat_collocation_and_validation(test_data, tmpdir):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = 's3a'
    varalias = 'Hs'
    twin = 30
    nID = 'cmems_L3_NRT'
    model = 'ww3_4km'
    # init satellite_object and check for polygon region
    sco = sc(sd=sd, ed=ed, nID=nID, name=name,
             varalias=varalias, twin=twin)
    # read data
    sco = sco.populate(reader='read_local_ncfiles',
                       path=str(test_data/"L3/s3a"))
    # crop to region
    sco = sco.crop_to_region(model)

    # collocate
    cco = cc(oco=sco, model=model, leadtime='best', distlim=6).populate()
    assert len(vars(cco).keys()) == 19
    assert len(cco.vars.keys()) == 9

    # validate


def test_insitu_collocation_and_validation(test_data, tmpdir):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    varalias = 'Hs'
    twin = 30
    model = 'ww3_4km'
    nID = 'D_Breisundet_wave'
    name = 'wavescan'

    # init insitu_object and check for polygon region
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias,
             name=name, twin=twin)

    # read data
    ico = ico.populate()

    # collocate
    cco = cc(oco=ico, model=model, leadtime='best', distlim=6).populate()
    assert len(vars(cco).keys()) == 19
    assert len(cco.vars.keys()) == 9

    # validate


def test_poi_collocation():
    # define poi dictionary for track
    dt = ["2023-7-1", "2023-7-2", "2023-7-3"]
    lats = [56.5, 59.3, 64.3]
    lons = [3.5, 1.8, 4.2]
    poi_dict = {'time': dt, 'lons': lons, 'lats': lats}

    # init poi_class
    pco = pc(poi_dict)

    # collocate
    cco = cc(oco=pco, model='ww3_4km', leadtime='best').populate()
    assert len(vars(cco).keys()) == 19
    assert len(cco.vars.keys()) == 9


#    # write to nc
#    cco.write_to_nc(pathtofile=tmpdir.join('test.nc'))
#    # test validation
#    cco.validate_collocated_values()
#
#def test_insitu_collocation_and_validation():
#    sd = "2021-8-2 01"
#    ed = "2021-8-2 03"
#    nID = 'D_Breisundet_wave'
#    sensor = 'wavescan'
#    ico = ic(nID,sd,ed,varalias=varalias,stwin=1,date_incr=1,sensor=sensor)
#    # collocate
#    cco = cc(model='mwam4',obs_obj_in=ico,distlim=6,
#             leadtime='best',date_incr=1)
#    # test validation
#    cco.validate_collocated_values()


def test_collocate_observations(test_data):
    from wavy.collocation_module import collocate_observations
    sd = '2023-07-04'
    ed = '2023-07-05'
    ico = ic(sd=sd, ed=ed, nID='MO_Draugen_monthly', name = 'Draugen').\
             populate(path=str(test_data/"insitu/monthly/Draugen/"))
    sco = sc(sd=sd, ed=ed, nID='cmems_L3_NRT', name='s3a').\
             populate(path=str(test_data/"L3/s3a"))

    ico_colloc, sco_colloc = collocate_observations(ico, sco)
    print(ico_colloc)
    print(len(ico_colloc.vars.keys()))
    assert len(ico_colloc.vars.keys()) == 3
    assert len(ico_colloc.vars.time.values) > 0
    print(sco_colloc)
    print(len(sco_colloc.vars.keys()))
    assert len(sco_colloc.vars.keys()) == 4
    assert len(sco_colloc.vars.time.values) > 0
