#import pytest
#
#from wavy.multiins import multiins_class as mi

#def test_consolidate_scos(test_data):
#    mio = mi(sdate = '2022-1-1',edate='2022-1-2',
#              tags=['E39','wave','pytest'],varalias='Hs')
#    assert len(list(vars(mio).keys())) == 14
#    assert len(mio.ocos) >= 1
