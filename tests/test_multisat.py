#import pytest
#
#from wavy.multisat import multisat_class as ms

#def test_multisat(test_data):
#    mso = ms(sdate="2020-11-1",edate="2020-11-3",region="NordicSeas",
#              mission=['s3a','s3b'], path_local=str(test_data/"L3"))
#    assert len(list(vars(mso).keys())) == 15
#    assert len(list(mso.ocos)) >= 1
