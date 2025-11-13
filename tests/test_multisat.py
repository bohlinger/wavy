import pytest
from wavy import ms

def test_multisat(test_data):
    sd = "2022-2-1 12"
    ed = "2022-2-1 12"
    name = ['s3a','s3b']
    varalias = 'Hs'

    # init multisat_object
    mso = ms(sd=sd,
         ed=ed, 	
         name=name,
         varalias = varalias, 
         path = [str(test_data/"L3/s3a"),
                 str(test_data/"L3/s3b")])
    # read data
    assert mso.__class__.__name__ == 'multisat_class'
    # compare number of available variables
    vlst = list(vars(mso).keys())
    assert len(vlst) == 17
    # compare number of available functions
    dlst = dir(mso)
    flst = [n for n in dlst if n not in vlst if '__' not in n]
    assert len(flst) >= 27
    assert type(mso.vars == 'xarray.core.dataset.Dataset')
    assert not 'error' in vars(mso).keys()
