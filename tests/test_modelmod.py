from wavy.modelmod import model_class as mc

def test_model_class_init():
    #get_model
    mco = mc(nID='ww3_4km', sd="2023-6-1", ed="2023-6-1 01")
    assert mco.__class__.__name__ == 'model_class'

def test_ww3_4km_reader():
    #get_model
    mco = mc(nID='ww3_4km', sd="2023-6-1", ed="2023-6-1 01")
    mco.populate()
    #assert len(vars(mco).keys()) >= 9
    #assert len(mco.vars.keys()) >= 7
