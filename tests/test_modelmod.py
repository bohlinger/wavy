from datetime import datetime

def test_model_class():
    from datetime import datetime
    from wavy.modelmod import model_class as mc
    model = "mwam4"
    fc_date = datetime(2020,10,1,12)
    #get_model
    mc_obj = mc(model=model,fc_date=fc_date)
    assert mc_obj.__class__.__name__ == 'model_class'
    assert len(vars(mc_obj).keys()) >= 9
    assert len(mc_obj.vars.keys()) >= 7
    assert not 'error' in vars(mc_obj).keys()
