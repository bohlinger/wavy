import sys
from datetime import datetime
from wavy import satmod

date = datetime(2017,8,1,12)
timewin = 30
distlim = 6
region = 'ARCMFC'

def test_get_model_cont():
    from datetime import datetime
    from satmod import sentinel_altimeter as sa
    from modelmod import get_model
    from modelmod import collocate
    expname = "OPEWAVE"
    model = "ARCMFCnew"
    date = datetime(2017,10,1,12)
    #get_sa_obj
    sa_obj = sa(date,edate=date,timewin=timewin,
                region=region,mode="ARCMFC")
    #get_model
    model_Hs,model_lats,model_lons,model_time,model_time_dt = \
            get_model(simmode="cont",model=model,sdate=date,
                    edate=date,expname=expname)
    #collocation
    print ("model_Hs.shape: ")
    print (model_Hs.shape)
    results_dict = collocate(model,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,date,distlim=distlim)
    assert model_Hs.__class__.__name__ == 'ndarray'
    assert len(model_time_dt) > 0
    assert results_dict.__class__.__name__ == 'dict'

def test_get_model_ARCMFC(self):
    from datetime import datetime
    from satmod import sentinel_altimeter as sa
    from modelmod import get_model
    from modelmod import collocate
    import numpy as np
    model = "ARCMFC"
    fc_date = datetime(2017,10,1,12)
    init_date = datetime(2017,10,1,12)
    #get_sa_obj
    sa_obj = sa(fc_date,edate=fc_date,timewin=self.timewin,
                region=self.region,mode="ARCMFC")
    simmode="fc"
    #get_model
    model_Hs,model_lats,model_lons,model_time,model_time_dt = \
        get_model(simmode="fc",model=model,fc_date=fc_date,
        init_date=init_date)
    #collocation
    results_dict = collocate(model,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=self.distlim)
    self.assertEqual(model_Hs.__class__.__name__,'ndarray')
    self.assertEqual(np.isnan(np.nanmean(model_Hs)),False)
    self.assertGreater(len(model_time_dt),0)
    self.assertEqual(results_dict.__class__.__name__,'dict')

def test_get_model_mwam4(self):
    from datetime import datetime
    from satmod import sentinel_altimeter as sa
    from modelmod import get_model
    from modelmod import collocate
    import numpy as np
    model = "mwam4"
    region = "mwam4"
    fc_date = datetime(2017,10,1,12)
    init_date = datetime(2017,10,1,12)
    leadtime = 0
    #get_sa_obj
    sa_obj = sa(fc_date,edate=fc_date,timewin=self.timewin,
                region=region)
    simmode="fc"
    #get_model
    model_Hs,model_lats,model_lons,model_time,model_time_dt = \
        get_model(simmode="fc",model=model,fc_date=fc_date,
        leadtime=leadtime)
    #collocation
    results_dict = collocate(model,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=self.distlim)
    self.assertEqual(model_Hs.__class__.__name__,'ndarray')
    self.assertEqual(np.isnan(np.nanmean(model_Hs)),False)
    self.assertGreater(len(model_time_dt),0)
    self.assertEqual(results_dict.__class__.__name__,'dict')

def test_get_model_mwam8(self):
    from datetime import datetime
    from satmod import sentinel_altimeter as sa
    from modelmod import get_model
    from modelmod import collocate
    import numpy as np
    model = "mwam8"
    region = "mwam8"
    fc_date = datetime(2018,10,1,6)
    init_date = datetime(2018,10,1,6)
    leadtime = 0
    #get_sa_obj
    sa_obj = sa(fc_date,edate=fc_date,timewin=self.timewin,
                region=region)
    simmode="fc"
    #get_model
    model_Hs,model_lats,model_lons,model_time,model_time_dt = \
        get_model(simmode="fc",model=model,fc_date=fc_date,
        leadtime=leadtime)
    #collocation
    results_dict = collocate(model,model_Hs,model_lats,
            model_lons,model_time_dt,sa_obj,fc_date,distlim=self.distlim)
    self.assertEqual(model_Hs.__class__.__name__,'ndarray')
    self.assertEqual(np.isnan(np.nanmean(model_Hs)),False)
    self.assertGreater(len(model_time_dt),0)
    self.assertEqual(results_dict.__class__.__name__,'dict')

def test_get_model_station(self):
    import numpy as np
    from datetime import datetime,timedelta
    from modelmod import get_model
    from stationmod import station_class as sc
    from stationmod import matchtime
    from stationmod import get_loc_idx
    model = "mwam4"
    region = "mwam4"
    fc_date = datetime(2018,10,1,12)
    init_date = datetime(2018,10,1,12)
    leadtime = 0
    sc_obj = sc("ekofiskL",fc_date,fc_date)
    self.assertEqual(sc_obj.__class__.__name__,'station_class')
    model_Hs,model_lats,model_lons,model_time,model_time_dt = \
        get_model(simmode="fc",model=model,fc_date=fc_date,
        leadtime=leadtime)
    mHs = model_Hs.squeeze()
    ctime, cidx = \
        matchtime(fc_date,fc_date,sc_obj.time,sc_obj.basedate)
    idx,idy,distM,mlats_new,mlons_new = \
        get_loc_idx(model_lats,model_lons,sc_obj.lat,sc_obj.lon,mHs)
    self.assertEqual(mHs.__class__.__name__,'ndarray')
    self.assertEqual(np.isnan(np.nanmean(mHs)),False)
    self.assertGreater(len(idx),0)
    self.assertEqual(np.isnan(idx[0]),False)



