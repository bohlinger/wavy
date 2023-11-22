import pytest
from datetime import datetime, timedelta
import yaml
import numpy as np
import os
from copy import deepcopy
from wavy.insitu_module import insitu_class as ic


def test_filter_runmean(test_data):
    varalias = 'Hs'  # default
    sd = "2023-8-20 00"
    ed = "2023-8-21 00"
    nID = 'MO_Draugen_daily'
    name = 'Draugen'
    ico = ic(nID=nID, sd=sd, ed=ed, varalias=varalias, name=name)
    print(ico)
    print(vars(ico).keys())

    ico = ico.populate(path=str(test_data/"insitu/daily/Draugen"))
    new = ico.filter_runmean(window=3,
                             chunk_min=3,
                             sampling_rate_Hz=1/600)
    print(new.vars.time)
    print(new.vars.Hs)
    assert len(new.vars.time) == 6
    assert not all(np.isnan(v) for v in ico.vars['Hs'])
    print(ico.vars.Hs[1:4])
    print(np.mean(ico.vars.Hs[1:4]))
    assert new.vars.Hs[2] == np.mean(ico.vars.Hs[1:4])
