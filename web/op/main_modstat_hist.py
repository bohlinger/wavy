'''
    Pogram to manage the workflow for station collocation
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/op')

import subprocess
import os

sd = '2019120100'
ed = '2020040100'

# --- stations --- #
### mwam4 ###
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Hs -mod mwam4 -sd ' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Mdir -mod mwam4' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Pdir -mod mwam4' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tm02 -mod mwam4' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tp -mod mwam4' + sd + ' -ed ' + ed
t = os.system(cmd)
### ecwam ###
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Hs -mod ecwam' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Mdir -mod ecwam' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tp -mod ecwam' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tm02 -mod ecwam' + sd + ' -ed ' + ed
t = os.system(cmd)
# wind from mwam3force
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var U -mod mwam3force' + sd + ' -ed ' + ed
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var V -mod mwam3force' + sd + ' -ed ' + ed
t = os.system(cmd)
