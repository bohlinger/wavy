'''
    Pogram to manage the workflow for station collocation
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/op')

import subprocess
import os

# --- stations --- #
### mwam4 ###
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Hs -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_Hs_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Mdir -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_Mdir_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Pdir -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_Pdir_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tm02 -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_Tm02_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tp -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_Tp_mwam4.out'
t = os.system(cmd)
### ecwam ###
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Hs -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_Hs_ecwam.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Mdir -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_Mdir_ecwam.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tp -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_Tp_ecwam.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all.py -var Tm02 -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_Tm02_ecwam.out'
t = os.system(cmd)
