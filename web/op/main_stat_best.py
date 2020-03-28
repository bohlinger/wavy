'''
    Pogram to manage the workflow for station collocation bestguess
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/op')

import subprocess
import os

# --- stations --- #
### mwam4 ###
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Hs_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Mdir -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Mdir_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Pdir -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Pdir_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Tm02 -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Tm02_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Tp -mod mwam4 >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Tp_mwam4.out'
t = os.system(cmd)
### ecwam ###
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Hs -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Hs_ecwam.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Mdir -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Mdir_ecwam.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Tp -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Tp_ecwam.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_model_at_location_all_bestguess.py -var Tm02 -mod ecwam >& /home/patrikb/wavy/logs/collect_model_at_location_all_bestguess_Tm02_ecwam.out'
t = os.system(cmd)
