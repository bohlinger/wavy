'''
    Pogram to manage the workflow for station data collection
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/op')

import subprocess
import os

# --- stations --- #
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var Hs_10min >& /home/patrikb/wavy/logs/collect_station_all_Hs_10min.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var Tm02_10min >& /home/patrikb/wavy/logs/collect_station_all_Tm02_10min.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var Tm01_10min >& /home/patrikb/wavy/logs/collect_station_all_Tm01_10min.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var Tm10_10min >& /home/patrikb/wavy/logs/collect_station_all_Tm10_10min.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var Mdir_10min >& /home/patrikb/wavy/logs/collect_station_all_Mdir_10min.out'
t = os.system(cmd)
#
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var Pdir_10min >& /home/patrikb/wavy/logs/collect_station_all_Pdir_10min.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var FF_10min_10m >& /home/patrikb/wavy/logs/collect_station_all_FF_10min_10m.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var FF_10min_sensor >& /home/patrikb/wavy/logs/collect_station_all_FF_10min_sensor.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/collect_station_all.py -var DD_10min_sensor >& /home/patrikb/wavy/logs/collect_station_all_DD_10min_sensor.out'
t = os.system(cmd)
