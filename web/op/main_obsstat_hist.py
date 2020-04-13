'''
    Pogram to manage the workflow for station data collection
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/op')

import subprocess
import os

sd = '2020010100'
ed = '2020040100'

# --- stations --- #
#cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py' 
#        + ' -var Hs_10min'
#        + ' -sd ' + sd
#        + ' -ed ' + ed)
#t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var Tm02_10min'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var Tm01_10min'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var Tm10_10min'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var Mdir_10min'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var Pdir_10min'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var FF_10min_10m'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var FF_10min_sensor'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
cmd = ( 'python /home/patrikb/wavy/op/collect_station_all.py'
        + ' -var DD_10min_sensor'
        + ' -sd ' + sd
        + ' -ed ' + ed)
t = os.system(cmd)
