'''
    Program to manage the workflow for the arcmfc product quality web page
    Executes in order:
        - arcmfc_collocate (collocation of obs and mod)
        - arcmfc_validate (compute arcmfc validation stats)
        - arcmfc_figures (make arcmfc validation figures)
        - webpage.sh (creates html webpage)
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
import subprocess
import os
"""
 Collocate
"""
sdstr = '2020040100'
edstr = '2020041400'
sat = 'cfo'

cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sat ' + sat + ' -mod ARCMFC3 -reg ARCMFC3 -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sat ' + sat + ' -mod ARCMFC3 -reg NordicSeas -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
"""
 Validate
"""
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sat ' + sat + ' -mod ARCMFC3 -reg ARCMFC3 -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sat ' + sat + ' -mod ARCMFC3 -reg NordicSeas -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
