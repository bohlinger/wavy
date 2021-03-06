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
sdstr = '2020070100'
edstr = '2020073118'
sat = 's3a'

cmd = 'python3 /home/patrikb/wavy/op/op_collocate.py -sat ' + sat + ' -mod ARCMFC3 -reg ARCMFC3 -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
cmd = 'python3 /home/patrikb/wavy/op/op_collocate.py -sat ' + sat + ' -mod ARCMFC3 -reg NordicSeas -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
"""
 Validate
"""
cmd = 'python3 /home/patrikb/wavy/op/op_validate.py -sat ' + sat + ' -mod ARCMFC3 -reg ARCMFC3 -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
cmd = 'python3 /home/patrikb/wavy/op/op_validate.py -sat ' + sat + ' -mod ARCMFC3 -reg NordicSeas -sd ' + sdstr + ' -ed ' + edstr
t = os.system(cmd)
"""
 Make figures
"""
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat ' + sat + ' -mod ARCMFC3 -reg ARCMFC3 -d 202007'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat ' + sat + ' -mod ARCMFC3 -reg NordicSeas -d 202007'
t = os.system(cmd)
"""
 Make homepage
"""
#cmd = 'sh /lustre/storeB/project/fou/om/waveverification/ARCMFC3/satellites/altimetry/s3a/WebPage/webpage.sh'
#t = os.system(cmd)
cmd = 'sh /lustre/storeB/project/fou/om/waveverification/ARCMFC3/satellites/altimetry/s3a/WebPage/webpage.sh ARCMFC3'
t = os.system(cmd)
cmd = 'sh /lustre/storeB/project/fou/om/waveverification/ARCMFC3/satellites/altimetry/s3a/WebPage/webpage.sh NordicSeas'
t = os.system(cmd)

