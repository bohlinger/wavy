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
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg ARCMFC3 &> arcmfc3_coll_test.out'
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sat s3a -mod ARCMFC3 -reg ARCMFC3'
t = os.system(cmd)
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg NordicSeas &> arcmfc3_coll_NordicSeas_test.out'
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sat s3a -mod ARCMFC3 -reg NordicSeas'
t = os.system(cmd)
"""
 Validate
"""
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg ARCMFC3 &> arcmfc3_val_test.out &'
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sat s3a -mod ARCMFC3 -reg ARCMFC3'
t = os.system(cmd)
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg NordicSeas &> arcmfc3_val_NordicSeas.out &'
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sat s3a -mod ARCMFC3 -reg NordicSeas'
t = os.system(cmd)
"""
 Make figures
"""
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg ARCMFC3 -d 202001'
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg ARCMFC3'
t = os.system(cmd)
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg NordicSeas -d 202001'
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg NordicSeas'
t = os.system(cmd)
"""
 Make homepage
"""
cmd = 'sh /lustre/storeB/project/fou/om/waveverification/ARCMFC3/satellites/altimetry/s3a/WebPage/webpage.sh'
t = os.system(cmd)
