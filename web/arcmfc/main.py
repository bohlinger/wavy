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
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate.py -sd 2019110100 -ed 2019110318'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate.py'
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg ARCMFC3 &> arcmfc3_coll_test.out'
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sat s3a -mod ARCMFC3 -reg ARCMFC3'
t = os.system(cmd)
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate_NordicSeas.py'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate_NordicSeas.py -sd 2019110100 -ed 2019110318'
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg NordicSeas &> arcmfc3_coll_NordicSeas_test.out'
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_collocate.py -sat s3a -mod ARCMFC3 -reg NordicSeas &> arcmfc3_coll_NordicSeas_test.out'
#t = os.system(cmd)

#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate.py -sd 2019110100 -ed 2019110318'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate.py'
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sd 2020020100 -ed 2020020818 -sat s3a -mod ARCMFC3 -reg ARCMFC3 &> arcmfc3_val_test.out &'
cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sat s3a -mod ARCMFC3 -reg ARCMFC3'
t = os.system(cmd)
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate_NordicSeas.py'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate_NordicSeas.py -sd 2019110100 -ed 2019110318'
#cmd = 'python3 /home/patrikb/wavy/arcmfc/arcmfc3_validate.py -sat s3a -mod ARCMFC3 -reg NordicSeas &> arcmfc3_val_test.out &'
#t = os.system(cmd)

#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_figures.py'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg ARCMFC3 -d 202001'
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg ARCMFC3'
t = os.system(cmd)
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_figures_NordicSeas.py'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg NordicSeas -d 202001'
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc3_figures.py -sat s3a -mod ARCMFC3 -reg NordicSeas'
#t = os.system(cmd)

#cmd = 'sh /home/patrikb/wavy/web/arcmfc/webpage.sh'
#cmd = 'sh /lustre/storeB/project/fou/om/ARCMFC/s3a/WebPage/webpage.sh'
cmd = 'sh /lustre/storeB/project/fou/om/waveverification/ARCMFC3/satellites/altimetry/s3a/WebPage/webpage.sh'
t = os.system(cmd)
