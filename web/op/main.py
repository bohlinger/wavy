'''
    Pogram to manage the workflow for the op-product quality web page
    Executes in order:
        - op_collocate (collocation of obs and mod)
        - op_validate (compute validation stats)
        - op_figures (make validation figures)
        - webpage.sh (creates html webpage)
'''
import sys
sys.path.append(r'/home/patrikb/wavy/wavy')
sys.path.append(r'/home/patrikb/wavy/wavy/op')

import subprocess
import os

#cmd = 'python /home/patrikb/wavy/op/op_collocate.py -m mwam4 -sd 2019020100 -ed 2019022123'
cmd = 'python /home/patrikb/wavy/op/op_collocate.py'
t = os.system(cmd)
#cmd = 'python /home/patrikb/wavy/op/op_validate.py -m mwam4 -sd 2019020100 -ed 2019022123'
cmd = 'python /home/patrikb/wavy/op/op_validate.py'
t = os.system(cmd)

cmd = 'python /home/patrikb/wavy/op/op_figures.py'
t = os.system(cmd)

cmd = 'sh /home/patrikb/wavy/web/op/webpage.sh'
t = os.system(cmd)

#cmd = 'cp index.html /lustre/storeB/project/fou/om/waveverification/s3a/'
#cmd = 'cp index.html /lustre/storeB/project/fou/om/waveverification/mwam4/s3a/'
#t = os.system(cmd)
