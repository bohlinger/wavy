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
sys.path.append(r'/home/patrikb/wavy/op')

import subprocess
import os

# --- satellites --- #
### mwam4 ###
#cmd = 'python /home/patrikb/wavy/op/op_collocate.py -mod mwam4 -sd 2019110100 -ed 2019110323'
cmd = 'python /home/patrikb/wavy/op/op_collocate.py'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/op_collocate.py -mod mwam4 -reg BalticSea -sat j3'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/op_collocate.py -mod mwam4 -reg BalticSea -sat al'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/op_collocate.py -mod mwam4 -reg BalticSea -sat c2'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/op_collocate.py -mod mwam4 -reg BalticSea -sat cfo'
t = os.system(cmd)
#cmd = 'python /home/patrikb/wavy/op/op_validate.py -mod mwam4 -sd 2019110100 -ed 2019110323'
cmd = 'python /home/patrikb/wavy/op/op_validate.py'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/op/op_figures.py'
t = os.system(cmd)
#cmd = 'sh /home/patrikb/wavy/web/op/webpage.sh'
#t = os.system(cmd)
#cmd = 'cp index.html /lustre/storeB/project/fou/om/waveverification/s3a/'
#cmd = 'cp index.html /lustre/storeB/project/fou/om/waveverification/mwam4/s3a/'
#t = os.system(cmd)

### ecwam ###
cmd = 'python /home/patrikb/wavy/web/op/op_collocate.py -mod ecwam -sat s3a -reg mwam4 &> /home/patrikb/wavy/logs/op_collocate_ecwam_s3a_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/web/op/op_collocate.py -mod ecwam -sat s3b -reg mwam4 &> /home/patrikb/wavy/logs/op_collocate_ecwam_s3b_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/web/op/op_collocate.py -mod ecwam -sat c2 -reg mwam4 &> /home/patrikb/wavy/logs/op_collocate_ecwam_c2_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/web/op/op_collocate.py -mod ecwam -sat j3 -reg mwam4 &> /home/patrikb/wavy/logs/op_collocate_ecwam_j3_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/web/op/op_collocate.py -mod ecwam -sat cfo -reg mwam4 &> /home/patrikb/wavy/logs/op_collocate_ecwam_cfo_mwam4.out'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/web/op/op_collocate.py -mod ecwam -sat al -reg mwam4 &> /home/patrikb/wavy/logs/op_collocate_ecwam_al_mwam4.out'
t = os.system(cmd)
