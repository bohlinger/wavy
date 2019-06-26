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
#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate.py -sd 2019030100 -ed 2019033118'
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate.py'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_collocate_NordicSeas.py'
t = os.system(cmd)
#p = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
#(output, err) = p.communicate()
#This lets the program wait until this process is finished
# before it proceeds with the next command
#p_status = p.wait()
#This will give you the output of the command being executed
#print("Command output: " + output)

#cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate.py -sd 2019030100 -ed 2019033118'
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate.py'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_validate_NordicSeas.py'
t = os.system(cmd)

cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_figures.py'
t = os.system(cmd)
cmd = 'python /home/patrikb/wavy/arcmfc/arcmfc_figures_NordicSeas.py'
t = os.system(cmd)

#cmd = 'sh /home/patrikb/wavy/web/arcmfc/webpage.sh'
cmd = 'sh /lustre/storeB/project/fou/om/ARCMFC/s3a/WebPage/webpage.sh'
t = os.system(cmd)
