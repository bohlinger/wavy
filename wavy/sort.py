# sort.py
"""
purpose: mv files to sub-folders of year and month
"""

# import libraries
import numpy as np
import os

def get_localfiles():
    """
    get list of local files
    """
    dirpath = os.getcwd()
    filelst = np.sort(os.listdir(dirpath))
    return dirpath,filelst

def sort_files(dirpath,filelst):
    for e in filelst:
        if os.path.isfile(dirpath + '/' + e):
            tmp = 'global_vavh_l3_rt_' + os.path.basename(dirpath) + '_'
            year, month = e[len(tmp):len(tmp)+4],e[len(tmp)+4:len(tmp)+6]
            folder = (dirpath + '/' + year + '/' + month)
            cmd = 'mkdir -p ' + folder
            os.system(cmd)
            cmd = 'mv ' + e + ' ' + folder
            os.system(cmd)

if __name__ == "__main__":
    dirpath, filelst = get_localfiles()
    sort_files(dirpath,filelst)
