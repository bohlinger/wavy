# sort.py
"""
purpose: mv files to sub-folders of year and month
"""

# import libraries
import numpy as np
import os

# flatten lists
flatten = lambda l: [item for sublist in l for item in sublist]

def get_localfiles():
    """
    get list of local files
    """
    dirpath = os.getcwd()
    filelst = np.sort(os.listdir(dirpath))
    return dirpath,filelst

def sort_files(dirpath,filelst):
    for e in filelst:
        tmp = e[22:30]
        year, month = tmp[0:4],tmp[4:6]
        folder = ('./' + year + '/' + month)
        os.system('mkdir -p ' + folder)
        os.system('mv ' + e + ' ' + folder)

if __name__ == "__main__":
    dirpath, filelst = get_localfiles()
    sort_files(dirpath,filelst)
