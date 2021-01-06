#!/usr/bin/env python

'''

PROCEDURE:
* get list of current file
* combine images from the same pointing into one directory
* I had already split by pointing and filter, but I want to run scamp and swarp in the directory that contains all images from a pointing


'''

import os
import shutil
from astropy.io import fits
import glob

homedir = os.getenv("HOME")
# define directory for all coadds
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
output_dir_coadds = os.path.join(working_dir,output_dir_coadds)
# overwrite output files if they exist
overwrite = True
flist1.sort()
# get list of pointings - cut off filter from dir name
pointings = []
for f in flist1:
    if f.find('pointing') > -1:
        pointings.append(f.split('-')[0])

# pointings will now have two listing for each filter b/c of r and Halpha
upointings = set(pointings) # unique pointings

for p in upointings:
    os.mkdir(p)
    os.system('cp '+p+'-r/* '+p+'/.')
    os.system('cp '+p+'-Halpha/* '+p+'/.')
    os.system('cp '+p+'-Ha6657/* '+p+'/.')    


