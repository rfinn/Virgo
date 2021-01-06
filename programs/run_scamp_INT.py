#!/usr/bin/env python

'''
get list of directories
move into each directory
run scamp and swarp


'''

import os
import shutil
from astropy.io import fits

homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
for subdir in flist1:
    print('WORKING ON DIRECTORY: ',f1)
    # loop through list
    if os.path.isdir(subdir) & (subdir.find('pointing') > -1):
        # move to subdirectory
        os.chdir(subdir)
        # run scamp
        scampflag=False
        try:
            os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --scamp --int --filestring WFC')
            scampflag=True
            os.system('ls WFC*PA.fits > swarp_inlist')
        except:
            print('##########################################')
            print('WARNING: problem running scamp for ',subdir)
            print('##########################################')
            
        # run swarp
        if scampflag:
            
            os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --filestring WFC')
