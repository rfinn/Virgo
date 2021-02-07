#!/usr/bin/env python

'''
GOAL: 
- go into each subdirectory and run getzp on coadded images
- delete any intermediate images to save space

Run this from, e.g. /home/rfinn/data/reduced/scratch-int-feb2019
- this directory has a subdirectory for each pointing


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib
matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
#print(flist1)
for subdir in flist1: # loop through list
    #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
    if os.path.isdir(subdir) & (subdir.find('pointing') > -1):
        print('##########################################')
        print('##########################################')        
        print('WORKING ON DIRECTORY: ',subdir)
        print('##########################################')
        print('##########################################')
            
        # move to subdirectory
        os.chdir(subdir)
        # get list of coadds
        # there will be two per directory per filter - one with background subtracted, and one without

        allcoadd = glob.glob(subdir+'*noback.coadd.fits')
        allcoadd.sort()
        for c in allcoadd:
            try:
                if c.find('Ha') > -1:
                    filter='ha'
                else:
                    filter='r'
                os.system('python ~/github/HalphaImaging/python3/getzp.py --instrument i --nexptime --image '+c+' --filter '+filter)
            except:
                print('##########################################')
                print('WARNING: problem running getzp for ',subdir, c)
                print('##########################################')

        os.chdir(working_dir)
        # just running on one directory for testing purposes
        #break



