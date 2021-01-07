#!/usr/bin/env python

'''
get list of directories
move into each directory
run scamp and swarp

This assumes that the files are already sorted by pointing, but not by filter.

'''

import os
import shutil
from astropy.io import fits


def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0


homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()

for subdir in flist1: # loop through list
    if os.path.isdir(subdir) & (subdir.find('pointing') > -1):
        print('##########################################')
        print('##########################################')        
        print('WORKING ON DIRECTORY: ',subdir)
        print('##########################################')
        print('##########################################')
            
        # move to subdirectory
        os.chdir(subdir)
        # run scamp
        scampflag=False
        try:
            os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --scamp --int --filestring WFC')
            scampflag=True
            os.system('ls WFC.r*PA.fits > '+subdir+'_r')
            os.system('ls WFC.Halpha*PA.fits > '+subdir+'_Halpha')
            os.system('ls WFC.Ha6657*PA.fits > '+subdir+'_Ha6657')            
        except:
            print('##########################################')
            print('WARNING: problem running scamp for ',subdir)
            print('##########################################')
            
        # run swarp
        if scampflag:
            # count lines r band file, run if more than 2 lines
            suffix = ['_r','_Halpha','_Ha6657']
            filelists = [subdir+i for i in suffix]
            #################################################3
            # make r-band mosaic
            #################################################3            
            nlines = count_lines(filelists[0])
            
            if nlines > 2:
                os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0])
                refimage = filelists[0]+'.coadd.fits'
                for i in [1,2]:
                    nlines = count_lines(filelists[i])
                    if nlines > 3:
                        os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[i]+' --refimage '+refimage)
            else:
                print('WARNING: No r-band mosaic, making remaining mosaics without alignment')
                for f in filelists:
                    nlines = count_lines(f)
                    if nlines > 2:
                        os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0])
                    
        os.chdir(working_dir)
