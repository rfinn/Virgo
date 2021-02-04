#!/usr/bin/env python

'''
get list of directories
move into each directory
run scamp and swarp

This assumes that the files are already sorted by pointing, but not by filter.

Run this from, e.g. /home/rfinn/data/reduced/scratch-int-feb2019
- this directory has a subdirectory for each pointing
- the subdirectory contains both the r and Halpha image

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
runscamp=False
runswarp=True
for subdir in flist1: # loop through list
    if os.path.isdir(subdir) & (subdir.find('pointing') > -1):# & (subdir.find('-') > -1):
        print('##########################################')
        print('##########################################')        
        print('WORKING ON DIRECTORY: ',subdir)
        print('##########################################')
        print('##########################################')
            
        # move to subdirectory
        os.chdir(subdir)
        # run scamp
        if runscamp:
            scampflag=False # keeps track if scamp finishes successfully
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
        if runswarp:
            # move short exposure times
            os.system('python ~/github/HalphaImaging/python3/move_short_exposures.py --filestring WFC')        
            # subtract median
            os.system('python ~/github/HalphaImaging/python3/subtract_median.py --filestring WFC')

            # gather files
            os.system('ls WFC.r*PA.fits > '+subdir+'_r')
            os.system('ls WFC.Halpha*PA.fits > '+subdir+'_Halpha')
            os.system('ls WFC.Ha6657*PA.fits > '+subdir+'_Ha6657')            
            
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
                        #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[i]+' --refimage '+refimage)
                        os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[i]+' --refimage '+refimage+' --noback')
                # run swarp again on the rband data, using the same refimage
                os.system('cp '+refimage+' refimage.fits')
                #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0]+' --refimage refimage.fits')
                os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0]+' --refimage refimage.fits --noback')
            else:
                print('WARNING: No r-band mosaic, making remaining mosaics without alignment')
                for f in filelists:
                    nlines = count_lines(f)
                    if nlines > 2:
                        #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+f)
                        os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+f+' --noback')
                    else:
                        print('WARNING: not enough images to make mosaic in ',f)
            # remove median subtracted images
            os.system('rm mWFC*.fits')
        os.chdir(working_dir)
        # just running on one directory for testing purposes
        #break
