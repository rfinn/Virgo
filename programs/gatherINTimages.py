#!/usr/bin/env python

'''

PROCEDURE:
* get list of current file

* run this from the base data directory, like ~/data/INT
  - this have subfolders arranged by date
* program will look in each subfolder, and then in each directory that has a "pointing" in the name
* it then looks in YYYYMMDD/pointingXXX-r/coadd-r for the coadd.fits file
* it will copy that to the output_dir_coadds directory specified below
* coadds will be renamed by RA, DEC, telescope, pointing, and filter

'''

import os
import shutil
from astropy.io import fits
import glob

homedir = os.getenv("HOME")
# define directory for all coadds
output_dir_coadds = 'allimages/'
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
output_dir_coadds = os.path.join(working_dir,output_dir_coadds)
# overwrite output files if they exist
overwrite = True
flist1.sort()
for f1 in flist1:
    # if item is a directory and starts with 20, then assume it is the date
    # save dir name as the date

    if os.path.isdir(f1) & (f1.find('pointing') > -1) & (f1.find('-') > -1):
        print('checking ',f1)
        flist2 = glob.glob(f1+'/WFC*.fits')

        #print(flist2)
        for f2 in flist2:
            print(f2,os.path.basename(f2))
            fname = os.path.basename(f2)
            imfile = os.path.join(working_dir,f2)
            outfile = os.path.join(output_dir_coadds,fname)
            if (not os.path.exists(outfile)) or overwrite:
                print('\t   mv ',imfile,' -> ',outfile)
                shutil.move(imfile,outfile)
            else:
                print('\t   '+outfile,' already exists. set overwrite if you want to copy anyway.')


