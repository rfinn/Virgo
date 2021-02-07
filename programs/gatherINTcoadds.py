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

homedir = os.getenv("HOME")
# define directory for all coadds
output_dir_coadds = homedir+'/data/virgo-coadds/'
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
for f1 in flist1:
    # if item is a directory and starts with 20, then assume it is the date
    # save dir name as the date
    if os.path.isdir(f1):
        dateobs = f1
        
        # get list of directory
        flist2 = os.listdir(f1)
        flist2.sort()
        print('WORKING ON DIRECTORY: ',f1)
        # loop through list
        for f2 in flist2:
            subdir = os.path.join(f1,f2) # e.g. 20190204/pointing001-r
            # if item is a directory and the name contains pointing, then assume it is a target
            if os.path.isdir(subdir) & (subdir.starts('pointing') ) & (subdir.find('-') > -1):

                # store pointing and filter
                long_pointing,filter = f2.split('-')
                pointing = long_pointing.replace('ointing','')

                #print('\t found a pointing',pointing,filter)
                # look for subdirectory named "coadd_"+filter
                coadd_path = os.path.join(subdir,'coadd_'+filter)
                # this is directory structure setup by theli
                # when I processed them myself, the coadds are in the main directory
                if os.path.exists(coadd_path):
                    print('\t found',coadd_path)
                    fff_file = 'coadd_'+filter+'/fffcoadd.fits'
                    ff_file = 'coadd_'+filter+'/fffcoadd.fits'
                    f_file = 'coadd_'+filter+'/fffcoadd.fits'
                    c_file = os.path.join(coadd_path,'coadd.fits')

                    if os.path.exists(c_file):
                        imfile = c_file
                    else:
                        print('\t WARNING: no coadd in ',dateobs,' ',f2)
                        continue

                    '''
                    if os.path.exists(fff_file):
                        imfile = fff_file
                    elif os.path.exists(ff_file):
                        imfile = ff_file
                    elif os.path.exists(f_file):
                        imfile = f_file
                    '''

                    weight_file = imfile.split('.fits')[0]+'.weight.fits'

                    h = fits.getheader(imfile)
                    ra = float(h['CRVAL1'])
                    dec = float(h['CRVAL2'])

                    # create string for output name
                    if float(dec) < 0:
                        outfile = output_dir_coadds+'VF-{:.4f}-{:.4f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter)
                    else:
                        outfile = output_dir_coadds+'VF-{:.4f}+{:.4f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter)

                    # copy imfile to outfile
                    out1 = outfile+'.fits'
                    out2 = outfile+'.weight.fits'
                    if (not os.path.exists(out1)) or overwrite:
                        print('\t   copy ',imfile,' -> ',out1)
                        shutil.copyfile(imfile,out1)
                    else:
                        print('\t   '+out1,' already exists. set overwrite if you want to copy anyway.')
                    if (not os.path.exists(out2)) or overwrite:
                        # copy weight file
                        shutil.copyfile(weight_file,out2)                    
                        print('\t   weight file = ',weight_file)
    


