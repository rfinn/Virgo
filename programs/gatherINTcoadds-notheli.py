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


RUNNING FOR data/reduced/scratch-int-feb2019/attempt2/
'''

import os
import shutil
from astropy.io import fits
import glob

homedir = os.getenv("HOME")
# define directory for all coadds
output_dir_coadds = homedir+'/data/reduced/virgo-coadds-feb2019-int/'
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
workingdir = os.getcwd()
for subdir in flist1:
    if os.path.isdir(subdir) & (subdir.find('pointing')>-1 ):# & (subdir.find('-') > -1):
        os.chdir(subdir)
        # store pointing and filter
        #long_pointing,filter = subdir.split('-')
        # redoing for when both filters are in the same directory
        long_pointing = subdir
        pointing = long_pointing.replace('ointing','')


        ## GRAB THE COADDS IN THIS DIRECTORY
                
        #print('\t found a pointing',pointing,filter)
        # look for subdirectory named "coadd_"+filter

        # this is directory structure setup by theli
        # when I processed them myself, the coadds are in the main directory
        filters = ['r','Halpha','Ha6657']
        for i,filter in enumerate(filters):
            if i == 0:
                imfile = 'fn'+long_pointing+'_r.noback.coadd.fits'
                fstring='r'
            elif i == 1:
                imfile = 'ffn'+long_pointing+'_Halpha.noback.coadd.fits'
                fstring='Halpha'
            elif i == 2:
                imfile = 'ffn'+long_pointing+'_Ha6657.noback.coadd.fits'
                fstring='Ha6657'
            if not os.path.exists(imfile):
                continue
            weight_file = imfile.strip('ffn').split('.fits')[0]+'.weight.fits'

            # grab date from one of the r-band files
            # get obs date for coadd file name
            rfiles = glob.glob('WFC.'+fstring+'*4PA.fits')            
            t = rfiles[0].split('T')
            dateobs = t[0].split(fstring+'.')[1].replace('-','')
            
            # read in one individual images to get airmass and obsdate to pass to coadd
            nmiddle = int(len(rfiles)/2)
            im1header = fits.getheader(rfiles[nmiddle])

            obs_date = im1header['DATE-OBS']
            airmass = im1header['AIRMASS']            

            # read in coadd so I can update the header to add date-obs and airmass
            hdu = fits.open(imfile)
            h = hdu[0].header
            ra = float(h['CRVAL1'])
            dec = float(h['CRVAL2'])
            hdu[0].header.set('DATE-OBS',value=obs_date)
            hdu[0].header.set('AIRMASS',value=airmass)
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
                # write coadd to new location with updated header
                hdu.writeto(out1,overwrite=True)
                hdu.close()
                #shutil.copyfile(imfile,out1)
            else:
                print('\t   '+out1,' already exists. set overwrite if you want to copy anyway.')
            if (not os.path.exists(out2)) or overwrite:
                # copy weight file
                shutil.copyfile(weight_file,out2)                    
                print('\t   weight file = ',weight_file)
                
        os.chdir(workingdir)
        #break


