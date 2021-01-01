#!/usr/bin/env python

'''
adding INSTRUMENT keyword to header for INT data so it can fit the geometric distortion for each chip separately.

USAGE:

python ~/github/Virgo/programs/INT_add_inst_scamp.py --ref r1442821.fit --image r1442823.fit --fixall

'''


from astropy.io import fits
from astropy import wcs
import glob
import os

matchstrings = ['WFC*1PA.fits','WFC*2PA.fits','WFC*3PA.fits','WFC*4PA.fits']
instruments = ['INTWFC1','INTWFC2','INTWFC3','INTWFC4']

for i in range(len(matchstrings)):
    files = glob.glob(matchstrings[i])
    print('chip ',i+1,' updating ',len(files),' files')
    for f in files:
        # read in image and header for the image that needs to be updated
        hdu = fits.open(f)
        hdu.verify()

        # trying again after implementing two commands that might fix the issue with CD1_1
        hdu[0].header.set('INSTRMNT',instruments[i])

        hdu.writeto(f,overwrite=True,output_verify='ignore')
