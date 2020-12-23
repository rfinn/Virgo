#!/usr/bin/env python

'''
Some INT images have missing header info

the most important fields are CRVAL1 and CRVAL2

going to write this so you can enter a ref image and then update the header of the
broken image using the coordinates in the ref image.

'''


from astropy.io import fits

import argparse
import os

parser = argparse.ArgumentParser(description ='Create a crazy big catalog from HL, AGC, NSA')
parser.add_argument('--image',dest = 'image', default=None,help='image with bad header')
parser.add_argument('--ref',dest = 'ref', default=None,help='reference image, to use RA and DEC from')
parser.add_argument('--fields',dest = 'fields', nargs='+',default=['CRVAL1','CRVAL2','AIRMASS'],help='list of header fields to update. default is CRVAL1 and CRVAL2')

# even for MEF files, the main header is 0, and this is what comes in 
parser.add_argument('--mef',dest = 'mef', action='store_true', default=False, help='set this if input is a multi-extension fits file')
parser.add_argument('--fixall',dest = 'fixall', action='store_true', default=False, help='set this to copy the entire header')

args = parser.parse_args()

if args.mef:
    args.fields = ['RA','DEC','EQUINOX']




href = fits.getheader(args.ref)

if args.fixall:
    h = fits.getheader(args.image)
    # find fields in reference header that are not in bad header
    # for INT data, the telescope telemetry fields are missing (rather than empty)
    # this may not work for other datasets...
    goodh = set(href)
    badh = set(h)
    fields = list(goodh.difference(badh))
    print('fields to update: ',fields)
else:
    fields = args.fields

for f in fields:
    #print(f)
    if (f.find('DUMMY') > -1) | (f.find('HISTORY') > -1):
        continue
    try: 
        newval = href[f]
    except:
        print('error with ',f)
        continue
    # getting an error when trying to write out the fits header for INT images
    # Card 'CD1_1' is not FITS standard (invalid value string: '-9.19444e-5').
    # so using sethead instead
    try:
        os.system('sethead '+args.image+' '+f+'='+str(newval))
    except KeyError:
        os.system('sethead -k '+args.image+' '+f+'='+str(newval))


