#!/usr/bin/env

"""
GOAL:

reproject legacy images to halpha pixel scale

PROCEDURE:
input: get directory name


"""
import sys
import os
from astropy.io import fits
import glob
from reproject import reproject_interp

def reproject_image(infile, reffile, outname):
    """reproject infile to reffile image"""
    if os.path.exists(outname):
        print("reprojected image exists - not redoing it")
        return
    
    hinfile = fits.open(infile)
    href = fits.open(reffile)
    # reproject input to referece image
    outim,footprint = reproject_interp(hinfile,href[0].header)

    fits.writeto(outname,outim,href[0].header,overwrite=True)
    hinfile.close()
    href.close()

if __name__ == '__main__':
    dirname = sys.argv[1]

    #get CS image
    reffile = os.path.join(dirname,dirname+'-CS.fits')

    if not os.path.exists(reffile):
        print("can't find CS image - exiting")
        sys.exit()
    #get legacy/*.fits
    legacy_images = glob.glob(os.path.join(dirname,'legacy/*.fits'))
    legacy_images.sort()

    for infile in legacy_images:
        outname = infile.replace('.fits','-ha.fits')
        
        # reproject onto halpha wcs
        reproject_image(infile,reffile,outname)
