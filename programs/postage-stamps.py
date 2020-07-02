#!/usr/bin/env python

'''
GOAL: create figure showing galaxy at many wavelengths
- legacy color - from legacy server
- UV
- WISE 12
- WISE 22
- R, Halpha



'''

def getlegacyimages(ra,dec):
    '''
    new function to download images in one fell swoop

    doing this so we can make cutouts of entire catalog, and then check each galaxy by hand

    This will need to be re-run each time the matching is altered, and the kitchen-sink catalog changes.
    
    '''
    for i in range(len(ra)):
        # name image files by ra and dec of galaxy
        gra = '%.5f'%(ra[i]) # accuracy is of order .1"
        gdec = '%.5f'%(dec[i])
        galnumber = gra+'-'+gdec
        rootname = 'cutouts/legacy-im-'+str(galnumber)+'-'+str(image_size)
        jpeg_name = rootname+'.jpg'

        fits_name = rootname+'.fits'
        # check if images already exist
        # if not download images
        if not(os.path.exists(jpeg_name)):
            print('retrieving ',jpeg_name)
            url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra[i])+'&dec='+str(dec[i])+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
            urlretrieve(url, jpeg_name)
        else:
            print('previously downloaded ',jpeg_name)
        if not(os.path.exists(fits_name)):
            print('retrieving ',fits_name)
            url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra[i])+'&dec='+str(dec[i])+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
            urlretrieve(url, fits_name)
        else:
            print('previously downloaded ',fits_name)
    pass


## need to decide if this is based on the halpha image
## could just get sizes from VFMAIN catalog

## let's start with an halpha based program
