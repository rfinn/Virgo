#!/usr/bin/env python

from astropy.io import fits
from virgoCommon import *

tablepath = gitpath+'Virgo/tables/'
cofile = 'CO-Masterfile-2017May15.fits'
codat = fits.getdata(tablepath+cofile)

ngcflag =  codat.filament == 'Filament1-Group5354'

# match with stellar mass, NSA, WISE


