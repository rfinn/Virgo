#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from virgoCommon import *

table_ext = gitpath+'Virgo/tables/'
#Line Match with Jablonka Catalog

vdat = fits.getdata(table_ext +'nsa.virgo.fits')
#jdat = fits.getdata(table_ext + 'CO-MasterFile-2017May15.fits')
#jdat = fits.getdata(table_ext + 'CO-MasterFile-2018Feb16.fits')
jdat = fits.getdata(table_ext + 'All-virgo-20feb18_env_H070-FITS.fits')


virgocat = SkyCoord(vdat.RA*u.degree,vdat.DEC*u.degree,frame='icrs')
jcat = SkyCoord(jdat.RA*u.degree,jdat.DEC*u.degree,frame='icrs')

index,dist2d,dist3d = virgocat.match_to_catalog_sky(jcat)

# only keep matches with matched RA and Dec w/in 1 arcsec
matchflag = dist2d.degree < 5./3600

# write out line-matched catalog
#outfile= table_ext + 'nsa_CO-HI.virgo.fits'
outfile= table_ext + 'nsa_CO-Gianluca.virgo.fits'
matchedarray1=np.zeros(len(vdat),dtype=jdat.dtype)
matchedarray1[matchflag] = jdat[index[matchflag]]

fits.writeto(outfile,matchedarray1,overwrite=True)
#t = fits.getdata(outfile)
#t.add_column()

## #Adding a column of whether a galaxy is in a filament
## dat = fits.open(table_ext + 'CO-HI_virgo.fits')[1].data

## filcol = np.zeros(len(matchedarray1))
## filcol[matchflag] = 1.
## newcol = []
## newcol = fits.Column(name ='IsInFilament', format = 'D', array = filcol)

## #dat.columns()
## origcols = dat.columns
## hdu = fits.BinTableHDU.from_columns(origcols + newcol)
## outfile= table_ext + 'CO-HI_virgo.fits'
## hdu.writeto(outfile, clobber = 'True')
