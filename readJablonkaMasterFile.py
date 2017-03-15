#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
infile = 'tables/Jablonka-MasterFile-Finn.list'

filaments = ['Filament1','Filament1-complement','N5353group','Filament3','LeoII-B','LeoII-A','LeoMinorFilament4','VirgoIIIcenter','VirgoIIIWest']

inf1 = open(infile,'r')
NEDname=[]
SDSS_ID = []
RA = []
DEC = []
CO = []
CO_DETECT = []
HI =[]
filament= []

i=0
for line in inf1:
    i+=1
    if i < 11: # skipp first 10 header lines
        continue
    if line.startswith('#'):
        j,current_filament = line.split('#')
        i += 1
        continue
    t=line.split()
    if len(t) == 6:
        NEDname.append(t[0])
        SDSS_ID.append(t[1])
        RA.append(t[2])
        DEC.append(t[3])
        CO.append(t[4])
        CO_DETECT.append(t[5])
        filament.append(current_filament)
        HI.append('-1')
    elif len(t) == 7:
        NEDname.append(t[0])
        SDSS_ID.append(t[1])
        RA.append(t[2])
        DEC.append(t[3])
        CO.append(t[4])
        CO_DETECT.append(t[5])
        HI.append(t[6])
        filament.append(current_filament)

coords = SkyCoord(RA,DEC, unit= (u.hour,u.deg))

# write out a fits file


col1 = fits.Column(name='NED_name',format='20A',array=np.array(NEDname))
col2 = fits.Column(name='SDSS_ID',format='20A',array=np.array(SDSS_ID))
col3 = fits.Column(name='RA',format='E',array=coords.ra.deg)
col4 = fits.Column(name='DEC',format='E',array=np.array(coords.dec.deg))
col5 = fits.Column(name='CO',format='2A',array=np.array(CO))
col6 = fits.Column(name='CO_DETECT',format='I',array=np.array(CO_DETECT))
col7 = fits.Column(name='HI',format='J',array=np.array(HI))
col8 = fits.Column(name='filament',format='20A',array=np.array(filament))

cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8])

tbhdu = fits.BinTableHDU.from_columns(cols)

tbhdu.writeto('tables/Jablonka-MasterFile.fits',clobber=True)

