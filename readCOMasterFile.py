#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import time
# I removed some extra lines from Pascale's file
infile = 'tables/Jablonka-MasterFile-Finn.list'

infile = 'tables/Master_file_fm.dat'
infile = 'tables/Combes-2017May09-masterfile.dat'
infile = 'tables/CO-masterfile-2017May09.dat'

'''
UGCA201 - changed 1? to 2 for CO 
NGC3595 - changed HI detect from xxxx to xxx
'''


filaments = ['Filament1','Filament1-complement','N5353group','Filament3','LeoII-B','LeoII-A','LeoMinorFilament4','VirgoIIIcenter','VirgoIIIWest']

inf1 = open(infile,'r')
NEDname=[]
SDSS_ID = []
RA = []
DEC = []
vr = []
CO = []
CO_DETECT = []
HI =[]
filament= []
HImass = []
H2mass = []
upperlimit=[]
i=0
for line in inf1:
    i+=1
    if i < 14: # skipp first 15 header lines - this may change with future versions - check!!
        continue
    if len(line) < 1: # skip blank lines
        continue
    if line.startswith('#'):
        t = line.split('#')
        j,current_filament  = t
        current_filament = "".join(current_filament.replace("&","-").split())
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
        CO.append(t[8])
        CO_DETECT.append(t[9])
        HI.append(t[6])
        filament.append(current_filament)
    elif len(t) == 14:
        NEDname.append(t[0])
        SDSS_ID.append(t[1])
        RA.append(t[2]+':'+t[3]+':'+t[4])
        DEC.append(t[5]+':'+t[6]+':'+t[7])
        vr.append(t[8])
        CO.append(t[9])
        CO_DETECT.append(t[10].replace('xxx','-1'))
        HI.append(t[11].replace('xxx','-1').replace('N','0'))
        HImass.append(t[12].replace('xxx','-1'))
        h2 = t[13].replace('xxx','-1')
        if h2.find('<') > -1:
            upperlimit.append('1')
            H2mass.append(h2.replace('<',''))  # needed to edit catalog- some have 4 x's
        else:
            upperlimit.append('0')
            H2mass.append(h2)
        filament.append(current_filament)

coords = SkyCoord(RA,DEC, unit= (u.hour,u.deg))

# write out a fits file


col1 = fits.Column(name='NED_name',format='20A',array=np.array(NEDname))
col2 = fits.Column(name='SDSS_ID',format='20A',array=np.array(SDSS_ID))
col3 = fits.Column(name='RA',format='E',array=coords.ra.deg)
col4 = fits.Column(name='DEC',format='E',array=np.array(coords.dec.deg))
col5 = fits.Column(name='CO',format='2A',array=np.array(CO))
col6 = fits.Column(name='CO_DETECT',format='J',array=np.array(CO_DETECT))
col7 = fits.Column(name='HI',format='J',array=np.array(HI))
col8 = fits.Column(name='HImass',format='E',array=np.array(HImass))
col9 = fits.Column(name='H2mass',format='F',array=np.array(H2mass))
col10 = fits.Column(name='H2mass_upperlimit',format='I',array=np.array(upperlimit))
col11 = fits.Column(name='filament',format='20A',array=np.array(filament))

cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11])

tbhdu = fits.BinTableHDU.from_columns(cols)

date_stamp = time.strftime("%Y%b%d")
tbhdu.writeto('tables/CO-MasterFile-'+date_stamp+'.fits',clobber=True)
#tbhdu.writeto('tables/Jablonka-MasterFile-'+date_stamp+'.fits',clobber=True)
#tbhdu.writeto('tables/Jablonka-MasterFile.fits',clobber=True)

