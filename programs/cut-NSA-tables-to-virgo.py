#!/usr/bin/env python

'''
GOAL
- read in catalogs
- write out fits tables for Virgo region only


CATALOGS
- NSA
- WISE
- John M.'s stellar masses
- Simard B/D decompositions
- Yang central/satellite
- Salim GALEX-SDSS-WISE Legacy Catalog


CO is matched in a separate program: match-CO-NSA.py

'''

from astropy.io import fits
import os
import numpy as np

mypath=os.getcwd()
if mypath.find('rfinn') > -1:
    print "Running on Rose's computer"
    #agcfile='/Users/rfinn/idl/programs/idl_alfa/agctotal.sav'
    gitpath='/Users/rfinn/github/'
    nsapath = '/Users/rfinn/research/NSA/'
    gswlpath = '/Users/rfinn/Dropbox/Research/GSWLC/'
    agcpath = '/Users/rfinn/research/AGC/'

    
nsafile=nsapath+'nsa_v0_1_2.fits'
nsa=fits.getdata(nsafile)
# read in WISE catalog


# select galaxies near Virgo
raflag = (nsa.RA > 120.) & (nsa.RA < 220.) 
decflag= (nsa.DEC > -10.) & (nsa.DEC < 70.) 
velflag = (nsa.Z*3.e5 > 1000.) & (nsa.Z*3.e5 < 3000.)
vflag = raflag & decflag & velflag
vra=187.69708
vdec=12.33694
rad_distance = np.sqrt((nsa.RA - vra)**2+(nsa.DEC-vdec)**2)


nsaagcfile=nsapath+'nsa_v0_1_2.fits'
nsa_agc=fits.getdata(nsaagcfile)

wisefile=nsapath+'nsa_v0_1_2_wise.fits'
wise=fits.getdata(wisefile)
# read in John's stellar masses
massfile=nsapath+'nsa_v1_2_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits'
jmass=fits.getdata(massfile)

infile=nsapath+'NSA-GSWLC.fits'
gswlc=fits.getdata(infile)

infile=nsapath+'Simard1ToNSA.fits'
simard1=fits.getdata(infile)




fits.writeto(gitpath+'Virgo/tables/nsa.virgo.fits',nsa[vflag],overwrite=True)

# just contains wise columns, line-matched to above table
fits.writeto(gitpath+'Virgo/tables/nsa_wise.virgo.fits',wise[vflag],overwrite=True)
fits.writeto(gitpath+'Virgo/tables/nsa_mstar.virgo.fits',jmass[vflag],overwrite=True)
fits.writeto(gitpath+'Virgo/tables/nsa_gswlc.virgo.fits',gswlc[vflag],overwrite=True)
fits.writeto(gitpath+'Virgo/tables/nsa_simard1.virgo.fits',gswlc[vflag],overwrite=True)

