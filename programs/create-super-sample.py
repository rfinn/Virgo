#!/usr/bin/env python

'''
GOAL:
- create a super list of all galaxies in

USAGE:
- here you go

'''

from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
## import argparse

## parser = argparse.ArgumentParser(description ='Match the Halpha observations with Virgo NSA catalog')

## parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
## parser.add_argument('--write-fits',dest = 'writefits', action='store_true',help='write out fits version of observing summary file?')
## parser.add_argument('--input',dest = 'input', default='Observing-Summary-Halpha-good-latest.csv',help='write out fits version of observing summary file?')
        
## args = parser.parse_args()



## read in HL catalog
gl = fits.getdata('/Users/rfinn/research/VirgoFilaments/Gianluca/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')

## read in NSA catalog

## read in AGC catalog
# got a new version from Martha on 11/19/19
agcfile = '/Users/rfinn/research/AGC/agcm1.sh191118.fits'
agc = fits.getdata(agcfile)

## start with HL - match to NSA and AGC
## consider all matches with offsets less than 5" to be the same galaxy


# match NSA to AGC
insa_n2a, d2d_n2a, d3d_n2a = acoord.match_to_catalog_sky(ncoord)
# first look at number with match w/in 10"
matchflag_n2a = d2d_n2a < 10./3600*u.deg
print('number of matches w/in 10 arcsec = ',sum(matchflag_n2a),'/',len(matchflag_n2a))

## track NSA and AGC names of matches

## add HL objects with closest match > 5"

## match remaining AGC and NSA

## consider all matches with offsets less than 5" to be the same galaxy

## list NSA name, track AGC name of match

## add remaining NSA galaxies with no HL and no AGC within 5"

## add remaining AGC galaxies with no HL and no NSA within 5"




vdat = fits.getdata(args.tablepath +'nsa.virgo.fits')
hdat = fits.getdata(args.tablepath + 'Observing-Summary-Halpha-good-latest.fits')

nsadict=dict((a,b) for a,b in zip(vdat.NSAID,np.arange(len(vdat.NSAID))))
# match by NSAID
index = np.zeros(len(hdat.nsa_id),'i')
matchflag = np.zeros(len(hdat.nsa_id),'bool')
for i in range(len(hdat.nsa_id)):
    try:
        index[i] = nsadict[int(hdat.nsa_id[i])]
        matchflag[i] = 1
    except:
        print 'could not match NSAID = ',hdat.nsa_id[i]
        #print 'recession velocity = ',vdat.Z[nsadict[hdat.nsa_id[i]]]*3.e5
    





