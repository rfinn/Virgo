#!/usr/bin/env python

'''
GOAL:
- create a super list of all galaxies in

USAGE:
- here you go

'''

from astropy.io import fits, ascii
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import os

homedir = os.getenv('HOME')
## import argparse

## parser = argparse.ArgumentParser(description ='Match the Halpha observations with Virgo NSA catalog')

## parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
## parser.add_argument('--write-fits',dest = 'writefits', action='store_true',help='write out fits version of observing summary file?')
## parser.add_argument('--input',dest = 'input', default='Observing-Summary-Halpha-good-latest.csv',help='write out fits version of observing summary file?')
        
## args = parser.parse_args()

class sample:
    def __init__(self):
        
        ## read in my HL catalog
        '''
        Hyperleda query:
        select
        objname,objtype,de2000,al2000,v,e_v,vopt,e_vopt,vrad,e_vrad,bt,e_bt,type,bar,ring,multiple,compactness,t,e_t,logd25,e_logd25,logr25,e_logr25,pa,incl,logdc,btc,itc,ubtc,bvtc,m21c,hic,mabs,agnclass,kt,e_kt,it,e_it,ut,vt,mfir,e_ut,e_vt
        where de2000 > -35 and de2000 < 75 and  al2000 < 280./360.*24.  and al2000 > 100./360.*24. and v < 3300 and objtype='G'
    
        not using Gialuca's catalog b/c I'm not sure if there were other cuts made already
        gl = fits.getdata('/Users/rfinn/research/VirgoFilaments/Gianluca/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')
    
        '''
    
        hlfile = homedir+'/github/Virgo/tables/hyperleda-finn-11nov2019-full.csv'
        self.hl = ascii.read(hlfile)
        self.hl = Table(self.hl)
    
        ## read in NSA catalog
        '''
        using new NSA catalog
        '''
        nfile = homedir+'/research/NSA/nsa_v1_0_1.fits'
        self.nsa = fits.getdata(nfile)
        ## read in AGC catalog
        # got a new version from Martha on 11/19/19
        agcfile = '/Users/rfinn/research/AGC/agcm1.sh191118.fits'
        self.agc = fits.getdata(agcfile)

        # flags to track nsa matches to HL and AGC
        self.hl_2_nsa_matchflag = np.zeros(len(self.hl['al2000']),'bool')
        self.hl_2_agc_matchflag = np.zeros(len(self.hl['de2000']),'bool')

        # flags to track nsa matches to HL and AGC
        self.nsa_2_hl_matchflag = np.zeros(len(self.nsa.RA),'bool')
        self.nsa_2_agc_matchflag = np.zeros(len(self.nsa.RA),'bool')

        # flags to track AGC matches to HL and NSA
        self.agc_2_hl_matchflag = np.zeros(len(self.agc.radeg),'bool')
        self.agc_2_nsa_matchflag = np.zeros(len(self.agc.radeg),'bool')
    def match_hl_2_nsa(self):
        pass

    def match_hl_2_agc(self):
        pass

    def identify_matches(self):
        '''
        HL catalog is the start of the super sample
        
        identify HL galaxies with AGC or NSA w/in 5"
         - track HL, NSA, and AGC name

        '''
        pass

    def match_nsa_2_agc(self):
        '''
        for NSA and AGC galaxies NOT already matched to HL, match NSA to AGC

        if offset is < 5", add galaxy to catalog

        track NSA and AGC names
        
        '''
        pass
    def add_remaining_NSA(self):
        '''
        add NSA galaxies with no counterpart in HL or AGC
        '''
        pass

    def add_remaining_AGC(self):
        '''
        add AGC galaxies with no counterpart in HL or NSA
        '''
        

    def 
if __name__ == '__main__':
        
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



    





