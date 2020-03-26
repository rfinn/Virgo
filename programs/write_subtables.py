#!/usr/bin/env python

'''
GOAL:
- read in final catalog vf_clean_sample.fits
- split table into line-matched tables containing
  - basic info: RA, DEC, vel, NEDname, HLname, A100name, NSAID, flag for each
  - HL
  - NSA
  - AGC
  - A100
  - unWISE
  
- write out several different views
  - one for matching with Leroy+2019 sample
  - one for Dustin Lang to get unWISE photometry
  - one for catalog paper


'''
from astropy.io import fits
import os

homedir = os.getenv("HOME")
masterfile = homedir+'/research/supersample/vf_clean_sample.fits'
outdir = homedir+'/research/tables/'
class catalog:
    def __init__(self,catalog):
        self.cat = fits.getdata(catalog)

    def catalog_for_z0MGS(self):
        '''
        need RA, DEC, and search radius
        in IPAC format table
        for matching with leroy+2019 Galaxy Synthesis WISE+GALEX
        table is served by IRSA
        '''
        search_radius = 10.*np.ones(len(self.cat)) # units are arcsec
        newtable = Table([self.cat['galnumber'],self.cat['RA'],self.cat['DEC'],search_radius],names=['galid','ra','dec','major'])
        newtable.write(outdir+'coords_for_z0MGS.txt',format='ipac',overwrite=True)
    def catalog_for_paper(self):
        # write out
        # ra, dec, velocity, HL, NSA id, A100 id
        # NED name
        colnames = ['galnumber','RA','DEC']
        newtable = Table([self.cat['galnumber'],self.cat['RA'],self.cat['DEC'],search_radius],names=['galid','ra','dec','major'])
        newtable.write(outdir+'coords_for_z0MGS.txt',format='ipac',overwrite=True)
        pass
    def write_table(self,colnames,outfile,format=None):
        if format is None:
            format = 'fits'
        else:
            format = format
        for c in colnames:
            mycolumns.append(self.cat[c])
        newtable = Table(mycolumns)
        newtable.write(outdir+outfile,format=format,overwrite=True)
if __name__ == '__main__':
    c = catalog(masterfile)
    
