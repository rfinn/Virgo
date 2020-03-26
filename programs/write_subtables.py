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
from astropy.table import Table
import os
import numpy as np

homedir = os.getenv("HOME")
masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample.fits'
outdir = homedir+'/research/Virgo/tables/'
class catalog:
    def __init__(self,catalog):
        self.cat = Table(fits.getdata(catalog))
    def catalog_for_z0MGS(self):
        '''
        need RA, DEC, and search radius
        in IPAC format table
        for matching with leroy+2019 Galaxy Synthesis WISE+GALEX
        table is served by IRSA
        '''
        search_radius = 10.*np.ones(len(self.cat)) # units are arcsec
        newtable = Table([self.cat['VFID'],self.cat['RA'],self.cat['DEC'],search_radius],names=['vfid','ra','dec','major'])
        newtable.write(outdir+'coords_for_z0MGS.txt',format='ipac',overwrite=True)
    def main_table(self):
        # write out
        # ra, dec, velocity, HL, NSA id, A100 id
        # NED name
        colnames = ['VFID','RA','DEC','vr','objname','NSAID','AGC','NEDname','HLflag','NSAflag','A100flag']
        self.write_table(colnames,'vf_main.fits',format='fits')
    def hyperleda_table(self):
        colnames = self.cat.colnames[0:43]
        self.write_table(colnames,'vf_hyperleda.fits',format='fits')
        # write out HL columns in line-matched table
        pass
    def nsa_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[90:200]
        # RA and DEC get changed to RA_2, DEC_2 during the matching process
        # changing them back to native NSA values
        newcolnames = colnames.copy()
        newcolnames[2] = 'RA'
        newcolnames[3] = 'DEC'
        self.write_table(colnames,'vf_nsa.fits',format='fits',names=newcolnames)
    def a100_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[205:230]
        self.write_table(colnames,'vf_a100.fits',format='fits')
    def a100_sdss_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[230:324]
        self.write_table(colnames,'vf_a100_sdssphot.fits',format='fits')
    def a100_unwise_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[324:389]
        # RA and DEC get changed to RA_2, DEC_2 during the matching process
        # changing them back to native NSA values
        newcolnames = colnames.copy()
        newcolnames[1] = 'ra'
        newcolnames[2] = 'dec'
        #print(colnames)
        self.write_table(colnames,'vf_a100_unwise.fits',format='fits',names=newcolnames)
    def write_table(self,colnames,outfile,format=None,names=None):
        if format is None:
            format = 'fits'
        else:
            format = format
        mycolumns = []
        for c in colnames:
            mycolumns.append(self.cat[c])
        if names is not None:
            newtable = Table(mycolumns,names=names)
        else:
            newtable = Table(mycolumns)
        newtable.write(outdir+outfile,format=format,overwrite=True)
if __name__ == '__main__':
    c = catalog(masterfile)    
    c.catalog_for_z0MGS()
    c.main_table()
    c.hyperleda_table()
    c.nsa_table()
    c.a100_table()
    c.a100_sdss_table()
    c.a100_unwise_table()
