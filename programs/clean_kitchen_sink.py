#!/usr/bin/env python

'''
GOAL:
* read in kitchen_sink
* read in results virgo_check_sample_by_eye
* edit kitchen sink to
  - remove bad sources
  - merge shredded galaxies

'''

from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
### INPUT FILES
kitchen_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink.fits'
byeye_classifications = '/home/rfinn/research/Virgo/supersample/virgo_check_sample_by_eye.csv'


###  CATALOG CLASS
class catalog:
    def __init__(self,kitchen_sink,byeye_classifications):
        self.kitchen = fits.getdata(kitchen_sink)
        self.kitchen = Table(self.kitchen)
        self.byeye = ascii.read(byeye_classifications, delimiter=',')

    def cut_catalog(self):
        self.cutflag = (self.byeye['class'] == 2) | (self.byeye['class'] == 4) | (self.byeye['class'] == 0)

        self.cleancat = self.byeye[~self.cutflag]
        #remove first column, which is a duplicate with second column
        n = self.cleancat.colnames
        self.cleancat.remove_column(n[0])
        self.cleancat.write('clean_sample.fits',format='fits',overwrite=True)
    def catalog_for_z0MGS(self):
        '''
        need RA, DEC, and search radius
        in IPAC format table
        for matching with leroy+2019 Galaxy Synthesis WISE+GALEX
        table is served by IRSA
        '''
        search_radius = 10.*np.ones(len(self.cleancat)) # units are arcsec
        newtable = Table([self.cleancat['galnumber'],self.cleancat['RA'],self.cleancat['DEC'],search_radius],names=['galid','ra','dec','major'])
        newtable.write('clean_sample.txt',format='ipac',overwrite=True)
if __name__ == '__main__':
    c = catalog(kitchen_sink,byeye_classifications)
    c.cut_catalog()
    c.catalog_for_z0MGS()
