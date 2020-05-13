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
import os
import numpy as np
import time

from astropy.io import fits
from astropy.table import Table, join, hstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.ned import Ned

from matplotlib import pyplot as plt
homedir = os.getenv("HOME")
#sys.path.append(homedir+'/github/appss/')
#from join_catalogs import make_new_cats, join_cats


masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample.fits'
outdir = homedir+'/research/Virgo/tables/'



# keep DEC > -1 galaxies only
NORTH_ONLY = False
if NORTH_ONLY:
    outdir = homedir+'/research/Virgo/tables-north/'
class catalog:
    def __init__(self,catalog):
        self.cat = Table(fits.getdata(catalog))

        if NORTH_ONLY:
            self.cat, self.keepnorth_flag = self.keep_north()
        self.basictable = self.cat['VFID','RA','DEC','NEDname']
        self.maintable = self.cat['VFID','RA','DEC','vr','objname','NSAID','AGC','NEDname','HLflag','NSAflag','A100flag']
        self.catcoord = SkyCoord(self.cat['RA'],self.cat['DEC'],frame='icrs',unit='deg')
        
    def keep_north(self, cat=None):
        # cut the catalog to keep Dec > -1 only
        if cat is None:
            cat = self.cat
        else:
            cat = cat
        keepnorth = cat['DEC'] > -1.3
        return cat[keepnorth], keepnorth
    
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
        self.maintable = self.cat[colnames]
        # make flags to denote if galaxy is in:
        # - CO sample
        self.get_CO()
        # - 2MASS
        # - z0MGS
        # - unWISE
        self.write_table(colnames,'vf_main.fits',format='fits')

    def get_CO(self,match_by_coords=False,match_by_name=True):
        # read in CO mastertable
        # match to main table
        cofile = homedir+'/github/Virgo/tables/CO-MasterFile-2018Feb16.fits'
        # this file has a new column with the exact NED names
        # can use this to match to mastertable NEDname column
        cofile = homedir+'/research/Virgo/tables/CO-MasterFile-2018Feb16-fixedNEDnames.fits'        
        self.co = Table(fits.getdata(cofile))


        

        cocoord = SkyCoord(self.co['RA'],self.co['DEC'],unit='deg',frame='icrs')
        if match_by_coords:
            # match co coords to mastertable 
            idx, d2d, d3d = self.catcoord.match_to_catalog_sky(cocoord)
            self.d2d = d2d
            self.idx = idx
            self.coflag = d2d < 15./3600*u.deg
            # match to 
            newco = Table(np.zeros(len(self.basictable),dtype=self.co.dtype))
            newco[self.coflag] = self.co[idx[self.coflag]]

            # join basic table and co

            self.cotable = hstack([self.basictable,newco])

        if match_by_name:
            #self.co.rename_column('NED_name','NEDname')
            #np.searchsorted(names1,names2)
            self.cotable = join(self.basictable,self.co,keys='NEDname',join_type='left')
            #self.coflag = len(self.co['CO']) > 0
            self.coflag = ~self.cotable['CO'].mask

            # also check to see which CO sources were not matched
            self.testtable = join(self.co,self.basictable,keys='NEDname',join_type='left')
            #self.coflag = len(self.co['CO']) > 0
            self.comatchflag = ~self.testtable['VFID'].mask
            print('CO sources with no match in mastertable:')
            print(self.testtable['NEDname','NED_name'][~self.comatchflag])
        plt.figure()
        plt.plot(self.testtable['RA_1'][~self.comatchflag],self.testtable['DEC_1'][~self.comatchflag],'bo')
        print('number of galaxies with CO matches = ',sum(self.coflag))
        self.cotable.add_column(Column(self.coflag),name='COFlag')
        self.cotable.write(outdir+'vf_co.fits',format='fits',overwrite=True)

        # print CO sources that are not in the table
        
    def get_2massflag(self,twomassfile=None):
        if twomassfile is None:
            print('need to provide the twomass file name')
            return
        else:
            twomass = ascii.read(twomassfile,format='ipac')
        # cut on declination


        #create a flag for mastertable
        self.twomassflag = np.zeros(len(self.cat),'bool')
        self.twomassflag[twomass['cntr_01']] = np.ones(len(twomass),'bool')

        # write out north version of file
    def get_z0MGS_flag(self,mgsfile=None):
        # get z0MGS

        # cut on declination

        # create a flag for mastertable

        # write out north file
        pass
    def match_steer17(self):
        # match to GL's steer catalog
        steercat = '/home/rfinn/research/Virgo/ancil-tables/Steer2017_cat_Virgo_field_H0_74_0.fits'
        self.steer = Table(fits.getdata(steercat))
        # GL suggests using
        # np.searchsorted(names1,names2)

        # I might otherwise do a loop
        # for each name in our catalog, look for match in steer catalog
        # probably astropy has a way to do this (topcat certainly does)

        self.basic_with_steer = join(basictable,self.steer,keys='NEDname',join_type='left')
        self.steerflag = self.basic_with_steer['Dmedian'] > 0.
        self.basic_with_steer.add_column(Column(self.steerflag),name='steerFlag')
        if NORTH_ONLY:
            outfile = outdir+'/vf_steer17_north.fits'
        else:
            outfile = outdir+'/vf_steer17.fits'            
        self.basic_with_steer.write(outfile, format='fits',overwrite=True)

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
    def write_main_table(self,colnames,outfile,format=None,names=None):
        if format is None:
            format = 'fits'
        else:
            format = format
        if names is not None:
            newtable = Table(self.cat[colnames],names=names)
        else:
            newtable = Table(self.cat[colnames])

        # get flag for CO sample
        
        # add flag for 2mass
        self.get_2massflag()
        c = Column(self.twomassflag,name='2MASS')
        # add flag for z0MGS

        # add flag for Steer+17

        # add flag for legacy survey
        
        newtable.write(outdir+outfile,format=format,overwrite=True)

if __name__ == '__main__':
    c = catalog(masterfile)
    #c.keep_north()
    #c.catalog_for_z0MGS()
    #c.main_table()
    #c.hyperleda_table()
    #c.nsa_table()
    #c.a100_table()
    #c.a100_sdss_table()
    #c.a100_unwise_table()
