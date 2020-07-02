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
import argparse

from astropy.io import fits
from astropy.table import Table, join, hstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.ned import Ned

from matplotlib import pyplot as plt
homedir = os.getenv("HOME")
#sys.path.append(homedir+'/github/appss/')
#from join_catalogs import make_new_cats, join_cats

from virgoCommon import *



masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample_wNEDname.fits'
outdir = homedir+'/research/Virgo/tables/v0/'
file_root = 'vf_v0'


parser = argparse.ArgumentParser(description ='write out subtables for virgo filaments catalog')

#parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
parser.add_argument('--north',dest = 'north', action='store_true',help='keep DEC > -1 galaxies')
     
args = parser.parse_args()


# keep DEC > -1 galaxies only
if (args.north):
    NORTH_ONLY = True
else:
    NORTH_ONLY = False
if NORTH_ONLY:
    outdir = homedir+'/research/Virgo/tables-north/v0/'
    file_root = 'vf_north_v0_'

def duplicates(table,column,flag=None):
    if flag is None:
        unique, counts = np.unique(table[column], return_counts=True)
    elif flag is not None:
        unique, counts = np.unique(table[column][flag], return_counts=True)
    print('number of duplicates = ',sum(counts > 1))
    #print('duplicates = ',unique[counts > 1])
    return unique, counts

class catalog:
    def __init__(self,catalog):
        self.cat = Table(fits.getdata(catalog))

        if NORTH_ONLY:
            self.cat, self.keepnorth_flag = self.keep_north()
        self.basictable = self.cat['VFID','RA','DEC','NEDname']
        self.maintable = self.cat['VFID','RA','DEC','vr','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag']
        self.maintable.rename_column('NSAID_2','NSAIDV0')
        self.maintable.rename_column('NSA0flag','NSAV0flag')        
        self.catcoord = SkyCoord(self.cat['RA'],self.cat['DEC'],frame='icrs',unit='deg')
    def get_unwise(self):
        self.unwise = Table.read(outdir+'vf_north_v0_main_unwise.fits')
        self.unwiseFlag = (self.unwise['x'] > 0)
    def keep_north(self, cat=None):
        # cut the catalog to keep Dec > -1 only
        if cat is None:
            cat = self.cat
        else:
            cat = cat
        keepnorth = cat['DEC'] > -1.3
        return cat[keepnorth], keepnorth
    def runall(self):
        self.get_unwise()
        self.get_z0MGS_flag()
        self.get_CO()
        self.get_halpha()        
        self.get_steer17()
        self.get_radius()
        self.main_table()
        self.print_stats()
        self.hyperleda_table()
        self.nsa_table()
        self.nsa_v0_table()        
        self.a100_table()
        self.a100_sdss_table()
        self.a100_unwise_table()
        self.get_size_for_JM()
        pass
    def print_stats(self):
        print('Number in sample = ',len(self.cat))
        print('Number with CO data = ',sum(self.coflag))
        print('Number with A100 data = %i (%.3f)'%(sum(self.cat['A100flag']),sum(self.cat['A100flag'])/len(self.cat['A100flag'])))
        print('Number with z0MGS matches = %i (%.3f)'%(sum(self.z0mgsFlag),sum(self.z0mgsFlag)/len(self.z0mgsFlag)))
        print('Number with steer17 matches = %i (%.3f)'%(sum(self.steerFlag),sum(self.steerFlag)/len(self.steerFlag)))
        f = self.unwiseFlag
        print('Number with unwise matches = %i (%.3f)'%(sum(f),sum(f)/len(f)))        
        f = self.unwiseFlag & (self.cat['NSAflag'] | self.cat['NSA0flag'])
        print("Number with unWISE and NSA = %i (%.3f)"%(sum(f),sum(f)/len(f)))

        print("CO SOURCES")
        nco = sum(self.coflag)
        f = self.coflag & self.z0mgsFlag
        print("\tNumber of CO sources in z0MGS = %i (%.2f)"%(sum(f),sum(f)/nco))
        f = self.coflag & self.steerFlag
        print("\tNumber of CO sources in Steer = %i (%.2f)"%(sum(f),sum(f)/nco))
        f = self.coflag & self.z0mgsFlag & self.steerFlag
        print("\tNumber of CO sources in z0MGS+Steer = %i (%.2f)"%(sum(f),sum(f)/nco))
        f = self.coflag & self.unwiseFlag & (self.cat['NSAflag'] | self.cat['NSA0flag'])
        print("\tNumber of CO sources with unWISE and NSA = %i (%.2f)"%(sum(f),sum(f)/nco))


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
    def get_radius(self):
        ## INCLUDE A RADIUS FOR EACH GALAXY
        ## USE D25 FROM HYPERLEDA IF IT EXISTS
        ## OTHERWISE USE NSA PETRO 90 IF IT EXISTS
        ## OTHERWISE USE AGC SDSS PETRO 90

        # galaxy has a valid value for D25 if the error is not equal to zero
        # and it is not equal to nan
        d25flag = (self.cat['e_logd25'] != 0) & (self.cat['e_logd25'] == self.cat['e_logd25'])
        print('number of galaxies with D25 measurement = {:d} ({:.2f}%)'.format(sum(d25flag),sum(d25flag)/len(d25flag)))
        self.radius = np.zeros(len(self.cat),'f')
        # convert to arcseconds
        # logD25 is in units of 0.1 arcmin
        # so add one, then raise to 10**,
        # then multiply by 60 to go from arcmin to arcsec
        # then divide by 2 to get radius
        self.radius[d25flag] = pow(10,(self.cat['logd25'][d25flag]-1))*60/2

        ## IF GAL DOESN'T HAVE D25
        ## THEN CHECK NSA V0 - AND USE PETROTH90 IF AVAILABLE
        ## NSA v1 = PETRO_TH90
        ## NSA v= = PETROTH90
        ## scale petro radius by 2 to get something closer to D25
        flag = ~d25flag & self.cat['NSAflag']
        self.radius[flag] = self.cat['PETRO_TH90'][flag]*1.3
        print('number of galaxies using NSA V1 Petro TH90 = {:d} ({:.2f}%)'.format(sum(flag),sum(flag)/len(flag)))
        flag = ~d25flag & ~self.cat['NSAflag'] & self.cat['NSA0flag']
        self.radius[flag] = self.cat['PETROTH90'][flag]*1.3
        print('number of galaxies using NSA V1 Petro TH90 = {:d} ({:.2f}%)'.format(sum(flag),sum(flag)/len(flag)))
        ## IF NOT IN HL OR NSA OR NSAV0,
        ## AND IF IN A100, THEN USE PETRO90 SIZE
        ## petroR90_r
        ## give an extra boost b/c they seem smaller than 
        flag = ~d25flag & ~self.cat['NSAflag'] & ~self.cat['NSA0flag'] & self.cat['A100flag'] & (self.cat['parentID'] != 0)& (self.cat['parentID'] != 999999)
        a100radius_flag = flag
        self.radius[flag] = self.cat['petroR90_r'][flag]*1.4
        print('number of galaxies using A100 sdss Petro TH90 = {:d} ({:.2f}%)'.format(sum(flag),sum(flag)/len(flag)))
        ## IF NONE OF THESE SIZES AVAILABLE, SET RADIUS TO 100 ARCSEC
        flag = ~d25flag & ~self.cat['NSAflag'] & ~self.cat['NSA0flag'] & ~a100radius_flag
        self.radius[flag] = 100*np.ones(sum(flag),'f')
        print('number of galaxies with no size measurement = {:d} ({:.2f}%)'.format(sum(flag),sum(flag)/len(flag)))
        radius_flag = ~flag
        c = Column(self.radius,name='radius',unit='arcsec')
        self.cat.add_column(c)
        c = Column(radius_flag,name='radius_flag',description='if False, then rad is set to 100 arcsec')
        self.cat.add_column(c)

        
    def main_table(self):
        # write out
        # ra, dec, velocity, HL, NSA id, A100 id
        # NED name
        # make flags to denote if galaxy is in:
        # - CO sample
        colnames = ['VFID','RA','DEC','vr','radius','radius_flag','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag']
        #self.maintable = self.cat['VFID','RA','DEC','vr','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag']
        
        self.maintable = self.cat[colnames]
        self.maintable.rename_column('NSAID_2','NSAIDV0')
        self.maintable.rename_column('NSA0flag','NSAV0flag')        


        c1 = Column(self.coflag,name='COflag')
        c2 = Column(self.z0mgsFlag,name='Z0MGSflag')
        c3 = Column(self.steerFlag,name='Steerflag')
        c4 = Column(self.unwiseFlag,name='unwiseflag')        
        self.maintable.add_columns([c1,c2,c3,c4])
        # - 2MASS
        # - z0MGS
        # - unWISE
        self.maintable.write(outdir+file_root+'main.fits',format='fits',overwrite=True)

    def get_CO(self,match_by_coords=False,match_by_name=True):
        # read in CO mastertable
        # match to main table
        cofile = homedir+'/github/Virgo/tables/CO-MasterFile-2018Feb16.fits'
        # this file has a new column with the exact NED names
        # can use this to match to mastertable NEDname column
        cofile = homedir+'/research/Virgo/tables/CO-MasterFile-2018Feb16-fixedNEDnames.fits'        
        self.co = Table(fits.getdata(cofile))

        # replace NED_name of SHOC206b with NED name for SHOC 206a
        # this is really the same galaxy, just has less common id in CO file
        # set to MCG +08-16-005;
        
        # fix case for UGC09348, until we run the full query again...
        flag = self.co['NED_name'] == 'SHOC206b'
        if sum(flag > 0):
            self.co['NEDname'][flag] = 'MCG +08-16-005'
        
        

        cocoord = SkyCoord(self.co['RA'],self.co['DEC'],unit='deg',frame='icrs')
        if match_by_coords:
            # match co coords to mastertable 
            idx, d2d, d3d = self.catcoord.match_to_catalog_sky(cocoord)
            self.d2d = d2d
            self.idx = idx
            self.coflag = d2d < 15./3600*u.deg
            # create a new, blank table with same # of lines as mastertable
            # but with columns like the co table
            newco = Table(np.zeros(len(self.basictable),dtype=self.co.dtype))
            # add co information into the new table for the galaxies with
            # a match to the CO sample
            # NOTE: this could match multiple galaxies to the same CO source
            newco[self.coflag] = self.co[idx[self.coflag]]

            # join basic table and co table

            self.cotable = hstack([self.basictable,newco])

        if match_by_name:
            #self.co.rename_column('NED_name','NEDname')
            #np.searchsorted(names1,names2)
            
            # match the basictable and the CO table by matching
            # entries by the NEDname colums
            self.cotable = myjoinleft(self.basictable,self.co,keys='NEDname')

            self.coflag = ~self.cotable['CO'].mask

            # also check to see which CO sources were not matched
            self.testtable = join(self.co,self.basictable,keys='NEDname',join_type='left')
            #self.coflag = len(self.co['CO']) > 0
            try:
                self.comatchflag = ~self.testtable['VFID'].mask
                print('CO sources with no match in mastertable:')
                print(self.testtable['NEDname','NED_name'][~self.comatchflag])
                ## plot the positions of CO galaxies that weren't matched to mastertable
                plt.figure()
                plt.plot(self.testtable['RA_1'][~self.comatchflag],self.testtable['DEC_1'][~self.comatchflag],'bo')


            except AttributeError:
                print('all CO sources have been matched. CONGRATULATIONS!!!!!!!')
                self.comatchflag = np.ones(len(self.testtable))
        ## print the CO galaxies with no matches in the mastertable 
        print('number of galaxies with CO matches = ',sum(self.coflag))

        ## look for CO galaxies that were matched to multiple galaxies in the
        ## mastertable
        unique, counts = duplicates(self.cotable,'NED_name')
        print("CO sources that are matched to multiple galaxies in the mastertable:")
        print(unique[counts>1])
        
        
        self.cotable.add_column(Column(self.coflag),name='COflag')
        self.cotable.write(outdir+file_root+'co.fits',format='fits',overwrite=True)
            

        # print CO sources that are not in the table
    def get_halpha(self,halphafile=None):
        # read in Halpha observing summary file
        if halphafile is None:
            #infile = '/home/rfinn/research/Virgo/Halpha/observing-summary-Halpha-latest.csv'
            infile = '/home/rfinn/research/Virgo/Halpha/observing-summary-Halpha-clean-04Jun2020.csv'            
        else:
            infile = halphafile
        self.ha = Table.read(infile,format='csv')
        self.ha.rename_column('NSA ID','NSAIDV0')

        # check for duplicates in the hafile
        unique, counts = duplicates(self.ha,'NSAIDV0')
        print("Halpha sources that are listed multiple times in the halpha file:")
        print(unique[counts>1])
        
        
        # match ha file to the base table using the NSA v0
        # match the basictable and the CO table by matching
        # entries by the NEDname colums
        #print(self.ha.colnames)
        #print(self.maintable.colnames)
        # note myjoinleft preserves the original ordering in the tables
        # rather than having the joined table sorted by according to the match key

        self.hatable = myjoinleft(self.maintable,self.ha,keys='NSAIDV0')
        print('Halpha table lengths')
        print(len(self.maintable),len(self.ha),len(self.hatable))        
        self.haflag = ~self.hatable['Date Obs'].mask
        self.hatable.add_column(Column(self.haflag),name='haflag')
        self.hatable.write(outdir+file_root+'ha.fits',format='fits',overwrite=True)
        
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
        #cat = Table.read('/home/rfinn/research/Virgo/tables/vf-z0MGS.tbl',format='ipac')
        cat = Table.read('/home/rfinn/research/Virgo/tables/vf_z0mgs_30arcsec_051920.tbl',format='ipac')
        # cut on declination
        if NORTH_ONLY:
            cat = cat[self.keepnorth_flag]
        # create a flag for mastertable
        self.z0mgsFlag = ~cat['pgc_name'].mask
        self.z0mgs_cat = cat
        c = Column(self.z0mgsFlag,name='Z0MGSflag')
        cat.add_column(c)
        # write out north file
        cat.write(outdir+file_root+'z0mgs.fits',format='fits',overwrite=True)
    def get_steer17(self):
        # match to GL's steer catalog
        steercat = '/home/rfinn/research/Virgo/ancil-tables/Steer2017_cat_Virgo_field_H0_74_0.fits'
        self.steer = Table(fits.getdata(steercat))
        # GL suggests using
        # np.searchsorted(names1,names2)

        # I might otherwise do a loop
        # for each name in our catalog, look for match in steer catalog
        # probably astropy has a way to do this (topcat certainly does)

        self.basic_with_steer = myjoinleft(self.basictable,self.steer,keys='NEDname')
        self.steerFlag = ~self.basic_with_steer['Dmedian'].mask
        self.basic_with_steer.add_column(Column(self.steerFlag),name='Steerflag')
        outfile = outdir+file_root+'steer17.fits'            
        self.basic_with_steer.write(outfile, format='fits',overwrite=True)

    def hyperleda_table(self):
        colnames = self.cat.colnames[0:43]
        
        self.write_table(colnames,outdir+file_root+'hyperleda.fits',format='fits')
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
        self.write_table(colnames,outdir+file_root+'nsa.fits',format='fits',names=newcolnames)
    def nsa_v0_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[204:348]
        # RA and DEC get changed to RA_2, DEC_2 during the matching process
        # changing them back to native NSA values
        newcolnames = colnames.copy()
        for i,n in enumerate(newcolnames):
            if n.find('_2') > -1:
                newcolnames[i] = n.strip('_2')
                
        newcolnames[2] = 'RA'
        newcolnames[3] = 'DEC'
        self.write_table(colnames,outdir+file_root+'nsa_v0.fits',format='fits',names=newcolnames)
    def a100_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[352:377]
        self.write_table(colnames,outdir+file_root+'a100.fits',format='fits')
    def a100_sdss_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[377:472]
        self.write_table(colnames,outdir+file_root+'a100_sdssphot.fits',format='fits')
    def a100_unwise_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[472:537]
        # RA and DEC get changed to RA_2, DEC_2 during the matching process
        # changing them back to native NSA values
        newcolnames = colnames.copy()
        newcolnames[1] = 'ra'
        newcolnames[2] = 'dec'
        #print(colnames)
        self.write_table(colnames,outdir+file_root+'a100_unwise.fits',format='fits',names=newcolnames)
    def write_table(self,colnames,outfile,format=None,names=None):
        if format is None:
            format = 'fits'
        else:
            format = format
        mycolumns = []
        for c in colnames:
            mycolumns.append(self.cat[c])
        if names is not None:
            subtable = Table(mycolumns,names=names)
        else:
            subtable = Table(mycolumns)
        newtable = hstack([self.basictable,subtable])
        newtable.write(outfile,format=format,overwrite=True)
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
        c1 = Column(self.twomassflag,name='2MASS')
        # add flag for z0MGS

        # add flag for Steer+17

        # add flag for legacy survey
        
        newtable.write(outfile,format=format,overwrite=True)
    def get_size_for_JM(self):
        # need to get John the RA, DEC, and size estimate for each galaxy
        # so he can run the legacy code.

        # use the hyperleda size if available
        # but need to convert from kpc to an angular size
        #
        # otherwise use the 
        pass
    def table_for_halphagui(self):
        # need RA, DEC, redshift, radius
        pass
if __name__ == '__main__':
    c = catalog(masterfile)
    c.runall()
