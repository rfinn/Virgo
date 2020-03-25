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
from astropy.table import Table, join, hstack, Column, MaskedColumn 
from astropy.coordinates import SkyCoord
import numpy as np
import sys
import os

homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/APPSS/')
from join_catalogs import make_new_cats, join_cats

### INPUT FILES
kitchen_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink.fits'
byeye_classifications = '/home/rfinn/research/Virgo/supersample/virgo_check_sample_by_eye.csv'


###  CATALOG CLASS
class catalog:
    def __init__(self,kitchen_sink,byeye_classifications):
        self.kitchen = fits.getdata(kitchen_sink)
        self.kitchen = Table(self.kitchen)
        self.byeye = ascii.read(byeye_classifications, delimiter=',')
        #self.byeye['parent'] = np.array(self.byeye['parent'],'i')
    def runall(self):
        self.merge_class4()
        #self.remove_agc_only()
        self.get_radec_first()
        self.cut_catalog()
        self.match_a100()
        self.check_new_a100()
    def merge_class4(self):
        # find class 4 objects
        class4 = (self.byeye['class'] == 4)
        # find parent object
        parents = self.byeye['parent'][class4]
        self.parents = parents
        galid = np.arange(len(self.byeye['parent']))[class4]        
        HLflag = ~(self.byeye['HL'].mask)
        NSAflag = (self.byeye['NSAID'] != 0)
        AGCflag = self.byeye['AGC'] != 0
        
        # figure out which survey data we have for 4, and which is in the parent
        HL4 = HLflag[class4]
        NSA4 = NSAflag[class4]
        AGC4 = AGCflag[class4]
        
        # merge any missing survey information that is associated with 4
        # with parent entry
        
        for i in range(len(HL4)):
            print(galid[i], 'parents = ',parents[i])
            parentflag = self.byeye['galnumber'] == int(float(parents[i]))
            parentid = np.arange(len(class4))[parentflag]            
            #print('len(parentid) = ',len(parentid))
            print(parents[i],parentid,galid[i])
            if len(parentid) == 0:
                print('error matching parent {} for galaxy {}'.format(parents[i],galid[i]))
                
            elif len(parentid) >= 1:

                #print('error matching parent {} for galaxy {}'.format(parents[i],galid[i]))
                # this will occur
                HLp = (~self.byeye['HL'].mask[parentflag])[0]
                NSAp = (self.byeye['NSAID'][parentflag] != 0)[0]
                AGCp = (self.byeye['AGC'][parentflag] != 0)[0]
                #print(HLp,HL4[i],parentid,sum(parentflag))
                if not(HLp) and HL4[i]:
                    self.merge_sources(parentid,galid[i],HL=True)
                    self.kitchen['HLflag'][parentid] = True                    
                if not(NSAp) and NSA4[i]:
                    self.merge_sources(parentid,galid[i],NSA=True)
                    self.kitchen['NSAflag'][parentid] = True                    
                if not(AGCp) and AGC4[i]:
                    self.merge_sources(parentid,galid[i],AGC=True)
                    self.kitchen['AGCflag'][parentid] = True
    def merge_sources(self,parentid,galid,HL=False,NSA=False,AGC=False):
        colnames = self.kitchen.colnames
        if HL:
            survey_columns = np.arange(0,45)
        if AGC:
            survey_columns = np.arange(44,83)
        if NSA:
            survey_columns = np.arange(90,200)
        for i in survey_columns:
            self.kitchen[colnames[i]][parentid] = self.kitchen[colnames[i]][galid]

    def get_radec_first(self):
        # ra is RA HL for those with HL data, then RA NSA 
        # same for DEC and recession velocity
        
        self.ra = self.kitchen['al2000']*15*self.kitchen['HLflag'] + self.kitchen['RA_2']*~self.kitchen['HLflag']
        self.dec = self.kitchen['de2000']*self.kitchen['HLflag'] + self.kitchen['DEC_2']*~self.kitchen['HLflag']
        self.vel = self.kitchen['v']*self.kitchen['HLflag'] + self.kitchen['Z']*3.e5*~self.kitchen['HLflag']        
        # for galaxies with class=16, use the RA and DEC in by-eye file
        class16 = self.byeye['class'] == 16
        self.ra[class16] = self.byeye['RA'][class16]
        self.dec[class16] = self.byeye['DEC'][class16]
        
        # if galaxy is class = 4, use the RA, DEC of parent
        # didn't implement this yet
        
        pass
    def cut_catalog(self):
        self.cutflag = (self.byeye['class'] == 2) | (self.byeye['class'] == 4) | (self.byeye['class'] == 0)

        # removing sources with AGC only
        # will match to A100 instead

        self.agconly = ~self.kitchen['HLflag'] & ~self.kitchen['NSAflag'] & self.kitchen['AGCflag']

        self.cutflag = self.cutflag | self.agconly

        self.cleancat = self.byeye[~self.cutflag]
        #remove first column, which is a duplicate with second column
        n = self.cleancat.colnames
        self.cleancat.remove_column(n[0])
        
        # append original galaxy id to new kitchen sink file
        c = Column(self.byeye['galnumber'],name='galnumber')
        self.kitchen.add_column(c)
        self.clean_kitchen = self.kitchen[~self.cutflag]
        
        # cut ra,dec,vel
        c1 = Column(self.ra[~self.cutflag],'RAtemp')
        c2 = Column(self.dec[~self.cutflag],'DECtemp')
        c3 = Column(self.vel[~self.cutflag],'vrtemp')
        self.clean_kitchen.add_columns([c1,c2,c3])
    def remove_agc_columns(self):
        pass
    def match_a100(self):
        # read in a100
        self.a100 = fits.getdata('/home/rfinn/research/Virgo/tables/a100-sdss-wise-virgo.fits')
        #a100coord = SkyCoord(a100['RAdeg_Use'],a100['DECdeg_Use'],frame='icrs',unit='deg')
        # define clean catalog coords
        #ccoord = SkyCoord(self.ra, self.dec, frame='icrs',unit='deg')

        ###############################################
        ## MATCH A100 TO  HYPERLEDA+NSA
        ###############################################    
        v1 = self.clean_kitchen['vrtemp']
        v2 = self.a100['Vhelio']
        veloffset=300.
        maxoffset = 15.
        hlnsa_2, hlnsa_matchflag, a100_2, a100_matchflag = make_new_cats(self.clean_kitchen, self.a100, RAkey1='RAtemp',DECkey1='DECtemp',RAkey2='RAdeg_Use',DECkey2='DECdeg_Use', velocity1=v1, velocity2=v2, maxveloffset = veloffset,maxoffset=maxoffset)
        # write out joined a100-sdss-nsa catalog
        joined_table2 = hstack([hlnsa_2,a100_2])
        c1 = Column(a100_matchflag,name='A100flag')
        joined_table2.add_column(c1)
        
        # for any new galaxies, set RA and DEC equal to A100 values
        ra = hlnsa_matchflag*joined_table2['RAtemp'] + ~hlnsa_matchflag*joined_table2['RAdeg_Use']
        dec = hlnsa_matchflag*joined_table2['DECtemp'] + ~hlnsa_matchflag*joined_table2['DECdeg_Use']
        vel = hlnsa_matchflag*joined_table2['vrtemp'] + ~hlnsa_matchflag*joined_table2['Vhelio']        
        joined_table2.remove_columns(['RAtemp','DECtemp','vrtemp'])
        c3 = Column(ra,name='RA',dtype='f')
        c4 = Column(dec,name='DEC',dtype='f')
        c5 = Column(vel,name='vr',dtype='f')

        joined_table2.add_columns([c3,c4,c5])

        ###############################################
        ## FIX BOOLEAN COLUMNS
        ###############################################
        # boolean columns are getting converted weird
        # when I write and then read the fits table
        '''
        try:
            joined_table2['NSAflag'] = (joined_table2['NSAflag'] == 84)
            joined_table2['HLflag'] = (joined_table2['HLflag'] == 84)
        except KeyError:
            print('trouble in paradise')
        '''
        self.clean_a100 = joined_table2
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)

    def check_new_a100(self):
        # inspect new a100 sources
        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']

        print(self.clean_a100['AGC'][a100flag])
        pass
        
    def get_NEDname(self):
        # look up NED name for each galaxy
        # https://astroquery.readthedocs.io/en/latest/ned/ned.html

        # if in HL, look up HL name
        
        pass
            
    def fix_8822(self):
        # add blue object, this is one with CO detection
        # the following information is from NED
        colid = ['NEDname', 'RA','DEC','vel']
        colvalues  = ['UGC8656 NOTES01', 205.129878,42.993819,2899]

        colnames = np.array(self.clean_kitchen.colnames)
        colindices = []
        # save column number for each dictionary key
        for i in range(len(colid)):
            colindices.append(np.arange(len(colnames))[colnames == colid[i]])
        colindices = np.array(colindices)
        new_row = np.zeros_like(self.clean_kitchen[0])

        for i in range(len(colindices)):
            new_row[colindices[i]] = colvalues[i] 
        self.clean_kitchen.add_row(new_row)
        
        self.cleancat.add_row()
        pass
    
    
    def write_clean_cat(self):
        self.cleancat.write('clean_sample.fits',format='fits',overwrite=True)
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)
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
    def catalog_for_paper(self):
        # write out
        # ra, dec, velocity, HL, NSA id, A100 id
        # NED name
        pass
        
if __name__ == '__main__':
    c = catalog(kitchen_sink,byeye_classifications)
    #c.merge_class4()
    #c.cut_catalog()
    #c.catalog_for_z0MGS()
    c.runall()
