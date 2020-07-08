#!/usr/bin/env python

"""
GOAL:
* read in the virgo filament tables 
* this can serve as a base class to other programs

"""
from astropy.table import Table

class readfulltables():
    '''
    Read in the Virgo filament tables.

    Args: 
    - tabledir = directory to virgo filament tables
    - tableprefix = prefix of the tables, e.g. vf_north_v0

    Returns:
    - class instance, with tables appended.
    - table names are listed in the methods below. 
    '''
    def __init__(self,tabledir='/home/rfinn/research/Virgo/tables-north/v0/',tableprefix='vf_north_v0_'):
        self.tabledir = tabledir
        self.tableprefix = tableprefix
        self.readall()
    def read_a_table(self,suffix):
        try:
            return Table.read(self.tabledir+self.tableprefix+suffix)
        except:
            print('trouble reading table ',suffix)
        #t= Table.read(self.tabledir+self.tableprefix+suffix)
        #return t        
    def read_maintable(self):
        '''Read in the main table as self.main.'''
        self.main = self.read_a_table('main.fits')
        
    def read_a100(self):
        '''Read in the a100 tables as self.a100, self.a100sdss,self.a100unwise.'''        
        self.a100 = self.read_a_table('a100.fits')
        self.a100sdss = self.read_a_table('a100_sdssphot.fits')
        self.a100unwise = self.read_a_table('a100_unwise.fits')                
        
    def read_co(self):
        '''Read in the CO mastertable table as self.co.'''                
        self.co = self.read_a_table('co.fits')

    def read_groups(self):
        '''Read in the JM group catalog table as self.groups.'''                
        self.groups = self.read_a_table('groups.fits')
            
    def read_ha(self):
        '''Read in the H-alpha table as self.ha.'''                        
        self.ha = self.read_a_table('ha.fits')

    def read_hyperleda(self):
        '''Read in the Hyperleda table as self.hl.'''        
        self.hl = self.read_a_table('hyperleda.fits')

    def read_environment(self):
        '''Read in the Gianluca's environment table as self.env.'''                
        self.env = self.read_a_table('main_env_prop_H0_74_0.fits')
            
    def read_filament_members(self):
        '''Read in the Gianluca's filament member table as self.filmemb.'''                        
        self.filmemb = self.read_a_table('main_filament_membership.fits')

    def read_unwise(self):
        '''Read in the unwise table as self.unwise.'''         
        self.unwise = self.read_a_table('main_unwise.fits')

    def read_mstar(self):
        '''Read in Greg's bell stellar mass tables as self.mstar and self.mstar0.''' 
        self.mstar = self.read_a_table('nsa_bellmasses.fits')
        self.mstar0 = self.read_a_table('nsa_v0_bellmasses.fits')
        
    def read_nsa(self):
        '''Read in nsa v1 table as self.nsa.'''  
        self.nsa = self.read_a_table('nsa.fits')

    def read_nsa0(self):
        '''Read in nsa v0 table as self.nsa0.'''            
        self.nsa0 = self.read_a_table('nsa_v0.fits')
            
    def read_sfr(self):
        '''Read in NUV+WISE Kennicutt&Evans sfr table as self.sfr.'''                
        self.sfr = self.read_a_table('sfr.fits')

    def read_steer(self):
        '''Read in Gianluca's Steer+17 table as self.steer.'''                
        self.steer = self.read_a_table('steer17.fits')

    def read_z0mgs(self):
        '''Read in Leroy+2019 z0MGS table as self.z0mgs.'''                
        self.z0mgs = self.read_a_table('z0mgs.fits')
    def readall(self):
        '''Read in all tables.'''
        self.read_maintable()
        self.read_a100()
        self.read_co()
        self.read_groups()
        self.read_ha()
        self.read_hyperleda()
        self.read_environment()
        self.read_filament_members()
        self.read_unwise()
        self.read_mstar()
        self.read_nsa()
        self.read_nsa0()
        self.read_sfr()
        self.read_steer()
        self.read_z0mgs()

class readCOtables():
    '''
    Read in the Virgo tables for CO sample.

    Args: 
    - tabledir = directory to Gianluca's CO tables
    - tableprefix = prefix of the tables, e.g. galaxy_sample_prop

    Returns:
    - class instance, with tables appended.
    - table names are listed in the methods below. 
    '''
    def __init__(self,tabledir='/home/rfinn/research/Virgo/tables-north/v0/',tableprefix='galaxy_sample_prop_'):
        self.tabledir = tabledir
        self.tableprefix = tableprefix
        self.readallco()
    def read_co_table(self,suffix):
        try:
            return Table.read(self.tabledir+self.tableprefix+suffix)
        except:
            print('trouble reading table ',suffix)
        
    def read_CO(self):
        '''Read in CO table as self.co.'''                        
        self.co = self.read_co_table('COliterature.fits')
    def read_general(self):
        '''Read in general table as self.gen.'''                                
        self.gen = self.read_co_table('general.fits')
    def read_HI(self):
        '''Read in HI table as self.HI.'''                                
        self.HI = self.read_co_table('HI.fits')
    def read_IRAM(self):
        '''Read in IRAM30m table as self.iram.'''                                
        self.iram = self.read_co_table('IRAM30m.fits')
    def readallco(self):
        '''Read all tables.'''
        self.read_CO()
        self.read_general()
        self.read_HI()
        self.read_IRAM()
        
if __name__ == '__main__':

    vf = readfulltables(tabledir='/home/rfinn/research/Virgo/tables-north/v0/',tableprefix='vf_north_v0_')
    vf.readall()

    co = readCOfiles(tabledir='/home/rfinn/research/Virgo/tables-north/v0/',tableprefix='galaxy_sample_prop_')
