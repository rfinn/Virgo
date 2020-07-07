#!/usr/bin/env python

"""
GOAL:
* read in the virgo filament tables 
* this can serve as a base class to other programs

"""
from astropy.table import Table

class readtables():
    '''
    Read in the Virgo filament tables.

    Args:

    Attributes:

    Returns:

    '''
    def __init__(self,tabledir='~/research/Virgo/tables-north/v0/',tableprefix='vf_north_v0_'):
        self.tabledir = tabledir
        self.tableprefix = tableprefix

    def read_a_table(self,suffix):
        try:
            return Table.read(self.tabledir+self.tableprefix+suffix)
        except:
            print('trouble reading table ',suffix)

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
        
    def read_maintable(self):
        '''Read in the main table.'''
        self.main = self.read_a_table('main.fits')
        
    def read_a100(self):        
        self.a100 = self.read_a_table('a100.fits')
        self.a100sdss = self.read_a_table('a100_sdssphot.fits')
        self.a100unwise = self.read_a_table('a100_unwise.fits')                
        
    def read_co(self):
        self.co = self.read_a_table('co.fits')

    def read_groups(self):
        self.groups = self.read_a_table('groups.fits')
            
    def read_ha(self):
        self.ha = self.read_a_table('ha.fits')

    def read_hyperleda(self):
        self.hl = self.read_a_table('hyperleda.fits')

    def read_environment(self):
        self.env = self.read_a_table('main_env_prop_H0_74_0.fits')
            
    def read_filament_members(self):
        self.filmemb = self.read_a_table('filament_membership.fits')

    def read_unwise(self):
        self.unwise = self.read_a_table('main_unwise.fits')

    def read_mstar(self):
        self.mstar = self.read_a_table('nsa_bellmasses.fits')
        self.mstar0 = self.read_a_table('nsa_v0_bellmasses.fits')
        
    def read_nsa(self):
        self.nsa = self.read_a_table('nsa.fits')

    def read_nsa0(self):
        self.nsa0 = self.read_a_table('nsa_v0.fits')
            
    def read_sfr(self):
        self.sfr = self.read_a_table('sfr.fits')

    def read_steer(self):
        self.steer = self.read_a_table('steer17.fits')

    def read_z0mgs(self):
        self.steer = self.read_a_table('steer17.fits')

class readCOfiles():
    def __init__(self,tabledir='~/research/Virgo/tables-north/v0/',tableprefix='galaxy_sample_prop_'):
        self.tabledir = tabledir
        self.tableprefix = tableprefix

    def read_a_table(self,suffix):
        try:
            return Table.read(self.tabledir+self.tableprefix+suffix)
        except:
            print('trouble reading table ',suffix)
    def readall(self):
        self.read_CO()
        self.read_general()
        self.read_HI()
        self.read_IRAM()
        
    def read_C0(self):
        self.co = self.read_a_table('COliterature.fits')
    def read_general(self):
        self.gen = self.read_a_table('general.fits')
    def read_HI(self):
        self.HI = self.read_a_table('HI.fits')
    def read_IRAM(self):
        self.iram = self.read_a_table('IRAM30m.fits')
        
if __name__ == '__main__':

    vf = readtables(tabledir='~/research/Virgo/tables-north/v0/',tableprefix='vf_north_v0_')
    vf.readall()

    co = readCOfiles(tabledir='~/research/Virgo/tables-north/v0/',tableprefix='galaxy_sample_prop_')
