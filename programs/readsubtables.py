#!/usr/bin/env python

'''
program to read in subtables for virgo

'''
import os
from astropy.table import Table,hstack

homedir = os.getenv("HOME")

class vtables:
    def __init__(self,tabledir,tableprefix):
        '''class containing all tables for virgo filament project  '''
        self.tabledir = tabledir
        self.tableprefix = tableprefix        
        pass
    def read_all(self):
        self.read_main()
        self.read_hyperleda()
        self.read_nsav0()
        self.read_nsav1()
        self.read_steer17()
        self.read_z0mgs()
        self.read_ned()
        self.read_env()
        self.read_a100()
        self.read_agc()
        self.read_co()
        self.read_halpha()                
    def read_main(self):
        ''' read in main table; store as self.main  '''
        self.main = Table.read(self.tabledir+self.tableprefix+'main.fits')
        pass
    def read_hyperleda(self):
        ''' read in Hyperleda table; store as self.hl  '''
        self.hl = Table.read(self.tabledir+self.tableprefix+'hyperleda.fits')                              
        pass
    def read_nsav0(self):
        ''' read in NSA v0 table; store as self.nsav0  '''
        self.nsav0 = Table.read(self.tabledir+self.tableprefix+'nsa_v0.fits')                              
        pass
    def read_nsav1(self):
        ''' read in NSA v1 table; store as self.nsav1  '''
        self.nsav1 = Table.read(self.tabledir+self.tableprefix+'nsa.fits')
        pass
    def read_steer17(self):
        ''' read in Steer+17 table; store as self.steer  '''
        self.steer = Table.read(self.tabledir+self.tableprefix+'steer17.fits')                               
        pass
    def read_z0mgs(self):
        ''' read in Z0MGS table (Leroy+2019); store as self.z0mgs  '''
        self.z0mgs = Table.read(self.tabledir+self.tableprefix+'z0mgs.fits')
        pass
    def read_ned(self):
        ''' read in NED query table; store as self.ned  '''
        self.ned = Table.read(self.tabledir+self.tableprefix+'nedquery.fits')                      
        pass
    def read_env(self):
        ''' read in GC's env and BV envsummary table; store as self.env  '''
        tab1 = Table.read(self.tabledir+self.tableprefix+'main_env_prop_H0_74_0.fits')
        tab2 = Table.read(self.tabledir+self.tableprefix+'main_envsummary.fits')
        self.env = hstack([tab1,tab2])
        pass
                              
    def read_a100(self):
        ''' read in ALFALFA 100 table; store as self.a100, self.a100sdss, self.a100unwise  '''
        self.a100 = Table.read(self.tabledir+self.tableprefix+'a100.fits')
        self.a100sdss = Table.read(self.tabledir+self.tableprefix+'a100_sdssphot.fits')
        self.a100unwise = Table.read(self.tabledir+self.tableprefix+'a100_unwise.fits')
        # calculate HI deficiency
        # use the R90 petro
        pass
    def read_agc(self):
        ''' read in AGC table; store as self.agc  '''
        self.agc = Table.read(self.tabledir+self.tableprefix+'agc.fits')                               
        pass
    def read_co(self):
        ''' read in CO table; store as self.co  '''
        self.co = Table.read(self.tabledir+self.tableprefix+'co.fits')                               
        pass
    def read_halpha(self):
        ''' read in halpha observations table; store as self.ha  '''
        self.ha = Table.read(self.tabledir+self.tableprefix+'ha.fits')                               
        pass
    def read_tempel(self):
        ''' read in Tempel table; store as self.te  '''
        tab1 = Table.read(self.tabledir+self.tableprefix+'Tempelgroups_infos.fits')
        tab2 = Table.read(self.tabledir+self.tableprefix+'Tempel_groupsinfo.fits')
        self.tempel = hstack([tab1,tab2])
        # should merge these, and maybe remove the columns from main

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v1/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_north_v1_', help = 'prefix for tables; default is vf_north_v1')                               
    args = parser.parse_args()
    
    v = vtables(args.tabledir,args.tableprefix) 
    v.read_all()
