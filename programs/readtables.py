#!/usr/bin/env python

'''
program to read in subtables for virgo



written for version 1 tables
11/17/20 by Rose Finn

'''
import os
from astropy.table import Table,hstack
import numpy as np
homedir = os.getenv("HOME")
H0 = 70.

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
        self.read_filaments()
        self.read_tempel()        
        self.read_a100()
        self.read_agc()
        self.read_co()
        self.read_unwise()        
        self.read_halpha()
        self.read_rphot()
        self.read_legacy()        
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
        tab1 = Table.read(self.tabledir+self.tableprefix+'main_env_prop_H0_74_0_Mr_max_-15_7.fits')
        #tab2 = Table.read(self.tabledir+self.tableprefix+'main_envsummary.fits')
        tab3 = Table.read(self.tabledir+self.tableprefix+'main_finalenvironments.fits')        
        #self.env = hstack([tab1,tab2,tab3])
        self.env = hstack([tab1,tab3])
        pass
    def read_filaments(self):
        ''' read in GC's filament_membership catalog  '''
        #self.fil = Table.read(self.tabledir+self.tableprefix+'main_filament_membership.fits')
        self.fil = Table.read(self.tabledir+self.tableprefix+'main_filament_membership_allgalaxies.fits')

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
    def read_unwise(self):
        ''' read in unWISE table; store as self.unwise  '''
        self.unwise = Table.read(self.tabledir+self.tableprefix+'unwise.fits')                               
        pass
    def read_halpha(self):
        ''' read in halpha observations table; store as self.ha  '''
        self.ha = Table.read(self.tabledir+self.tableprefix+'ha.fits')                               
        pass
    def read_tempel(self):
        ''' read in Tempel table; store as self.tempel  '''
        tab1 = Table.read(self.tabledir+self.tableprefix+'main_Tempelgroups_infos.fits')
        tab2 = Table.read(self.tabledir+self.tableprefix+'matchTempel_groupinfo.fits')
        self.tempel = hstack([tab1,tab2])
        # should merge these, and maybe remove the columns from main
        
    def read_rphot(self):
        ''' read in rband phot  '''
        self.rphot = Table.read(self.tabledir+self.tableprefix+'r_photometry.fits')
    def read_legacy(self):
        ''' read in rband phot  '''
        dr9 = Table.read(self.tabledir+self.tableprefix+'legacy_dr9.fits')
        g = 22.5 - 2.5*np.log10(dr9['FLUX_G'])
        r = 22.5 - 2.5*np.log10(dr9['FLUX_R'])
        z = 22.5 - 2.5*np.log10(dr9['FLUX_Z'])
        d_pc = self.env['Vcosmic']/H0*1.e6
        const = 5*np.log10(d_pc) - 5
        MG = g - const
        MR = r - const
        MZ = z - const
        newtab = Table([g,r,z,MG,MR,MZ],names=['g','r','z','Mg','Mr','Mz'])
        self.dr9 = hstack([dr9,newtab])
                                 
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v1/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_north_v1_', help = 'prefix for tables; default is vf_north_v1')                               
    args = parser.parse_args()
    
    v = vtables(args.tabledir,args.tableprefix) 
    v.read_all()
