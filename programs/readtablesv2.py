#!/usr/bin/env python

'''
program to read in subtables for virgo



written for version 2 tables
10/02/21 by Rose Finn

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
        self.read_a100()
        self.read_agc()        
        self.read_co()
        self.read_halpha()
        self.read_hyperleda()
        self.read_ned()
        self.read_nsav0()
        self.read_nsav1()
        self.read_steer17()
        self.read_z0mgs()
        self.read_unwise()        
        self.read_env()
        self.read_filaments()
        #self.read_magphys()
        self.read_paper1()
        #except FileNotFoundError:
        #    print("WARNING: magphys file not found (this is probably ok)")
        try:
            self.read_magphys()
        except FileNotFoundError:
            print("WARNING: magphys file not found (this is probably ok)")
        #self.read_tempel()        
        #self.read_rphot()
        self.read_legacy()
        self.read_ephot()
        self.read_galfit()        
        #self.read_extinction()
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
        self.nsav1 = Table.read(self.tabledir+self.tableprefix+'nsa_v1.fits')
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
        #tab1 = Table.read(self.tabledir+self.tableprefix+'main_env_prop_H0_74_0_Mr_max_-15_7.fits')
        #tab2 = Table.read(self.tabledir+self.tableprefix+'main_envsummary.fits')
        #tab3 = Table.read(self.tabledir+self.tableprefix+'main_environment.fits')        
        #self.env = hstack([tab1,tab2,tab3])
        self.env = Table.read(self.tabledir+self.tableprefix+'environment.fits')    
        pass
    def read_filaments(self):
        ''' read in GC's filament_membership catalog  '''
        #self.fil = Table.read(self.tabledir+self.tableprefix+'main_filament_membership.fits')
        #self.fil = Table.read(self.tabledir+self.tableprefix+'main_filament_membership_allgalaxies.fits')
        self.fil = Table.read(self.tabledir+self.tableprefix+'filament_distances.fits')

        pass
    def read_paper1(self):
        ''' read in GC's table from paper1  '''
        #self.fil = Table.read(self.tabledir+self.tableprefix+'main_filament_membership.fits')
        #self.fil = Table.read(self.tabledir+self.tableprefix+'main_filament_membership_allgalaxies.fits')
        self.paper1 = Table.read(self.tabledir+self.tableprefix+'paper1.fits')

        pass
                              
    def read_a100(self):
        ''' read in ALFALFA 100 table; store as self.a100  '''
        self.a100 = Table.read(self.tabledir+self.tableprefix+'a100.fits')
        #self.a100sdss = Table.read(self.tabledir+self.tableprefix+'a100_sdssphot.fits')
        #self.a100unwise = Table.read(self.tabledir+self.tableprefix+'a100_unwise.fits')
        # calculate HI deficiency
        # use the R90 petro
        pass
    def read_agc(self):
        ''' read in AGC table; store as self.agc  '''
        self.agc = Table.read(self.tabledir+self.tableprefix+'agc.fits')                               
        pass
    def read_co(self):
        ''' read in CO table; store as self.co  '''
        #self.co = Table.read(self.tabledir+self.tableprefix+'co.fits')
        self.co = Table.read(self.tabledir+self.tableprefix+'CO_HI.fits')
        pass
    def read_unwise(self):
        ''' read in unWISE table; store as self.unwise  '''
        self.unwise = Table.read(self.tabledir+self.tableprefix+'unwise.fits')                               
        pass
    def read_galfit(self):
        ''' read in galfit tables of single-component Sersic fits; store as self.galfit_{} [g,r,z,W1,W2,W3,W4]'''
        
        self.galfit_g = Table.read(self.tabledir+self.tableprefix+'galfit_g.fits')
        self.galfit_r = Table.read(self.tabledir+self.tableprefix+'galfit_r.fits')
        self.galfit_z = Table.read(self.tabledir+self.tableprefix+'galfit_z.fits')

        self.galfit_W1 = Table.read(self.tabledir+self.tableprefix+'galfit_W1.fits')
        self.galfit_W2 = Table.read(self.tabledir+self.tableprefix+'galfit_W2.fits')        
        self.galfit_W3 = Table.read(self.tabledir+self.tableprefix+'galfit_W3.fits')
        self.galfit_W4 = Table.read(self.tabledir+self.tableprefix+'galfit_W4.fits')        

    def read_halpha(self):
        ''' read in halpha observations table; store as self.halpha; table generated from web coadds is self.haobs  '''


        self.halpha = Table.read(self.tabledir+self.tableprefix+'halpha.fits')
        # read in table that is created from making web pages for gui
        self.haobs = Table.read(self.tabledir+self.tableprefix+'halpha_obs.fits')
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
        ''' read in legacy dr9 photometry file; store as self.dr9  '''
        dr9 = Table.read(self.tabledir+self.tableprefix+'legacy_dr9.fits')
        self.dr9 = dr9
        # the following quantities are already in legacy_dr9.fits
        #g = 22.5 - 2.5*np.log10(dr9['FLUX_G'])
        #r = 22.5 - 2.5*np.log10(dr9['FLUX_R'])
        #z = 22.5 - 2.5*np.log10(dr9['FLUX_Z'])
        #d_pc = self.env['Vcosmic']/H0*1.e6
        #const = 5*np.log10(d_pc) - 5
        #MG = g - const
        #MR = r - const
        #MZ = z - const
        #newtab = Table([g,r,z,MG,MR,MZ],names=['g','r','z','Mg','Mr','Mz'])
        #self.dr9 = hstack([dr9,newtab])

        
    def read_magphys(self):
        """
        self.magphys_lext = read in magphys table for legacy extinction with z-band
        self.magphys_noz_lext = read in magphys table for legacy extinction, without z-band
        self.magphys = no zband in the North, and zband in the South
        """
        #tab1 = Table.read(self.tabledir+self.tableprefix+'main_env_prop_H0_74_0_Mr_max_-15_7.fits')
        #tab2 = Table.read(self.tabledir+self.tableprefix+'main_envsummary.fits')
        #tab3 = Table.read(self.tabledir+self.tableprefix+'main_environment.fits')        
        #self.env = hstack([tab1,tab2,tab3])

        # not keeping all of these variations - just the Legacy Extinction, with and w/out zband
        #self.magphys_noext = Table.read(self.tabledir+self.tableprefix+'magphys_10-Jul-2023.fits')


        # these are using version V3b on John Moustakas's tables
        self.magphys_lext = Table.read(self.tabledir+self.tableprefix+'magphys_legacyExt_16-Feb-2024.fits')
        self.magphys_noz_lext = Table.read(self.tabledir+self.tableprefix+'magphys_nozband_legacyExt_16-Feb-2024.fits')

        
        #self.magphys_sext = Table.read(self.tabledir+self.tableprefix+'magphys_salimExt_11-Jul-2023.fits')


        outtab = self.tabledir+self.tableprefix+'magphys_legacyExt_final.fits'
        self.magphys = Table.read(outtab)

    def read_extinction(self):
        ''' read in extinction table '''
        self.extinct = Table.read(self.tabledir+self.tableprefix+'extinction.fits')        
        pass
    #added by GHR on 17.April.2023
    def read_ephot(self):
        ''' read in elliptical aperture photometry from John '''
        self.ephot = Table.read(self.tabledir+self.tableprefix+'legacy_ephot.fits')        
        pass
        

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v2/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_v2_', help = 'prefix for tables; default is vf_v2')                               
    args = parser.parse_args()

    if args.tabledir.startswith('/home/rfinn/'):
        homedir = os.getenv("HOME")
        args.tabledir = args.tabledir.replace('/home/rfinn',homedir)
    v = vtables(args.tabledir,args.tableprefix) 
    v.read_all()
