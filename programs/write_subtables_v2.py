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

ARGS:

OUTPUT:
- write out several different views
  - one for matching with Leroy+2019 sample
  - one for Dustin Lang to get unWISE photometry
  - one for catalog paper

NOTES:
* This program is being updated for v2.  

* v2 has different number of galaxies as v1.  
  * we had one missing galaxies that was supposed to be merged from two input sources, but we got rid of both. (We decided NOT to include this b/c it would open a can of worms to go back to step one...)
  * JM found one source that was part of UGC04499
  * we then matched the catalog to itself, looking for sources w/in 50" of each other. 
    * we found 5 stars in this comparison, as well as 3 galaxies that need their centers adjusted
  * we then looked for matches with Dustin Lang's catalog of Tycho and Gaia sources.  found 140 matches w/in 5".  inspected these.  most are the saturated centers of bright galaxies (brings up issue with photometry for these gals..). Found 2 stars that we already had identified using the internal match.  we found another two sources that need to have coords updated.
  * v2 has 6780 sources, whereas v1 had 6797
  * the sources that are eliminated are listed in remove_stars
* In addition to this difference in 17 sources that have been removed, we have cleaned the tables and eliminated duplicate columns.

* The v2 catalogs have different VFIDs.  I have included the v1 ids as well.

***********************************************************************************
* Aug 2022 - data specialists at NED contacted us and said some of the NED names were wrong on our tables.
Gianluca and I met on 8/18/2022 to look through catalogs and see what was going on.

Conclusion: if a
galaxy did not have a HL name, I would query NED using the NSA id.  Unfortunately, I used the NSA V1 ids, rather
than the NSA V0 ids.  NED does not knowt the V1 id, so NED interprets the name as a NSA v0 id.  This of course
leads to the wrong galaxy...

Fix:
we searched again using the NSAID v0 for the NED name for the 219 galaxies for
which we originally searched using the NSAID v1.  The NED name is stored in

/home/rfinn/research/Virgo/tables-north/vf_v2_correct_NEDname_for_NSAsourcegals.fits

I will update the main table and the nedquery table
*************************************************************************************
'''
import os
import shutil
import sys
import numpy as np
#import time
import argparse

from astropy.io import fits
from astropy.table import Table, QTable, hstack#, join, vstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy import constants
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo

H0 = cosmo.H0.value
#from astroquery.ned import Ned

#from matplotlib import pyplot as plt
homedir = os.getenv("HOME")
#sys.path.append(homedir+'/github/appss/')
#from join_catalogs import make_new_cats, join_cats

from virgoCommon import *



###################################################################
#### CONSTANTS
###################################################################
sdsspixelscale=0.396127#conversion for isophotal radii from pixels to arcseconds
H0 = 70.
h = H0/100 #H0/100






###################################################################
#### FUNCTIONS
###################################################################
def duplicates(table,column,flag=None):
    if flag is None:
        unique, counts = np.unique(table[column], return_counts=True)
    elif flag is not None:
        unique, counts = np.unique(table[column][flag], return_counts=True)
    print('number of duplicates = ',sum(counts > 1))
    #print('duplicates = ',unique[counts > 1])
    return unique, counts

###################################################################
#### CLASSES
###################################################################
class catalog:
    def __init__(self,catalog,hav1=False,version='v1'):
        self.cat = Table(fits.getdata(catalog))

        if NORTH_ONLY:
            self.cat, self.keepnorth_flag = self.keep_north()
        # V2 change
        self.hav1flag = hav1
        self.version = version


    def def_basictable(self):
        # removing RA, DEC and NEDname as part of basic table
        self.basictable = self.cat['VFID','RA','DEC','NEDname']
        # keeping NEDname for now b/c that is how I am matching some catalogs
        #self.basictable_wNED = self.cat['VFID','NEDname']
        #self.basictable = self.cat['VFID']
        self.vfidtable = self.cat['VFID']

    def def_maintable(self):
        self.maintable = self.cat['VFID','RA','DEC','vr','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag','VFID_V1']
        self.maintable.rename_column('NSAID_2','NSAIDV0')
        self.maintable.rename_column('NSA0flag','NSAV0flag')
        
        self.maintable.rename_column('NSAID','NSAIDV1')
        self.maintable.rename_column('NSAflag','NSAV1flag')

        #c1 = Column(self.vfids_v1

    def keep_north(self, cat=None):
        # cut the catalog to keep Dec > -1 only
        if cat is None:
            cat = self.cat
        else:
            cat = cat
        keepnorth = cat['DEC'] > -1.3
        return cat[keepnorth], keepnorth

    def recenter_coords(self):
        ''' recenter coords of some galaxies from v1 tables '''
        new_coords={357:(245.2100,63.1211),\
                    3765:(170.9218,19.2717),\
                    4275:(182.3882,14.4823),\
                    1818:(261.0407,43.4527),\
                    5199:(145.9635,9.6550)}
        for id in new_coords.keys():
            self.cat['RA'][id] = new_coords[id][0]
            self.cat['DEC'][id] = new_coords[id][1]
                    
        pass
        
    def remove_stars(self):
        ''' remove stars that were indentified in the v1 tables '''
        # did a 90" internal match to v1 catalog and then
        # reviewed matches using the legacy viewer.
        
        # stars and globular clusters
        stars = [591,2998,4299,4582,4662, 4668,\
                 4669,4705,4725,4728,4726,4730,\
                 4735,4832,4839,5438]
        # parts of galaxies or extremely questionable
        shredded = [1241]
        bad_ids = stars+shredded

        # cutout bad sources
        keepflag = np.ones(len(self.cat),'bool')
        for id in bad_ids:
            keepflag[id]=False

        # we are hereby changing the catalog length,
        # so that v2 ids are not the same as v1 ids
        self.cat = self.cat[keepflag]
        self.v2keepflag = keepflag

        pass
    
    def redefine_vfid(self):
        ''' renumber the VFIDs after the stars are removed '''
        self.cat.rename_column('VFID','VFID_V1')
        self.vfid_v2 = ['VFID{:04d}'.format(i) for i in np.arange(len(self.cat))]
        c = Column(self.vfid_v2,name="VFID")
        self.cat.add_column(c,index=0)
                    
        pass
    

    def get_unwise(self):        
        self.unwise = Table.read(unwise_cat)
        self.unwiseFlag = (self.unwise['x'] > 0)
        # clean catalog
        unwise_keepcols = self.unwise.colnames[22:75]
        self.unwise = self.unwise[unwise_keepcols]
        self.unwise.rename_column('ra_2','ra_unwise')
        self.unwise.rename_column('dec_2','dec_unwise')

        # cut to length of v2

        self.unwise = self.unwise[self.v2keepflag]
        self.unwiseFlag = self.unwiseFlag[self.v2keepflag]
        
        # insert the VFID into the first column
        self.unwise.add_column(self.basictable['VFID'],index=0)
        self.unwise.write(outdir+file_root+'unwise.fits',format='fits',overwrite=True)

    def get_pgcname(self):
        ''' get PGC name to include in table 1 '''
        pgcfile = Table.read(homedir+'/research/Virgo/ancil-tables/vf_north_v1_hyperleda_pgc.fits')
        pgcdict = dict((a,b) for a,b in zip(pgcfile['VFID'],pgcfile['pgc']))
        self.pgcname = np.zeros(len(self.cat),'i')
        for i,vf in enumerate(self.cat['VFID_V1']):
            try:
                self.pgcname[i] = pgcdict[vf]
            except KeyError:
                print('no pgc match for ',vf)

    def fix_NED_duplicates(self):
        '''
        GOAL: there are 7 pairs of galaxies matched to same NED source

        - in most cases, this is the wrong match for one gal in the pair
        - going to remove NED name for this galaxy
        - do this before matching to CO table to avoid duplicate CO entries
        - one is a merger, so leaving this

        ALSO August 2022 Update
        - 219 galaxies were matched to NED using the v1 NSAID, but NED thinks this is v0
        - fixing the NED names for these

        '''
        VFID_wrong_NEDname = ['VFID0356','VFID3406','VFID5564',\
                              'VFID4586','VFID6053','VFID3503']

        for vfid in VFID_wrong_NEDname:
            self.cat['NEDname'][self.cat['VFID'] == vfid]=''


        # DONE - TODO fix NED names for ones that were mismatched with NSA v1 ids
        # how we do this depends on whether we have cut the catalog yet
        # the table with corrected NED names
        #
        # /home/rfinn/research/Virgo/tables-north/v2_v2_correct_NEDname_for_NSAsourcegals.fits
        #
        # has 6780 lines.  Need to make sure self.cat does as well

        # read in table
        fixedNED = Table.read(homedir+'/research/Virgo/tables-north/vf_v2_corrected_NEDname_for_NSAsourcegals.fits')
        nindices = np.arange(len(fixedNED))[fixedNED['NSAsourceflag']]

        # check to make sure self.cat and fixedNED are same length
        if len(self.cat) == len(fixedNED):
            # fix self.cat['NEDinput', 'NEDname', 'NEDra', 'NEDdec']

            for i in nindices:
                self.cat['NEDinput'][i] = self.cat['NSAID'][i]
                self.cat['NEDname'][i] = fixedNED['NEDname_NSAV0_input'][i]
                self.cat['NEDra'][i] = fixedNED['NSAra'][i]
                self.cat['NEDdec'][i] = fixedNED['NSAdec'][i]
        else:
            print('ERROR: NED table and self.cat are not the same length')
            sys.exit()
        pass
        
    def runall(self):
        #self.catalog_for_z0MGS()

        
        # updates for V2
        self.recenter_coords()
        self.remove_stars()
        self.redefine_vfid()

        # remove galaxies matched to wrong NED name
        self.fix_NED_duplicates()
        
        # back to v1 stuff
        self.def_basictable()
        #self.def_maintable()

        self.catcoord = SkyCoord(self.cat['RA'],self.cat['DEC'],frame='icrs',unit='deg')
        


        # GL is making the steer table now
        self.get_steer17()
        self.get_radius()

        # skipping SFR calculations in v2
        # these will be replace with SED fitting from legacy photometry
        self.get_sfr()
        
        self.get_HIdef()

        
        # add pgc name
        self.get_pgcname()
        # rematching to z0MGS with new v2 length
        self.get_z0MGS_flag()


        # routines that match to catalog
        self.get_CO()
        self.get_CO_paper1()
        # updating halpha routine to match to VFID_V1
        self.get_halpha()
        self.get_unwise()        
        self.main_table()
        
        self.hyperleda_table()
        self.nsa_v1_table()
        self.nsa_v0_table()
        self.agc_table()        
        self.a100_table()
        #self.a100_sdss_table()

        # V2 change
        # not including a100_unwise table
        #self.a100_unwise_table()
        
        self.get_size_for_JM()
        
        # convert John's legacy ellip phot to line-matched for v2
        self.get_JMphot_table()

        # combine magphys tables, no zband N of 32 dec, zband in south
        self.get_magphys()
        
        self.ned_table() # NED input, ra, dec, and NEDname
        self.print_stats()
        self.convert_rphot_v1_2_v2()        
        self.combine_env_v1_2_v2()

        # make names of filament spine files to be consistent with v2
        # naming conventions
        self.update_spine_filenames_v1_2_v2()

        # fix table name and convert flag column from 0/1 to boolean
        self.convert_kourkchi_v1_2_v2()
        
        #self.convert_v1_2_v2()
        self.merge_legacy_phot_tables()
        self.make_legacy_viewer_table()

        self.get_extinction()

        # write the combined magphys table
        
    def print_stats(self):
        print('Number in sample = ',len(self.cat))
        print('Number with CO data = ',sum(self.coflag))
        print('Number with A100 data = %i (%.3f)'%(sum(self.cat['A100flag']),sum(self.cat['A100flag'])/len(self.cat['A100flag'])))
        print('Number with z0MGS matches = %i (%.3f)'%(sum(self.z0mgsFlag),sum(self.z0mgsFlag)/len(self.z0mgsFlag)))
        try:
            print('Number with steer17 matches = %i (%.3f)'%(sum(self.steerFlag),sum(self.steerFlag)/len(self.steerFlag)))
        except AttributeError:
            print("problem finding steer flag")
        try:
            f = self.unwiseFlag
            print('Number with unwise matches = %i (%.3f)'%(sum(f.data),sum(f.data)/len(f)))        
            f = self.unwiseFlag & (self.cat['NSAflag'] | self.cat['NSA0flag'])
            print("Number with unWISE and NSA = %i (%.3f)"%(sum(f.data),sum(f.data)/len(f)))
        except AttributeError:
            pass
        print("CO SOURCES")
        nco = sum(self.coflag)
        f = self.coflag & self.z0mgsFlag
        print("\tNumber of CO sources in z0MGS = %i (%.2f)"%(sum(f.data),sum(f.data)/nco))
        try:
            f = self.coflag & self.steerFlag
            print("\tNumber of CO sources in Steer = %i (%.2f)"%(sum(f.data),sum(f.data)/nco))
            f = self.coflag & self.z0mgsFlag & self.steerFlag
            print("\tNumber of CO sources in z0MGS+Steer = %i (%.2f)"%(sum(f.data),sum(f.data)/nco))
        except AttributeError:
            print("no steer flag")
        try:
            f = self.coflag & self.unwiseFlag & (self.cat['NSAflag'] | self.cat['NSA0flag'])
            print("\tNumber of CO sources with unWISE and NSA = %i (%.2f)"%(sum(f.data),sum(f.data)/nco))
        except AttributeError:
            print('FYI: no unwise flag')


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
        ''' get estimate of radius of objects, using D25 or a scaling from another radius measurement '''
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

        # use scaling relation from
        ## https://www.aanda.org/articles/aa/full_html/2013/12/aa21326-13/T1.html
        ## M. Argudo-Fernández1, 2013, A&A, 560, A9
        ## Answers vary between 1.23 and 1.58 depending on sample
        ##
        ## 1.58 for isolated galaxies
        ## 
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
        ## IF NONE OF THESE SIZES AVAILABLE, SET RADIUS TO 30 ARCSEC
        flag = ~d25flag & ~self.cat['NSAflag'] & ~self.cat['NSA0flag'] & ~a100radius_flag
        self.radius[flag] = 30*np.ones(sum(flag),'f')
        print('number of galaxies with no size measurement = {:d} ({:.2f}%)'.format(sum(flag),sum(flag)/len(flag)))
        radius_flag = ~flag

        ## READ IN SPREADSHEET FROM MATCHING VF CATALOG WITH JM'S SMA CATALOG
        ## CHECK TO SEE IF GALAXY HAS A FIX_RAD=1
        ## IF SO, UPDATE self.radius but DON'T update radius_flag b/c this will be used to calculate HI def

        vfsheet = Table.read(homedir+'/research/Virgo/google-tables/VF-notin-SGA-10arcsecmatch-bestmatch-symmetric.csv',format='ascii')
        fixrad_flag = (vfsheet['fix_rad'] == 1) & (vfsheet['keep?'] == 1)
        print('adjusting radius for ',sum(fixrad_flag),' galaxies based on legacy images')
        
        vfsheet = vfsheet[fixrad_flag]
        print('number of galaxies with updated values of radius = ',len(vfsheet))
        # FIND MATCHING GALAXY BASED ON objname or NSAID or NSAIDV0
        catindex = np.arange(len(self.cat))
        for i in range(len(vfsheet)):
            matchflag=False
            newradius = vfsheet['radius_new'][i]
            # match to objname if found
            # ran into some trouble with empty values - they are "masked"
            # but the whole column isn't masked.  had to cluge together a solution
            #print(i)
            try:
                t = len(vfsheet['objname'][i])
                gmatch = (vfsheet['objname'][i] == self.cat['objname'])
                #print(i,vfsheet['objname'][i],gmatch)
                if sum(gmatch) == 1:
                    matchflag = True
                    matchindex = catindex[gmatch]
            except TypeError: # will go here if there is no value for objname
                try:
                    t=len(vfsheet['NSAID'][i])
                    # match to NSAID if found
                    gmatch = (vfsheet['NSAID'][i] == self.cat['NSAID'])
                    if sum(gmatch) == 1:
                        matchflag = True
                        matchindex = catindex[gmatch]
                except TypeError: # go here for NSAIDV0
                    try:
                        t = len(vfsheet['NSAIDV0'][i])
                        # match to NSAIDV0 if found
                        gmatch = (vfsheet['NSAIDV0'][i] == self.cat['NSAIDV0'])
                        if sum(gmatch) == 1:
                            matchflag = True
                            matchindex = catindex[gmatch]

                    except TypeError: # will do here if NSAIDV0 name is zero
                        # try matching by RA and DEC
                        d = np.sqrt((self.cat['RA']-vfsheet['RA'][i])**2+(self.cat['DEC']-vfsheet['DEC'][i])**2)
                        if min(d) < 3./3600.: # require match to be within 3 arcsec
                            gmatch = d == min(d)
                            matchindex = catindex[gmatch]
                            matchflag=True
                        pass
            if matchflag:
                self.radius[matchindex] = newradius
                #print('matched a galaxy to update radius!')
            else:
                print('could not match galaxy ')
                print(vfsheet[i])
        
        c = Column(self.radius,name='radius',unit='arcsec')
        self.cat.add_column(c)
        c = Column(radius_flag,name='radius_flag',description='if False, then rad is set to 30 arcsec')
        self.cat.add_column(c)

    def get_sfr(self):
        print('CALCULATING SFRS')
        ### calculate the UV + IR SFR

        # estimate distance from recession velocity
        distance = self.cat['vr']/70 
        # NUV is 230 nm, according to Kennicutt & Evans
        wavelength_NUV = 230.e-9*u.m
        freq_NUV = constants.c/wavelength_NUV
        
        # convert NSA NUV mag to nuLnu_NUV
        # NSA magnitudes are already corrected for galactic extinction
        nuv_mag = np.zeros(len(self.cat),'d')
        fnu_nuv = np.zeros(len(self.cat),'d')*u.Jy
        flagNUV = self.cat['SERSIC_NMGY'][:,1] > 0.
        nuv_mag[flagNUV] = 22.5 - np.log10(self.cat['SERSIC_NMGY'][:,1][flagNUV])
        fnu_nuv[flagNUV] = 3631*10**(-1*nuv_mag[flagNUV]/2.5)*u.Jy
        nuLnu_NUV = fnu_nuv*4*np.pi*(distance*u.Mpc)**2*freq_NUV
        self.nuLnu_NUV = nuLnu_NUV

        
        # GET IR VALUES
        wavelength_22 = 22*u.micron
        freq_22 = constants.c/wavelength_22
        # need to convert W4 flux from vega magnitude to Jansky
        # AB to Vega conversion is about 6 mag for W4
        w4_ab_mag = self.cat['w4_mag']+6.620

        flag22 = self.cat['w4_mag'] > 0
        # flux zp in AB mag is 3630 Jy
        fluxzp_22_jy = 3631.
        Fnu22 = np.zeros(len(flag22),'d')*u.Jy
        Fnu22[flag22] = fluxzp_22_jy*10**(-1*w4_ab_mag[flag22]/2.5)*u.Jy

        # caculate nuFnu22
        nuFnu22 = Fnu22*freq_22
        # then calculate nu L_nu, using distance
        nuLnu22_ZDIST = nuFnu22 * 4 * np.pi * (distance*u.Mpc)**2                

        self.logSFR_IR_KE = np.zeros(len(self.cat),'d')
        self.logSFR_IR_KE[flag22] = np.log10(nuLnu22_ZDIST.cgs.value[flag22])-42.69
        print(max(self.logSFR_IR_KE))
        self.Fnu22 = Fnu22
        self.nuFnu22 = nuFnu22
        self.nuLnu22 = nuLnu22_ZDIST                
        self.distance = distance
        c0 = Column(self.logSFR_IR_KE,name='logSFR22_KE')
        c0a = Column(flag22,name='W4flag')        

        
        # correct NUV luminosity by IR flux
        myunit = nuLnu_NUV.unit
        nuLnu_NUV_cor = np.zeros(len(nuLnu_NUV),'d')*myunit
        flag = self.cat['w4_mag'] > 0.
        #self.nuLnu_NUV_cor = self.nuLnu_NUV

        nuLnu_NUV_cor = nuLnu_NUV + 2.26*nuLnu22_ZDIST
        # need relation for calculating SFR from UV only
        #
        # eqn 12
        # log SFR(Msun/yr) = log Lx - log Cx
        # NUV - log Cx = 43.17
        # 24um - logCx = 42.69
        # Halpha - log Cx = 41.27

        ### ZERO OUT MEANINGLESS NUMBERS
        self.logSFR_NUV_KE = np.zeros(len(self.cat),'d')
        self.logSFR_NUVIR_KE = np.zeros(len(self.cat),'d')
        self.logSFR_NUV_KE[flagNUV] = np.log10(nuLnu_NUV.cgs.value[flagNUV]) - 43.17
        flag = flagNUV & flag22
        self.logSFR_NUVIR_KE[flag] = np.log10(nuLnu_NUV_cor.cgs.value[flag]) - 43.17
    
        # write columns out to table
        #c0 = MaskedColumn(self.nuLnu_NUV,name='nuLnu_NUV')


        c1 = Column(self.logSFR_NUV_KE,name='logSFR_NUV_KE')
        c2 = Column(self.logSFR_NUVIR_KE,name='logSFR_NUVIR_KE')
        c2a = Column(flagNUV,name='NUVflag')                
        sfrtable = Table([c0,c1,c2,c0a,c2a])
        sfrtab = hstack([self.basictable,sfrtable])
        ### write out to separate table
        sfrtab.write(outdir+file_root+'sfr.fits',format='fits',overwrite=True)

    def get_sfr_JMphot(self):
        """ updating SFR calculations to use legacy photometry from JM  """
        print('CALCULATING SFRS from JM phot')
        ### calculate the UV + IR SFR

        # estimate distance from recession velocity
        distance = self.cat['vr']/70 
        # NUV is 230 nm, according to Kennicutt & Evans
        wavelength_NUV = 230.e-9*u.m
        freq_NUV = constants.c/wavelength_NUV
        
        # convert NSA NUV mag to nuLnu_NUV
        # NSA magnitudes are already corrected for galactic extinction
        nuv_mag = np.zeros(len(self.cat),'d')
        fnu_nuv = np.zeros(len(self.cat),'d')*u.Jy
        flagNUV = self.cat['SERSIC_NMGY'][:,1] > 0.
        nuv_mag[flagNUV] = 22.5 - np.log10(self.cat['SERSIC_NMGY'][:,1][flagNUV])
        fnu_nuv[flagNUV] = 3631*10**(-1*nuv_mag[flagNUV]/2.5)*u.Jy
        nuLnu_NUV = fnu_nuv*4*np.pi*(distance*u.Mpc)**2*freq_NUV
        self.nuLnu_NUV = nuLnu_NUV

        
        # GET IR VALUES
        wavelength_22 = 22*u.micron
        freq_22 = constants.c/wavelength_22
        # need to convert W4 flux from vega magnitude to Jansky
        # AB to Vega conversion is about 6 mag for W4
        w4_ab_mag = self.cat['w4_mag']+6.620

        flag22 = self.cat['w4_mag'] > 0
        # flux zp in AB mag is 3630 Jy
        fluxzp_22_jy = 3631.
        Fnu22 = np.zeros(len(flag22),'d')*u.Jy
        Fnu22[flag22] = fluxzp_22_jy*10**(-1*w4_ab_mag[flag22]/2.5)*u.Jy

        # caculate nuFnu22
        nuFnu22 = Fnu22*freq_22
        # then calculate nu L_nu, using distance
        nuLnu22_ZDIST = nuFnu22 * 4 * np.pi * (distance*u.Mpc)**2                

        self.logSFR_IR_KE = np.zeros(len(self.cat),'d')
        self.logSFR_IR_KE[flag22] = np.log10(nuLnu22_ZDIST.cgs.value[flag22])-42.69
        print(max(self.logSFR_IR_KE))
        self.Fnu22 = Fnu22
        self.nuFnu22 = nuFnu22
        self.nuLnu22 = nuLnu22_ZDIST                
        self.distance = distance
        c0 = Column(self.logSFR_IR_KE,name='logSFR22_KE')
        c0a = Column(flag22,name='W4flag')        

        
        # correct NUV luminosity by IR flux
        myunit = nuLnu_NUV.unit
        nuLnu_NUV_cor = np.zeros(len(nuLnu_NUV),'d')*myunit
        flag = self.cat['w4_mag'] > 0.
        #self.nuLnu_NUV_cor = self.nuLnu_NUV

        nuLnu_NUV_cor = nuLnu_NUV + 2.26*nuLnu22_ZDIST
        # need relation for calculating SFR from UV only
        #
        # eqn 12
        # log SFR(Msun/yr) = log Lx - log Cx
        # NUV - log Cx = 43.17
        # 24um - logCx = 42.69
        # Halpha - log Cx = 41.27

        ### ZERO OUT MEANINGLESS NUMBERS
        self.logSFR_NUV_KE = np.zeros(len(self.cat),'d')
        self.logSFR_NUVIR_KE = np.zeros(len(self.cat),'d')
        self.logSFR_NUV_KE[flagNUV] = np.log10(nuLnu_NUV.cgs.value[flagNUV]) - 43.17
        flag = flagNUV & flag22
        self.logSFR_NUVIR_KE[flag] = np.log10(nuLnu_NUV_cor.cgs.value[flag]) - 43.17
    
        # write columns out to table
        #c0 = MaskedColumn(self.nuLnu_NUV,name='nuLnu_NUV')


        c1 = Column(self.logSFR_NUV_KE,name='logSFR_NUV_KE')
        c2 = Column(self.logSFR_NUVIR_KE,name='logSFR_NUVIR_KE')
        c2a = Column(flagNUV,name='NUVflag')                
        sfrtable = Table([c0,c1,c2,c0a,c2a])
        sfrtab = hstack([self.basictable,sfrtable])
        ### write out to separate table
        sfrtab.write(outdir+file_root+'sfr.fits',format='fits',overwrite=True)


    def get_HIdef(self):
        '''
        use relationship from Boselli+Gavazzi 2009
        https://www.aanda.org/articles/aa/full_html/2009/46/aa12658-09/aa12658-09.html

        HI def = log MHI_ref - log MHI_obs
        h^2 MHI_ref = c + d log(hdiam)^2

        Type 	c 	d 	Ref.
        E-S0a 	6.88 	0.89 	HG84
        Sa-Sab 	7.75 	0.59 	S96
        Sb 	7.82 	0.62 	S96
        Sbc 	7.84 	0.61 	S96
        Sc 	7.16 	0.87 	S96
        Scd-Im-BCD 	7.45 	0.70 	G10 

        email from Martha:

        First, Carmen Toribio looked at a subset of isolated galaxies in
        a40 to derive scaling relations
        https://ui.adsabs.harvard.edu/abs/2011ApJ...732...93T
        https://ui.adsabs.harvard.edu/abs/2015ApJ...802...72T/     (erratum to
        above)
        
        More recently, Mike led an effort to establish relations for the AMIGA
        sample (of isolated galaxies):
        https://ui.adsabs.harvard.edu/abs/2018A%26A...609A..17J
        
        They use isophotal radii, and refer to a discussion of deriving D25 from
        SDSS in the
        earlier paper by Argudo-Fernandez 2013:
        https://ui.adsabs.harvard.edu/abs/2013A%26A...560A...9A
        
        I think there are some others too, but these are what I am most familiar
        with.
        
        Martha
        
        '''
        # get distance from GL env cata
        #env = Table.read('temp')
        #distance = self.env['Vcosmic']/H0
        distance = self.cat['vr']/cosmo.H0
        # calculate MHI (even though A100 catalog has this - need to use consistent distance
        
        self.logMHI = np.log10(2.3566e5*self.cat['HIflux'])+2*np.log10(self.cat['Dist']) #Dist^2
        
        # calculate HI deficiency using Toribio et al 2011 results
        # their relation is
        # log(M_HI/Msun) = 8.72 + 1.25 log(D_25,r/kpc)
        # and
        # log D_25 = log D_25(obs) + beta log(b/a), where beta = 0.35 in r-band
        # NOTE: SDSS isophotal radii are given in pixels!!!!
        # convert from arcsec to kpc with self.AngDistance (which is in units of kpc/arcsec)
        # multiply by 2 to convert from radius to diameter
        # multiply by sdss pixel scale (0.39) b/c isophotal radii are given in pixels

        # use the radius measurements that I collated for John's group catalog
        # these should approximate D25
        # only calculate for galaxies with radiusflag = True

        # returns DA in Mpc/radians
        DA=cosmo.angular_diameter_distance(self.cat['Vhelio']/3.e5)

        # try using distance in A100
        vr = self.cat['Dist']*cosmo.H0
        DA = cosmo.angular_diameter_distance(vr/3.e5)
        D25obskpc=2.*self.radius/206264*DA.to('kpc')

        
        # apply correction from toribio et al 2011 
        logD25kpc=np.log10(D25obskpc.value) + 0.35*np.log10(self.cat['expAB_r'])
        self.logD25kpc = logD25kpc
        # use toribio et al relation to predict the expected HI mass,
        # including factor of 2 correction
        logHImassExpected = 8.72 + 1.25*(logD25kpc-np.log10(2.))
        
        # use jones+2018, A&A, 609, A17 (AMIGA sample
        # relation to predict the expected HI mass
        # using Maximum Liklihood Estimator for Detections
        logHImassExpected_jones = 7.32 + 0.86*(2*logD25kpc)
        # relation for ALL galaxies is similar
        # intrinsic scatter is 0.21
        #
        # fits by morphological type
        # <3 	1.04 ± 0.21 	6.44 ± 0.59 	0.27 ± 0.08 	0.67
        # 3–5 	0.93 ± 0.06 	7.14 ± 0.18 	0.16 ± 0.02 	0.74
        # >5 	0.81 ± 0.09 	7.53 ± 0.24 	0.17 ± 0.03 	0.73

        # calculate deficiency as log expected - log observed
        self.HIdef = np.zeros(len(self.cat),'f')
        self.HIdef_jones = np.zeros(len(self.cat),'f')
        self.HIdef_jones_bytype = np.zeros(len(self.cat),'f')        

        # report values if a100 flag and radius flag, meaning there was some measurement of radius
        flag = self.cat['A100flag'] & self.cat['radius_flag']
        self.HIdef[flag] = logHImassExpected[flag] - self.logMHI[flag]
        self.HIdef_jones[flag] = logHImassExpected_jones[flag] - self.logMHI[flag]        
        self.HIdef_flag = flag

        flag = ~np.isnan(self.cat['t']) & (self.cat['t'] < 3)
        self.HIdef_jones_bytype[flag] = 6.44 + 1.04*(2*logD25kpc[flag])- self.logMHI[flag]        
        flag = ~np.isnan(self.cat['t']) & (self.cat['t'] > 5)
        self.HIdef_jones_bytype[flag] = 7.53 + .81*(2*logD25kpc[flag]) - self.logMHI[flag]        
        flag = ~np.isnan(self.cat['t']) & (self.cat['t'] <= 5) & (self.cat['t'] >= 3)
        self.HIdef_jones_bytype[flag] = 7.14 + .93*(2*logD25kpc[flag]) - self.logMHI[flag]        

        # calculate HI def by type for hyperleda galaxies with Ttype
        
        # add HIdef to a100 catalog

        # boselli & gavazzi prescription
        # use value for Sb
        c = 7.82
        d = 0.62
        logh2MHIref = (c + 2*d*logD25kpc)
        self.HIdef_bos = logh2MHIref - self.logMHI
        pass

    def main_table(self):
        # write out
        # ra, dec, velocity, HL, NSA id, A100 id
        # NED name
        # make flags to denote if galaxy is in:
        # - CO sample
        colnames = ['VFID','RA','DEC','vr','radius','radius_flag','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag','VFID_V1']
        #self.maintable = self.cat['VFID','RA','DEC','vr','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag']
        
        self.maintable = self.cat[colnames]
        self.maintable.rename_column('NSAID_2','NSAIDV0')
        self.maintable.rename_column('NSA0flag','NSAV0flag')

        # V2 change
        self.maintable.rename_column('NSAflag','NSAV1flag')                
        self.maintable.rename_column('NSAID','NSAIDV1')        

        c1 = Column(self.coflag,name='COflag')
        c1a = Column((self.hatable['HAflag']), name='HAflag')
        c1b = Column(self.hatable['HAobsflag'],name='HAobsflag')
        c2 = Column(self.z0mgsFlag,name='Z0MGSflag')
        try:
            c3 = Column(self.steerFlag,name='Steerflag')
        except AttributeError:
            print("creating dummy column for steer flag")
            c3 = Column(np.ones(len(self.z0mgsFlag),'bool'),name='Steerflag')
        # removing unwise for now - need to update code
        
        c4 = Column(self.unwiseFlag,name='unwiseflag')        
        self.maintable.add_columns([c1,c1a,c1b,c2,c3,c4])

        c5 = Column(self.pgcname,name='PGC')
        self.maintable.add_column(c5,index=7)
        #self.maintable.add_columns([c1,c1a,c1b,c2,c3])

        
        nedname=[]
        for i in range(len(self.maintable)):
            nedname.append(self.maintable['VFID'][i]+'-'+self.maintable['NEDname'][i].replace(" ","").replace("[","").replace("]","").replace("/",""))
        
        c0= Column(nedname,name='prefix')
        self.maintable.add_column(c0)

        # V2 change
        # update AGC names for galaxies in AGC but not in ALFALFA
        

        
        # add a column called "name" for use with the legacy image server
        # this will be the VFID

        # made a separate file instead for the legacy viewer
        #c0= Column(self.cat['VFID'],name='name')
        #self.maintable.add_column(c0)

        
        # - 2MASS
        # - z0MGS
        # - unWISE
        print('writing main table :',outdir+file_root+'main.fits')
        self.maintable.write(outdir+file_root+'main.fits',format='fits',overwrite=True)

        # write out table in csv format
        self.maintable.write(outdir+file_root+'main.csv',format='ascii',delimiter=',',overwrite=True)
        
        # write out sample into smaller tables
        # 500 lines per table
        # to use when uploading tables to the legacy imager viewer
        #
        # we will use these to inspect our galaxy coordinates

        '''
        nlines = 500
        ntotal = len(self.maintable)
        ntables = (ntotal/nlines)
        if (ntables - np.floor(ntables)) > 0:
            ntables += 1
        ntables = int(ntables)
        fileroot = outdir+file_root+'main_'
        for i in range(ntables):
            vfid_min = i*nlines
            vfid_max = (i+1)*nlines
            if vfid_max > ntotal:
                vfid_max = ntotal
            fname = fileroot+'%04d_%04d.fits'%(vfid_min,vfid_max-1)
            self.maintable[vfid_min:vfid_max].write(fname,format='fits',overwrite=True)
        '''


    def get_CO(self):
        """
        I made the first version of this using topcat to do the matching while at Trieste in Apr 2022 with BV.  
        We found 389 matches within 15"

        We should find the same number of matches when running automatically.

        """
        match_by_coords=True
        print('\n######################')
        print('Running get_CO')
        print('######################\n')

        # read in CO mastertable
        # topcat can read this fine, but astropy table can't.  I get some weird result with one row
        # going to write as a fits basic table using topcat
        #
        #
        cofile = homedir+'/research/Virgo/gianluca_input_tables/for_vf_gastable.fits'
        cofile = homedir+'/research/Virgo/gianluca_input_tables/for_vf_gastable_fitsbasic.fits'        
        self.co = Table(fits.getdata(cofile))


        cocoord = SkyCoord(self.co['RA'],self.co['DEC'],unit='deg',frame='icrs')
        if match_by_coords:
            # match co coords to mastertable 
            idx, d2d, d3d = self.catcoord.match_to_catalog_sky(cocoord)
            self.d2d = d2d
            self.idx = idx
            self.coflag = d2d < 15./3600*u.deg
        
            # create a new, blank table with same # of lines as mastertable
            # but with columns like the co table
            newco = Table(np.zeros(len(self.cat),dtype=self.co.dtype))
            # add co information into the new table for the galaxies with
            # a match to the CO sample
            # NOTE: this could match multiple galaxies to the same CO source
            newco[self.coflag] = self.co[idx[self.coflag]]

            # join basic table and co table

            self.cotable = hstack([self.basictable,newco])

        ## print the CO galaxies with no matches in the mastertable 
        print(f'number of galaxies with CO matches = {np.sum(self.coflag)} (should be 389)')


        ## look for CO galaxies that were matched to multiple galaxies in the
        ## mastertable
        
        self.cotable.add_column(Column(self.coflag),name='COflag')
        self.cotable.write(outdir+file_root+'co.fits',format='fits',overwrite=True)
            
    def get_CO_old(self,match_by_coords=False,match_by_name=True):

        print('\n######################')
        print('Running get_CO')
        print('######################\n')

        # read in CO mastertable
        # match to main table
        cofile = homedir+'/github/Virgo/tables/CO-MasterFile-2018Feb16.fits'
        # this file has a new column with the exact NED names
        # can use this to match to mastertable NEDname column
        cofile = homedir+'/research/Virgo/tables/CO-MasterFile-2018Feb16-fixedNEDnames.fits'
        ned_column = 'NED_name'
        # CO file 245 sources
        # from June 2019
        cofile = homedir+'/research/Virgo/tables/All-virgo_Master_file_19Jun2019-fixedNEDnames.fits'
        cofile = homedir+'/research/Virgo/tables/galaxy_sample_prop_general_2020Oct28-fixedNEDnames.fits'
        ned_column='NEDname'
        self.co = Table(fits.getdata(cofile))

        # replace NED_name of SHOC206b with NED name for SHOC 206a
        # this is really the same galaxy, just has less common id in CO file
        # set to MCG +08-16-005;
        
        # fix case for UGC09348, until we run the full query again...
        flag = self.co[ned_column] == 'SHOC206b'
        if sum(flag > 0):
            self.co[ned_column][flag] = 'MCG +08-16-005'
        
        

        cocoord = SkyCoord(self.co['RA'],self.co['DEC'],unit='deg',frame='icrs')
        if match_by_coords:
            # match co coords to mastertable 
            idx, d2d, d3d = self.catcoord.match_to_catalog_sky(cocoord)
            self.d2d = d2d
            self.idx = idx
            self.coflag = d2d < 15./3600*u.deg
        
            # create a new, blank table with same # of lines as mastertable
            # but with columns like the co table
            newco = Table(np.zeros(len(self.cat),dtype=self.co.dtype))
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

            self.coflag = self.cotable['alphaCO'] > .1
            # having issues with this being a masked array, so
            # saving mask as the flag...
            self.coflag = ~self.coflag.mask
            #self.coflag = self.cotable['COflag']

            # also check to see which CO sources were not matched
            self.testtable = join(self.co,self.basictable,keys='NEDname',join_type='left')
            #self.coflag = len(self.co['CO']) > 0
            try:
                self.comatchflag = ~self.testtable['VFID'].mask
                print('CO sources with no match in mastertable:')
                print(self.testtable['NEDname','source_name','VH'][~self.comatchflag])
                ## plot the positions of CO galaxies that weren't matched to mastertable
                # V2 changes - RA_1 not valid...
                # don't need the plot so I am commenting this out
                #plt.figure()
                #plt.plot(self.testtable['RA_1'][~self.comatchflag],self.testtable['DEC_1'][~self.comatchflag],'bo')


            except AttributeError:
                print('all CO sources have been matched. CONGRATULATIONS!!!!!!!')
                self.comatchflag = np.ones(len(self.testtable))
        ## print the CO galaxies with no matches in the mastertable 
        print('number of galaxies with CO matches = ',sum(self.coflag))

        ## look for CO galaxies that were matched to multiple galaxies in the
        ## mastertable
        unique, counts = duplicates(self.cotable,'NEDname')
        print("CO sources that are matched to multiple galaxies in the mastertable:")
        print(unique[counts>1])
        
        
        self.cotable.add_column(Column(self.coflag),name='COflag')
        self.cotable.write(outdir+file_root+'co.fits',format='fits',overwrite=True)
            

        # print CO sources that are not in the table
    def get_CO_paper1(self,match_by_coords=False,match_by_name=True):


        print('\n######################')
        print('Running get_CO_paper1')
        print('######################\n')
        
        # CO file 245 sources
        # from Gianluca
        cofile = homedir+'/research/Virgo/gianluca_input_tables/galaxy_sample_prop_general_mod_fixedNEDnames.fits'


        # match by name
        ned_column='NEDname'
        self.co = Table(fits.getdata(cofile))

        # replace NED_name of SHOC206b with NED name for SHOC 206a
        # this is really the same galaxy, just has less common id in CO file
        # set to MCG +08-16-005;
        
        # fix case for UGC09348, until we run the full query again...
        flag = self.co[ned_column] == 'SHOC206b'
        if sum(flag > 0):
            self.co[ned_column][flag] = 'MCG +08-16-005'
        
        

        cocoord = SkyCoord(self.co['RA'],self.co['DEC'],unit='deg',frame='icrs')
        if match_by_coords:
            # match co coords to mastertable 
            idx, d2d, d3d = self.catcoord.match_to_catalog_sky(cocoord)
            self.d2d = d2d
            self.idx = idx
            self.coflag = d2d < 15./3600*u.deg
        
            # create a new, blank table with same # of lines as mastertable
            # but with columns like the co table
            newco = Table(np.zeros(len(self.cat),dtype=self.co.dtype))
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

            self.coflag = self.cotable['alphaCO'] > .1
            # having issues with this being a masked array, so
            # saving mask as the flag...
            self.coflag = ~self.coflag.mask
            #self.coflag = self.cotable['COflag']

            # also check to see which CO sources were not matched
            self.testtable = join(self.co,self.basictable,keys='NEDname',join_type='left')
            #self.coflag = len(self.co['CO']) > 0
            try:
                self.comatchflag = ~self.testtable['VFID'].mask
                print('Paper 1 CO sources with no match in mastertable:')
                print(self.testtable['NEDname','source_name','Vcosmic'][~self.comatchflag])
                ## plot the positions of CO galaxies that weren't matched to mastertable
                # V2 changes - RA_1 not valid...
                # don't need the plot so I am commenting this out
                #plt.figure()
                #plt.plot(self.testtable['RA_1'][~self.comatchflag],self.testtable['DEC_1'][~self.comatchflag],'bo')


            except AttributeError:
                print('all paper 1 CO sources have been matched. CONGRATULATIONS!!!!!!!')
                self.comatchflag = np.ones(len(self.testtable))
        ## print the CO galaxies with no matches in the mastertable 
        print('number of galaxies with Paper 1 CO matches = ',sum(self.coflag))

        ## look for CO galaxies that were matched to multiple galaxies in the
        ## mastertable
        unique, counts = duplicates(self.cotable,'NEDname')
        print("Paper 1 CO sources that are matched to multiple galaxies in the mastertable:")
        print(unique[counts>1])


        ## cut the columns 
        ## remove:
        ## env variables: n5th,n5th_2D, nCl08, filament_dist
        remove_columns = ['n5th','n5th_err','n5th_2D','n5th_2D_err','nCI08','nCI08_err',\
                          'filament_dist_3D','filament_dist_2D',\
                          'SGX','SGY','SGZ','Vcosmic',
                          'NEDname', 'source_name', 'objname_HL','id_nsa']
                          
                          
        self.cotable.remove_columns(remove_columns)
        ## HL
        # TODO - remove RA_1 and DEC_1 - these are the VFID coords
        
        self.cotable.add_column(Column(self.coflag),name='COflag')
        self.cotable.write(outdir+file_root+'paper1.fits',format='fits',overwrite=True)
            

        
    def get_halpha_old(self,halphafile=None):
        # read in Halpha observing summary file
        if halphafile is None:
            #infile = homedir+'/research/Virgo/Halpha/observing-summary-Halpha-latest.csv'
            infile = homedir+'/research/Virgo/Halpha/observing-summary-Halpha-clean-04Jun2020.csv'            
        else:
            infile = halphafile
        self.ha = Table.read(infile,format='csv')
        self.ha.rename_column('NSA ID','NSAIDV0')

        # check for duplicates in the hafile
        unique, counts = duplicates(self.ha,'NSAIDV0')
        print("Halpha sources that are listed multiple times in the halpha file:")
        #print(unique[counts>1])
        
        
        # match ha file to the base table using the NSA v0
        # match the basictable and the CO table by matching
        # entries by the NEDname colums
        #print(self.ha.colnames)
        #print(self.maintable.colnames)
        # note myjoinleft preserves the original ordering in the tables
        # rather than having the joined table sorted by according to the match key

        self.hatable = myjoinleft(self.maintable,self.ha,keys='NSAIDV0')
        print('Halpha table lengths')
        #print(len(self.maintable),len(self.ha),len(self.hatable))        
        self.haflag = ~self.hatable['Date Obs'].mask
        self.hatable.add_column(Column(self.haflag),name='haflag')
        self.hatable.write(outdir+file_root+'halpha.fits',format='fits',overwrite=True)

    def get_halpha(self,halphafile=None):
        print('\n######################')
        print('Running get_halpha')
        print('######################\n')
        
        # read in Halpha observing summary file
        if halphafile is None:
            #infile = homedir+'/research/Virgo/Halpha/observing-summary-Halpha-latest.csv'
            #infile = homedir+'/research/Virgo/Halpha/observing-summary-Halpha-clean-04Jun2020.csv'
            infile = homedir+'/research/Virgo/halpha-tables/halpha-05Jul2020.fits'
            infile = homedir+'/research/Virgo/halpha-tables-20210311/halphagui-output-combined-2021-Mar-25.fits'
            # ids are v2 VFIDs
            infile = '/home/rfinn/research/Virgo/halpha-tables/2023-July-11/halphagui-output-combined-2023-Jul-11.fits'
            infile = '/home/rfinn/research/Virgo/halpha-tables/halphagui-output-combined-2023-Jul-13.fits'
            infile = homedir+'/research/Virgo/halpha-tables/2023-July-11/halphagui-output-combined-2023-Jul-11.fits'
            infile = homedir+'/research/Virgo/halpha-tables/halphagui-output-combined-2023-Jul-13.fits'
            infile = homedir+'/research/Virgo/halpha-tables/halphagui-output-combined-2023-Aug-19.fits'
            # in the notebook duplicates, I select the best option of the duplicates
            # and write out a file with no duplicates
            infile = homedir+'/research/Virgo/halpha-tables/halphagui-output-combined-2023-Aug-27.noduplicates.fits'
            infile = homedir+'/research/Virgo/halpha-tables/halphagui-output-combined-2024-Mar-17.fits'
            # in the notebook duplicates, I select the best option of the duplicates
            # and write out a file with no duplicates
            
            infile = homedir+'/research/Virgo/halpha-tables/halphagui-output-combined-2023-Aug-27.noduplicates.fits'
            infile = homedir+'/research/Virgo/halpha-tables/halphagui-output-combined-2024-Mar-22.noduplicates.fits'            
        else:
            infile = halphafile
        self.ha = Table.read(infile,format='fits')            

        ##
        # should add some code to clean the halpha table?
        # like a cut on filter correction, or if M16 is zero?
        # maybe not, just want to decide how to handle the duplicates
        ##
        unique, counts = duplicates(self.ha,'VFID')
        print("Halpha sources that are listed multiple times in the halpha file:")
        print("\t1 observation : ",len(unique[(counts== 1)]))
        print("\t2 observations: ",len(unique[(counts>1) & (counts < 3)]))
        print("\t3 observations: ",len(unique[counts == 3]))        

        ##
        # NEED TO REMOVE DUPLICATE ROWS
        # could give preference to galaxies that have
        # - more measurements reported: GALFIT, ELLIP, SMORPH
        # - lower filter correction
        # - better seeing
        # - lower sky noise
        # - better coverage on CCD (like not near edge) - but I don't really track this...
        #   - maybe a square image vs rectangle?
        ##

        vfid_haduplicates = unique[counts > 1]

        # need to 
        ### skipping this for now.  not sure if program will crash...
        #VFID = ['VFID2080','VFID2154','VFID2145','VFID2136','VFID2060',\
        #        'VFID2157','VFID2076','VFID2089','VFID2095','VFID2095','VFID2141']
        #PID = ['v20p35','v17p12','v17p12','v17p12','v20p35',\
        #        'v17p12','v18p54','v20p35','v20p35','v18p54','v17p12']
        #for i in range(len(VFID)):
        #    idel = np.where((self.ha['VFID'] == VFID[i]) & (self.ha['POINTING'] == PID[i]))
        #    #print(idel)
        #    self.ha.remove_rows(idel)
        #print('number of lines in ha file after deleting rows = ',len(self.ha))

        ### DELETE SOME UNNECESSARY COLUMNS
        dcolnames = ['HA_FLAG','GAL_HRA','GAL_HDEC','GAL_2SERSIC','GAL_2SERSIC_ERR',\
                     'GAL_2SERSIC_ERROR','GAL_2SERSIC_CHISQ',\
                     'GAL_HXC', 'GAL_HXC_ERR','GAL_HYC', 'GAL_HYC_ERR',\
                     'GAL_HMAG', 'GAL_HMAG_ERR','GAL_HRE', 'GAL_HRE_ERR',\
                     'GAL_HN', 'GAL_HN_ERR','GAL_HBA', 'GAL_HBA_ERR',\
                     'GAL_HPA', 'GAL_HPA_ERR','GAL_HSKY', 'GAL_HCHISQ',\
                     'GAL_H2SERSIC', 'GAL_H2SERSIC_ERR', \
                     'GAL_H2SERSIC_ERROR','GAL_H2SERSIC_CHISQ',\
                     'GAL_HSERSASYM', 'GAL_HSERSASYM_ERR', \
                     'GAL_HSERSASYM_ERROR','GAL_HSERSASYM_CHISQ',\
                     'GAL_HSERSASYM_RA','GAL_HSERSASYM_DEC',\
                     'CONTSUB_FLAG','MERGER_FLAG','SCATLIGHT_FLAG',\
                     'OVERSTAR_FLAG','OVERGAL_FLAG',\
                     'PARTIAL_FLAG','EDGEON_FLAG','NUC_HA']

        self.ha.remove_columns(dcolnames)
        print('number of columns after removing cols = ',len(self.ha[0]))
        c0 = Column(np.ones(len(self.ha),'bool'),name='HAobsflag',description='observed in halpha')
        self.ha.add_column(c0)
        
        # match ha file to the base table using the NSA v0
        # match the basictable and the CO table by matching
        # entries by the NEDname colums
        #print(self.ha.colnames)
        #print(self.maintable.colnames)
        # note myjoinleft preserves the original ordering in the tables
        # rather than having the joined table sorted by according to the match key

        ### MAKE A TABLE - ALL ZEROS, WITH DATA TYPE LIKE HALPHA TABLE
        ### AND LENGTH OF BASICTABLE

        self.hatable = QTable(np.zeros(len(self.basictable),dtype=self.ha.dtype))


        ##
        # Duplicates are removed in duplicates.ipynb
        # TODO - this is probably not ideal - should convert notebook cells into script
        ##

        ### ADD VFID FOR ALL
        self.hatable['VFID'] = self.basictable['VFID']

        # V2 changes - removing RA and DEC from halpha table
        #### ADD RA FOR ALL
        self.hatable['RA'] = self.basictable['RA']
        #### ADD DEC FOR ALL
        self.hatable['DEC'] = self.basictable['DEC']

        ##
        # TODO - REMOVE ADDITIONAL COLUMNS THAT ARE DUPLICATED IN OTHER TABLES
        ##
        


        
        self.hatable.add_column(Column(self.cat['VFID_V1'],name='VFID_V1'))

        ### ADD HALPHA SOURCES TO THEIR CORRESPONDING ROWS


        ##
        # why am I matching by RA and DEC, when both catalogs have the VFID?
        ##

        
        ## for halpha measurements that reference the input catalogs, match by RA and DEC
        ## this will make the matching independent of the version number of the VFID 
        ha_input_coord = SkyCoord(self.ha['RA'],self.ha['DEC'],unit='deg',frame='icrs')

        ##
        # what happens with the duplicate observations here?
        ##
        
        # match steer to vf catalog
        idx, d2d, d3d = self.catcoord.match_to_catalog_sky(ha_input_coord)
        # consider as matches any galaxies within 30 arcsec
        matchflag = d2d < 20./3600*u.deg
        print('number of Halpha matches = ',sum(matchflag))

        #print(matchflag)
        #print(idx)
        self.hatable[matchflag] = self.ha[idx[matchflag]]
        dcolnames = ['RA','DEC','vr','NEDname','REDSHIFT','ZDIST'] 
        self.hatable.remove_columns(dcolnames)
        
        ############################################################################
        ### READ IN HALPHA ANALYSIS FILES THAT WERE CREATED BEFORE SWITCH TO V2
        ### match these by VFID
        ############################################################################
        if self.hav1flag:
            if halphafile is None:
                #infile = homedir+'/research/Virgo/Halpha/observing-summary-Halpha-latest.csv'
                #infile = homedir+'/research/Virgo/Halpha/observing-summary-Halpha-clean-04Jun2020.csv'
                infile = homedir+'/research/Virgo/halpha-tables-v1/halpha-10Feb2021.fits'
                infile = homedir+'/research/Virgo/halpha-tables-v1/halpha-10Feb2021.fits'
                infile = homedir+'/research/Virgo/halpha-tables-v1/halpha-26Feb2021.fits'                

            else:
                infile = halphafile
            print('reading in halpha v1 table ',infile)
            self.hav1 = Table.read(infile,format='fits')            
            ### DELETE SOME UNNECESSARY COLUMNS
            self.hav1.remove_columns(dcolnames)
            print('number of columns after removing cols = ',len(self.hav1[0]))
            c0 = Column(np.ones(len(self.hav1),'bool'),name='HAobsflag',description='observed in halpha')
            self.hav1.add_column(c0)

            # add in VFID_V1
            c0 = Column(self.hav1['VFID'],name='VFID_V1',description='VFID V1')
            self.hav1.add_column(c0)
        
            for i in range(len(self.hav1)):
                flag =(self.hatable['VFID_V1'] == self.hav1['VFID'][i].rstrip())
                if sum(flag) > 0:
                    #print('found halpha v1 match ',self.hav1['VFID'][i])
                    VFID_V2 = self.hatable['VFID'][i]
                    #print(np.sum(flag),len(self.hatable[flag]),len(self.hav1[i]))
                    self.hatable[flag] = self.hav1[i]
                    self.hatable['VFID'][i] = VFID_V2
        
        #self.hatable = join(self.basictable,self.ha,keys='VFID',join_type='left')
        #print('Halpha table lengths')
        #rename_cols = ['RA','DEC','vr']
        #print(len(self.maintable),len(self.ha),len(self.hatable))
        print('calculating snr')
        self.haflag = ((self.hatable['HF_R24']/self.hatable['HF_R24_ERR']) > 1.)
        print('setting HAflag')
        self.hatable['HAflag'] = self.haflag
        #self.hatable['HAflag'].description = 'halpha snr > 1'

        ## adding haobs flag for galaxies that were observed with
        ## bok in 2021A.  they haven't been processed yet, but we
        ## want to remove them from the list of targets to be observed

        # read in file with list of VFIDs that were observed
        # set HAobsflag to true for these galaxies

        use_old_ha_tables = False
        if use_old_ha_tables:
            self.use_old_hafiles(True)
        print('writing hafile')
        #fits.writeto(outdir+file_root+'ha.fits',np.array(self.hatable),overwrite=True)
        self.hatable.write(outdir+file_root+'halpha.fits',format='fits',overwrite=True)
        print('finished writing hafile')
        
    def use_old_hafiles(self,use_old_ha_tables):
        """
        obsolete,this uses hand-kept lists that I created 
        to keep track of what was observed while observations 
        were still being carried out and before data were reduced.
        keeping for posterity.
        """

        if use_old_ha_tables:
            htab = Table.read(homedir+'/github/Virgo/halpha/ha-observed-bok-21A.txt',format='ascii')
            for i,v in enumerate(htab['VFID']):
                flag =(self.hatable['VFID_V1'] == v.rstrip())
                if sum(flag) > 0:
                    #print('found match to bok obs target')
                    self.hatable['HAobsflag'][flag] = True
                    #print(self.hatable['VFID'][flag],self.hatable['HAobsflag'][flag])
                else:
                    print('no match to bok obs target ',v)
        

            print('after matching to spring 2021, number of halpha obs = {}'.format(np.sum(self.hatable['HAobsflag'])))

            # read in the list of targets from BOK 2022 Apr run
            htab = Table.read(homedir+'/github/Virgo/halpha/ha-observed-bok-22A.csv',format='csv')
            print("BOK 22A table colnames = ",htab.colnames)

            for i,v in enumerate(htab['VFID_V1']):
                #idnumber = v.replace('VFID','')
                #print(idnumber)
                flag =(self.hatable['VFID_V1'] == v)
                if sum(flag) > 0:
                    #print('found match to bok obs target')
                    self.hatable['HAobsflag'][flag] = True
                    #print(self.hatable['VFID'][flag],self.hatable['HAobsflag'][flag])
                else:
                    print('no match to bok obs target ',v)


            print('after matching to spring 2022, number of halpha obs = {}'.format(np.sum(self.hatable['HAobsflag'])))


            # read in the list of targets from BOK 2022 Apr run
            htab = Table.read(homedir+'/github/Virgo/halpha/ha-observed-int-may-22A.csv',format='csv')
            print("INT MAY 22A table colnames = ",htab.colnames)

            for i,v in enumerate(htab['VFID_V1']):
                #idnumber = v.replace('VFID','')
                #print(idnumber)
                flag =(self.hatable['VFID_V1'] == v)
                if sum(flag) > 0:
                    #print('found match to bok obs target')
                    self.hatable['HAobsflag'][flag] = True
                    #print(self.hatable['VFID'][flag],self.hatable['HAobsflag'][flag])
                else:
                    print('no match to bok obs target ',v)


            print('after matching to spring 2022, number of halpha obs = {}'.format(np.sum(self.hatable['HAobsflag'])))

        
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
        #cat = Table.read(homedir+'/research/Virgo/tables/vf-z0MGS.tbl',format='ipac')
        #cat = Table.read(homedir+'/research/Virgo/tables/vf_z0mgs_30arcsec_051920.tbl',format='ipac')
        #cat = Table.read(homedir+'/research/Virgo/tables/vf_v1_z0mgs_30arcsec_102820.tbl',format='ipac')
        #cat = Table.read(homedir+'/research/Virgo/tables/vf_v1_z0mgs_10arcsec_102920_9153gal.tbl',format='ipac')
        cat = Table.read(z0mgs_cat,format='ipac')
        print('number of lines in z0MGS cat = ',len(cat))
        print('length of keepnorth_flag = ',len(self.keepnorth_flag))
        # cut on declination
        if NORTH_ONLY:
            cat = cat[self.keepnorth_flag]
        # cut for v2
        cat = cat[self.v2keepflag]
        # create a flag for mastertable
        self.z0mgsFlag = ~cat['pgc_name'].mask
        self.z0mgs_cat = cat
        c = Column(self.z0mgsFlag,name='Z0MGSflag')
        #cat.add_column(c)
        cat.rename_column('ra','ra_z0mgs')
        cat.rename_column('dec','dec_z0mgs')
        cat.rename_column('cntr_01','cntr')
        cat.rename_column('galid_01','galid')
        cat.rename_column('major_01','major')         
        cat.remove_columns(['ra_01','dec_01'])
        cat.add_column(self.basictable['VFID'],index=0)
        # write out north file
        cat.write(outdir+file_root+'z0mgs.fits',format='fits',overwrite=True)
    def get_steer17(self):
        # match to GL's steer catalog
        steercat = homedir+'/research/Virgo/ancil-tables/Steer2017_cat_Virgo_field_H0_74_0.fits'
        self.steer = Table(fits.getdata(steercat))
        # GL suggests using
        # np.searchsorted(names1,names2)

        # I might otherwise do a loop
        # for each name in our catalog, look for match in steer catalog
        # probably astropy has a way to do this (topcat certainly does)
        
        # need to change this to match by coordinates
        #self.basic_with_steer = myjoinleft(self.basictable,self.steer,keys='NEDname')

        # create a Skycoord for each galaxy
        # self.catcoord is the main table
        scat =  SkyCoord(self.steer['raNED'],self.steer['decNED'],frame='icrs',unit='deg')
        # match steer to vf catalog
        idx, d2d, d3d = self.catcoord.match_to_catalog_sky(scat)

        # consider as matches any galaxies within 30 arcsec
        matchflag = d2d < 30./3600*u.deg

        
        # create null table with datatype like steer but length of vf
        self.temp = Table(np.empty_like(self.steer,dtype=self.steer.dtype,shape=(len(self.basictable))))

        self.temp[matchflag] = self.steer[idx[matchflag]]

        c = Column(self.cat['VFID'],name='VFID')
        self.temp.add_column(c,index=0)
        sorted_colnames = ['VFID','Type_NED','raNED','decNED','zNED',\
                           'Dall','Dall_err','Num. D Measures','Dmin','Dmax',\
                           'Dmedian','Dmean','Dmedian_Steer','Dmean_Steer']
        output_colnames = ['VFID','Type_NED','raNED','decNED','zNED',\
                           'Dall','Dall_err','number_D','Dmin','Dmax',\
                           'Dmedian','Dmean','Dmedian_Steer','Dmean_Steer']
        
        temp_sorted = Table(self.temp[sorted_colnames],names=output_colnames)
        outfile = outdir+file_root+'steer17.fits'            
        #self.basic_with_steer.write(outfile, format='fits',overwrite=True)
        temp_sorted.write(outfile, format='fits',overwrite=True)

    def hyperleda_table(self):
        # V2 change
        # removing objname from hyperleda table b/c this is in the main table
        # colnames = self.cat.colnames[0:43]
        self.cat.rename_column('vopt_1','vopt')
        self.cat.rename_column('type_1','type')        
        colnames = self.cat.colnames[1:44]

        self.write_table(colnames,outdir+file_root+'hyperleda.fits',format='fits')
        # write out HL columns in line-matched table
        pass
    def nsa_v1_table(self):
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[91:200]
        # RA and DEC get changed to RA_2, DEC_2 during the matching process
        # changing them back to native NSA values
        newcolnames = colnames.copy()
        newcolnames[2] = 'RA'
        newcolnames[3] = 'DEC'
        self.write_table(colnames,outdir+file_root+'nsa_v1.fits',format='fits',names=newcolnames)
    def nsa_v0_table(self):
        # write out NSA columns in line-matched table

        # V2 change
        # remove NSAflag at the end of this file
        # colnames = self.cat.colnames[204:348]
        colnames = self.cat.colnames[204:347]
        # removing some joint velocity
        colnames = self.cat.colnames[205:348]        
        # RA and DEC get changed to RA_2, DEC_2 during the matching process
        # changing them back to native NSA values
        newcolnames = colnames.copy()
        for i,n in enumerate(newcolnames):
            if n.find('_2') > -1:
                newcolnames[i] = n.strip('_2')
                
        newcolnames[2] = 'RA'
        newcolnames[3] = 'DEC'
        self.write_table(colnames,outdir+file_root+'nsa_v0.fits',format='fits',names=newcolnames)
    def agc_table(self):
        '''write out columns from agc catalog'''
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[46:83]
        self.write_table(colnames,outdir+file_root+'agc.fits',format='fits')

    def a100_table(self):
        '''write out columns from a100 catalog, plus HI deficiency'''
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[364:389]
        subtable = Table(self.cat[colnames])
        newtable = hstack([self.vfidtable,subtable])
        c0 = Column(self.HIdef,name='HIdef')
        c1 = Column(self.HIdef_flag,name='HIdef_flag')        
        c2 = Column(self.HIdef_bos,name='HIdef_bos')
        c3 = Column(self.HIdef_jones,name='HIdef_jon')
        c4 = Column(self.HIdef_jones_bytype,name='HIdef_bytype')        
        c5 = Column(self.cat['logMstarTaylor'],name='logMstarTaylor')
        newtable.add_columns([c0,c1,c2,c3,c4,c5])
        sdsscolnames = self.cat.colnames[394:482]
        subtable = Table(self.cat[sdsscolnames])
        newtable = hstack([newtable,subtable])

        # columns to remove
        remove_cols = ['lnLExp_r','lnLDeV_r','fracDev_g','fracDev_r','fracDev_i',\
                       'agcnum','icode','pcode','ipcode',\
                       'flags_u','flags_g','flags_r','flags_i','flags_z','flags',\
                       'G_Shao','I_Shao','gmi_Shao',\
                       'gamma_g','gamma_i']
        for c in remove_cols:
            newtable.remove_column(c)
        
        #newtable.add_columns(self.cat[sdsscolnames])
        newtable.write(outdir+file_root+'a100.fits',format='fits',overwrite=True)

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
    def ned_table(self):
        '''write out columns from ned query'''


        # make a column with NEDname, but no spaces in the name
        NEDname_nospace = []
        for i in range(len(self.cat)):
            NEDname_nospace.append(self.cat['NEDname'][i].replace(' ',''))
        
        colnames = ['NEDinput','NEDra','NEDdec','NEDname']
        subtable = Table(self.cat[colnames])
        newtable = hstack([self.basictable,subtable])
        c1 = Column(NEDname_nospace,name='NEDname_nospace')
        newtable.add_column(c1)
        newtable.write(outdir+file_root+'nedquery.fits',format='fits',overwrite=True)
    def update_spine_filenames_v1_2_v2(self):
        ''' 
        * read in v1 spine files and update filament names to remove Virgo 
        * spines are kept in a subfolder of v2 tables

        '''
        inputfiles = glob.glob(homedir+'/research/Virgo/tables-north/spines/filament*.fits')
        outdir = homedir+'/research/Virgo/tables-north/v2/spines/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        # remove Virgo from filename
        for f in inputfiles:
            spine_file = os.path.basename(f)
            out_file = os.path.join(outdir,spine_file.replace("_Filament","").replace("Virgo_",""))
            print('just checking, spine file = ',out_file)
            shutil.copy(f,out_file)

        # fix up all_filament_spines.fits as well
        infile =homedir+'/research/Virgo/tables-north/spines/all_filament_spines.fits'
        outfile = homedir+'/research/Virgo/tables-north/v2/spines/all_filament_spines.fits'
        alltab = Table.read(infile)

        # replace filament names
        newfilname = []
        for fname in alltab['filament']:
            newfilname.append(fname.replace('_Filament','').replace('Virgo_',''))
        alltab.remove_column('filament')

        c = Column(newfilname,name='filament')
        alltab.add_column(c)
        alltab.write(outfile,format='fits',overwrite=True)
        
        pass


    def write_table_v1(self,colnames,outfile,format=None,names=None):
        '''function for writing out a subset of columns into a new table.'''
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

        # remove RA, DEC and NEDname b/c these are duplicates from main table
        try:
            newtable.remove_columns(['RA'])
        except KeyError:
            print('no RA column to remove')
            pass
        try:
            newtable.remove_columns(['DEC'])
        except KeyError:
            print('no DEC column to remove')            
            pass
        try:
            newtable.remove_columns(['NEDname'])
        except KeyError:
            print('no NEDname column to remove')            
            pass
        # for NSA tables, their RA went to RA_2
        try:
            newtable.remove_columns(['RA_1'])
            newtable.rename_column('RA_2','RA')            
        except KeyError:
            print('no RA_1 column to remove')
            pass
        try:
            newtable.remove_columns(['DEC_1'])
            newtable.rename_column('DEC_2','DEC')            
        except KeyError:
            print('no DEC_1 column to remove')            
            pass
        
        newtable.write(outfile,format=format,overwrite=True)
    def write_table(self,colnames,outfile,format=None,names=None):
        '''function for writing out a subset of columns into a new table.  does not include basictable'''
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
        c = Column(self.cat['VFID'],name='VFID')
        subtable.add_column(c,index=0)
        subtable.write(outfile,format=format,overwrite=True)
        
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

    def merge_legacy_phot_tables(self):
        tabledir = os.path.join(homedir,'research','Virgo','tables-moustakas',"")
        vfdir = os.path.join(homedir,'research','Virgo','tables-north',"v1","")
        northcat = Table.read(tabledir+'vf_north_v1_main_dr9north.fits')
        southcat = Table.read(tabledir+'vf_north_v1_main_dr9south.fits')
        legcat = Table(np.empty(len(northcat),dtype=northcat.dtype))
        
        vmain = Table.read(vfdir+"vf_north_v1_main.fits")
        nflag = vmain['DEC'] >= 32.375
        legcat[nflag] = northcat[nflag]
        legcat[~nflag] = southcat[~nflag]
        legcatv2 = legcat[self.v2keepflag]
        legcatv2.rename_column('VFID','VFID_V1')
        c = Column(self.cat['VFID'],name='VFID')
        legcatv2.add_column(c,index=0)

        g = 22.5 - 2.5*np.log10(legcatv2['FLUX_G'])
        r = 22.5 - 2.5*np.log10(legcatv2['FLUX_R'])
        z = 22.5 - 2.5*np.log10(legcatv2['FLUX_Z'])
        d_pc = self.env['Vcosmic']/H0*1.e6
        const = 5*np.log10(d_pc) - 5
        MG = g - const
        MR = r - const
        MZ = z - const
        newtab = Table([g,r,z,MG,MR,MZ],names=['g','r','z','Mg','Mr','Mz'])
        legcatv2 = hstack([legcatv2,newtab])
        legcatv2.write(outdir+file_root+'legacy_dr9.fits',overwrite=True)
        
    def make_legacy_viewer_table(self):
        legacy_table = Table(self.cat['VFID','RA','DEC','radius'],names=['name','RA','DEC','radius'])
        legacy_table.write(outdir+file_root+'legacy_viewer.fits',format='fits',overwrite=True)

    def get_JMphot_table(self):
        '''convert John's custom photometry file to line-matched version   '''

        photfile = homedir+'/research/Virgo/legacy-phot/virgofilaments-v2-legacyphot.fits'

        ##
        # 2023-07-10 : Updating to use John's catalog from 4/29/2023
        ## 
        photfile = homedir+'/research/Virgo/legacy-phot/virgofilaments-v3-legacyphot.fits'

        # updating Feb 2024 to use John's catalog after correcting for the phot bug
        photfile = homedir+'/research/Virgo/legacy-phot/virgofilaments-v3b-legacyphot.fits'
        mef_table = fits.open(photfile)
        # changing to extension 2 for file that john sent on Aug 14, 2021
        #ephot = Table.read(photfile,1) # first hdu is the elliptical photometry
        ephot = Table(mef_table['ELLIPSE'].data) # first hdu is the elliptical photometry
        ephot1 = Table(mef_table['PARENT'].data) # this one has RA and DEC, which we need to separate N and S

        ephot_tractor = Table(mef_table['TRACTOR'].data) # this one has RA and DEC, which we need to separate N and S

        

        mef_table.close()
        
        # keep SGA_ID from parent table for now
        ephot.add_column(ephot1['SGA_ID'],index=0)

        # add info about parent/group status from PARENT table        
        # adding columns from parent catalog - columns associated with group
        #ephot.add_column(ephot1['GALAXY']) # this is already in the ellipse table
        ephot.add_column(ephot1['GROUP_NAME'])
        ephot.add_column(ephot1['GROUP_MULT'])
        ephot.add_column(ephot1['GROUP_PRIMARY'])

        # adding information from Tractor catalog
        # https://www.legacysurvey.org/dr9/bitmasks/#maskbits
        ephot.add_column(ephot_tractor['MASKBITS'])

        # calculate flags to see if galaxy is saturated in g,r,z
        # also track presence of nearby star
        # then add columns to ephot tables

        bitflag = [2,3,4,8,9,11]
        flagname= ['GSATURATE','RSATURATE','ZSATURATE','WISEM1','WISEM2','NEARBYSTAR']
        for i,b in enumerate(bitflag):
            maskflag = (ephot_tractor['MASKBITS'] & 2**b) == 2**b
            ephot.add_column(maskflag,name=flagname[i])


        # set up flag to indicate presence of photometry
        ephot.add_column(np.zeros(len(ephot),'bool'),index=0,name='photFlag')

        ### MAKE A TABLE - ALL ZEROS, WITH DATA TYPE LIKE EPHOT TABLE
        ### AND LENGTH OF BASICTABLE

        phottable = QTable(np.empty(len(self.basictable),dtype=ephot.dtype))
        phottable[ephot['VF_ID']] = ephot
        phottable['photFlag'][ephot['VF_ID']] = np.ones(len(ephot),'bool')        
        phottable.add_column(self.cat['VFID'],index=0)
        
        out_prefix = homedir+'/research/Virgo/tables-north/v2/vf_v2_'
        phottable.write(out_prefix+'legacy_ephot.fits',format='fits',overwrite=True)
        pass

    def get_magphys(self):
        tabledir = homedir+'/research/Virgo/tables-north/v2/vf_v2_'
        #self.magphys_lext = Table.read(tabledir+'magphys_legacyExt_17-Feb-2024.fits')
        #self.magphys_noz_lext = Table.read(tabledir+'magphys_nozband_legacyExt_17-Feb-2024.fits')
        #self.magphys_lext = Table.read(tabledir+'magphys_legacyExt_23-Mar-2024.fits')
        #self.magphys_noz_lext = Table.read(tabledir+'magphys_nozband_legacyExt_23-Mar-2024.fits')

        ################################################################################
        # fixed error in gatherMagphys.py in how I was parsing the percentile lines 
        ################################################################################       
        self.magphys_lext = Table.read(tabledir+'magphys_legacyExt_25-Apr-2024.fits')
        self.magphys_noz_lext = Table.read(tabledir+'magphys_nozband_legacyExt_25-Apr-2024.fits')

        # testing to see if the number of weird galaxies changes
        #self.magphys_lext = Table.read(tabledir+'magphys_legacyExt_18-Apr-2024.fits')
        #self.magphys_noz_lext = Table.read(tabledir+'magphys_nozband_legacyExt_18-Apr-2024.fits')

        
        #self.magphys_sext = Table.read(self.tabledir+self.tableprefix+'magphys_salimExt_11-Jul-2023.fits')


        outtab = tabledir+'magphys_legacyExt_mergedNS.fits'
        self.magphys = Table.read(tabledir+'magphys_legacyExt_17-Feb-2024.fits')
        self.magphys = Table.read(tabledir+'magphys_legacyExt_23-Mar-2024.fits')
        #self.magphys = Table.read(tabledir+'magphys_legacyExt_25-Apr-2024.fits')
        #self.magphys = Table.read(tabledir+'magphys_legacyExt_18-Apr-2024.fits') 
        Nflag = (self.maintable['DEC'] >= 32.375)
        Sflag = (self.maintable['DEC'] < 32.375)

        # combine results
        self.magphys[Nflag] = self.magphys_noz_lext[Nflag]
        
        # write out this table as magphys_final.fits
            
        self.magphys.write(outtab,format='fits',overwrite=True)

        
    def convert_kourkchi_v1_2_v2(self):
        # read in the v1 environment tables, remove bad sources, and save for v2
        prefix = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_'
        out_prefix = homedir+'/research/Virgo/tables-north/v2/vf_v2_'
        
        alltables = ['main_Tempelgroups_infos.fits','matchTempel_groupinfo.fits']

        
        for tname in alltables:
            mytab = Table.read(prefix+tname)
            # cut table
            mytab = mytab[self.v2keepflag]
            # remove v1 VFIDs
            mytab.remove_column('VFID')
            # add v2 VFIDs in first columns
            mytab.add_column(self.cat['VFID'],index=0)
            mytab.write(out_prefix+tname.replace('main_',''),format='fits',overwrite=True)

        # BV Kourchi table
        # change name to include v2, and change last column from 0/1 to boolean
        # I downloaded the input table from google drive v2 table
        infile = homedir+'/research/Virgo/tables-north/BV-kourchi-tables/vf_kourkchi_galaxies.fits'
        ktab = Table.read(infile)
        # keep v2 members
        ktab = ktab[self.v2keepflag]
        ktab.remove_column('VFID')
        # add v2 VFIDs in first columns
        ktab.add_column(self.cat['VFID'],index=0)        
        # convert Kflag to boolean
        newflag = np.array(ktab['Kflag'],'bool')
        ktab.remove_column('Kflag')
        ktab.add_column(newflag,name='Kflag')
        outfile = homedir+'/research/Virgo/tables-north/v2/vf_v2_kourkchi_galaxies.fits'
        ktab.write(outfile,format='fits',overwrite=True)
    def convert_rphot_v1_2_v2(self):
        # read in the v1 environment tables, remove bad sources, and save for v2
        prefix = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_'
        out_prefix = homedir+'/research/Virgo/tables-north/v2/vf_v2_'
        
        alltables = ['r_photometry.fits']

        for tname in alltables:
            mytab = Table.read(prefix+tname)
            # cut table
            mytab = mytab[self.v2keepflag]
            # remove v1 VFIDs
            mytab.remove_column('VFID')
            # add v2 VFIDs in first columns
            mytab.add_column(self.cat['VFID'],index=0)
            
            mytab.write(out_prefix+tname.replace('main_',''),format='fits',overwrite=True)
    def combine_env_v1_2_v2(self):
        # read in the v1 environment tables, remove bad sources, and save for v2
        # also need to pare down the columns
        prefix = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_'
        out_prefix = homedir+'/research/Virgo/tables-north/v2/vf_v2_'
        
        alltables = ['main_env_prop_H0_74_0_Mr_max_-15_7.fits','main_finalenvironments.fits',
                     'main_filament_membership_allgalaxies.fits','main_cluster_membership.fits']


                     
        # make a sub table of each table, and then do an hstack
        env = Table.read(prefix+alltables[0])
        keepcols = ['VFID','DM','SGX','SGY','SGZ',\
                    'nCI08','nCI08_err',\
                    'distSGX_Virgo','distSGY_Virgo','distSGZ_Virgo',\
                    'n5th_2D','n5th_2D_err','n5th','n5th_err',\
                    'Vcosmic','Vmodel']
        envcut = env[keepcols]

        vmedian = Column(env['Dmedian']*74,name='Vmedian')
        envcut.add_column(vmedian)
        
        finalenv = Table.read(prefix+alltables[1])
        poor_group_memb = finalenv['flag_gr'] == 1
        rich_group_memb = finalenv['flag_gr'] == 2
        pure_field_memb = finalenv['flag_pf'] == 1
        finalenvcut = Table([poor_group_memb,rich_group_memb,pure_field_memb],names=['poor_group_memb','rich_group_memb','pure_field'])

        # rename columns in this one
        filmemb = Table.read(prefix+alltables[2])
        filmemb.rename_column('filament_dist_2D','nearest_filament_dist_2D')
        filmemb.rename_column('filament_dist_3D','nearest_filament_dist_3D')        

        keepcols = ['nearest_filament_dist_2D','nearest_filament_dist_3D','filament',\
                    'filament_PA','orientation_wrt_filament','filament_member']
        filmembcut = filmemb[keepcols]

        # rename filament
        for i in range(len(filmembcut)):
            filmembcut['filament'][i] = filmembcut['filament'][i].replace('_Filament','').replace('Virgo_','')

        
        clustermemb = Table.read(prefix+alltables[3])                        
        clustermemb.remove_column('VFID')
        clustermemb.remove_column('vr')
        clustermemb.rename_column('mem','cluster_member')

        envtable = hstack([envcut,finalenvcut,filmembcut,clustermemb])
        # this table has v1 IDs - need to replace with v2
        # but this table hasn't been cut to v2 length yet
        envtable.remove_column('VFID')
        # cut table
        envtable = envtable[self.v2keepflag]
        # add in v2 ids
        envtable.add_column(self.cat['VFID'],index=0)
        envtable.write(out_prefix+'environment.fits',format='fits',overwrite=True)
        self.env = envtable
        # write out the filament distances in a separate table
        keepcols = []
        # list of filaments is in virgoCommon file
        for f in filaments:
            keepcols.append('dist_2D_'+f)
            keepcols.append('dist_3D_'+f)            

        colnames = []
        for n in keepcols:
            colnames.append(n.replace('_Filament','').replace('Virgo_',''))
        fildist = Table(filmemb[keepcols],names=colnames)
        # this table has v1 IDs - need to replace with v2
        # but this table hasn't been cut to v2 length yet
        
        #fildist.remove_column('VFID') # apparently this doesn't have VFID

        
        # cut table
        fildist = fildist[self.v2keepflag]
        # add in v2 ids
        fildist.add_column(self.cat['VFID'],index=0)
        
        fildist.write(out_prefix+'filament_distances.fits',format='fits',overwrite=True)
    def get_extinction(self):
        '''
        GOAL:
        * get extinction coefficients for each galaxy

        PROCEDURE:
        * use extinction coefficients from Legacy (g-W4)
        * use extinction coefficients from Wyder+2007 for GALEX
        * get E(B-V) from irsa extiction
        '''

        ### NOTE: THESE TABLES ARE SAVED TO A SEPARATE DIRECTORY B/C
        ### THEY ARE NOT OF GENERAL INTEREST TO VF COLLABORATION
        out_prefix = homedir+'/research/Virgo/tables-north/v2/vf_v2_'        
        size = 2*np.ones(len(self.cat),'f')
        
        irsa_input = Table([self.cat['RA'],self.cat['DEC'],size],names=['RA','DEC','size'])
        self.irsa_input = irsa_input
        irsa_input.write(homedir+'/research/Virgo/tables/v2/vf_v2_irsa_input.tbl',format='ipac',overwrite=True)

        # upload this to https://irsa.ipac.caltech.edu/applications/DUST/
        #
        # this gives results from Schlegel, Finkbeiner and Davis, 1998, ApJ, 500, 525
        # and
        # Schlafly and Finkbeiner, 2011, ApJ, 737, 103
        #
        # download results
        # once that is uploaded, download the results to:
        infile = homedir+'/research/Virgo/tables/v2/vf_v2_irsa_extinction.tbl'
        etab = Table.read(infile,format='ipac')
        self.etab = etab
        print('length of etab = ',len(etab),len(self.cat['VFID']))
        etab.add_column(self.cat['VFID'],index=0)


        # these must be values from legacy survey?
        # galex coefficients are from Wyder+07, referenced in legacy code below
        # https://github.com/moustakas/impro/blob/master/pro/galaxy/im_galex_to_maggies.pro#L64-L81
        #
        # some details at: https://www.legacysurvey.org/dr9/catalogs/#galactic-extinction-coefficients
        #
        filters = ['FUV','NUV','G','R','Z','W1','W2','W3','W4']
        
        Rvalues = np.array([8.24, 8.20, 3.214, 2.165, 1.211, 0.184, 0.113, 0.0241, 0.00910],'d')

        EBV_colnames = ['E_B_V_SandF','E_B_V_SFD']
        for cname in EBV_colnames:
            for i,f in enumerate(filters):
                alambda = np.array(etab[cname])*Rvalues[i]
                columnname = f"A({f})_{cname.replace('E_B_V_','')}"
                newcol = Column(alambda,name=columnname,unit='mag')
                etab.add_column(newcol)

        # add Salim +2016 values
        # https://iopscience.iop.org/article/10.3847/0067-0049/227/1/2#apjsaa4425s3
        # section 3.3
        # they use Peek & Schiminovich (2013) for UV and Yuan+2013 for optical
        # they use SFD extinction
        RvaluesS16 = np.array([8.24, 8.20, 3.30, 2.31, 1.29, 0.184, 0.113, 0.0241, 0.00910],'d')
        ebv = np.array(etab['E_B_V_SFD'])
        for i,f in enumerate(filters):

            if i == 0:
                alambda = 10.47 + 8.59*ebv -82.2*ebv**2
                alambda =  np.minimum(alambda, 0.2)
            elif i == 1:
                alambda = 9.36 + 14.3*ebv -82.2*ebv**2
                alambda =  np.minimum(alambda, 0.2)

            else:
                alambda = ebv*Rvalues[i]                
            columnname = f"A({f})_S16"
            newcol = Column(alambda,name=columnname,unit='mag')
            etab.add_column(newcol)
        

        
        etab.write(out_prefix+'extinction.fits',format='fits',overwrite=True)
        
        
if __name__ == '__main__':
    ###################################################################
    #### SET UP ARGPARSE
    #### ADD OPTION FOR NORTH ONLY WHEN MAKING TABLES
    #### (NO GUARANTEE THAT THIS WORKS FOR THE FULL TABLES ANYMORE!)
    ###################################################################

    parser = argparse.ArgumentParser(description ='write out subtables for virgo filaments catalog')
    parser.add_argument('--version',dest = 'version', default='v2',help='version of tables. default is v2')
    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--north',dest = 'north', action='store_true',help='keep DEC > -1 galaxies')
    parser.add_argument('--hav1',dest = 'hav1', action='store_true',default=False,help='set this if there is halpha data post transition to virgo v1 catalogs')
    
    args = parser.parse_args()

    # keep DEC > -1 galaxies only
    if (args.north):
        NORTH_ONLY = True
    else:
        NORTH_ONLY = False

    ###################################################################
    #### SET VERSION NUMBER
    ###################################################################
    version = args.version
    outfile_suffix = '_'+args.version
    ###################################################################
    #### INPUT FILES
    ###################################################################
    masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample_wNEDname_'+version+'.fits'
    masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample_wNEDname_v1_evcc.fits'
    z0mgs_cat = homedir+'/research/Virgo/ancil-tables/vf_v1_z0mgs_30arcsec_120220.tbl'
    unwise_cat = homedir+'/research/Virgo/ancil-tables/vf_north_v1_unwise_dl_30arcsec_20201205.fits'
    ###################################################################
    #### SET UP TABLE DIRECTORIES
    ###################################################################
    outdir = homedir+'/research/Virgo/tables/'+version+'/'
    if NORTH_ONLY:
        outdir = homedir+'/research/Virgo/tables-north/'+version+'/'

    ###################################################################
    #### SET TABLE PREFIX
    ###################################################################
    file_root = 'vf_'+version+'_'
    if NORTH_ONLY:
        file_root = 'vf_north_'+version+'_'
        # dropping the north notation b/c we don't plan to expand to south
        file_root = 'vf_'+version+'_'        
        
    ###################################################################
    #### CREATE TABLE DIRECTORY IF IT DOESN'T EXIST
    ###################################################################
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    print('outdir = ',outdir)
    c = catalog(masterfile,hav1=args.hav1,version=args.version)
    c.runall()
