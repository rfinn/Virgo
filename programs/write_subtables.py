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
* when making new version, need to download z0mgs

* when making new version, need to rematch to unwise

'''
import os
import numpy as np
import time
import argparse

from astropy.io import fits
from astropy.table import Table, QTable, join, hstack, vstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy import constants
from astropy import units as u
from astropy.cosmology import WMAP9 as cosmo

from astroquery.ned import Ned

from matplotlib import pyplot as plt
homedir = os.getenv("HOME")
#sys.path.append(homedir+'/github/appss/')
#from join_catalogs import make_new_cats, join_cats

from virgoCommon import *

###################################################################
#### SET UP ARGPARSE
#### ADD OPTION FOR NORTH ONLY WHEN MAKING TABLES
#### (NO GUARANTEE THAT THIS WORKS FOR THE FULL TABLES ANYMORE!)
###################################################################

parser = argparse.ArgumentParser(description ='write out subtables for virgo filaments catalog')
parser.add_argument('--version',dest = 'version', default='v1',help='version of tables. default is v1')
#parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
parser.add_argument('--north',dest = 'north', action='store_true',help='keep DEC > -1 galaxies')

args = parser.parse_args()
outfile_suffix = '_'+args.version


###################################################################
#### CONSTANTS
###################################################################
sdsspixelscale=0.396127#conversion for isophotal radii from pixels to arcseconds
H0 = 70.
h = H0/100 #H0/100

###################################################################
#### SET VERSION NUMBER
###################################################################
version = args.version

###################################################################
#### INPUT FILES
###################################################################
masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample_wNEDname_'+version+'.fits'
masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample_wNEDname_'+version+'_evcc.fits'
z0mgs_cat = '/home/rfinn/research/Virgo/ancil-tables/vf_v1_z0mgs_30arcsec_120220.tbl'
unwise_cat = '/home/rfinn/research/Virgo/ancil-tables/vf_north_v1_unwise_dl_30arcsec_20201205.fits'
###################################################################
#### SET UP DIRECTORIES
###################################################################
outdir = homedir+'/research/Virgo/tables/'+version+'/'
file_root = 'vf_'+version+'_'




# keep DEC > -1 galaxies only
if (args.north):
    NORTH_ONLY = True
else:
    NORTH_ONLY = False
if NORTH_ONLY:
    outdir = homedir+'/research/Virgo/tables-north/'+version+'/'
    file_root = 'vf_north_'+version+'_'

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
        self.unwise = Table.read(unwise_cat)
        self.unwiseFlag = (self.unwise['x'] > 0)
        # clean catalog
        unwise_keepcols = self.unwise.colnames[22:75]
        self.unwise = self.unwise[unwise_keepcols]
        self.unwise.write(outdir+file_root+'unwise.fits',format='fits',overwrite=True)
    def keep_north(self, cat=None):
        # cut the catalog to keep Dec > -1 only
        if cat is None:
            cat = self.cat
        else:
            cat = cat
        keepnorth = cat['DEC'] > -1.3
        return cat[keepnorth], keepnorth
    def runall(self):
        #self.catalog_for_z0MGS()
        self.catalog_for_z0MGS()
        self.get_unwise()        
        self.get_CO()
        self.get_halpha()        
        self.get_steer17()
        self.get_radius()
        self.get_sfr()
        self.get_HIdef()
        self.get_z0MGS_flag()        
        self.main_table()

        self.hyperleda_table()
        self.nsa_table()
        self.nsa_v0_table()
        self.agc_table()        
        self.a100_table()
        self.a100_sdss_table()
        self.a100_unwise_table()

        self.get_size_for_JM()
        self.ned_table() # NED input, ra, dec, and NEDname
        self.print_stats()        
        

        pass
    def print_stats(self):
        print('Number in sample = ',len(self.cat))
        print('Number with CO data = ',sum(self.coflag))
        print('Number with A100 data = %i (%.3f)'%(sum(self.cat['A100flag']),sum(self.cat['A100flag'])/len(self.cat['A100flag'])))
        print('Number with z0MGS matches = %i (%.3f)'%(sum(self.z0mgsFlag),sum(self.z0mgsFlag)/len(self.z0mgsFlag)))
        print('Number with steer17 matches = %i (%.3f)'%(sum(self.steerFlag),sum(self.steerFlag)/len(self.steerFlag)))
        try:
            f = self.unwiseFlag
            print('Number with unwise matches = %i (%.3f)'%(sum(f),sum(f)/len(f)))        
            f = self.unwiseFlag & (self.cat['NSAflag'] | self.cat['NSA0flag'])
            print("Number with unWISE and NSA = %i (%.3f)"%(sum(f),sum(f)/len(f)))
        except AttributeError:
            pass
        print("CO SOURCES")
        nco = sum(self.coflag)
        f = self.coflag & self.z0mgsFlag
        print("\tNumber of CO sources in z0MGS = %i (%.2f)"%(sum(f),sum(f)/nco))
        f = self.coflag & self.steerFlag
        print("\tNumber of CO sources in Steer = %i (%.2f)"%(sum(f),sum(f)/nco))
        f = self.coflag & self.z0mgsFlag & self.steerFlag
        print("\tNumber of CO sources in z0MGS+Steer = %i (%.2f)"%(sum(f),sum(f)/nco))
        try:
            f = self.coflag & self.unwiseFlag & (self.cat['NSAflag'] | self.cat['NSA0flag'])
            print("\tNumber of CO sources with unWISE and NSA = %i (%.2f)"%(sum(f),sum(f)/nco))
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
        colnames = ['VFID','RA','DEC','vr','radius','radius_flag','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag']
        #self.maintable = self.cat['VFID','RA','DEC','vr','objname','NSAID','NSAID_2','AGC','NEDname','HLflag','NSAflag','NSA0flag','A100flag']
        
        self.maintable = self.cat[colnames]
        self.maintable.rename_column('NSAID_2','NSAIDV0')
        self.maintable.rename_column('NSA0flag','NSAV0flag')        
        

        c1 = Column(self.coflag,name='COflag')
        c1a = Column((self.hatable['HAflag']), name='HAflag')
        c1b = Column(self.hatable['HAobsflag'],name='HAobsflag')
        c2 = Column(self.z0mgsFlag,name='Z0MGSflag')
        c3 = Column(self.steerFlag,name='Steerflag')

        # removing unwise for now - need to update code
        
        c4 = Column(self.unwiseFlag,name='unwiseflag')        
        self.maintable.add_columns([c1,c1a,c1b,c2,c3,c4])
        #self.maintable.add_columns([c1,c1a,c1b,c2,c3])

        
        nedname=[]
        for i in range(len(self.maintable)):
            nedname.append(self.maintable['VFID'][i]+'-'+self.maintable['NEDname'][i].replace(" ","").replace("[","").replace("]","").replace("/",""))
        
        c0= Column(nedname,name='prefix')
        self.maintable.add_column(c0)

        # add a column called "name" for use with the legacy image server
        # this will be the VFID
        c0= Column(self.cat['VFID'],name='name')
        self.maintable.add_column(c0)

        
        # - 2MASS
        # - z0MGS
        # - unWISE
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


    def get_CO(self,match_by_coords=False,match_by_name=True):
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
                plt.figure()
                plt.plot(self.testtable['RA_1'][~self.comatchflag],self.testtable['DEC_1'][~self.comatchflag],'bo')


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
    def get_halpha_old(self,halphafile=None):
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
        #print(len(self.maintable),len(self.ha),len(self.hatable))        
        self.haflag = ~self.hatable['Date Obs'].mask
        self.hatable.add_column(Column(self.haflag),name='haflag')
        self.hatable.write(outdir+file_root+'ha.fits',format='fits',overwrite=True)

    def get_halpha(self,halphafile=None):
        # read in Halpha observing summary file
        if halphafile is None:
            #infile = '/home/rfinn/research/Virgo/Halpha/observing-summary-Halpha-latest.csv'
            #infile = '/home/rfinn/research/Virgo/Halpha/observing-summary-Halpha-clean-04Jun2020.csv'
            infile = '/home/rfinn/research/Virgo/halpha-tables/halpha-05Jul2020.fits'

        else:
            infile = halphafile
        self.ha = Table.read(infile,format='fits')            


        unique, counts = duplicates(self.ha,'VFID')
        print("Halpha sources that are listed multiple times in the halpha file:")
        print(unique[counts>1])

        ### REMOVE DUPLICATE ROWS FOR NOW
        ### getting rid of second in each pair for now...
        VFID = ['VFID2080','VFID2154','VFID2145','VFID2136','VFID2060',\
                'VFID2157','VFID2076','VFID2089','VFID2095','VFID2095','VFID2141']
        PID = ['v20p35','v17p12','v17p12','v17p12','v20p35',\
                'v17p12','v18p54','v20p35','v20p35','v18p54','v17p12']
        for i in range(len(VFID)):
            idel = np.where((self.ha['VFID'] == VFID[i]) & (self.ha['POINTING'] == PID[i]))
            #print(idel)
            self.ha.remove_rows(idel)
        print('number of lines in ha file after deleting rows = ',len(self.ha))

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
                     'ASYMR_FLAG','ASYMHA_FLAG','OVERSTAR_FLAG','OVERGAL_FLAG',\
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

        ### ADD VFID FOR ALL
        self.hatable['VFID'] = self.basictable['VFID']        
        ### ADD RA FOR ALL
        self.hatable['RA'] = self.basictable['RA']
        ### ADD DEC FOR ALL
        self.hatable['DEC'] = self.basictable['DEC']        

        ### ADD HALPHA SOURCES TO THEIR CORRESPONDING ROWS

        for i in range(len(self.ha)):
            #print(i,self.ha['VFID'][i])

            ### RAF, 12/14/20
            # will need to update this to use V1 tables, or else
            # I will have to redo all of the Halpha analysis, which would be painful
            flag =(self.hatable['VFID'] == self.ha['VFID'][i])
            self.hatable[flag] = self.ha[i]

        
        #self.hatable = join(self.basictable,self.ha,keys='VFID',join_type='left')
        #print('Halpha table lengths')
        #rename_cols = ['RA','DEC','vr']
        #print(len(self.maintable),len(self.ha),len(self.hatable))
        print('calculating snr')
        self.haflag = ((self.hatable['HF_R24']/self.hatable['HF_R24_ERR']) > 1.)
        print('setting HAflag')
        self.hatable['HAflag'] = self.haflag
        #self.hatable['HAflag'].description = 'halpha snr > 1'
        print('writing hafile')
        #fits.writeto(outdir+file_root+'ha.fits',np.array(self.hatable),overwrite=True)
        self.hatable.write(outdir+file_root+'ha.fits',format='fits',overwrite=True)
        print('finished writing hafile')                
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
        #cat = Table.read('/home/rfinn/research/Virgo/tables/vf_z0mgs_30arcsec_051920.tbl',format='ipac')
        #cat = Table.read('/home/rfinn/research/Virgo/tables/vf_v1_z0mgs_30arcsec_102820.tbl',format='ipac')
        #cat = Table.read('/home/rfinn/research/Virgo/tables/vf_v1_z0mgs_10arcsec_102920_9153gal.tbl',format='ipac')
        cat = Table.read(z0mgs_cat,format='ipac')
        print('number of lines in z0MGS cat = ',len(cat))
        print('length of keepnorth_flag = ',len(self.keepnorth_flag))
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

        # join with basic

        self.basic_with_steer = hstack([self.basictable,self.temp])

        # some of the NED names change over time, so that's not a reliable way to match
        # all of the galaxies from Benedetta's low vr catalog should 
        self.steerFlag = self.basic_with_steer['Num. D Measures'] > 0
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
        '''write out columns from a100 catalog, plus HI deficiency'''
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[352:377]
        subtable = Table(self.cat[colnames])
        newtable = hstack([self.basictable,subtable])
        c0 = Column(self.HIdef,name='HIdef')
        c1 = Column(self.HIdef_flag,name='HIdef_flag')        
        c2 = Column(self.HIdef_bos,name='HIdef_bos')
        c3 = Column(self.HIdef_jones,name='HIdef_jon')
        c4 = Column(self.HIdef_jones_bytype,name='HIdef_bytype')        

        newtable.add_columns([c0,c1,c2,c3,c4])
        newtable.write(outdir+file_root+'a100.fits',format='fits',overwrite=True)
    def agc_table(self):
        '''write out columns from agc catalog'''
        # write out NSA columns in line-matched table
        colnames = self.cat.colnames[44:83]
        self.write_table(colnames,outdir+file_root+'agc.fits',format='fits')

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

    def write_table(self,colnames,outfile,format=None,names=None):
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
    #c.runall()
