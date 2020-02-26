#!/usr/bin/env python

'''
GOAL:
- create a super list of all galaxies in

USAGE:
- example will print out when you run

- to create new mastertable
    s = sample()
    s.get_smart()

- to create new cutouts, run from supersample directory
- make sure there is a cutouts and plots subdirectory
    t = fulltable()
    t.plot_all()
    



'''
import numpy as np
import os
import sys

from astropy.io import fits, ascii
from astropy.table import Table, join, hstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.visualization import simple_norm

from astroquery.skyview import SkyView

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

from urllib.parse import urlencode
from urllib.request import urlretrieve

import pandas as pd

homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/APPSS/')
from join_catalogs import make_new_cats, join_cats


## import argparse

## parser = argparse.ArgumentParser(description ='Match the Halpha observations with Virgo NSA catalog')

## parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
## parser.add_argument('--write-fits',dest = 'writefits', action='store_true',help='write out fits version of observing summary file?')
## parser.add_argument('--input',dest = 'input', default='Observing-Summary-Halpha-good-latest.csv',help='write out fits version of observing summary file?')
        
## args = parser.parse_args()


# max offset in arcsec b/w sources from different catalogs
# and still be considered the same source
max_match_offset = 5.


## BOUNDARIES OF SURVEY REGION
decmin = -1.2
decmin = -35
decmax = 75 
ramax = 280.
ramin = 100. 
vmax = 3300.
vmin = 500.

## LEGACY SURVEY
legacy_pixel_scale = 0.262 # arcsec/pixel, don't actually use this
image_size = 60 # image size to download from legacy survey, in pixels

def duplicates(table,column,flag=None):
    if flag is None:
        unique, counts = np.unique(table[column], return_counts=True)
    elif flag is not None:
        unique, counts = np.unique(table[column][flag], return_counts=True)
    print('number of duplicates = ',sum(counts > 1))
    #print('duplicates = ',unique[counts > 1])
    return unique, counts

def getlegacyimages(ra,dec):
    '''
    new function to download images in one fell swoop

    doing this so we can make cutouts of entire catalog, and then check each galaxy by hand

    This will need to be re-run each time the matching is altered, and the kitchen-sink catalog changes.
    
    '''
    for i in range(len(ra)):
        # name image files by ra and dec of galaxy
        gra = '%.5f'%(ra[i]) # accuracy is of order .1"
        gdec = '%.5f'%(dec[i])
        galnumber = gra+'-'+gdec
        rootname = 'cutouts/legacy-im-'+str(galnumber)+'-'+str(image_size)
        jpeg_name = rootname+'.jpg'

        fits_name = rootname+'.fits'
        # check if images already exist
        # if not download images
        if not(os.path.exists(jpeg_name)):
            print('retrieving ',jpeg_name)
            url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra[i])+'&dec='+str(dec[i])+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
            urlretrieve(url, jpeg_name)
        else:
            print('previously downloaded ',jpeg_name)
        if not(os.path.exists(fits_name)):
            print('retrieving ',fits_name)
            url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra[i])+'&dec='+str(dec[i])+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
            urlretrieve(url, fits_name)
        else:
            print('previously downloaded ',fits_name)
    pass

def getlegacy(ra1,dec1,ra2=None,dec2=None, ra3=None,dec3=None,agcflag=False,onlyflag=False,jpeg=True):
    gra = '%.5f'%(ra1) # accuracy is of order .1"
    gdec = '%.5f'%(dec1)
    galnumber = gra+'-'+gdec
    rootname = 'cutouts/legacy-im-'+str(galnumber)+'-'+str(image_size)
    jpeg_name = rootname+'.jpg'

    fits_name = rootname+'.fits'

    # check if images already exist
    # if not download images
    if not(os.path.exists(jpeg_name)):
        print('downloading image ',jpeg_name)
        url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
        urlretrieve(url, jpeg_name)
    #else:
    #    pass
    #    print('found image ',jpeg_name)
    if not(os.path.exists(fits_name)):
        print('downloading image ',fits_name)
        url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
        urlretrieve(url, fits_name)
    #else:
    #    pass
    #    print('found image ',fits_name)
            
    try:
        t,h = fits.getdata(fits_name,header=True)
    except IndexError:
        print('problem accessing image')
        print(fits_name)
        url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
        print(url)
        return None
    '''
        try: # try redownloading it
            print('downloading image ',jpeg_name)
            url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
            urlretrieve(url, jpeg_name)
            print('downloading image ',fits_name)
            url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
            urlretrieve(url, fits_name)
            t,h = fits.getdata(fits_name,header=True)
        except IndexError:
    '''
    
    # write out r-band image
    # nevermind - John M figured out how to use MEF with WCS
    #fits.writeto('r-test.fits',t[1],header=h,overwrite=True)
    if np.mean(t[1]) == 0:
        return None
    norm = simple_norm(t[1],stretch='asinh',percent=99.5)
    if jpeg:
        t = Image.open(jpeg_name)
        plt.imshow(t,origin='lower')
    else:
        plt.imshow(t[1],origin='upper',cmap='gray_r', norm=norm)

    # update this to plot the jpeg color image
    
    dx=20
    if agcflag:
        colors = ['cyan','blue','red']
    else:
        colors = ['red','blue','cyan']
    if onlyflag:
        colors = ['k','k','k']
    if (ra2 is not(None)) & (ra3 is not(None)):
        ra = np.array([ra1,ra2, ra3])
        dec = np.array([dec1,dec2, dec3])   
        #w = WCS('r-test.fits')
        #px,py = w.wcs_world2pix(ra,dec)
        w = WCS(fits_name,naxis=2)
        px,py = w.wcs_world2pix(ra,dec,1)
        #print(px,py)
        r1 = Rectangle((px[0]-dx/2, py[0]-dx/2), dx, dx, edgecolor=colors[0], facecolor='none')
        dx=17.5
        #r2 = Rectangle((px[1]-dx/2, py[1]-dx/2), dx, dx, edgecolor=colors[1], facecolor='none')
        #dx=15
        #r3 = Rectangle((px[2]-dx/2, py[2]-dx/2), dx, dx, edgecolor=colors[2], facecolor='none')
        plt.gca().add_patch(r1)
        #plt.gca().add_patch(r2)
        #plt.gca().add_patch(r3)
        return w


class sample:
    def __init__(self, max_match_offset=7.5):
        
        ## read in my HL catalog
        '''
        Hyperleda query:  http://leda.univ-lyon1.fr/fullsql.html

        parameters described here: http://leda.univ-lyon1.fr/leda/meandata.html

        SQL QUERY:
        
        select
        objname,objtype,de2000,al2000,v,e_v,vopt,e_vopt,vrad,e_vrad,bt,e_bt,type,bar,ring,multiple,compactness,t,e_t,logd25,e_logd25,logr25,e_logr25,pa,incl,logdc,btc,itc,ubtc,bvtc,m21c,hic,mabs,agnclass,kt,e_kt,it,e_it,ut,vt,mfir,e_ut,e_vt, modz, e_modz, mod0, e_mod0,vmaxg, e_vmaxg, vmaxs,e_vmaxs,vdis,e_vdis
        
        where

        de2000 > -35 and de2000 < 75 and  al2000 < 280./360.*24.  and al2000 > 100./360.*24. and v < 3300 and v > 500 and objtype='G'

        - output as csv, separator is other, ,


        
        not using Gialuca's catalog b/c I'm not sure if there were other cuts made already
        gl = fits.getdata('/Users/rfinn/research/VirgoFilaments/Gianluca/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')
    
        '''
        self.max_match_offset = max_match_offset

        ################################################################
        ## READ IN HYPERLEDA CATALO
        ################################################################
        hlfile = homedir+'/github/Virgo/tables/hyperleda-finn-09dec2019-full.csv'
        hlfile = homedir+'/github/Virgo/tables/hyperleda-finn-05Feb20.csv'
        hlfile = homedir+'/github/Virgo/tables/hyperleda-finn-24Feb20.csv'        
        self.hl = ascii.read(hlfile)
        self.hl = Table(self.hl)
        self.cull_hl()
        c1 = Column(self.hl['al2000']*15,name='RAdeg')
        self.hl.add_column(c1)
        self.hl.write('temp.fits',format='fits',overwrite=True)

        self.hl = fits.getdata('temp.fits')

        os.remove('temp.fits')

        ################################################################
        ## READ IN NED CATALOG
        ################################################################        
        #  downloaded in 10-dec-2019
        '''
        Search by by Parameters
        decmin = -1.2
        decmax = 75
        ramax = 280.
        ramin = 100.
        vmax = 3300.
        vmin = 500.
        object = Galaxy
        
        RA has to be in hours
        ramin = 6.6666666667
        ramax = 18.666666667

        output = text, ascii, bar separated
        velocity lower limit = -99

        Dowloaded text file - there was garbage at top of file (unnecessary info) that I deleted.

        Saved as

        /Users/rfinn/github/Virgo/tables/ned-noprolog-10dec2019.txt
        '''

        '''
        nedfile = homedir+'/github/Virgo/tables/ned-noprolog-10dec2019.txt'
        self.ned = ascii.read(nedfile,delimiter='|')

        # having issues with ned, maybe because of masked array?
        # going to write out fits and read it back in
        #
        self.ned.write('temp.fits',format='fits',overwrite=True)
        self.ned = fits.getdata('temp.fits')
        os.remove('temp.fits')
        self.cull_ned()

        '''
        ################################################################
        ## read in NSA catalog
        #  using newest version of NSA
        ################################################################
        '''
        using new NSA catalog
        '''
        nfile = homedir+'/research/NSA/nsa_v1_0_1.fits'
        self.nsa = fits.getdata(nfile)
        self.cull_nsa()
        #self.nsa = Table(self.nsa)

        ################################################################        
        ## read in AGC catalog
        ################################################################
        # got a new version from Martha on 11/19/19
        #agcfile = '/Users/rfinn/research/AGC/agcm1.sh191118.fits'
        #self.agc = fits.getdata(agcfile)
        #self.cull_agc()

        # using full agc so I can remove HIonly
        agcfile = '/Users/rfinn/research/AGC/agcnorthminus1.2019Sep24.fits'
        agcfile = homedir+'/research/AGC/agcnorthminus1.fits'        
        self.agc = fits.getdata(agcfile)
        self.agc = Table(self.agc)
        self.cull_agc_full()
        #self.cull_agc()
        #self.agc = Table(self.agc)
        
        # flags to track nsa matches to HL and AGC
        self.hl_2_nsa_matchflag = np.zeros(len(self.hl['al2000']),'bool')
        self.hl_2_agc_matchflag = np.zeros(len(self.hl['de2000']),'bool')

        # flags to track nsa matches to HL and AGC
        self.nsa_2_hl_matchflag = np.zeros(len(self.nsa['RA']),'bool')
        self.nsa_2_agc_matchflag = np.zeros(len(self.nsa['RA']),'bool')

        # flags to track AGC matches to HL and NSA
        self.agc_2_hl_matchflag = np.zeros(len(self.agc[self.agc_ra_key]),'bool')
        self.agc_2_nsa_matchflag = np.zeros(len(self.agc[self.agc_dec_key]),'bool')

        self.hcoord = SkyCoord(self.hl['al2000']*u.hr,self.hl['de2000']*u.deg,frame='icrs')
        self.ncoord = SkyCoord(self.nsa['RA']*u.deg,self.nsa['DEC']*u.deg,frame='icrs')
        self.acoord = SkyCoord(self.agc[self.agc_ra_key]*u.deg,self.agc[self.agc_dec_key]*u.deg, frame='icrs')
        #self.nedcoord = SkyCoord(self.ned['RA']*u.deg,self.ned['DEC']*u.deg, frame='icrs')

        self.hvel = self.hl['v'] # mean of optical and radio velocities
        self.nvel = self.nsa['Z']*3.e5
        self.avel = self.agc_vbest
        #self.nedvel = self.ned['Velocity']
        
    def run_it(self, maxoffset=10):
        self.max_match_offset = maxoffset
        self.match_nsa_2_hl()
        self.match_agc_2_hl()
        self.match_nsa_2_agc()
        self.count_sample()
        
    def cull_hl(self):
        vbest = self.hl['v']
        vflag = (vbest > vmin) & (vbest < vmax)
        raflag = (self.hl['al2000']*15 > ramin) & (self.hl['al2000']*15 < ramax) 
        decflag = (self.hl['de2000'] < decmax) & (self.hl['de2000']> decmin)
        overlap =  raflag & decflag
        self.hl = self.hl[overlap]

    def cull_nsa(self):
        vbest = self.nsa['Z']*3.e5
        vflag = (vbest > vmin) & (vbest < vmax)
        raflag = (self.nsa['RA'] > ramin) & (self.nsa['RA'] < ramax) 
        decflag = (self.nsa['DEC'] < decmax) & (self.nsa['DEC'] > decmin)
        overlap = vflag & raflag & decflag
        self.nsa = self.nsa[overlap]

    def cull_ned(self):
        vbest = self.ned['Velocity']
        vflag = (vbest > vmin) & (vbest < vmax)
        raflag = (self.ned['RA'] > ramin) & (self.ned['RA'] < ramax) 
        decflag = (self.ned['DEC'] < decmax) & (self.ned['DEC'] > decmin)
        # only keep objects with spectroscopic redshifts
        # https://ned.ipac.caltech.edu/help/faq5.html#5f
        
        speczflag = (self.ned['Redshift Flag'] == 'SPEC') | ((self.ned['Redshift Flag'] == 'N/A') & (self.ned['Redshift Points'] > 2.1))
        print('ned speczflag = ',sum(speczflag))
        #speczflag =  (self.ned['Redshift Flag'] == 'N/A') 
        overlap = vflag & raflag & decflag & speczflag
        self.ned = self.ned[overlap]

    def cull_agc(self):
        # create velocity that is V21 if present, and VOPT otherwise
        flag = self.agc['V21'] > 1.
        vbest = ~flag*self.agc['VOPT'] + flag*self.agc['V21']

        # or create velocit that is vopt if present, and V21 otherwise
        flag = self.agc['OPT'] > 1.
        vbest = flag*self.agc['VOPT'] + ~flag*self.agc['V21']

        # newest version of agc has vhelagc
        # using this as velocity
        vbest = self.agc['vhelagc']
        #avflag1 = (agc['VOPT'] > vmin) & (agc['VOPT'] < vmax)
        #avflag2 = (agc['V21'] > vmin) & (agc['V21'] < vmax)
        avflag = (vbest > vmin) & (vbest < vmax)
        raflag = (self.agc['radeg'] > ramin) & (self.agc['radeg'] < ramax) 
        decflag = (self.agc['decdeg'] < decmax) & (self.agc['decdeg'] > decmin)

        # cut based on iposition
        # keep iposition > 7

        # cut based on description - remove HIonly sources
        overlap = avflag & raflag & decflag
        self.agc = self.agc[overlap]
        self.agc_vbest = vbest[overlap]
        self.agc_ra_key = 'radeg'
        self.agc_dec_key = 'decdeg'
        
    def cull_agc_full(self):
        # create velocity that is V21 if present, and VOPT otherwise
        v21flag = self.agc['v21'] > 1.
        vbest = ~v21flag*self.agc['vopt'] + v21flag*self.agc['v21']

        # newest version of agc has vhelagc
        # using this as velocity
        #vbest = self.agc['vhelagc']
        #avflag1 = (agc['VOPT'] > vmin) & (agc['VOPT'] < vmax)
        #avflag2 = (agc['V21'] > vmin) & (agc['V21'] < vmax)
        avflag = (vbest > vmin) & (vbest < vmax)
        #raflag = (self.agc['RA'] > ramin) & (self.agc['RA'] < ramax) 
        #decflag = (self.agc['DEC'] < decmax) & (self.agc['DEC'] > decmin)
        self.agc_ra_key = 'RA'
        self.agc_dec_key = 'DEC'
        raflag = (self.agc['radeg'] > ramin) & (self.agc['radeg'] < ramax) 
        decflag = (self.agc['decdeg'] < decmax) & (self.agc['decdeg'] > decmin)
        self.agc_ra_key = 'radeg'
        self.agc_dec_key = 'decdeg'

        # cut based on iposition
        # keep iposition > 7
        ipflag = self.agc['iposition'] > 7
        # cut based on description - remove HIonly sources
        description_flag = self.agc['description'] != 'HIonly'
        overlap = avflag & raflag & decflag & ipflag & description_flag
        self.agc = self.agc[overlap]
        self.agc_vbest = vbest[overlap]
        c = Column(self.agc_vbest,name='vhelagc')
        self.agc.add_column(c)
        c1 = Column(self.agc['radeg'],name='RA')
        c2 = Column(self.agc['decdeg'],name='DEC')
        self.agc.add_columns([c1,c2])
    def match_nsa_2_hl(self):
        '''
        HL catalog is the start of the super sample
        
        identify HL galaxies with AGC or NSA w/in 5"
         - track HL, NSA, and AGC name

        '''
        # match NSA to Hlleda
        self.insa, d2d, d3d = self.hcoord.match_to_catalog_sky(self.ncoord)
        # first look at number with match w/in 10"
        self.n2h_matchflag = d2d < self.max_match_offset/3600*u.deg
        print('MATCHING NSA TO HYPERLEDA')
        print('number of matches w/in ',str(self.max_match_offset),' arcsec = ',sum(self.n2h_matchflag),'/',len(self.n2h_matchflag))

        # need to also keep track of which AGC galaxy was matched to HL source
        # and remove these from further searches

        self.nsa_matched2_hl = np.zeros(len(self.ncoord.ra),'bool')
        self.nsa_matched2_hl[self.insa[self.n2h_matchflag]] = np.ones(sum(self.n2h_matchflag),'bool') 

    def match_agc_2_hl(self):
        '''
        HL catalog is the start of the super sample
        
        identify HL galaxies with AGC or NSA w/in 5"
         - track HL, NSA, and AGC name

        '''
        # match NSA to Hlleda
        self.iagc, agcd2d, agcd3d = self.hcoord.match_to_catalog_sky(self.acoord)
        # first look at number with match w/in 10"
        self.a2h_matchflag = agcd2d < self.max_match_offset/3600*u.deg
        print('MATCHING AGC TO HYPERLEDA')
        print('number of matches w/in ',str(self.max_match_offset),' arcsec = ',sum(self.a2h_matchflag),'/',len(self.a2h_matchflag))

        # need to also keep track of which AGC galaxy was matched to HL source
        # and remove these from further searches

        self.agc_matched2_hl = np.zeros(len(self.acoord.ra),'bool')
        self.agc_matched2_hl[self.iagc[self.a2h_matchflag]] = np.ones(sum(self.a2h_matchflag),'bool') 

    def match_nsa_2_agc(self):
        '''
        for NSA and AGC galaxies NOT already matched to HL, match NSA to AGC

        if offset is < 5", add galaxy to catalog

        track NSA and AGC names
        
        '''
        # match NSA to AGC
        self.insa_n2a, d2d_n2a, d3d_n2a = self.acoord.match_to_catalog_sky(self.ncoord)
        # first look at number with match w/in 10"
        self.n2a_matchflag = d2d_n2a < self.max_match_offset/3600*u.deg
        print('MATCHING NSA TO AGC')
        print('number of matches w/in ',str(self.max_match_offset),' arcsec = ',sum(self.n2a_matchflag),'/',len(self.n2a_matchflag))
        print('')
        print('number of AGC sources matched to either HL or AGC =  ',sum(self.n2a_matchflag | self.agc_matched2_hl))
        # keep track of NSA galaxies that are not matched to HL or AGC

        self.nsa_matched2_agc = np.zeros(len(self.ncoord.ra),'bool')
        self.nsa_matched2_agc[self.insa_n2a[self.n2a_matchflag]] = np.ones(sum(self.n2a_matchflag),'bool') 
        
    def count_sample(self):
        # start with HL catalog
        n1 = len(self.hl)

        # add nsa to agc galaxy matches
        # that weren't matched to HL
        n2 = sum(self.n2a_matchflag & ~self.agc_matched2_hl)

        # add agc that weren't matched to either
        n3 = sum(~self.n2a_matchflag & ~self.agc_matched2_hl)

        # add nsa that weren't matched to agc or hl

        n4 = sum(~self.nsa_matched2_hl & ~self.nsa_matched2_agc)
        print(n1,n2,n3,n4)
        print('total number of galaxies in sample = ',n1+n2+n3+n4)
        ntotal = n1+n2+n3+n4

        # BUILD SAMPLE
        
        hlname =  np.zeros(ntotal, dtype=self.hl['objname'].dtype)
        hlra = np.zeros(ntotal,dtype=self.hl['al2000'].dtype)
        hldec = np.zeros(ntotal,dtype=self.hl['de2000'].dtype)
        hvel = np.zeros(ntotal,dtype=self.hvel.dtype)

        aname =  np.zeros(ntotal, dtype=self.agc['AGCnr'].dtype)
        ara = np.zeros(ntotal,dtype=self.agc[self.agc_ra_key].dtype)
        adec = np.zeros(ntotal,dtype=self.agc[self.agc_dec_key].dtype)
        avel = np.zeros(ntotal,dtype=self.avel.dtype)

        nname =  np.zeros(ntotal, dtype=self.nsa['NSAID'].dtype)
        nra = np.zeros(ntotal,dtype=self.nsa['RA'].dtype)
        ndec = np.zeros(ntotal,dtype=self.nsa['DEC'].dtype)
        nvel = np.zeros(ntotal,dtype=self.nvel.dtype)

        hflag = np.zeros(ntotal,'bool')
        aflag = np.zeros(ntotal,'bool')
        nflag = np.zeros(ntotal,'bool')

        
        # first section includes HL with matches to AGC and NSA

        out_columns = [hlname,hlra,hldec,hvel,\
                       aname,ara,adec,avel,\
                       nname,nra,ndec,nvel]
        data_columns = [self.hl['objname'],15.*self.hl['al2000'],self.hl['de2000'],self.hvel,\
                        self.agc['AGCnr'],self.agc[self.agc_ra_key],self.agc[self.agc_dec_key],self.avel,\
                        self.nsa['NSAID'],self.nsa['RA'],self.nsa['DEC'],self.nvel]

        for i in range(4):
            out_columns[i][0:n1] = data_columns[i][0:n1]

        for i in range(4,8):
            out_columns[i][0:n1][self.a2h_matchflag] = data_columns[i][self.iagc[self.a2h_matchflag]]
            
        for i in range(8,12):
            out_columns[i][0:n1][self.n2h_matchflag] = data_columns[i][self.insa[self.n2h_matchflag]]
            
        hflag[0:n1] = np.ones(n1,'bool')
        aflag[0:n1][self.a2h_matchflag] = np.ones(sum(self.a2h_matchflag),'bool')
        nflag[0:n1][self.n2h_matchflag] = np.ones(sum(self.n2h_matchflag),'bool')

        # add in additional NSA matches to AGC
        
        iagc = np.arange(len(self.agc))[(self.n2a_matchflag & ~self.agc_matched2_hl)]
        insa = self.insa_n2a[(self.n2a_matchflag & ~self.agc_matched2_hl)]

        for i in range(4,8):
            out_columns[i][n1:n1+n2] = data_columns[i][iagc]
            
        for i in range(8,12):
            out_columns[i][n1:n1+n2] = data_columns[i][insa]
            
        aflag[n1:n1+n2] = np.ones(len(iagc),'bool')
        nflag[n1:n1+n2] = np.ones(len(insa),'bool')

        # add remainder of AGC
        iagc = np.arange(len(self.agc))[(~self.n2a_matchflag & ~self.agc_matched2_hl)]
        for i in range(4,8):
            out_columns[i][n1+n2:n1+n2+n3] = data_columns[i][iagc]

        aflag[n1+n2:n1+n2+n3] = np.ones(len(iagc),'bool')


        # add remainder of NSA
        insa = np.arange(len(self.nsa))[(~self.nsa_matched2_hl & ~self.nsa_matched2_agc)]
        for i in range(8,12):
            out_columns[i][n1+n2+n3:n1+n2+n3+n4] = data_columns[i][insa]
        nflag[n1+n2+n3:n1+n2+n3+n4] = np.ones(len(insa),'bool')        
          
        out_column_names = ['HL_name','hlra','hldec','hvel','HLflag',\
                            'AGC_name','ara','adec','avel','AGCflag',\
                            'NSA_name','nra','ndec','nvel','NSAflag']

        self.sample_table = Table([hlname,hlra,hldec,hvel,hflag,aname,ara,adec,avel,aflag,nname,nra,ndec,nvel,nflag],\
                                  names = out_column_names)
        self.sample_table.write('kitchen_sink.fits',format='fits',overwrite=True)
        self.check_duplicates_table1()
        
    def check_duplicates_table1(self):
        print('METHOD 1')
        fields = ['HL','AGC','NSA']
        for n in fields:
            print('checking HL ',n,' name')
            duplicates(self.sample_table,n+'_name',flag=self.sample_table[n+'flag'])
        
        
    def get_smart(self,maxoffset=10,veloffset=300.):
        '''
        matching offset in arcsec

        veloffset in km/s
        '''
        # use code I already wrote to match catalogs for A100+SDSS!!!
        self.max_match_offset = maxoffset

        ###############################################
        ## FIRST MATCH AGC AND HYPERLEDA
        ###############################################
        velocity1 = self.hvel
        velocity2 = self.avel

        # don't use velocity matching for now        
        hl_2, hl_matchflag, agc_2, agc_matchflag = make_new_cats(self.hl, self.agc,RAkey1='RAdeg',DECkey1='de2000',RAkey2=self.agc_ra_key,DECkey2=self.agc_dec_key, velocity1=velocity1, velocity2=velocity2, maxveloffset = veloffset,maxoffset=maxoffset)

        # getting data coercion error
        # testing by matching agc and nsa - match
        # still get the same problem - odd

        # now trying without converting fits tables to Table
        # this worked fine!!!
        # so need to convert Hyperleda table to fits table
        # will write this out and read it back in in the __init__ function...
        # hl_2, hl_matchflag, agc_2, agc__matchflag = make_new_cats(self.nsa, self.agc,RAkey1='RA',DECkey1='DEC',RAkey2='radeg',DECkey2='decdeg', velocity1=None, velocity2=None, maxveloffset = voffset,maxoffset=max_match_offset)
        
        ###############################################
        # join HL and AGC into one table
        ###############################################
        joined_table = hstack([hl_2,agc_2])
        
        ###############################################
        # add columns that track if galaxy is in agc and in nsa
        ###############################################
        c1 = Column(hl_matchflag,name='HLflag')
        c2 = Column(agc_matchflag,name='AGCflag')
        ra = np.zeros(len(hl_matchflag),'f')
        dec = np.zeros(len(hl_matchflag),'f')
        ra = hl_matchflag*joined_table['RAdeg'] + ~hl_matchflag*joined_table[self.agc_ra_key]
        dec = hl_matchflag*joined_table['de2000'] + ~hl_matchflag*joined_table[self.agc_dec_key]

        c3 = Column(ra,name='RA-HL-AGC',dtype='f')
        c4 = Column(dec,name='DEC-HL-AGC',dtype='f')

        # adding another column to track velocity
        # uses HL velocity for all objects in HL catalog
        # uses AGC velocity for any objects in AGC but NOT in HL
        vel = hl_matchflag*joined_table['v'] + ~hl_matchflag*joined_table['vhelagc']
        c5 = Column(vel,name='HL-AGC-VEL',dtype='f')
        joined_table.add_columns([c1,c2,c3,c4,c5])

        joined_table.write('temp.fits',format='fits',overwrite=True)
        joined_table = fits.getdata('temp.fits')

        self.table1 = joined_table
        print('METHOD 2: AFTER FIRST MERGE')
        columns=['objname','AGCnr']
        fields = ['HL','AGC']
        for i,n in enumerate(fields):
            print('checking HL ',n,' name')
            
        ###############################################
        ## SECOND MATCH NSA TO  AGC+HYPERLEDA
        ###############################################    
        # now repeat - join NSA to HL+AGC table
        v1 = joined_table['HL-AGC-VEL']
        v2 = self.nsa['Z']*3.e5
        hlagc_2, hlagc_matchflag, nsa_2, nsa_matchflag = make_new_cats(joined_table, self.nsa, RAkey1='RA-HL-AGC',DECkey1='DEC-HL-AGC',RAkey2='RA',DECkey2='DEC', velocity1=v1, velocity2=v2, maxveloffset = veloffset,maxoffset=maxoffset)
        # write out joined a100-sdss-nsa catalog
        joined_table2 = hstack([hlagc_2,nsa_2])
        c1 = Column(nsa_matchflag,name='NSAflag')
        joined_table2.add_column(c1)

        ra = hlagc_matchflag*joined_table2['RA-HL-AGC'] + ~hlagc_matchflag*joined_table2['RA_2']
        dec = hlagc_matchflag*joined_table2['DEC-HL-AGC'] + ~hlagc_matchflag*joined_table2['DEC_2']

        #c3 = Column(ra,name='RA-HL-AGC-NSA',dtype='f')
        #c4 = Column(dec,name='DEC-HL-AGC-NSA',dtype='f')
        c3 = Column(ra,name='RA-COMBINED',dtype='f')
        c4 = Column(dec,name='DEC-COMBINED',dtype='f')

        # adding another column to track velocity
        # uses HL velocity for all objects in HL catalog
        # uses AGC velocity for any objects in AGC but NOT in HL
        vel = hlagc_matchflag*joined_table2['HL-AGC-VEL'] + ~hlagc_matchflag*joined_table2['Z']*3.e5
        c5 = Column(vel,name='HL-AGC-NSA-VEL',dtype='f')
        joined_table2.add_columns([c3,c4,c5])

        ###############################################
        ## FIX BOOLEAN COLUMNS
        ###############################################
        # boolean columns are getting converted weird
        # when I write and then read the fits table

        try:
            joined_table2['AGCflag'] = (joined_table2['AGCflag'] == 84)
            joined_table2['HLflag'] = (joined_table2['HLflag'] == 84)
        except KeyError:
            print('trouble in paradise')

        # trying this to see if it avoids KeyError that I am getting
        #
        # KeyError was just me being stupid
        # this does fix the data coercion error though...
        #
        # need to figure out how to convert the table to fits format
        # without writing it out and reading it back in
        joined_table2.write('smart_kitchen_sink.fits',format='fits',overwrite=True)
        self.table2 = joined_table2
        #self.check_duplicates_t1()
        #joined_table2.write('temp.fits',format='fits',overwrite=True)
        #joined_table2 = fits.getdata('temp.fits')

        #####  SKIPPING NED MATCH FOR NOW ##############

        ## REDO THIS TO MATCH TO GIANLUCA'S CATALOG
        
        ## ###########################################
        ## ## THIRD MATCH NED TO  AGC+HYPERLEDA+NSA
        ## ###############################################
        ## self.test = joined_table2
        ## v1 = joined_table2['HL-AGC-NSA-VEL']
        ## v2 = self.ned['Velocity']

        ## hlagc_3, hlagc_matchflag3, ned_2, ned_matchflag = make_new_cats(joined_table2, self.ned, RAkey1='RA_1',DECkey1='DEC_1',RAkey2='RA',DECkey2='DEC', velocity1=v1, velocity2=v2, maxveloffset = veloffset,maxoffset=maxoffset)
        ## # write out joined a100-sdss-nsa catalog
        ## joined_table3 = hstack([hlagc_3,ned_2])
        ## c1 = Column(ned_matchflag,name='NEDflag')
        ## joined_table3.add_column(c1)

        ## ra = hlagc_matchflag3*joined_table3['RA-HL-AGC-NSA'] + ~hlagc_matchflag3*joined_table3['RA']
        ## dec = hlagc_matchflag3*joined_table3['DEC-HL-AGC-NSA'] + ~hlagc_matchflag3*joined_table3['DEC']

        ## c3 = Column(ra,name='RA-COMBINED',dtype='f')
        ## c4 = Column(dec,name='DEC-COMBINED',dtype='f')

        ## # adding another column to track velocity
        ## # uses HL velocity for all objects in HL catalog
        ## # uses AGC velocity for any objects in AGC but NOT in HL
        ## vel = hlagc_matchflag3*joined_table3['HL-AGC-NSA-VEL'] + ~hlagc_matchflag3*joined_table3['Velocity']
        ## c5 = Column(vel,name='VEL-COMBINED',dtype='f')
        ## joined_table3.add_columns([c3,c4,c5])
        
        ## ###############################################
        ## ## FIX BOOLEAN COLUMNS...AGAIN
        ## ###############################################
        ## # boolean columns are getting converted weird
        ## # when I write and then read the fits table
        ## try:
        ##     joined_table3['AGCflag'] = (joined_table3['AGCflag'] == 84)
        ##     joined_table3['HLflag'] = (joined_table3['HLflag'] == 84)
        ##     joined_table3['NSAflag'] = (joined_table3['NSAflag'] == 84)
        ## except KeyError:
        ##     print('trouble in paradise')
            
        ## joined_table3.write('smart_kitchen_sink.fits',format='fits',overwrite=True)
        ## self.table2 = joined_table3
        ## self.check_duplicates_t2()
    def check_duplicates_t2(self):
        print('METHOD 2')
        columns=['objname','AGCnr','NSAID','Object Name']
        fields = ['HL','AGC','NSA','NED']
        for i,n in enumerate(fields):
            print('checking HL ',n,' name')
            duplicates(self.table2,columns[i],flag=self.table2[n+'flag'])
    def check_vel(self, table2flag=False):
        # compare recession velocities of "matches"
        plt.figure(figsize=(8,6))
        if table2flag:
            mytable = self.table2
            fields=['v','vhelagc','Z']
        else:
            mytable = self.sample_table
            fields=['hvel','avel','nvel']
        plt.plot(mytable[fields[0]],mytable[fields[1]],'bo',label='AGC',alpha=.5)
        if table2flag:
            plt.plot(mytable[fields[0]],mytable[fields[2]]*3.e5,'ro',label='NSA',alpha=.5)
        else:
            plt.plot(mytable[fields[0]],mytable[fields[2]],'ro',label='NSA',alpha=.5)
        xmin,xmax = plt.xlim()
        xl = np.linspace(xmin,xmax,100)
        plt.plot(xl,xl,'k-')
        plt.plot(xl,xl-300,'k--')
        plt.plot(xl,xl+300,'k--')
        plt.xlabel('Hyperleda v_r')
        plt.ylabel('v_r of matched galaxy')
        plt.legend()
        plt.savefig('dv_of_matches.pdf')

class panel_plots:
    def plotimages(self,flag, outfile_string='test',agcflag=False,nsaflag=False,nedflag=False,onlyflag=False,startindex=0):
        nsaindex = self.t.NSAID[flag]
        hra1 = self.hcoord.ra.deg[flag]
        hdec1 = self.hcoord.dec.deg[flag]
        nra2 = self.ncoord.ra.deg[flag]
        ndec2 = self.ncoord.dec.deg[flag]
        ara3 = self.acoord.ra.deg[flag]
        adec3 = self.acoord.dec.deg[flag]
        #nedra = self.nedcoord.ra.deg[flag]
        #neddec = self.nedcoord.dec.deg[flag]
        w21 = self.t.width[flag]
        hlname = self.t.objname[flag]
        nsaid = self.t.NSAID[flag]
        agcnumber = self.t.AGCnr[flag]
        #nedname = self.t['Object Name'][flag]
        galnumber = np.arange(len(self.t.NSAID))[flag]
        plt.figure(figsize=(12,10))
        plt.subplots_adjust(bottom=.05,left=.05,top=.9,right=.95,hspace=.5)
        # plots suddenly stopped working for AGC
        # tryting to see if starting at a different index will help
        if agcflag:
            startindex=10

        i=0 + startindex
        nsubplot = 1
        nrow=4
        ncol=5
        if sum(flag) == 0:
            print('no duplicates to plot')
            return
        while nsubplot < nrow*ncol+1:#for i in range(9):
            plt.subplot(nrow,ncol,nsubplot)
            #print('flag index = ',i)
            #try:
            if agcflag:
                print('agcflag is set',i,nsubplot)
                w = getlegacy(ara3[i],adec3[i],ra2=nra2[i],dec2=ndec2[i],ra3=hra1[i],dec3=hdec1[i],agcflag=agcflag,onlyflag=onlyflag)
            elif nsaflag:
                w = getlegacy(nra2[i], ndec2[i],ra2=hra1[i],dec2=hdec1[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
            elif nedflag:
                w = getlegacy(nedra[i], neddec[i],ra2=hra1[i],dec2=hdec1[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
            else:
                w = getlegacy(hra1[i], hdec1[i],ra2=nra2[i],dec2=ndec2[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
            '''
            except:
                i = i + 1
                print('trouble in paradise',i)
                print('maybe coords are outside Legacy Survey?')
                if agcflag:
                    print(ara3[i],adec3[i])
                elif nsaflag:
                    print(nra2[i],ndec2[i])
                else:
                    print(hra1[i],hdec1[i])
                continue
            '''
            #plt.axis([50,200,50,200])
            #plt.axis([75,175,75,175])
            #plt.title(str(hlname[i])+'\n'+nedname[i]+'\n NSA '+str(nsaid[i])+' / AGC '+str(agcnumber[i]),fontsize=8)
            plt.title(str(hlname[i])+'\n'+'\n NSA '+str(nsaid[i])+' / AGC '+str(agcnumber[i]),fontsize=8)
            if nsubplot == 1:
                plt.text(80, 205,str(outfile_string),fontsize=16,horizontalalignment='left')
            self.add_allgals(w, agcflag=agcflag)
            if w21[i] > .1:
                plt.text(.05, .05,'W21='+str(w21[i]),fontsize=8,c='.7', transform=plt.gca().transAxes)
            plt.text(.05,.85,'Gal '+str(galnumber[i]),fontsize=8,c='.7', transform=plt.gca().transAxes)
            if self.HLflag[i]:
                gname = self.t['objname'][i]
            elif self.AGCflag[i]:
                gname = 'AGC'+self.t['AGCnr'][i]
            elif self.NSAflag[i]:
                gname = 'NSAID'+self.t['NSAID'][i]
            plt.text(.05,.1,'Gal '+str(gname),fontsize=8,c='.7', transform=plt.gca().transAxes)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            i = i + 1
            nsubplot += 1
        plt.savefig('../plots/AGC-HL-NSA-'+outfile_string+'.png')
    def densearray(self,flag, outfile_string='test',agcflag=False,nsaflag=False,nedflag=False,onlyflag=False,startindex=0,endindex=None):
        nsaindex = self.t.NSAID[flag]
        hra1 = self.hcoord.ra.deg[flag]
        hdec1 = self.hcoord.dec.deg[flag]
        nra2 = self.ncoord.ra.deg[flag]
        ndec2 = self.ncoord.dec.deg[flag]
        ara3 = self.acoord.ra.deg[flag]
        adec3 = self.acoord.dec.deg[flag]
        #nedra = self.nedcoord.ra.deg[flag]
        #neddec = self.nedcoord.dec.deg[flag]
        super_ra = self.t['RA-COMBINED']
        super_dec = self.t['DEC-COMBINED']
        
        w21 = self.t.width[flag]
        hlname = self.t.objname[flag]
        nsaid = self.t.NSAID[flag]
        agcnumber = self.t.AGCnr[flag]
        #nedname = self.t['Object Name'][flag]
        galnumber = np.arange(len(hra1))[flag]
        plt.figure(figsize=(12,7))
        plt.subplots_adjust(bottom=.05,left=.05,top=.9,right=.95,hspace=.01,wspace=.01)
        # plots suddenly stopped working for AGC
        # tryting to see if starting at a different index will help
        if agcflag:
            startindex=10

        i=0 + startindex
        nsubplot = 1
        nrow=5
        ncol=10
        galids_in_fov = []
        if sum(flag) == 0:
            print('no duplicates to plot')
            return
        if endindex is not None:
            maxcount = endindex-startindex+1
        else:
            maxcount = nrow*ncol+1
        while nsubplot < maxcount:
            #print(i,nsubplot,maxcount)
            plt.subplot(nrow,ncol,nsubplot)
            #print('flag index = ',i)
            #try:
            massflag=False
            w = getlegacy(super_ra[i], super_dec[i],ra2=nra2[i],dec2=ndec2[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
            jpegflag=True
            if w is None:
                jpegflag=False
                print('trouble in paradise',i)
                print('maybe coords are outside Legacy Survey?')
                print(super_ra[i],super_dec[i])
                # try to get 2MASS J image
                # check to see if 2MASS image exists
                gra = '%.5f'%(super_ra[i]) # accuracy is of order .1"
                gdec = '%.5f'%(super_dec[i])
                galpos = gra+'-'+gdec
                rootname = 'cutouts/2MASS-J-'+str(galpos)
                rootname = 'cutouts/DSS2-'+str(galpos)+'-'+str(image_size)+'-1arcsecpix'     

                fits_name = rootname+'.fits'
                if not(os.path.exists(fits_name)):
                    #print('downloading 2MASS J image ')
                    print('downloading DSS2 Image ')                    
                #
                    c = SkyCoord(ra=super_ra[i]*u.deg,dec=super_dec[i]*u.deg)
                    x = SkyView.get_images(position=c,survey=['DSS2 Red'],pixels=[60,60])
                    # save fits image
                    fits.writeto(fits_name, x[0][0].data, header=x[0][0].header)
                else:
                    print('using 2mass image ',fits_name)
                im, h = fits.getdata(fits_name,header=True)
                w = WCS(h)
                norm = simple_norm(im,stretch='asinh',percent=99.5)
                plt.imshow(im,origin='upper',cmap='gray_r', norm=norm)
                # pixel scale is 1 arcsec
                # therefore, to show a 60x60 arcsec image, want to set boundary to center-30:center+30
                im_nrow,im_ncol=im.shape

                massflag=True
            #plt.axis([50,200,50,200])
            #plt.axis([75,175,75,175])

            #plt.title(str(hlname[i])+'\n'+nedname[i]+'\n NSA '+str(nsaid[i])+' / AGC '+str(agcnumber[i]),fontsize=8)
            #if nsubplot == 1:
            #    plt.text(10, 205,str(outfile_string), dtype='i'),fontsize=16,horizontalalignment='left')
            ids = self.add_allgals(w, agcflag=agcflag,jpegflag=jpegflag)
            galids_in_fov.append(ids)
            if massflag:
                text_color='k'
            else:
                text_color='0.7'
            if w21[i] > .1:
                plt.text(.05, .05,'W21='+str(w21[i]),fontsize=8,c=text_color, transform=plt.gca().transAxes)
            plt.text(.05,.85,'Gal '+str(galnumber[i]),fontsize=8,c=text_color, transform=plt.gca().transAxes)
            # remove ticks for internal images
            #print(nsubplot,np.mod(nsubplot,ncol))
            # adjust ticksize of outer left and bottom images
            if massflag:
                plt.axis([int(im_nrow/2-image_size/2),int(im_nrow/2+image_size/2),int(im_ncol/2-image_size/2),int(im_ncol/2+image_size/2)])
            else:
                plt.xticks(np.arange(0,image_size,20),fontsize=8)
                plt.yticks(np.arange(0,image_size,20),fontsize=8)

            #plt.axis([20,80,20,80])
            if (nsubplot < (nrow-1)*(ncol)):
                plt.xticks([],[])
            if (np.mod(nsubplot,ncol) > 1) | (np.mod(nsubplot,ncol) == 0) :
                #print('no y labels')
                plt.yticks([],[])
            i = i + 1
            nsubplot += 1
        #plt.savefig('../plots/densearray-'+outfile_string+'.png')
        return galids_in_fov

    def one_gal(self,i,dssflag=False):
        plt.figure(figsize=(4,4))
        flag = np.ones_like(self.AGCflag, dtype='bool')
        agcflag=False
        onlyflag=False
        nsaindex = self.t.NSAID[flag]
        hra1 = self.hcoord.ra.deg[flag]
        hdec1 = self.hcoord.dec.deg[flag]
        nra2 = self.ncoord.ra.deg[flag]
        ndec2 = self.ncoord.dec.deg[flag]
        ara3 = self.acoord.ra.deg[flag]
        adec3 = self.acoord.dec.deg[flag]
        #nedra = self.nedcoord.ra.deg[flag]
        #neddec = self.nedcoord.dec.deg[flag]
        super_ra = self.t['RA-COMBINED']
        super_dec = self.t['DEC-COMBINED']
        w21 = self.t.width[flag]
        if dssflag:
            w = None
            jpegflag = False
        else:
            w = getlegacy(super_ra[i], super_dec[i],ra2=nra2[i],dec2=ndec2[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
            jpegflag = True
        if w is None:
            jpegflag = False
            print('trouble in paradise',i)
            print('maybe coords are outside Legacy Survey?')
            print(super_ra[i],super_dec[i])
            # try to get 2MASS J image
            # check to see if 2MASS image exists
            gra = '%.5f'%(super_ra[i]) # accuracy is of order .1"
            gdec = '%.5f'%(super_dec[i])
            galpos = gra+'-'+gdec
            rootname = 'cutouts/2MASS-J-'+str(galpos)
            rootname = 'cutouts/DSS2-'+str(galpos)+'-'+str(image_size)+'-1arcsecpix'     

            fits_name = rootname+'.fits'
            if not(os.path.exists(fits_name)):
                #print('downloading 2MASS J image ')
                print('downloading DSS2 Image ')                    
                #
                c = SkyCoord(ra=super_ra[i]*u.deg,dec=super_dec[i]*u.deg)
                x = SkyView.get_images(position=c,survey=['DSS2 Red'],pixels=[60,60])
                # save fits image
                fits.writeto(fits_name, x[0][0].data, header=x[0][0].header)
            else:
                print('using 2mass image ',fits_name)
            im, h = fits.getdata(fits_name,header=True)
            w = WCS(h)
            norm = simple_norm(im,stretch='asinh',percent=99.5)
            plt.imshow(im,origin='upper',cmap='gray_r', norm=norm)

        ids = self.add_allgals(w, agcflag=agcflag,jpegflag = jpegflag)
        text_color='0.7'
        if w21[i] > .1:
            plt.text(.05, .05,'W21='+str(w21[i]),fontsize=8,c=text_color, transform=plt.gca().transAxes)
            plt.text(.05,.85,'Gal '+str(galnumber[i]),fontsize=8,c=text_color, transform=plt.gca().transAxes)
            # remove ticks for internal images
            #print(nsubplot,np.mod(nsubplot,ncol))
            # adjust ticksize of outer left and bottom images
        return ids

    def add_allgals(self,w,agcflag=False,twomass_flag=False,jpegflag=False):
        
        cats = [self.acoord, self.ncoord, self.hcoord,self.glcoord]
        symbols=['co','b*','r+']
        edgecolors = ['c','w','r']
        symbols=['co','r^','yD','gs']
        edgecolors = ['c','b','r','xkcd:goldenrod', 'g']
        edgecolors = ['c','r','y', 'g']
        #edgecolors = ['c0','c1','c2','c3','c4']
        if agcflag:
            facecolors = ['None','None','None','None','None']
        else:
            facecolors = ['None','None','None','None','None']#['c','b','r','None','g']
        sizes = [14,14,14,16,18]
        text_offsets = [(10,14),(10,7),(10,0),(10,-7),(10,-14)]
        gals_in_fov = []
        #w = WCS('hyper-nsa-test.fits',naxis=2)
        for i,c in enumerate(cats):
            px,py = w.wcs_world2pix(c.ra.deg,c.dec.deg,1)
            galnumber = np.arange(len(c.ra.deg))
            #print('number of galaxies in catalog = ',len(c.ra.deg))
            # only keep objects on image
            if twomass_flag:
                #assume image is 300 pixels(default), pix
                keepflag = (px > 0) & (py > 0) & (px < image_size) & (py < image_size)
            else:
                keepflag = (px > 0) & (py > 0) & (px < image_size) & (py < image_size)
            if jpegflag:
                plt.plot(px[keepflag],image_size - py[keepflag],symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i])
            else:
                plt.plot(px[keepflag],py[keepflag],symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i])
            # label points
            #print('number of galaxies in FOV = ',sum(keepflag))
            gnumbers = galnumber[keepflag]
            j = 0
            if i < (len(edgecolors)-1): # by stopping at 3, I am not labeling Gianluca's catalog id
                gals_in_fov.append(gnumbers.tolist())
                for x,y in zip(px[keepflag],py[keepflag]):
                    if jpegflag:
                        y = image_size - y
                    label = str(gnumbers[j])

                    plt.annotate(label, # this is the text
                             (x,y), # this is the point to label
                             textcoords="offset points", # how to position the text
                             xytext=text_offsets[i], # distance from text to points (x,y)
                             ha='center', # horizontal alignment can be left, right or center
                             color=edgecolors[i],fontsize=8)
                    j += 1
        flattened_gnumbers = [val for sublist in gals_in_fov for val in sublist]
        gals_in_fov = list(sorted(set(flattened_gnumbers)))
        return gals_in_fov

class fulltable(panel_plots):
    def __init__(self,nedflag=False):
        # this file contains all columns from HL, AGC and NSA
        self.t = fits.getdata('smart_kitchen_sink.fits')
        #self.t = fits.getdata('smart_kitchen_sink_05feb2020.fits')
        self.hcoord = SkyCoord(self.t['al2000']*u.hr,self.t['de2000']*u.deg,frame='icrs')
        self.ncoord = SkyCoord(self.t['RA_2']*u.deg,self.t['DEC_2']*u.deg,frame='icrs')
        self.acoord = SkyCoord(self.t['RA_1']*u.deg,self.t['DEC_1']*u.deg, frame='icrs')


        self.hvel = self.t['v'] # mean of optical and radio velocities
        self.nvel = self.t['Z']*3.e5
        self.avel = self.t['vhelagc']

        

        # Boolean columns are not preserved correctly
        # topcat recognizes them ok, but they come in as integers (84=True, 70-something= False)
        # resetting flags here to boolean arrays
        self.AGCflag = self.t.AGCflag == 84
        self.HLflag = self.t.HLflag == 84
        self.NSAflag = self.t.NSAflag == 84

        

        self.hl_only_flag = ~self.AGCflag & self.HLflag & ~self.NSAflag 
        self.agc_only_flag = self.AGCflag & ~self.HLflag & ~self.NSAflag
        self.nsa_only_flag = ~self.AGCflag & ~self.HLflag & self.NSAflag


        if nedflag:
            self.nedcoord = SkyCoord(self.t['RA']*u.deg,self.t['DEC']*u.deg, frame='icrs')
            self.NEDflag = self.t.NEDflag == 84
            self.nedvel = self.t['Velocity']
            self.ned_only_flag = ~self.AGCflag & ~self.HLflag & ~self.NSAflag & self.NEDflag
            self.hl_only_flag = ~self.AGCflag & self.HLflag & ~self.NSAflag & ~self.NEDflag
            self.agc_only_flag = self.AGCflag & ~self.HLflag & ~self.NSAflag & ~self.NEDflag
            self.nsa_only_flag = ~self.AGCflag & ~self.HLflag & self.NSAflag & ~self.NEDflag

        self.nedflag = nedflag
        self.fields = ['HL','AGC','NSA','NED']
        self.fields = ['HL','AGC','NSA','GL']
        #self.temp_fix_for_ned()
        self.gl = fits.getdata(homedir+'/research/VirgoFilaments/Gianluca/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.18Dec2020.fits')
        keepflag = self.gl['v_HL'] > 500.
        self.gl = self.gl[keepflag]
        self.glcoord = SkyCoord(self.gl['RA']*u.deg,self.gl['DEC']*u.deg, frame='icrs')
        self.gvel = self.gl['v_HL']
        self.galids_in_fov = len(self.AGCflag)*[None]
        print('LENGTH GALIDS_IN_FOV = ',len(self.galids_in_fov))
    def summary_statistics(self):
        print('number in combined mastertable = ',len(self.AGCflag))
        if self.nedflag:
            print('number in all cats = ',sum(self.AGCflag & self.HLflag & self.NSAflag & self.NEDflag))
            print('number in all cats except NED = ',sum(self.AGCflag & self.HLflag & self.NSAflag & ~self.NEDflag))
            totalflag = np.array(self.AGCflag,'i') + np.array(self.HLflag,'i') + np.array(self.NSAflag,'i') + np.array(self.NEDflag,'i')
            print('number in 3 cats = ',sum(totalflag > 2))
            print('number in 2 cats = ',sum(totalflag > 1))
            print('number in HL only = ',sum(~self.AGCflag & self.HLflag & ~self.NSAflag & ~self.NEDflag))
            print('number in AGC only = ',sum(self.AGCflag & ~self.HLflag & ~self.NSAflag & ~self.NEDflag))
            print('number in NSA only = ',sum(~self.AGCflag & ~self.HLflag & self.NSAflag & ~self.NEDflag))
            print('number in NED only = ',sum(~self.AGCflag & ~self.HLflag & ~self.NSAflag & self.NEDflag))        
        else:
            print('number in all cats = ',sum(self.AGCflag & self.HLflag & self.NSAflag ))
            totalflag = np.array(self.AGCflag,'i') + np.array(self.HLflag,'i') + np.array(self.NSAflag,'i') 
            print('number in 3 cats = ',sum(totalflag > 2))
            print('number in 2 cats = ',sum(totalflag > 1))
            print('number in HL only = ',sum(~self.AGCflag & self.HLflag & ~self.NSAflag ))
            print('number in AGC only = ',sum(self.AGCflag & ~self.HLflag & ~self.NSAflag))
            print('number in NSA only = ',sum(~self.AGCflag & ~self.HLflag & self.NSAflag))
              
    def temp_fix_for_ned(self):
        # as a first test, I want to plot positions of NED sources on image cutouts
        # I haven't encorporated the NED columns into smart_kitchen_sink.fits
        # so just going to read in catalog separately so that it's available to plot in plotimage
                
        nedfile = homedir+'/github/Virgo/tables/ned-noprolog-10dec2019.txt'
        self.ned = ascii.read(nedfile,delimiter='|')
        self.nedcoord = SkyCoord(self.ned['RA']*u.deg,self.ned['DEC']*u.deg, frame='icrs')
        # updating allgals to include NED
    def download_images(self):
        getlegacyimages(self.t['RA-COMBINED'],self.t['DEC-COMBINED'])
    def all_cats(self):
        # keep galaxies that are in AGC and NOT in (HL and NSA)
        flag = self.AGCflag & self.HLflag & self.NSAflag & self.NEDflag
        self.plotimages(flag,outfile_string='All Catalogs',agcflag=False,onlyflag=False)
        #self.plotimages(flag,outfile_string='AGConly',onlyflag=True)
        plt.savefig('All-catalogs.png')
        pass
    def agc_only(self):
        # keep galaxies that are in AGC and NOT in (HL and NSA)
        if self.nedflag:
            flag = self.AGCflag & ~self.HLflag & ~self.NSAflag & ~self.NEDflag
        else:
            flag = self.AGCflag & ~self.HLflag & ~self.NSAflag 
        self.plotimages(flag,outfile_string='AGConly',agcflag=True,onlyflag=True)
        #self.plotimages(flag,outfile_string='AGConly',onlyflag=True)
        self.agc_only_flag = flag
        plt.savefig('AGConly.png')
        pass
    def hl_only(self):
        if self.nedflag:
            flag = ~self.AGCflag & self.HLflag & ~self.NSAflag & ~self.NEDflag
        else:
            flag = ~self.AGCflag & self.HLflag & ~self.NSAflag 
        self.plotimages(flag,outfile_string='HLonly',onlyflag=True)
        self.hl_only_flag = flag
        plt.savefig('HLonly.png')
        pass
    def nsa_only(self):
        if self.nedflag:
            flag = ~self.AGCflag & ~self.HLflag & self.NSAflag & ~self.NEDflag
        else:
            flag = ~self.AGCflag & ~self.HLflag & self.NSAflag 
        self.plotimages(flag,outfile_string='NSAonly',nsaflag=True,onlyflag=True)
        self.nsa_only_flag = flag
        plt.savefig('NSAonly.png')
        pass
    def ned_only(self):
        flag = ~self.AGCflag & ~self.HLflag & ~self.NSAflag & self.NEDflag
        self.plotimages(flag,outfile_string='NEDonly',nedflag=True,onlyflag=True)
        self.ned_only_flag = flag
        plt.savefig('NEDonly.png')
        pass
    def plot_only(self):
        self.agc_only()
        self.hl_only()
        self.nsa_only()
        if not(self.nedflag):
            self.ned_only()
    def plot_duplicates(self):

        columns=['objname','AGCnr','NSAID','Object Name']
        if self.nedflag:
            fields = ['HL','AGC','NSA','NED']
            flags = [self.HLflag,self.AGCflag,self.NSAflag,self.NEDflag]
        else:
            fields = ['HL','AGC','NSA']
            flags = [self.HLflag,self.AGCflag,self.NSAflag]
        for i,n in enumerate(fields):
            print('checking HL ',n,' name')
            unique, counts = duplicates(self.t,columns[i],flag=flags[i])
            plotids = unique[counts > 1]
            plotflag = np.zeros(len(flags[i]),'bool')
            for id in plotids:
                matchflag = id == self.t[columns[i]]
                plotflag[matchflag] = np.ones(sum(matchflag),'bool')


            if i == 1:
                agcflag=True
            else:
                agcflag=False
            if i == 2:
                nsaflag=True
            else:
                nsaflag=False
            if i == 3:
                nedflag=True
            else:
                nedflag=False
            self.plotimages(plotflag,outfile_string=n+'-duplicates',nsaflag=nsaflag,agcflag=agcflag,nedflag=nedflag,onlyflag=True)
            print(i)
            plt.savefig(n+'-duplicates.png')
    def plot_special_cases(self):
        self.plot_only()
        self.plot_duplicates()
    def plot_all(self,startgal=None):
        plt.close('all')
        flag = np.ones_like(self.AGCflag, dtype='bool')
        #print('LENGTH OF GALIDS IN FOV = ',len(self.galids_in_fov))
        #self.plotimages(flag,outfile_string='All Galaxies',agcflag=False,onlyflag=True)
        ngal = len(self.AGCflag)
        ngalperplot = 50
        nplots = np.floor(ngal/ngalperplot)
        #galids_in_fov = []
        if (ngal/ngalperplot - nplots) > 0:
            nplots += 1
        nplots = int(nplots)
        endindex = None
        if startgal is None:
            allplots = [i for i in range(nplots)]
        else:
            first_plot = int(np.floor(startgal/ngalperplot))
            allplots = [i for i in range(first_plot,nplots)]
        for i in allplots:
        #for i in range(1):
            plt.close('all')
            startindex = i*ngalperplot
            s1 = '%04d'%(startindex)
            n2 = startindex+49
            if n2 > (ngal-1):
                n2 = ngal-1
                endindex=n2
                print('MAKING LAST PLOT')
            s2 = '%04d'%(n2)
            print(s1,s2)

            galids =  self.densearray(flag,outfile_string='All-Galaxies',agcflag=False,onlyflag=True,startindex = startindex, endindex=endindex)

            #self.plotimages(flag,outfile_string='AGConly',onlyflag=True)
            #print('LENGTH OF GALIDS IN FOV = ',len(self.galids_in_fov))
            #print('startindex, n2+1 = ',startindex,n2+1)
            if endindex is None:     
                self.galids_in_fov[startindex:n2+1] = galids
            else:
                self.galids_in_fov[startindex:n2] = galids
            #print('LENGTH OF GALIDS IN FOV = ',len(self.galids_in_fov))
            # include range of galaxy ids in name of pdf file
            plt.savefig('plots/gcutouts-'+s1+'-'+s2+'.pdf')

            self.write_spreadsheet()
        pass
    def write_spreadsheet(self):
        '''
        spreadsheet to use in reviewing each galaxy in the sample.

        * col 1: the galaxy id (row number in catalog)
        * col 2: our classification, where:
        - 1 = galaxy is fine (this will be the default setting)
        - 0 = galaxy should be removed (nothing there in the image)
        - 2 = galaxy should be merged (provide id of parent galaxy in col 3)
        - 3 = can't tell, more detailed followup required
        * col 3: parent id, if galaxy should be merged with another source (default is blank)
        * col 4: list of other galaxy ids that appear within the FOV of the cutout (sometimes the numbers are hard to see)
        * col 6-9: list of galaxy name in HL, NED, NSA, AGC surveys

        using pandas, based on this tutorial

        https://xlsxwriter.readthedocs.io/example_pandas_multiple.html
        '''
        galnumber = np.arange(len(self.AGCflag))
        flag = np.ones_like(self.AGCflag, dtype='i')
        parent = np.full(len(self.AGCflag),'')
        #other_gals = np.full(len(self.AGCflag),'')
        other_gals = self.galids_in_fov

        # with NED
        #t = Table([galnumber,flag,parent,other_gals,\
        #           self.t['objname'],self.t['Object Name'],self.t['NSAID'],self.t['AGCnr'],self.t['RA-COMBINED'],self.t['DEC-COMBINED'],\
        #           self.t['v'],self.t['vhelagc'],self.t['Z']*3.e5,self.t['Velocity'],self.t['Redshift Flag'],self.t['Redshift Points']],
        #           names=['galnumber','class','parent','objid_in_fov','HL','NED','NSAID','AGC','RA','DEC',\
        #                  'vHL','vAGC','vNSA','vNED','NEDzflag','NEDzpoints'])
        t = Table([galnumber,flag,parent,other_gals,\
                   self.t['objname'],self.t['NSAID'],self.t['AGCnr'],self.t['RA-COMBINED'],self.t['DEC-COMBINED'],\
                   self.t['v'],self.t['vhelagc'],self.t['Z']*3.e5],
                   names=['galnumber','class','parent','objid_in_fov','HL','NSAID','AGC','RA','DEC',\
                          'vHL','vAGC','vNSA'])
        #cols = [galnumber,flag,parent,other_gals,self.t['objname'],self.t['Object Name'],self.t['NSAID'],self.t['AGCnr'],self.t['RA-COMBINED'],self.t['DEC-COMBINED']]
        #for a in cols:
        #    print(len(a))
        #t = Table([galnumber,flag,parent,other_gals,\
        #           self.t['objname'],self.t['Object Name'],self.t['NSAID'],self.t['AGCnr'],self.t['RA-COMBINED'],self.t['DEC-COMBINED']])
                   #names=['galnumber','class','parent','objid_in_fov','HL','NED','NSAID','AGC','RA','DEC'])
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter('virgo_check_sample_by_eye.xlsx', engine='xlsxwriter')
        gals_per_sheet=1000
        ngal = len(galnumber)
        nsheets = int(ngal/gals_per_sheet)
        if (ngal/gals_per_sheet - nsheets) > 0:
            nsheets += 1
        for i in range(nsheets):
            start_index = i*1000
            end_index = start_index + 1000
            if end_index > len(self.AGCflag):
                pdt = t[start_index:len(self.AGCflag)].to_pandas()
            else:
                pdt = t[start_index:end_index].to_pandas()

            sheet_name = 'Sheet'+str(i)
            # Write each dataframe to a different worksheet.
            pdt.to_excel(writer, sheet_name=sheet_name)

        # Close the Pandas Excel writer and output the Excel file.
        writer.save()
    def velhist(self):
        plt.figure()
        if self.nedflag:
            arrays = [self.hvel,self.avel,self.nvel,self.nedvel]
            names = ['HL','AGC','NSA','NED']
            flags = [self.HLflag,self.AGCflag,self.NSAflag,self.NEDflag]
        else:
            arrays = [self.hvel,self.avel,self.nvel]
            names = ['HL','AGC','NSA']
            flags = [self.HLflag,self.AGCflag,self.NSAflag]

        colors=['r','c','b','#ff7f0e']#'xkcd:goldenrod']
        mybins = np.linspace(450,3500,100)
        for i in range(len(names)):
            t = plt.hist(arrays[i][flags[i]],label=names[i],histtype='step',bins=mybins,color=colors[i])
        plt.legend(loc='upper right')
        plt.xlabel('Recession Velocity')
        plt.savefig('velhist.png')
    def positions_only(self):
        plt.figure(figsize=(10,8))
        self.nedonly = False
        if self.nedonly:
            ra = [self.hcoord.ra,self.acoord.ra,self.ncoord.ra,self.nedcoord.ra]
            dec = [self.hcoord.dec,self.acoord.dec,self.ncoord.dec,self.nedcoord.dec]
            vel = [self.hvel,self.avel,self.nvel,self.nedvel]
            flags = [self.hl_only_flag,self.agc_only_flag,self.nsa_only_flag,self.ned_only_flag]
            colors=['r','c','b','xkcd:goldenrod']
            symbols=['r+','co','b*','kD']
            edgecolors = ['r','c','w','xkcd:goldenrod']
            facecolors = ['r','c','b','None']
            sizes = [14,10,14,8]
        else:
            ra = [self.glcoord.ra,self.hcoord.ra,self.acoord.ra,self.ncoord.ra]
            dec = [self.glcoord.dec,self.hcoord.dec,self.acoord.dec,self.ncoord.dec]
            vel = [self.gvel,self.hvel,self.avel,self.nvel]
            flags = [np.ones(len(self.glcoord.ra),'bool'),self.hl_only_flag,self.agc_only_flag,self.nsa_only_flag]
            #colors=['r','c','b','xkcd:goldenrod']
            #symbols=['r+','co','b*','kD']
            #edgecolors = ['r','c','w','xkcd:goldenrod']
            #facecolors = ['r','c','b','None']
            symbols=['gs','yD','co','r^']
            facecolors = ['None','None','None','None']
            edgecolors = ['g','y','c','r']
            sizes = [2,4,4,3]
    
        
        for i,r in enumerate(ra):
            #plt.scatter(r,dec[i],c=colors[i],label=self.fields[i]+' only')
            plt.plot(r[flags[i]],dec[i][flags[i]],symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i],label=self.fields[i]+' only')
        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')
        plt.legend()
        plt.savefig('positions-only.png')
    def positions(self):
        plt.figure(figsize=(10,8))
        self.nedonly = False
        ra = [self.glcoord.ra,self.hcoord.ra,self.acoord.ra,self.ncoord.ra]
        dec = [self.glcoord.dec,self.hcoord.dec,self.acoord.dec,self.ncoord.dec]
        vel = [self.gvel,self.hvel,self.avel,self.nvel]
        flags = [np.ones(len(self.glcoord.ra),'bool'),self.HLflag,self.AGCflag,self.NSAflag]
        symbols=['k.','yD','co','ro']
        facecolors = ['None','None','None','None']
        edgecolors = ['k','y','c','r']
        sizes = [1,4,6,7]
        fields=['GL','HL','AGC','NSA']
    
        
        for i in np.array([3,2,1,0],'i'):
            #plt.scatter(r,dec[i],c=colors[i],label=self.fields[i]+' only')
            plt.plot(ra[i][flags[i]],dec[i][flags[i]],symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i],label=fields[i])
        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')
        plt.legend()
        plt.savefig('positions-only.png')
                   
        pass
if __name__ == '__main__':
        
    ## start with HL - match to NSA and AGC
    ## consider all matches with offsets less than 5" to be the same galaxy
    

    #s = sample()
    #s = sample(max_match_offset=7.)
    print('Welcome!')
    print('')
    print('To build catalogs, try: \n\n \t s=sample()\n \t s.get_smart() \n \n OR\n\t s.run_it()')
    print('\n\nTo read table and plot images, try: \n\n \t t=fulltable()\n \t t.agc_only() )')
    ## track NSA and AGC names of matches

    ## add HL objects with closest match > 5"

    ## match remaining AGC and NSA

    ## consider all matches with offsets less than 5" to be the same galaxy

    ## list NSA name, track AGC name of match

    ## add remaining NSA galaxies with no HL and no AGC within 5"

    ## add remaining AGC galaxies with no HL and no NSA within 5"

    #t = fulltable()
    #t.plot_all(startgal=6850)
    #t.plot_all(startgal=8950)
    #t.plot_all()
    





