#!/usr/bin/env python

'''
GOAL:
- create a super list of all galaxies in

USAGE:
- here you go

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

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

from urllib.parse import urlencode
from urllib.request import urlretrieve

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
decmax = 75 
ramax = 280.
ramin = 100. 
vmax = 3300.
vmin = 500.

def duplicates(table,column,flag=None):
    if flag is None:
        unique, counts = np.unique(table[column], return_counts=True)
    elif flag is not None:
        unique, counts = np.unique(table[column][flag], return_counts=True)
    print('number of duplicates = ',sum(counts > 1))
    #print('duplicates = ',unique[counts > 1])
    return unique, counts

def getlegacy(ra1,dec1,ra2=None,dec2=None, ra3=None,dec3=None,agcflag=False,onlyflag=False):
    #url='http://legacysurvey.org/viewer/fits-cutout?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=mzls+bass-dr6&pixscale=0.27&bands=r'
    url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&pixscale=1.00'
    #http://legacysurvey.org/viewer/cutout.fits?ra=156.2778&dec=28.0920&layer=dr8&pixscale=1.00
    fname = 'hyper-nsa-test.fits'
    urlretrieve(url, fname)
    t,h = fits.getdata(fname,header=True)
    # write out r-band image
    # nevermind - John M figured out how to use MEF with WCS
    #fits.writeto('r-test.fits',t[1],header=h,overwrite=True)
    norm = simple_norm(t[1],stretch='asinh',percent=99.5)
    plt.imshow(t[1],origin='upper',cmap='gray_r', norm=norm)
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
        w = WCS('hyper-nsa-test.fits',naxis=2)
        px,py = w.wcs_world2pix(ra,dec,1)
        #print(px,py)
        r1 = Rectangle((px[0]-dx/2, py[0]-dx/2), dx, dx, edgecolor=colors[0], facecolor='none')
        dx=17.5
        r2 = Rectangle((px[1]-dx/2, py[1]-dx/2), dx, dx, edgecolor=colors[1], facecolor='none')
        dx=15
        r3 = Rectangle((px[2]-dx/2, py[2]-dx/2), dx, dx, edgecolor=colors[2], facecolor='none')
        plt.gca().add_patch(r1)
        plt.gca().add_patch(r2)
        plt.gca().add_patch(r3)
        return w


class sample:
    def __init__(self, max_match_offset=5.):
        
        ## read in my HL catalog
        '''
        Hyperleda query:  http://leda.univ-lyon1.fr/fullsql.html

        parameters described here: http://leda.univ-lyon1.fr/leda/meandata.html

        SQL QUERY:
        
        select
        objname,objtype,de2000,al2000,v,e_v,vopt,e_vopt,vrad,e_vrad,bt,e_bt,type,bar,ring,multiple,compactness,t,e_t,logd25,e_logd25,logr25,e_logr25,pa,incl,logdc,btc,itc,ubtc,bvtc,m21c,hic,mabs,agnclass,kt,e_kt,it,e_it,ut,vt,mfir,e_ut,e_vt, modz, e_modz, mod0, e_mod0,vmaxg, e_vmaxg, vmaxs,e_vmaxs,vdis,e_vdis
        
        where

        de2000 > -35 and de2000 < 75 and  al2000 < 280./360.*24.  and al2000 > 100./360.*24. and v < 3300 and v > 500 and objtype='G'

        - output as csv


        
        not using Gialuca's catalog b/c I'm not sure if there were other cuts made already
        gl = fits.getdata('/Users/rfinn/research/VirgoFilaments/Gianluca/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')
    
        '''
        self.max_match_offset = max_match_offset

        ################################################################
        ## READ IN HYPERLEDA CATALO
        ################################################################
        hlfile = homedir+'/github/Virgo/tables/hyperleda-finn-09dec2019-full.csv'
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

        nedfile = homedir+'/github/Virgo/tables/ned-noprolog-10dec2019.txt'
        self.ned = ascii.read(nedfile,delimiter='|')

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
        agcfile = '/Users/rfinn/research/AGC/agcm1.sh191118.fits'
        self.agc = fits.getdata(agcfile)
        self.cull_agc()
        #self.agc = Table(self.agc)
        
        # flags to track nsa matches to HL and AGC
        self.hl_2_nsa_matchflag = np.zeros(len(self.hl['al2000']),'bool')
        self.hl_2_agc_matchflag = np.zeros(len(self.hl['de2000']),'bool')

        # flags to track nsa matches to HL and AGC
        self.nsa_2_hl_matchflag = np.zeros(len(self.nsa['RA']),'bool')
        self.nsa_2_agc_matchflag = np.zeros(len(self.nsa['RA']),'bool')

        # flags to track AGC matches to HL and NSA
        self.agc_2_hl_matchflag = np.zeros(len(self.agc['radeg']),'bool')
        self.agc_2_nsa_matchflag = np.zeros(len(self.agc['radeg']),'bool')

        self.hcoord = SkyCoord(self.hl['al2000']*u.hr,self.hl['de2000']*u.deg,frame='icrs')
        self.ncoord = SkyCoord(self.nsa['RA']*u.deg,self.nsa['DEC']*u.deg,frame='icrs')
        self.acoord = SkyCoord(self.agc['radeg']*u.deg,self.agc['decdeg']*u.deg, frame='icrs')
        self.nedcoord = SkyCoord(self.ned['RA']*u.deg,self.ned['DEC']*u.deg, frame='icrs')

        self.hvel = self.hl['v'] # mean of optical and radio velocities
        self.nvel = self.nsa['Z']*3.e5
        self.avel = self.agc_vbest
        self.nedvel = self.ned['Velocity']
        
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

    def cull_agc(self):
        # create velocity that is VOPT if present, and V21 otherwise
        voptflag = self.agc['VOPT'] > 1.
        vbest = voptflag*self.agc['VOPT'] + ~voptflag*self.agc['V21']
        vbest = self.agc['vhelagc']
        #avflag1 = (agc['VOPT'] > vmin) & (agc['VOPT'] < vmax)
        #avflag2 = (agc['V21'] > vmin) & (agc['V21'] < vmax)
        avflag = (vbest > vmin) & (vbest < vmax)
        raflag = (self.agc['radeg'] > ramin) & (self.agc['radeg'] < ramax) 
        decflag = (self.agc['decdeg'] < decmax) & (self.agc['decdeg'] > decmin)
        overlap = avflag & raflag & decflag
        self.agc = self.agc[overlap]
        self.agc_vbest = vbest[overlap]
        
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
        ara = np.zeros(ntotal,dtype=self.agc['radeg'].dtype)
        adec = np.zeros(ntotal,dtype=self.agc['decdeg'].dtype)
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
                        self.agc['AGCnr'],self.agc['radeg'],self.agc['decdeg'],self.avel,\
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
        
        
    def get_smart(self,maxoffset=10.):
        self.max_match_offset = maxoffset
        # and use code I already wrote to match catalogs!!!

        ## FIRST MATCH AGC AND HYPERLEDA

        velocity1 = self.avel
        velocity2 = self.nvel
        voffset = 300.
        voffset = None
        
        hl_2, hl_matchflag, agc_2, agc_matchflag = make_new_cats(self.hl, self.agc,RAkey1='RAdeg',DECkey1='de2000',RAkey2='radeg',DECkey2='decdeg', velocity1=None, velocity2=None, maxveloffset = voffset,maxoffset=self.max_match_offset)

        # getting data coercion error
        # testing by matching agc and nsa - match
        # still get the same problem - odd

        # now trying without converting fits tables to Table
        # this worked fine!!!
        # so need to convert Hyperleda table to fits table
        # will write this out and read it back in in the __init__ function...
        # hl_2, hl_matchflag, agc_2, agc__matchflag = make_new_cats(self.nsa, self.agc,RAkey1='RA',DECkey1='DEC',RAkey2='radeg',DECkey2='decdeg', velocity1=None, velocity2=None, maxveloffset = voffset,maxoffset=max_match_offset)
        
        # join HL and AGC into one table
        joined_table = hstack([hl_2,agc_2])
        # add columns that track if galaxy is in agc and in nsa
        c1 = Column(hl_matchflag,name='HLflag')
        c2 = Column(agc_matchflag,name='AGCflag')
        ra = np.zeros(len(hl_matchflag),'f')
        dec = np.zeros(len(hl_matchflag),'f')
        ra = hl_matchflag*joined_table['RAdeg'] + ~hl_matchflag*joined_table['radeg']
        dec = hl_matchflag*joined_table['de2000'] + ~hl_matchflag*joined_table['decdeg']
        c3 = Column(ra,name='RA',dtype='f')
        c4 = Column(dec,name='DEC',dtype='f')
        joined_table.add_columns([c1,c2,c3,c4])

        joined_table.write('temp.fits',format='fits',overwrite=True)
        joined_table = fits.getdata('temp.fits')

        self.table1 = joined_table
        print('METHOD 2: AFTER FIRST MERGE')
        columns=['objname','AGCnr']
        fields = ['HL','AGC']
        for i,n in enumerate(fields):
            print('checking HL ',n,' name')
        ## SECOND MATCH NSA TO  AGC+HYPERLEDA        
        # now repeat - join NSA to HL+AGC table
        hlagc_2, hlagc_matchflag, nsa_2, nsa_matchflag = make_new_cats(joined_table, self.nsa, RAkey1='RA',DECkey1='DEC',RAkey2='RA',DECkey2='DEC', velocity1=None, velocity2=None, maxveloffset = voffset,maxoffset=self.max_match_offset)
        # write out joined a100-sdss-nsa catalog
        joined_table2 = hstack([hlagc_2,nsa_2])
        c1 = Column(nsa_matchflag,name='NSAflag')
        joined_table2.add_column(c1)
        self.table2 = joined_table2
        
        # boolean columns are getting converted weird
        try:
            joined_table2['AGCflag'] = (joined_table2['AGCflag'] == 84)
            joined_table2['HLflag'] = (joined_table2['HLflag'] == 84)
        except KeyError:
            print('trouble in paradise')

        joined_table2.write('smart_kitchen_sink.fits',format='fits',overwrite=True)


        self.check_duplicates_t2()
    def check_duplicates_t2(self):
        print('METHOD 2')
        columns=['objname','AGCnr','NSAID']
        fields = ['HL','AGC','NSA']
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
    def plotimages(self,flag, outfile_string='test',agcflag=False,nsaflag=False,onlyflag=False,startindex=0):
        
        nsaindex = self.t.NSAID[flag]
        hra1 = self.hcoord.ra.deg[flag]
        hdec1 = self.hcoord.dec.deg[flag]
        nra2 = self.ncoord.ra.deg[flag]
        ndec2 = self.ncoord.dec.deg[flag]
        ara3 = self.acoord.ra.deg[flag]
        adec3 = self.acoord.dec.deg[flag]
        w21 = self.t.width[flag]
        hlname = self.t.objname[flag]
        nsaid = self.t.NSAID[flag]
        agcnumber = self.t.AGCnr[flag]
        
        plt.figure(figsize=(12,10))
        plt.subplots_adjust(bottom=.05,left=.05,top=.9,right=.95,hspace=.35)
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
            print('flag index = ',i)
            try:
                if agcflag:
                    print('agcflag is set',i,nsubplot)
                    w = getlegacy(ara3[i],adec3[i],ra2=nra2[i],dec2=ndec2[i],ra3=hra1[i],dec3=hdec1[i],agcflag=agcflag,onlyflag=onlyflag)
                elif nsaflag:
                    w = getlegacy(nra2[i], ndec2[i],ra2=hra1[i],dec2=hdec1[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
                else:
                    w = getlegacy(hra1[i], hdec1[i],ra2=nra2[i],dec2=ndec2[i],ra3=ara3[i],dec3=adec3[i],agcflag=agcflag,onlyflag=onlyflag)
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
            #plt.axis([50,200,50,200])
            plt.axis([75,175,75,175])
            plt.title(str(hlname[i])+'\n NSA '+str(nsaid[i])+' / AGC '+str(agcnumber[i]))
            if nsubplot == 1:
                plt.text(80, 205,str(outfile_string),fontsize=16,horizontalalignment='left')
            self.add_allgals(w, agcflag=agcflag)
            if w21[i] > .1:
                plt.text(80, 80,'W21='+str(w21[i]),fontsize=10)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            i = i + 1
            nsubplot += 1
        plt.savefig('../plots/AGC-HL-NSA-'+outfile_string+'.png')
    def add_allgals(self,w,agcflag=False):
        cats = [self.acoord, self.ncoord, self.hcoord,self.nedcoord]
        symbols=['co','b*','r+']
        edgecolors = ['c','w','r']
        symbols=['co','b*','r+','kD']
        edgecolors = ['c','w','r','xkcd:goldenrod']
        if agcflag:
            facecolors = ['None','b','r','None']
        else:
            facecolors = ['c','b','r','None']
        sizes = [10,14,14,20]
    
        w = WCS('hyper-nsa-test.fits',naxis=2)
        for i,c in enumerate(cats):
            px,py = w.wcs_world2pix(c.ra.deg,c.dec.deg,1)
            plt.plot(px,py,symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i])
        
class fulltable(panel_plots):
    def __init__(self):
        # this file contains all columns from HL, AGC and NSA
        self.t = fits.getdata('smart_kitchen_sink.fits')
        self.hcoord = SkyCoord(self.t['al2000']*u.hr,self.t['de2000']*u.deg,frame='icrs')
        self.ncoord = SkyCoord(self.t['RA_2']*u.deg,self.t['DEC_2']*u.deg,frame='icrs')
        self.acoord = SkyCoord(self.t['radeg']*u.deg,self.t['decdeg']*u.deg, frame='icrs')

        self.hvel = self.t['v'] # mean of optical and radio velocities
        self.nvel = self.t['Z']*3.e5
        self.avel = self.t['vhelagc']

        # Boolean columns are not preserved correctly
        # topcat recognizes them ok, but they come in as integers (84=True, 70-something= False)
        # resetting flags here to boolean arrays
        self.AGCflag = self.t.AGCflag == 84
        self.HLflag = self.t.HLflag == 84
        self.NSAflag = self.t.NSAflag == 84
        self.temp_fix_for_ned()
    def temp_fix_for_ned(self):
        # as a first test, I want to plot positions of NED sources on image cutouts
        # I haven't encorporated the NED columns into smart_kitchen_sink.fits
        # so just going to read in catalog separately so that it's available to plot in plotimage
                
        nedfile = homedir+'/github/Virgo/tables/ned-noprolog-10dec2019.txt'
        self.ned = ascii.read(nedfile,delimiter='|')
        self.nedcoord = SkyCoord(self.ned['RA']*u.deg,self.ned['DEC']*u.deg, frame='icrs')
        # updating allgals to include NED
    def agc_only(self):
        # keep galaxies that are in AGC and NOT in (HL and NSA)
        flag = self.AGCflag & ~self.HLflag & ~self.NSAflag
        self.plotimages(flag,outfile_string='AGConly',agcflag=True,onlyflag=True)
        #self.plotimages(flag,outfile_string='AGConly',onlyflag=True)
        self.agc_only_flag = flag
        plt.savefig('AGConly.png')
        pass
    def hl_only(self):
        flag = ~self.AGCflag & self.HLflag & ~self.NSAflag
        self.plotimages(flag,outfile_string='HLonly',onlyflag=True)
        self.hl_only_flag = flag
        plt.savefig('HLonly.png')
        pass
    def nsa_only(self):
        flag = ~self.AGCflag & ~self.HLflag & self.NSAflag
        self.plotimages(flag,outfile_string='NSAonly',nsaflag=True,onlyflag=True)
        self.nsa_only_flag = flag
        plt.savefig('NSAonly.png')
        pass
    def plot_only(self):
        self.agc_only()
        self.hl_only()
        self.nsa_only()
    def plot_duplicates(self):

        columns=['objname','AGCnr','NSAID']
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
            self.plotimages(plotflag,outfile_string=n+'-duplicates',nsaflag=nsaflag,agcflag=agcflag,onlyflag=True)
            print(i)
            plt.savefig(n+'-duplicates.png')
    def plotall(self):
        self.plot_only()
        self.plot_duplicates()
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



    





