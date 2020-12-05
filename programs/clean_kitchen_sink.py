#!/usr/bin/env python

'''
GOAL:
* read in kitchen_sink
* read in results virgo_check_sample_by_eye
* edit kitchen sink to
  - remove bad sources
  - merge shredded galaxies

INPUT:
* smart_kitchen_sink_v2.fits
* 

OUTPUT:
* output_catalog, as specified below (e.g. vf_clean_sample.v2.fits)
  - this should be updated for new version releases
* ipac table suitable for matching on irsa
  - specified below (e.g. clean_sample.v2.txt)
'''
import numpy as np
import sys
import os
import time

from astropy.io import fits, ascii
from astropy.table import Table, join, hstack, vstack, Column, MaskedColumn 
from astropy.coordinates import SkyCoord
import astropy.units as u

from astroquery.ned import Ned
from matplotlib import pyplot as plt
from mksupersample import getlegacy, getlegacyimages

homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/APPSS/')
from join_catalogs import make_new_cats, join_cats


###################################################################
#### SET UP ARGPARSE
###################################################################
import argparse
parser = argparse.ArgumentParser(description ='Create a crazy big catalog from HL, AGC, NSA')
parser.add_argument('--version',dest = 'version', default='v1',help='version of tables. default is v1')
parser.add_argument('--evcc',dest = 'evcc', default=False,action='store_true',help='include evcc catalog containing galaxies not in our original table')
        
args = parser.parse_args()
if args.evcc:
    outfile_suffix = '_'+args.version+'_evcc'
else:
    outfile_suffix = '_'+args.version

###################################################################
#### SELECT VERSION OF NSA
###################################################################
## CATALOG VERSION NUMBER
## V1 = USE NEWER VERSION OF NSA (this is what we used for all of visual classifications)
## V2 = USE ORIGINAL VERSION OF NSA
VERSION = 2

### INPUT FILES
if VERSION == 1:
    kitchen_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink.fits'
elif VERSION == 2:
    kitchen_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink_v2.fits'
if args.evcc:
    evcc_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink_v2'+outfile_suffix+'.fits'
byeye_classifications = '/home/rfinn/research/Virgo/google-tables/virgo_check_sample_by_eye_v1.csv'

### OUTPUT FILE
output_catalog = 'vf_clean_sample_'+args.version+'.fits'
output_clean = 'clean_sample_'+args.version+'.fits' # I'm not sure why I have two output files, just a holdover I think
ipac_table = 'clean_sample_'+args.version+'.txt'
## BOUNDARIES OF SURVEY REGION
decmin = -35
decmax = 75 
ramax = 280.
ramin = 100. 
vmax = 3300.
vmin = 500.

image_size = 60 # image size to download from legacy survey, in pixels
default_image_size = image_size

def duplicates(table,column,flag=None):
    if flag is None:
        unique, counts = np.unique(table[column], return_counts=True)
    elif flag is not None:
        unique, counts = np.unique(table[column][flag], return_counts=True)
    print('number of duplicates = ',sum(counts > 1))
    #print('duplicates = ',unique[counts > 1])
    return unique, counts

class cutouts:
    def densearray(self,ra=None,dec=None,names=None,ra2=None,dec2=None,outfile_string='test',agcflag=False,onlyflag=True,startindex=None,endindex=None,indices=None,image_size=60):
        plt.figure(figsize=(12,7))
        plt.subplots_adjust(bottom=.05,left=.05,top=.9,right=.95,hspace=.01,wspace=.01)
        if ra is None:
            print('need an ra and dec')
            return
        nsubplot = 1
        nrow=5
        ncol=10
        if endindex is not None:
            maxcount = endindex-startindex+1
        else:
            maxcount = nrow*ncol+1
        if startindex is not None:
            i = startindex
        else:
            i = 0
        while nsubplot < maxcount:
            jpgflag=True
            #print(i,nsubplot,maxcount)
            plt.subplot(nrow,ncol,nsubplot)
            #print('flag index = ',i)
            #try:
            massflag=False
            # get ra and dec

            w = getlegacy(ra[i], dec[i],jpeg=jpgflag,imsize=image_size)

            if w is None:
                jpgflag=False
                print('trouble in paradise',i)
                print('maybe coords are outside Legacy Survey?')
                print(ra[i],dec[i])
                # try to get 2MASS J image
                # check to see if 2MASS image exists
                gra = '%.5f'%(ra[i]) # accuracy is of order .1"
                gdec = '%.5f'%(dec[i])
                galpos = gra+'-'+gdec
                rootname = 'cutouts/DSS2-'+str(galpos)+'-'+str(image_size)+'-1arcsecpix'     
                
                fits_name = rootname+'.fits'
                if not(os.path.exists(fits_name)):
                    print('downloading DSS2 Image ')                    
                    #
                    c = SkyCoord(ra=ra[i]*u.deg,dec=dec[i]*u.deg)
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

            if massflag:
                text_color='k'
            else:
                text_color='0.7'
            plt.text(.05,.85,'AGC '+str(names[i]),fontsize=8,c=text_color, transform=plt.gca().transAxes)
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

            #print('jpegflag = ',jpgflag)
            gfov = self.addgals(w,ra,dec,ra2,dec2,jpegflag=jpgflag,image_size=image_size)
            if indices is not None:
                print(indices[i],'AGC ',names[i],': ',gfov)
            else:
                print(i,'AGC ',names[i],': ',gfov)
            i = i + 1
            nsubplot += 1
    def addgals(self,w,ra1,dec1,ra2,dec2,jpegflag=True,image_size=60):
        c1 = SkyCoord(ra1*u.deg,dec1*u.deg,frame='icrs')
        c2 = SkyCoord(ra2*u.deg,dec2*u.deg,frame='icrs')
        cats = [c1,c2]
        symbols=['co','b*','r+']
        edgecolors = ['c','w','r']
        symbols=['co','r^','yD','gs']
        edgecolors = ['c','b','r','xkcd:goldenrod', 'g']
        edgecolors = ['c','r','y', 'g']
        facecolors = ['None','None','None','None','None']
        sizes = [14,14,14,16,18]
        # this is for labeling different sources in the FOV from parent catalog
        text_offsets = [(0,0),(5,0),(5,0),(5,0),(5,0)]
        allgals = []
        for i,c in enumerate(cats):
            px,py = w.wcs_world2pix(c.ra.deg,c.dec.deg,1)
            galnumber = np.arange(len(c.ra.deg))
            #print('number of galaxies in catalog = ',len(c.ra.deg))
            # only keep objects on image
            keepflag = (px > 0) & (py > 0) & (px < image_size) & (py < image_size)
            if jpegflag:
                plt.plot(px[keepflag],image_size - py[keepflag],symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i])
            else:
                plt.plot(px[keepflag],py[keepflag],symbols[i],mec=edgecolors[i],mfc=facecolors[i],markersize=sizes[i])
            # label points
            #print('number of galaxies in FOV = ',sum(keepflag))
            gnumbers = galnumber[keepflag]

            if i > 0:
                x = px[keepflag]
                y = py[keepflag]

                for j in range(len(gnumbers)):
                    # added image_size - y on 10/26/2020
                    # figured this out in mksupersample but didn't propagate fix to this program
                    plt.text(x[j]+text_offsets[0][0],image_size-y[j]+text_offsets[0][0],gnumbers[j],color=edgecolors[i],fontsize=8)
                    print('%i, %i, x=%.1f, y=%.1f, gnumb=%i'%(i,j,x[j],y[j],gnumbers[j]))
        #print(allgals)
        allgals = gnumbers.tolist()
        return allgals

    def plot_all(self,ra=None,dec=None,flag=None,names=None,startgal=None,ra2=None,dec2=None,image_size=60):
        plt.close('all')

        #print('LENGTH OF GALIDS IN FOV = ',len(self.galids_in_fov))
        #self.plotimages(flag,outfile_string='All Galaxies',agcflag=False,onlyflag=True)
        
        if ra is not None:
            ra = ra[flag]
            dec = dec[flag]
            names = names[flag]
            indices = np.arange(len(flag))[flag]
        else:
            print('must give an array of ra to plot_all')
        ngal = sum(flag)
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
            if (n2 > (ngal)):
                n2 = ngal
                endindex=n2
                print('MAKING LAST PLOT')
            s2 = '%04d'%(n2)
            print(s1,s2)

            self.densearray(ra=ra,dec=dec,names=names,outfile_string='All-Galaxies',agcflag=False,onlyflag=True,startindex = startindex, endindex=endindex,ra2=ra2,dec2=dec2,indices=indices,image_size=image_size )

            plt.savefig('plots/gcutouts-'+s1+'-'+s2+'-'+str(image_size)+'.pdf')
            plt.savefig('plots/gcutouts-'+s1+'-'+s2+'-'+str(image_size)+'.png')

###  CATALOG CLASS
###  This is the main class
class catalog(cutouts):
    def __init__(self,kitchen_sink,byeye_classifications):
        self.kitchen = fits.getdata(kitchen_sink)

        self.kitchen = Table(self.kitchen)
        self.byeye = ascii.read(byeye_classifications, delimiter=',')
        #self.byeye['parent'] = np.array(self.byeye['parent'],'i')

        if args.evcc:
            self.evcc = Table(fits.getdata(evcc_sink))
            
    def runall(self):
        self.merge_class4()
        #self.remove_agc_only()

        self.fix_bad_HL_name()

        ## add a column that has the galid that corresponds to bye-eye cutouts
        self.add_byeye_galid()
        
        if args.evcc:
            # do the same for the evcc table, but fill with void numbers
            self.add_byeye_galid_evcc()
        ## remove sources that were flagged in bye-eye classifications
        ## creates clean_kitchen table
        self.cut_catalogs_byeye()

        
        # stack with evcc
        if args.evcc:
            self.stack_cleankitchen_evcc() # hereafter, array is clean_kitchen


        self.get_super_radec_vr_name()                
        ## match the catalog to the a100 catalog
        ## joined table in clean_a100
        self.match_a100()
        
        ## remove AGC galaxies with multiple entries in table
        self.remove_a100_duplicates()

        ## check the a100 sources that weren't matched to an existing
        ## entry in the table
        #self.check_new_a100(plotflag=True)

        ## remove/merge bad/offset AGC sources
        self.clean_new_a100()
        

        


        ## this galaxy is not in the A100,
        ## so ok to add this after matching with A100
        ## this is the smaller galaxy that was observed in CO
        print('length of clean_a100 before adding UGC gal = ',len(self.clean_a100))        
        self.fix_8822()
        print('length of clean_a100 after adding UGC gal = ',len(self.clean_a100))
        ## match to CO mastertable
        ## this is the 229 galaxies in the file from Francoise and Pascale
        ## moved this to write_subtables.py
        #self.get_CO_sources()


        
        
        ## sort catalog by declination
        self.sort_by_dec()

        ## add a unique VF id for each galaxy
        self.add_vf_galid()

        ## write the cleaned table
        self.write_clean()

        ## write out ipac catalog
        self.catalog_for_z0MGS()
        
    def fix_bad_HL_name(self):
        # replace HL name for UGC 09348
        # this has a HL name = NGC 5658, but NED says these are not the same galaxy
        #
        # replace NGC5658 with UGC09348
        flag = self.kitchen['objname'] == 'NGC5658'
        self.kitchen['objname'][flag] = 'UGC09348'
        
    def fix_8822(self,cat=None):
        # add blue object, this is one with CO detection
        # the following information is from NED
        if cat is None:
            cat = self.clean_a100
        else:
            cat = cat
        colid = ['superName', 'RA','DEC','vr']
        colvalues  = ['UGC08656 NOTES01', 205.129878,42.993819,2899]

        # add an empty row
        cat.add_row()

        new_index = len(self.clean_a100)-1

        for i in range(len(colid)):
            cat[colid[i]][new_index] = colvalues[i]

        self.write_clean(cat=cat)
        self.clean_a100 = cat
    def sort_by_dec(self):
        # sort catalog by declination
        sort_index = np.argsort(self.clean_a100['DEC'])[::-1]
        self.clean_a100 = self.clean_a100[sort_index]
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
        NSA0flag = self.kitchen['NSA0flag']
        # figure out which survey data we have for 4, and which is in the parent
        HL4 = HLflag[class4]
        NSA4 = NSAflag[class4]
        AGC4 = AGCflag[class4]
        # none of the NSA v0 sources need to merge
        # so I don't need to update this part of the code after
        # adding NSA v0_1_2 as the 4th catalog
        NSA04 = NSA0flag[class4]
        
        # merge any missing survey information that is associated with 4
        # with parent entry
        
        for i in range(len(HL4)):
            #print(galid[i], 'parents = ',parents[i])
            parentflag = self.byeye['galnumber'] == int(float(parents[i]))
            parentid = np.arange(len(class4))[parentflag]            
            #print('len(parentid) = ',len(parentid))
            #print(parents[i],parentid,galid[i])
            if len(parentid) == 0:
                print('error matching parent {} for galaxy {}'.format(parents[i],galid[i]))
                
            elif len(parentid) >= 1:

                #print('error matching parent {} for galaxy {}'.format(parents[i],galid[i]))
                # this will occur
                HLp = (~self.byeye['HL'].mask[parentflag])[0]
                NSAp = (self.byeye['NSAID'][parentflag] != 0)[0]
                AGCp = (self.byeye['AGC'][parentflag] != 0)[0]
                NSA0p = (self.kitchen['NSA0flag'][parentflag] != 0)[0]                
                #print(HLp,HL4[i],parentid,sum(parentflag))
                if not(HLp) and HL4[i]:
                    self.merge_sources(parentid,galid[i],HL=True,cat=self.kitchen)
                    self.kitchen['HLflag'][parentid] = True                    
                if not(NSAp) and NSA4[i]:
                    self.merge_sources(parentid,galid[i],NSA=True,cat=self.kitchen)
                    self.kitchen['NSAflag'][parentid] = True                    
                if not(AGCp) and AGC4[i]:
                    self.merge_sources(parentid,galid[i],AGC=True,cat=self.kitchen)
                    self.kitchen['AGCflag'][parentid] = True
                if not(NSA0p) and NSA04[i]:
                    self.merge_sources(parentid,galid[i],NSA0=True,cat=self.kitchen)
                    self.kitchen['NSAflag'][parentid] = True                    

    def merge_sources(self,parentid,galid,cat=None,HL=False,NSA=False,AGC=False,A100=False,NSA0=False):
        if cat is None:
            cat = self.kitchen
        else:
            cat = cat
            
        colnames = cat.colnames
        if HL:
            survey_columns = np.arange(0,44)
        if AGC:
            survey_columns = np.arange(44,83)
        if NSA:
            survey_columns = np.arange(90,200)
        if A100:
            #survey_columns = np.arange(205,389)
            survey_columns = np.arange(352,536)            
        if NSA0:
            survey_columns = np.arange(204,347)
        for i in survey_columns:
            cat[colnames[i]][parentid] = cat[colnames[i]][galid]

    def get_super_radec_vr_name(self):
        # ra is RA HL for those with HL data, then RA NSA 
        # same for DEC and recession velocity
        # for class 16 objects, RA and DEC are overwritten by by-eye values

        cat = self.clean_kitchen
        self.ra = np.zeros(len(cat),'f8')
        self.dec = np.zeros(len(cat),'f8')
        self.vel = np.zeros(len(cat),'f8')
        self.objectname = np.zeros(len(cat),'|S26')

        # use HL if it's available
        flag1 = cat['HLflag']

        # if HL not available, next choice is NSA v0 if available
        flag2 = ~flag1 & cat['NSA0flag']


        flag3 = ~flag1 & cat['NSAflag'] & ~cat['NSA0flag']
        # what about a100 sources???
        flags = [flag1, flag2, flag3]
        # for hyperleda sources
        self.ra[flag1] =  cat['al2000'][flag1]*15
        self.dec[flag1] = cat['de2000'][flag1]
        self.vel[flag1] = cat['v'][flag1]
        self.objectname[flag1] = cat['objname'][flag1]        

        # for NSA v0 sources
        self.ra[flag2] =  cat['RA_NSA0'][flag2]
        self.dec[flag2] = cat['DEC_NSA0'][flag2]
        self.vel[flag2] = cat['Z_2'][flag2]*3.e5
        matchindices = np.arange(len(flag2))[flag3]
        for i in matchindices:
            self.objectname[i] = 'NSA '+str(cat['NSAID_2'][i])


        # for objects not in HL and NSA v0, use NSA v1 sources
        self.ra[flag3] =  cat['RA_2'][flag3]
        self.dec[flag3] = cat['DEC_2'][flag3]
        self.vel[flag3] = cat['Z'][flag3]*3.e5
        matchindices = np.arange(len(flag3))[flag3]
        for i in matchindices:
            self.objectname[i] = 'NSA '+str(cat['NSAID'][i])
            

        
        # for galaxies with class=16, use the RA and DEC in by-eye file

        class16 = self.cleancat['class'] == 16
        aindex = np.arange(len(self.cleancat))[class16]        
        self.ra[aindex] = self.cleancat['RA'][class16]
        self.dec[aindex] = self.cleancat['DEC'][class16]
        
        # if galaxy is class = 4, use the RA, DEC of parent
        # didn't implement this yet
        
        # cut ra,dec,vel
        c1 = Column(self.ra,'RAtemp',dtype=np.float64)
        c2 = Column(self.dec,'DECtemp',dtype=np.float64)
        c3 = Column(self.vel,'vrtemp')
        c4 = Column(self.objectname,'superName')        
        self.clean_kitchen.add_columns([c1,c2,c3,c4])
        
    def add_byeye_galid(self):
        # append original galaxy id to new kitchen sink file
        galnumber = self.byeye['galnumber']
        c = Column(galnumber,name='galnumber')
        self.kitchen.add_column(c)
    def add_byeye_galid_evcc(self):
        # append original galaxy id to new kitchen sink file
        galnumber = -999*np.ones(len(self.evcc))
        c = Column(galnumber,name='galnumber')
        self.evcc.add_column(c)

    def cut_catalogs_byeye(self):
        self.cutflag = (self.byeye['class'] == 2) | (self.byeye['class'] == 4) | (self.byeye['class'] == 0)

        # removing sources with AGC only
        # will match to A100 instead

        self.agconly = ~self.kitchen['HLflag'] & ~self.kitchen['NSAflag'] & self.kitchen['AGCflag'] & ~self.kitchen['NSA0flag']

        self.cutflag = self.cutflag | self.agconly
        
        self.cleancat = self.byeye[~self.cutflag]
        #remove first column, which is a duplicate with second column
        n = self.cleancat.colnames
        self.cleancat.remove_column(n[0])
        self.clean_kitchen = self.kitchen[~self.cutflag]
    def stack_cleankitchen_evcc(self):
        self.clean_kitchen = vstack([self.clean_kitchen,self.evcc])
    def remove_agc_columns(self):
        pass
    def match_a100(self):
        # read in a100
        #self.a100 = fits.getdata('/home/rfinn/research/Virgo/tables/a100-sdss-wise-virgo.fits')
        self.a100 = fits.getdata('/home/rfinn/research/Virgo/ancil-tables/a100-sdss-wise-virgo.fits')
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

        superNames = np.zeros(len(hlnsa_matchflag),'|S26')
        superNames[hlnsa_matchflag] = joined_table2['superName'][hlnsa_matchflag]
        agcindices = np.arange(len(hlnsa_matchflag))[~hlnsa_matchflag]
        agcnames = []
        for i in agcindices:
            agcnames.append('AGC'+str(joined_table2['AGC'][i]))
        superNames[~hlnsa_matchflag] = agcnames
        
        joined_table2.remove_columns(['RAtemp','DECtemp','vrtemp','superName'])
        c3 = Column(ra,name='RA',dtype=np.float64)
        c4 = Column(dec,name='DEC',dtype=np.float64)
        c5 = Column(vel,name='vr',dtype='f8')
        c6 = Column(superNames,name='superName',dtype='S26')        

        joined_table2.add_columns([c3,c4,c5,c6])

        self.clean_a100 = joined_table2
        # inspect new a100 sources
        self.a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']
        print('number of A100-only before cleaning = ',sum(self.a100flag))
    def remove_a100_duplicates(self):
        # look for any duplicate entries from A100 (probably from WISE matching)
        keepflag = np.ones(len(self.clean_a100['AGC']),'bool')
        agcnames = self.clean_a100['AGC']

        ## not sure that this loop works
        ## doesn't it set both duplicate entries to False???

        # get names with more than one entry
        names, counts = duplicates(self.clean_a100[self.a100flag],'AGC')
        print('AGC-only galaxies with duplicate entries')
        print(names[counts > 1])
        for i,a in enumerate(agcnames):
            if self.a100flag[i]:
                # if this is an a100-only galaxy,
                # and it's already matched to another galaxy,
                # then remove the a100-only galaxy
                matches = np.sum(self.clean_a100['AGC'] == a)
                if matches > 1:
                    print('\tdouble entry for AGC ',self.clean_a100['AGC'][i])
                    keepflag[i] = False
        self.clean_a100 = self.clean_a100[keepflag]


        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']
        print('number of A100-only after removing duplicates = ',sum(a100flag))
        print('A100-only galaxies are: ',self.clean_a100['AGC'][a100flag])
        self.write_clean()
    def write_clean(self,cat=None):
        if cat is None:
            self.clean_a100.write(output_catalog,format='fits',overwrite=True)
        else:
            cat.write(output_catalog,format='fits',overwrite=True)

    def check_new_a100(self,plotflag=False):
                           
        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']
        print('inside check_new_a100: ',sum(a100flag))
        if plotflag:
            self.plot_all(ra=self.clean_a100['RAdeg_Use'],dec=self.clean_a100['DECdeg_Use'],flag=a100flag,names=self.clean_a100['AGC'],ra2=self.clean_a100['RA'],dec2=self.clean_a100['DEC'],image_size=90)
        else:
            pass
        # results from visual inspection
        # the following sources should be merged
        plt.savefig('check_a100_only.pdf')
    def clean_new_a100(self):
        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']        
        keepflag = np.ones(len(self.clean_a100['RA']),'bool')

        # decided to delete AGC 258535 after visual inspection - didn't match with SGA, nothing visible in legacy images, deleting from v1
        
        # AGC 223361 is crosslisted on NED with NSA 170605, PGC 039340, NED name is SDSS J121639.61+085004.4

        # AGC 208736
        
        # for 8949, 9151 pair, AGC 208736, use AGC coordinates
        # delete rows corresponding to children
        
        #child = np.array([9206, 9207, 9209, 9213],'i')
        #parent = np.array([7423, 6575, 6638, 8949],'i')
        child = np.array([7423, 9207, 9209, 8952],'i')
        parent = np.array([9206, 6575, 6638, 9213],'i')

        # for version 1 tables
        child = np.array([7331, 9032, 9034, 8851],'i')
        parent = np.array([9031, 6483, 6546, 9038],'i')

        # for version 1 tables, after adding evcc galaxies
        child = np.array([7331, 9146, 9148, 8851],'i')
        parent = np.array([9145, 6483, 6546, 9152],'i')
        
        # for version 1 tables, after adding evcc galaxies and removing two duplicates
        child = np.array([7330, 9144, 9146, 8850],'i')
        parent = np.array([9143, 6482, 6545, 9150],'i')

        for i in range(len(child)):
            print('merging {} with {}'.format(child[i],parent[i]))
            self.merge_sources(parent[i],child[i],cat=self.clean_a100,HL=False,NSA=False,AGC=False,A100=True)
            self.clean_a100['A100flag'][parent[i]] = True
        keepflag[child] = np.zeros(len(child),'bool')
        # remove AGC 258535        
        delflag = self.clean_a100['AGC'] == 258535
        keepflag[delflag] = False

        self.clean_a100 = self.clean_a100[keepflag]
    def read_ned(self):
        nedfile = homedir+'/research/Virgo/supersample/ned-noprolog-25mar2020.txt'
        self.ned = ascii.read(nedfile,delimiter='|')

        # having issues with ned, maybe because of masked array?
        # going to write out fits and read it back in
        #
        self.ned.write(homedir+'/research/Virgo/tables/ned-noprolog-10dec2019.fits',format='fits',overwrite=True)
        self.ned = Table(fits.getdata(homedir+'/research/Virgo/tables/ned-noprolog-10dec2019.fits',format='fits',overwrite=True))
    def cull_ned(self):
        vbest = self.ned['Velocity']
        vflag = (vbest > vmin) & (vbest < vmax)
        raflag = (self.ned['RA'] > ramin) & (self.ned['RA'] < ramax) 
        decflag = (self.ned['DEC'] < decmax) & (self.ned['DEC'] > decmin)
        # only keep objects with spectroscopic redshifts
        # https://ned.ipac.caltech.edu/help/faq5.html#5f

        # skipping this step for now b/c just using catalog to find NEDname
        # adding this back in because I am having trouble with some matches,
        # e.g. getting an SDSS name instead of a NGC name.  maybe this will remove
        # the yucky sources
        speczflag = (self.ned['Redshift Flag'] == 'SPEC') | ((self.ned['Redshift Flag'] == 'N/A') & (self.ned['Redshift Points'] > 2.1))
        print('ned speczflag = ',sum(speczflag))
        #speczflag =  (self.ned['Redshift Flag'] == 'N/A') 
        overlap = vflag & raflag & decflag & speczflag
        self.ned = self.ned[overlap]

    def add_vf_galid(self, cat=None):
        if cat is None:
            cat = self.clean_a100
        else:
            cat = cat
        # assign ids by descending declination
        
        vfid = ['VFID{0:04d}'.format(i) for i in np.arange(len(cat))]
        c = Column(vfid,name='VFID')
        cat.add_column(c)
        if cat is None:
            self.clean_a100 = cat
        else:
            return cat
    def read_clean_a100(self):
        self.clean_a100 = Table(fits.getdata(output_catalog))
    
    def write_clean_cat(self):
        self.cleancat.write(output_clean,format='fits',overwrite=True)
        self.clean_a100.write(output_catalog,format='fits',overwrite=True)

    def catalog_for_z0MGS(self,cat=None):
        '''
        need RA, DEC, and search radius
        in IPAC format table
        for matching with leroy+2019 Galaxy Synthesis WISE+GALEX
        table is served by IRSA
        '''
        if cat is None:
            cat = self.clean_a100
        else:
            cat = cat
        search_radius = 10.*np.ones(len(cat)) # units are arcsec
        newtable = Table([cat['galnumber'],cat['RA'],cat['DEC'],search_radius],names=['galid','ra','dec','major'])
        newtable.write(ipac_table,format='ipac',overwrite=True)
        
if __name__ == '__main__':
    c = catalog(kitchen_sink,byeye_classifications)
    #c.merge_class4()
    #c.cut_catalog()
    #c.catalog_for_z0MGS()
    c.runall()
