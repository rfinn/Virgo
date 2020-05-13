#!/usr/bin/env python

'''
GOAL:
* read in kitchen_sink
* read in results virgo_check_sample_by_eye
* edit kitchen sink to
  - remove bad sources
  - merge shredded galaxies

'''
import numpy as np
import sys
import os
import time

from astropy.io import fits, ascii
from astropy.table import Table, join, hstack, Column, MaskedColumn 
from astropy.coordinates import SkyCoord
import astropy.units as u

from astroquery.ned import Ned

from matplotlib import pyplot as plt

from mksupersample import getlegacy, getlegacyimages

homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/APPSS/')
from join_catalogs import make_new_cats, join_cats

## CATALOG VERSION NUMBER
## V1 = USE NEWER VERSION OF NSA (this is what we used for all of visual classifications)
## V2 = USE ORIGINAL VERSION OF NSA
VERSION = 2

### INPUT FILES
if VERSION == 1:
    kitchen_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink.fits'
elif VERSION == 2:
    kitchen_sink = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink_v2.fits'    
byeye_classifications = '/home/rfinn/research/Virgo/supersample/virgo_check_sample_by_eye.csv'

## BOUNDARIES OF SURVEY REGION
decmin = -35
decmax = 75 
ramax = 280.
ramin = 100. 
vmax = 3300.
vmin = 500.

image_size = 60 # image size to download from legacy survey, in pixels
default_image_size = image_size


class cutouts:
    def densearray(self,ra=None,dec=None,names=None,ra2=None,dec2=None,outfile_string='test',agcflag=False,onlyflag=True,startindex=None,endindex=None):
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
            gfov = self.addgals(w,ra,dec,ra2,dec2,jpegflag=jpgflag)
            print('AGC ',names[i],': ',gfov)
            i = i + 1
            nsubplot += 1
    def addgals(self,w,ra1,dec1,ra2,dec2,jpegflag=True):
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
        text_offsets = [(10,14),(10,7),(10,0),(10,-7),(10,-14)]
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
            
        #print(allgals)
        allgals = gnumbers.tolist()
        return allgals

    def plot_all(self,ra=None,dec=None,flag=None,names=None,startgal=None,ra2=None,dec2=None):
        plt.close('all')

        #print('LENGTH OF GALIDS IN FOV = ',len(self.galids_in_fov))
        #self.plotimages(flag,outfile_string='All Galaxies',agcflag=False,onlyflag=True)
        if ra is not None:
            ra = ra[flag]
            dec = dec[flag]
            names = names[flag]            
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
            if n2 > (ngal-1):
                n2 = ngal-1
                endindex=n2
                print('MAKING LAST PLOT')
            s2 = '%04d'%(n2)
            print(s1,s2)

            self.densearray(ra=ra,dec=dec,names=names,outfile_string='All-Galaxies',agcflag=False,onlyflag=True,startindex = startindex, endindex=endindex,ra2=ra2,dec2=dec2)

            plt.savefig('plots/gcutouts-'+s1+'-'+s2+'.pdf')
            plt.savefig('plots/gcutouts-'+s1+'-'+s2+'.png')

###  CATALOG CLASS
###  This is the main class
class catalog(cutouts):
    def __init__(self,kitchen_sink,byeye_classifications):
        self.kitchen = fits.getdata(kitchen_sink)

        self.kitchen = Table(self.kitchen)
        self.byeye = ascii.read(byeye_classifications, delimiter=',')
        #self.byeye['parent'] = np.array(self.byeye['parent'],'i')
    def runall(self):
        self.merge_class4()
        #self.remove_agc_only()
        self.get_radec_first()
        self.cut_catalog()
        # match_a100 writes RA, DEC, and VEL to table
        self.match_a100()
        self.check_new_a100()
        #self.get_NEDname()

        # skip for now until other parts are working
        # if file with NED names already exists from a previous query,
        # read the file
        if os.path.exists('ned_names.fits'):
            print('found file ned_names.fits.\nUsing this instead of querying NED')
            nednames = fits.getdata('ned_names.fits')
            c = Column(nednames['NEDname'],name='NEDname')
            
            self.clean_a100.add_column(c)
            self.write_clean()
        else:
            #sys.exit() # for testing purposes
            self.get_NEDname_query()

        # this galaxy is not in the A100,
        # so ok to add this after matching with A100
        self.fix_8822()
        self.sort_by_dec()
        #self.write_clean()        
        self.add_galid()
        self.write_clean()        
        #self.read_clean_a100()
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

    def get_radec_first(self):
        # ra is RA HL for those with HL data, then RA NSA 
        # same for DEC and recession velocity

        self.ra = np.zeros(len(self.kitchen),'f')
        self.dec = np.zeros(len(self.kitchen),'f')
        self.vel = np.zeros(len(self.kitchen),'f')
        flag1 = self.kitchen['HLflag']
        flag2 = ~flag1 & self.kitchen['NSAflag']
        
        flag3 = ~flag1 & ~self.kitchen['NSAflag'] & self.kitchen['NSA0flag']
        # what about a100 sources???
        flags = [flag1, flag2, flag3]
        # for hyperleda sources
        self.ra[flag1] = self.kitchen['al2000'][flag1]*15
        self.dec[flag1] = self.kitchen['de2000'][flag1]
        self.vel[flag1] = self.kitchen['v'][flag1]

        # for NSA v2 sources
        self.ra[flag2] = self.kitchen['RA_2'][flag2]
        self.dec[flag2] = self.kitchen['DEC_2'][flag2]
        self.vel[flag2] = self.kitchen['Z'][flag2]*3.e5

        self.ra[flag3] = self.kitchen['RA_NSA0'][flag3]
        self.dec[flag3] = self.kitchen['DEC_NSA0'][flag3]
        self.vel[flag3] = self.kitchen['Z_2'][flag3]*3.e5
        
        # for galaxies with class=16, use the RA and DEC in by-eye file
        class16 = self.byeye['class'] == 16
        self.ra[class16] = self.byeye['RA'][class16]
        self.dec[class16] = self.byeye['DEC'][class16]
        
        # if galaxy is class = 4, use the RA, DEC of parent
        # didn't implement this yet

    def cut_catalog(self):
        self.cutflag = (self.byeye['class'] == 2) | (self.byeye['class'] == 4) | (self.byeye['class'] == 0)

        # removing sources with AGC only
        # will match to A100 instead

        self.agconly = ~self.kitchen['HLflag'] & ~self.kitchen['NSAflag'] & self.kitchen['AGCflag'] & ~self.kitchen['NSA0flag']

        self.cutflag = self.cutflag | self.agconly

        self.cleancat = self.byeye[~self.cutflag]
        #remove first column, which is a duplicate with second column
        n = self.cleancat.colnames
        self.cleancat.remove_column(n[0])
        
        # append original galaxy id to new kitchen sink file
        c = Column(self.byeye['galnumber'],name='galnumber')
        self.kitchen.add_column(c)
        self.clean_kitchen = self.kitchen[~self.cutflag]
        
        # cut ra,dec,vel
        c1 = Column(self.ra[~self.cutflag],'RAtemp')
        c2 = Column(self.dec[~self.cutflag],'DECtemp')
        c3 = Column(self.vel[~self.cutflag],'vrtemp')
        self.clean_kitchen.add_columns([c1,c2,c3])
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
        joined_table2.remove_columns(['RAtemp','DECtemp','vrtemp'])
        c3 = Column(ra,name='RA',dtype='f')
        c4 = Column(dec,name='DEC',dtype='f')
        c5 = Column(vel,name='vr',dtype='f')

        joined_table2.add_columns([c3,c4,c5])

        ###############################################
        ## FIX BOOLEAN COLUMNS
        ###############################################
        # boolean columns are getting converted weird
        # when I write and then read the fits table
        '''
        try:
            joined_table2['NSAflag'] = (joined_table2['NSAflag'] == 84)
            joined_table2['HLflag'] = (joined_table2['HLflag'] == 84)
        except KeyError:
            print('trouble in paradise')
        '''
        self.clean_a100 = joined_table2
        # inspect new a100 sources
        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']
        print('number of A100-only before cleaning = ',sum(a100flag))
        # look for any duplicate entries from A100 (probably from WISE matching)
        keepflag = np.ones(len(self.clean_a100['AGC']),'bool')
        agcnames = self.clean_a100['AGC']
        
        for i,a in enumerate(agcnames):
            if a100flag[i]:
                matches = np.sum(self.clean_a100['AGC'] == a)
                if matches > 1:
                    keepflag[i] = False
        self.clean_a100 = self.clean_a100[keepflag]


        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']
        print('number of A100-only after cleaning = ',sum(a100flag))

        self.write_clean()
    def write_clean(self):
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)

    def check_new_a100(self,plotflag=False):
                           
        a100flag = ~self.clean_a100['HLflag'] & ~self.clean_a100['NSAflag'] & self.clean_a100['A100flag']
        if plotflag:
            self.plot_all(ra=self.clean_a100['RAdeg_Use'],dec=self.clean_a100['DECdeg_Use'],flag=a100flag,names=self.clean_a100['AGC'],ra2=self.clean_a100['RA'],dec2=self.clean_a100['DEC'])

        # results from visual inspection
        # the following sources should be merged

        keepflag = np.ones(len(self.clean_a100['RA']),'bool')

        child = np.array([9144, 9145, 9147, 9151],'i')
        parent = np.array([7423, 6575, 6638, 8949],'i')
        for i in range(len(child)):
            print('merging {} with {}'.format(child[i],parent[i]))
            self.merge_sources(parent[i],child[i],cat=self.clean_a100,HL=False,NSA=False,AGC=False,A100=True)
            self.clean_a100['A100flag'][parent[i]] = True
        # for 8949, 9151 pair, AGC 208736, use AGC coordinates
        # delete rows corresponding to children
        keepflag[child] = np.zeros(len(child),'bool')
        self.clean_a100 = self.clean_a100[keepflag]
    def get_NEDname_query(self,startindex=0):
        ## GOT BLOCKED BY NED FOR TOO MANY QUERIES
        ## TRYING ANOTHER APPROACH - TO MATCH TO CATALOG I DOWNLOADED FROM DEC
        
        # look up NED name for each galaxy
        # https://astroquery.readthedocs.io/en/latest/ned/ned.html


        # if in HL, look up HL name
        NEDid = []
        NEDra = []
        NEDdec = []        
        for i in range(startindex,len(self.clean_a100['objname'])):
            # check if HL id exists, look up NED name
            # if yes, break
            # if no NED comes back
            # check if NSA ID exists for that galaxy. if yes, look up NED name
            # if NED name is found, break
            # if no NED name comes back
            # check if A100 ID exists for that galaxy. if yes, look up NED name
            # if no Ned name comes back, then no NED name!
            foundit=False
            if self.clean_a100['HLflag'][i]:
                time.sleep(1)
                try:
                    t = Ned.query_object(self.clean_a100['objname'][i])
                    NEDid.append(t['Object Name'][0])
                    foundit=True
                    continue
                except IndexError:
                    pass
                except:
                    print(i,'2 NED did not like ',self.clean_a100['objname'][i])
                    pass
                

            if self.clean_a100['NSAflag'][i]:
                if not(foundit):
                    time.sleep(1)
                    try:
                        t = Ned.query_object('NSA '+str(self.clean_a100['NSAID'][i]))
                        NEDid.append(t['Object Name'][0])
                        foundit=True
                        continue
                    except IndexError:
                        pass
                    except:
                        print(i,'2 NED did not like ','NSA '+str(self.clean_a100['NSAID'][i]))
                        pass
            if self.clean_a100['A100flag'][i]:
                if not(foundit):
                    time.sleep(1)
                    try:
                        #print('AGC'+str(self.clean_a100['AGC'][i]))
                        t = Ned.query_object('AGC'+str(self.clean_a100['AGC'][i]))
                        NEDid.append(t['Object Name'][0])
                        NEDra.append(t['RA'][0])
                        NEDdec.append(t['DEC'][0])                        
                        foundit=True
                        continue
                    except IndexError:
                        pass
                    except:
                        print(i,'2 NED did not like ','AGC'+str(self.clean_a100['AGC'][i]))
                        pass
            if self.clean_a100['NSA0flag'][i]:
                if not(foundit):
                    time.sleep(1)
                    try:
                        t = Ned.query_object('NSA '+str(self.clean_a100['NSAID_2'][i]))
                        NEDid.append(t['Object Name'][0])
                        NEDra.append(t['RA'][0])
                        NEDdec.append(t['DEC'][0])                        
                        foundit=True
                        continue
                    except IndexError:
                        pass
                    except:
                        print(i,'2 NED did not like ','NSA '+str(self.clean_a100['NSAID'][i]))
                        pass
            if not(foundit):
                # search by coordinates
                time.sleep(1)
                coord = SkyCoord(self.clean_a100['RA'][i],self.clean_a100['DEC'][i],unit=(u.deg,u.deg), frame='icrs')
                print(i,self.clean_a100['RA'][i],self.clean_a100['DEC'][i])

                try:
                    t = Ned.query_region(coord,radius=10./3600*u.deg, equinox='J2000')
                    print(t)
                    #t = Ned.query_object('NSA '+str(self.clean_a100['NSAID_2'][i]))
                    NEDid.append(t['Object Name'][0])
                    NEDra.append(t['RA'][0])
                    NEDdec.append(t['DEC'][0])                        
                    foundit=True
                    continue
                except IndexError:
                    pass
                except:
                    print(i,'NED did not like search by coordinate for this object')

            if not(foundit):
                print("oh no - could not find NED name! so sorry...",i)
                NEDid.append('')
                NEDra.append(-999)
                NEDdec.append(-999)                        
        c1 = Column(NEDid,name='NEDname')
        c2 = Column(np.array(NEDra,'f'),name='NEDra')
        c3 = Column(np.array(NEDdec,'f'),name='NEDdec')        
        try:
            nedtable = Table([c1,c2,c3]).write('ned_names.fits',format='fits',overwrite=True)
        except:
            print("couldn't write ned names")
            
        self.clean_a100.add_columns([c1,c2,c3])
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)
    def get_NEDname(self):
        # matching to catalog I downloaded on March 25, 2020
        # read in NED file
        '''
        Search by by Parameters
        decmin = -35  ** this is different from first download
        decmax = 75
        ramax = 280.
        ramin = 100.
        vmax = 3300. (units km/s)
        vmin = 500.
        object = Galaxy
        
        RA has to be in hours
        ramin = 6.6666666667
        ramax = 18.666666667

        output = text, ascii, bar separated
        velocity lower limit = -99

        Had to run this on virgo b/c I got locked out of NED.
        can't access the website from any machine in my house!

        Dowloaded text file - there was garbage at top of file (unnecessary info) that I deleted.

        Saved as

        /Users/rfinn/github/Virgo/tables/ned-noprolog-25mar2020.txt
        '''
        
        self.read_ned()
        self.cull_ned()

        # create empty array of names
        NEDid = np.zeros(len(self.clean_a100['RA']),dtype='|S25')

        
        # set up SkyCoord for clean cat and NED
        n = SkyCoord(self.ned['RA'],self.ned['DEC'],frame='icrs',unit='deg')
        c = SkyCoord(self.clean_a100['RA'],self.clean_a100['DEC'],frame='icrs',unit='deg')        
        # find matches
        idn, d2d, d3d = c.match_to_catalog_sky(n) # matches NED to VF catalog
        #print(d2d)

        nedmatchflag = (d2d < 10.*u.arcsec) # keep matches within 15arcsec

        indices = np.arange(len(idn))
        #print(indices[0:50])
        #print(len(indices),len(NEDid),len(idn))
        #NEDid[indices[nedmatchflag]] = self.ned['Object Name'][idn[nedmatchflag]]
        # transpose matched NED names into clean_a100
        
        # for each galaxy in cleaned catalog, find ned sources within 10 arcsec
        self.NEDmultiples = np.zeros(len(self.clean_a100),'bool')

        multiples = 0
        '''
        procedure:
        - calc offset between VF galaxy and NED galaxies
        - keep those w/in 15" offset
        - if only one source, assign VF galaxy this NED id
        - if multiples
          - try to match the HL name to the matching NED sources, if HL name exists
            - set multiple flag to False if match is found
          - look for source with name starting with NGC, UGC, IC.  use this if found
          - otherwise set to first NED source in the list
        '''
        
        for i in range(len(self.clean_a100)):
            d = np.sqrt((self.clean_a100['RA'][i] - self.ned['RA'])**2 \
                        + (self.clean_a100['DEC'][i] - self.ned['DEC'])**2)
            matchflag = d < 15./3600
            if np.sum(matchflag) > 1:
                self.NEDmultiples[i] = True
                #if multiples < 30:
                #    print('number of matches = ',np.sum(matchflag),multiples)
                #    print(self.ned['Object Name', 'Redshift Flag','Redshift Points'][matchflag],'\nHLname = ',self.clean_a100['objname'][i],'\n')
                multiples += 1

                
                # if we get multiple matches, then look at names and find one with
                # closest match to objname if it exist.
                indexmultiples = np.arange(len(self.ned))[matchflag]
                if len(self.clean_a100['objname'][i]) > 1: # has a HL name
                    n1 = self.clean_a100['objname'][i] # HL name
                    foundmatch = False
                    # look for match with HL name
                    for j in indexmultiples:
                        n2 = self.ned['Object Name'][j] # NED name
                        if (n2[0:2] == n1[0:2]) & (n2[-2:] == n1[-2:]):
                            # start and end of names match
                            # so this is the correct NED name
                            NEDid[i] = n2
                            foundmatch = True
                            self.NEDmultiples[i] = False
                            print('found HL name: ',n1,n2)
                            break
                if (foundmatch == False): # continue looking for the best match
                    # if not, look for more common names like NGC, UGC, IC
                    for j in indexmultiples:
                        if (self.ned['Object Name'][j].startswith('NGC')) | \
                           (self.ned['Object Name'][j].startswith('UGC')) | \
                           (self.ned['Object Name'][j].startswith('IC')):
                            # NED name
                            # start and end of names match
                            NEDid[i] = self.ned['Object Name'][j]
                            foundmatch = True
                            print('found a nice name: ',self.clean_a100['objname'][i],NEDid[i])
                            break
                if (foundmatch == False):
                    # assign closest match
                    # set nedquery flag to true, and query these objects individually
                    closest_match = np.arange(len(self.clean_a100))[np.where(d == min(d))[0][0]]
                    NEDid[i] = self.ned['Object Name'][closest_match]
                    print('taking closest match',self.clean_a100['objname'][i],NEDid[i])
                    # could also look at number of redshifts
            elif np.sum(matchflag) == 1:
                NEDid[i] = self.ned['Object Name'][np.where(matchflag)[0][0]]
                
                        
        print('number with multiple NED matches = ',multiples)
        print('number of NED matches = {}/{}'.format(sum(nedmatchflag),len(nedmatchflag)))
        
        c1 = Column(NEDid,name='NEDname')
        c2 = Column(self.NEDmultiples,name='NEDmultiples')# check these by hand?
        self.clean_a100.add_columns([c1,c2,c3])
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)
    def get_GL_NEDname(self):
        # for galaxies with no NED match, use GL's catalog to match
        # HL name to NEDname (he did a position match for those that didn't return a NED name
        # when searching by name
        ref = fits.getdata(homedir+'/github/Virgo/tables/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')
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

    def add_galid(self, cat=None):
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
        self.clean_a100 = Table(fits.getdata('vf_clean_sample.fits'))
    
    def fix_8822(self,cat=None):
        # add blue object, this is one with CO detection
        # the following information is from NED
        if cat is None:
            cat = self.clean_a100
        else:
            cat = cat
        colid = ['NEDname', 'RA','DEC','vr']
        colvalues  = ['UGC8656 NOTES01', 205.129878,42.993819,2899]

        # add an empty row
        cat.add_row()

        new_index = len(self.clean_a100)-1

        for i in range(len(colid)):
            cat[colid[i]][new_index] = colvalues[i]
        


    def write_clean_cat(self):
        self.cleancat.write('clean_sample.fits',format='fits',overwrite=True)
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)

    def catalog_for_z0MGS(self):
        '''
        need RA, DEC, and search radius
        in IPAC format table
        for matching with leroy+2019 Galaxy Synthesis WISE+GALEX
        table is served by IRSA
        '''
        search_radius = 10.*np.ones(len(self.cleancat)) # units are arcsec
        newtable = Table([self.cleancat['galnumber'],self.cleancat['RA'],self.cleancat['DEC'],search_radius],names=['galid','ra','dec','major'])
        newtable.write('clean_sample.txt',format='ipac',overwrite=True)
        
if __name__ == '__main__':
    c = catalog(kitchen_sink,byeye_classifications)
    #c.merge_class4()
    #c.cut_catalog()
    #c.catalog_for_z0MGS()
    c.runall()
