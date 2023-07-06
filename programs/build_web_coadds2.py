#!/usr/bin/env python

'''
GOAL:
* create web page to inspect the coadds, zp calibration, and psf
* run this from the html directory

/home/rfinn/research/Virgo-dev/html-dev

or more recently

/data-pool/Halpha/html_dev/


NOTES:
* rewriting after I have adopted a uniform naming convention - 2023-May-18

* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)



'''
# DONE TODO - use radius from main table for cutout size - maybe 2.5x bigger than boxes shown on full image
# DONE TODO - show cutout from CS-ZP if available
# DINE TODO - add legacy jpg image to cutouts
# DONE TODO - add r and Halpha to cutouts
# TODO - check color correction for 90prime filters with Matteo - color term persists
# eg https://facultyweb.siena.edu/~rfinn/virgo/coadds/VF-135.196+45.266-BOK-20210315-VFID1728/VF-135.196+45.266-BOK-20210315-VFID1728.html
# TODO - figure out why scale between BOK jpg and r-band/halpha images is wrong - must have wrong pixel scale for bok
# TODO - 

import os
import numpy as np
import glob
import sys

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table, Column
from astropy import wcs
from astropy.coordinates import SkyCoord

from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip
from astropy.time import Time
from scipy.stats import scoreatpercentile

from urllib.parse import urlencode
from urllib.request import urlretrieve

import multiprocessing as mp
#import pathos.multiprocessing as mp
#from concurrent.futures import ProcessPoolExecutor
#mp.set_start_method('spawn')

import argparse

from PIL import Image

homedir = os.getenv("HOME")
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'

import build_web_common as buildweb

sys.path.append(homedir+'/github/halphagui/')
import filter_transmission as ft

OVERWRITE = False
VERBOSE = False
###########################################################
####  FUNCTIONS
###########################################################
def buildone(rimages,i,coadd_dir,psfdir,zpdir,fratiodir):
    """ code to build webpage for one coadd """    
    rimage = rimages[i]
    #try:

    
    #print()
    #print('###################################')
    print(f'r-band image: {rimage} ({i}/{len(rfiles)})')
    #print('###################################')        
    #print()
    # find matching ha4 coadd
    #if rimage.find('shifted.fits') > -1:
    #    # not sure what this section is doing
    #    h1files = glob.glob(coadd_dir+'VF*-Halpha.fits')
    #    try:
    #        haimage = h1files[0]
    #    except IndexError:

    #        h2files = glob.glob(coadd_dir+'VF*-Ha6657.fits')
    #        
    #        haimage = h2files[0]
    #        print(haimage)
    #    if not os.path.exists(haimage):
    #        print('WHAT IS HAPPENING???')
    #        continue
    if i < 0:
        print("just kidding...")
    else:
        rheader = fits.getheader(rimage)
        try:
            haimage = os.path.join(coadd_dir,rheader['HAIMAGE'])
        except KeyError:
            print("couldn't find the halpha image name in header of ", rimage)                
            return
            
        if not os.path.exists(haimage):
            print("couldn't find the halpha image ",haimage, rimage)                
            return



        csimage = haimage.replace('.fits','-CS.fits')
        if not os.path.exists(csimage):
            print("couldn't find the CS halpha image ",csimage)                
            return

    #print('###  Halpha image = ',haimage)
    # define previous gal for html links
    if i > 0:
        previous = os.path.basename(rfiles[i-1]).replace('-R.fits','').replace('-shifted','').replace('-r.fits','').replace('.fits','').replace('-R','').replace('-r','')
        #print('previous = ',previous)
    else:
        previous = None
    if i < len(rfiles)-1:
        next = os.path.basename(rfiles[i+1]).replace('-R.fits','').replace('-shifted','').replace('-r.fits','').replace('.fits','').replace('-R','').replace('-r','')
        #print('next = ',next)
    else:
        next = None
    # define pointing name - remove fits and filter information
    pname = os.path.basename(rimage).replace('-R.fits','').replace('-shifted','').replace('-r.fits','').replace('.fits','').replace('-r','').replace('-R','')
    # create a d
    poutdir = os.path.join(outdir,pname)
    #print(poutdir)
    p = pointing(rimage=rimage,haimage=haimage,psfdir=psfdir,zpdir=zpdir,fratiodir = fratiodir, outdir=poutdir)
    h = build_html_pointing(p,outdir=poutdir,next=next,previous=previous)

    #try:
    #     p = pointing(rimage=rimage,haimage=haimage,psfdir=psfdir,zpdir=zpdir,outdir=poutdir)
    #     h = build_html_pointing(p,outdir=poutdir,next=next,previous=previous)
    #    
    #except:
    #    print("")
    #    print('WE HAVE A PROBLEM!!!',rimage)
    #    print("")            
    plt.close('all')
    #except KeyError:
    #    print("WARNING: could not build page for ",rimages[i])


image_results = []
def collect_results(result):

    global results
    image_results.append(result)

def get_legacy_jpg(ra,dec,galid='VFID0',pixscale=1,imsize='60',subfolder=None):
    """
    Download legacy image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001); used for naming the image files
    * imsize = size of cutout in pixels
    * pixscale = pixel scale of cutout in arcsec; native is 0.262 for legacy
    * subfolder = default is None; you can specify a name of a subfolder to 
                  save the data in, e.g., subfolder='legacy-images'
    Returns:
    * jpeg_name = jpeg image name
    """
    imsize = int(imsize)

    # make output image names
    if subfolder is not None:
        # check if subfolder exists. if not, make it.
        if not os.path.exists(subfolder):
            os.mkdir(subfolder)
        rootname = subfolder+'/'+str(galid)+'-legacy-'+str(imsize)
    else:
        rootname = str(galid)+'-legacy-'+str(imsize)        
    jpeg_name = rootname+'.jpg'

    if VERBOSE:
        print('legacy image name = ',jpeg_name)
        print('legacy imsize = ',imsize)
    
    # check if images already exist
    # if not download images
    if not(os.path.exists(jpeg_name)):
        #print('retrieving ',jpeg_name)
        url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)
        #print('legacy url = ',url)
        urlretrieve(url, jpeg_name)
    else:
        
        print('previously downloaded ',jpeg_name)


    # return the name of the fits images and jpeg image
    return jpeg_name

def display_image(image,percent=99.5,lowrange=False,mask=None,sigclip=True,csimage=False):
    lowrange=False
    # use inner 80% of image
    xdim,ydim = image.shape
    if xdim > 1000:
        xmin = int(.1*xdim)
        xmax = int(.9*xdim)    
        ymin = int(.1*ydim)
        ymax = int(.9*ydim)
    else:
        xmin = 1
        xmax = xdim
        ymin = 1
        ymax = ydim
    if sigclip:
        clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=2,sigma_upper=3)#,grow=10)
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]
    
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        try:
            norm = simple_norm(clipped_data, stretch='asinh',percent=percent)
        except IndexError:
            norm = None

    if norm is not None:
        plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    else:
        v1,v2=scoreatpercentile(image,[.5,99.5])            
        plt.imshow(image, cmap='gray_r',vmin=v1,vmax=v2,origin='lower')    



def write_coadd_prop_table(html,filter,zp,fwhm_arcsec):
    html.write('<h3>Image Characteristics</h3>\n')
    html.write('<table>\n')
    html.write('<tr>')
    html.write('<th ">Filter</th>\n')
    html.write('<th ">ZP<br />mag</th>\n')
    html.write('<th ">PSF FWHM <br />arcsec</th>\n')
    html.write('</tr>\n')
    html.write('<tr><td>{}</td><td>{:.2f}</td><td>{:.2f}</td>\n'.format(filter,zp,fwhm_arcsec))
    html.write('</tr>\n')
    html.write('</table>\n')


###########################################################
####  CLASSES
###########################################################

class coadd_image():

    def __init__(self,imagename,psfimage=None,plotdir=None,cat=None,zpdir=None,filter=None):
        self.imagename = imagename
        self.psf_allstars_png = None
        self.psf_png = None                
        if psfimage is not None:
            self.psf_flag = True
            self.psf_image = psfimage
            self.psfdir = os.path.split(psfimage)[0]
        else:
            self.psf_flag = False
            self.fwhm_arcsec = None
            self.sefwhm_arcsec = None            

        if plotdir is None:
            self.plotdir = os.getcwd()
        else:
            self.plotdir = plotdir
        if cat is None:
            self.cat = fits.getdata(VFMAIN_PATH)
        self.filter = filter
        self.plotprefix = os.path.join(self.plotdir,filter+'-')
        self.zpdir = zpdir
        # might need to comment this out
        temp = self.imagename.replace('-shifted','').replace('.fits','')
        pointing = temp.split('-')[-2].replace('p','pointing')
        self.intprefix = "{}*_{}".format(pointing,temp[-1])
        #print('INT plot prefix = ',self.intprefix)
    def generate_plots(self):
        self.get_image()
        self.make_coadd_png()
        if self.psf_flag:
            self.get_psf_image()
            if self.found_psf:
                self.make_psf_png()
                self.get_psf_allstars()

        #try:
        #    self.get_zpplot_firstpass()
        #except:
        #    print('WARNING: problem getting zp calibration images ',self.imagename)
        self.zp_flag = True
        #self.get_zpimage_firstpass()
        if self.filter != 'CS':
            self.get_zpplot_firstpass()            
            self.get_zp_magcomp_firstpass()        
        #self.get_zpplot_secondpass()                
        #self.get_zpimage_secondpass()
        #self.get_zp_magcomp_secondpass()


        if self.filter == 'ha':
            self.get_gredshift_filter_curve()
    def get_image(self):
        '''  read in image, save data and header '''
        self.imdata,self.imheader = fits.getdata(self.imagename,header=True)
        self.wcs = wcs.WCS(self.imheader)
        self.xdim,self.ydim = self.imdata.shape
        self.racenter,self.deccenter = self.wcs.wcs_pix2world(self.xdim/2,self.ydim/2,1)
        try:
            self.zp = self.imheader['PHOTZP']
        except KeyError:
            self.zp = -1
        try:
            self.pscale = np.abs(self.imheader['PIXSCAL1'])
        except KeyError:
            try:
                self.pscale = np.abs(self.imheader['CD1_1'])*3600
            except KeyError:
                self.pscale = np.abs(self.imheader['CDELT1'])*3600
        try:
            self.exptime = self.imheader['ORIGEXPT']
        except KeyError:
            self.exptime = self.imheader['EXPTIME']
        try:
            self.sefwhm_arcsec = self.imheader['SEFWHM']
        except KeyError:
            try:
                # look at the unshifted image - because I made a boo boo
                imdata,imheader = fits.getdata(self.imagename.replace('-shifted',''),header=True)            
                self.sefwhm_arcsec = imheader['SEFWHM']
                self.imheader.set('SEFWHM',self.sefwhm_arcsec)
            except KeyError:
                try:
                    # look at the unshifted image - because I made a boo boo
                    imdata,imheader = fits.getdata(self.imagename,header=True)            
                    self.sefwhm_arcsec = imheader['SEFWHM']
                    self.imheader.set('SEFWHM',self.sefwhm_arcsec)
                except KeyError:
                    self.sefwhm_arcsec = None
        try:
            t = self.imheader['DATE-OBS']
            t = Time(t,format='isot')
            self.dateobs = t.iso.split()[0]
            self.utobs = t.iso.split()[1]
            
        except KeyError:
            try:
                t = self.imheader['EPOCH']
                # convert to year, month,day
                t = Time(t,format='decimalyear')
                self.dateobs = t.iso.split()[0]
                self.utobs = t.iso.split()[1]
            except:
                self.dateobs = None
                self.utobs = None

    def make_coadd_png(self):
        ''' display image, and mark position of galaxies '''
        self.coadd_png = self.plotprefix+'coadd.png'

        imx,imy,keepflag = buildweb.get_galaxies_fov(self.imagename,self.cat['RA'],self.cat['DEC'])
        self.galfov_imx = imx[keepflag]
        self.galfov_imy = imy[keepflag]
        # where are we cutting based on filter redshift?        
        # need to cut to keep the galaxies within the right filter
        self.keepflag = keepflag        
        if os.path.exists(self.coadd_png) and not OVERWRITE:
            print('Found {}.  not remaking this.'.format(self.coadd_png))
        else:
            plt.figure(figsize=(8,8))
            ax = plt.subplot(projection=wcs.WCS(self.imheader))
            plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
            if self.filter == 'CS':
                display_image(self.imdata,csimage=True,sigclip=True)
            else:
                display_image(self.imdata)
                
            galsize=60/self.pscale
            buildweb.plot_vf_gals(imx,imy,keepflag,self.cat,ax,galsize=galsize)
            ax.set_xlabel('RA (deg)',fontsize=16)
            ax.set_ylabel('DEC (deg)',fontsize=16)        
            plt.savefig(self.coadd_png)
            plt.close()
    def get_psf_image(self):
        ''' get psf image, store FWHM '''
        if os.path.exists(self.psf_image):
            self.psfdata,self.psfheader = fits.getdata(self.psf_image,header=True)
            self.fwhm_pix = self.psfheader['FWHM']
            self.fwhm_arcsec = self.fwhm_pix*self.pscale
            self.found_psf = True
        else:
            self.found_psf = False
            self.fwhm_arcsec = -1
            self.fwhm_pix = -1
        pass
    
    def make_psf_png(self):
        ''' read in psf image, save data and header '''
        # check that png file exists
        # display png file
        norm = simple_norm(self.psfdata, 'log', percent=99.)
        plt.figure(figsize=(8,8))
        plt.subplots_adjust(right=.9,top=.95,left=.1,bottom=.05)
        plt.imshow(self.psfdata, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar(fraction=.046,pad=.04)
        #plt.show()
        self.psf_png = self.plotprefix+'psf.png'
        plt.savefig(self.psf_png)        
        plt.close()
    def get_psf_allstars(self):
        ''' display psf image mosaic of 100 stars '''
        # check that png file exists
        # display png file
        imagebase = os.path.basename(self.imagename).replace('.fits','').replace('-shifted','')
        #print(imagebase)
        #print('plotdir = ',self.plotdir)
        zpsurf = os.path.join(self.psfdir,'plots',imagebase+"-allstars.png")
        if not os.path.exists(zpsurf):
            imagebase = os.path.basename(self.imagename).replace('.fits','')
            zpsurf = os.path.join(self.psfdir,'plots',imagebase+"-allstars.png")
        print('allstars source = ',zpsurf)
        self.psf_allstars_png = os.path.join(self.plotdir,imagebase+"-allstars.png")
        #print('allstars destination = ',self.psf_allstars_png)
        os.system('cp '+zpsurf+' '+self.psf_allstars_png)



    def get_zpplot_firstpass_old(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        #print(self.imagename)

        #cluge b/c Becky's filenames threw a wrench in naming conventions
        
        #print(imagebase)
        #print('plotdir = ',self.plotdir)
        try:
            imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','')
            zpsurf = glob.glob(os.path.join(self.zpdir,imagebase+"*getzp-xyresidual-fitted.png"))[0]
        except IndexError: # the above won't work for INT data b/c my naming conventions are a mess
            # and because I ran the zp calibration from diff directory
            #print('search path = ',os.path.join(self.zpdir,self.intprefix+"*getzp-xyresidual-fitted.png"))
            zpsurf = glob.glob(os.path.join(self.zpdir,self.intprefix+"*getzp*fitted*.png"))[0]
        zpplot_png = os.path.join(self.plotdir,'n'+imagebase+"-getzp-xyresidual-fitted.png")
        if os.path.exists(self.zpplot_png):
            os.system('cp '+zpsurf+' '+self.zpplot_png)
        else:
            if os.path.exists(zpplot_png):
                self.pzplot_png = zpplot_png
        pass
    def get_zpplot_firstpass(self):
        ''' get the zp image, first pass'''
        imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','')
        print(self.zpdir,imagebase,"-getzp-xyresidual-fitted.png")
        zpsurf = os.path.join(self.zpdir,imagebase+"-getzp-xyresidual-fitted.png")
        self.zpplot_png = os.path.join(self.plotdir,imagebase+"-getzp-xyresidual-fitted.png")
        os.system('cp '+zpsurf+' '+self.zpplot_png)

        
        pass
    def get_zpplot_secondpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        #print(self.imagename)
        try:
            imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
            #print(imagebase)
            #print('plotdir = ',self.plotdir)
            zpplot2 = glob.glob(os.path.join(self.zpdir,"f"+imagebase+"*getzp-xyresidual-fitted.png"))[0]
            self.zpplot2_png = os.path.join(self.plotdir,"f"+imagebase+"-getzp-xyresidual-fitted.png")
            os.system('cp '+zpplot2+' '+self.zpplot2_png)
        except:
            self.zpplot2_png = None
        pass
    def get_zpimage_firstpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        #print(self.imagename)
        try:
            imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
            #print(imagebase)
            #print('plotdir = ',self.plotdir)
            #print("Looking for : ",os.path.join(self.zpdir,imagebase+"*imsurfit-2*.png"))
            zpsurf = glob.glob(os.path.join(self.zpdir,imagebase+"*imsurfit-2*.png"))[0]
        except IndexError:
            zpsurf = glob.glob(os.path.join(self.zpdir,self.intprefix+"*imsurfit-2*.png"))[0]
        try:
            self.zpsurf_png = os.path.join(self.plotdir,imagebase+"-imsurfit-2.png")
            os.system('cp '+zpsurf+' '+self.zpsurf_png)
        except:
            self.zpsurf_png = None
        pass
    
    def get_zpimage_secondpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        try:
            imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
            #print(imagebase)
            #print('plotdir = ',self.plotdir)
            zpsurf = glob.glob(os.path.join(self.zpdir,"f"+imagebase+"*imsurfit-2-round2.png"))[0]
        except IndexError:
            zpsurf = glob.glob(os.path.join(self.zpdir,"f"+self.intprefix+"*imsurfit-2-round2.png"))[0]
        try:
            self.zpsurf2_png = os.path.join(self.plotdir,"f"+imagebase+"-imsurfit-2-round2.png")
            os.system('cp '+zpsurf+' '+self.zpsurf2_png)
        except:
            self.zpsurf2_png = None
        pass
        

    def get_zp_magcomp_firstpass(self):
        ''' get the final plot of inst mag vs panstarrs mag'''
        # check that png file exists
        # display png file
        #print('imagename = ',self.imagename)

        
        imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','')
        ##print('imagebase = ',imagebase)        
        #pancomp = glob.glob(self.zpdir+'/'+imagebase+"*se-pan-flux.png")[0]
        #self.pancomp_png = self.plotdir+'/'+imagebase+"-se-pan-flux.png"
        ##print('pancomp_png = ',self.pancomp_png)
        #os.system('cp '+pancomp+' '+self.pancomp_png)

        # flux comp
        fluxcomp = os.path.join(self.zpdir,imagebase+"-se-pan-flux.png")
        self.pancomp_png = os.path.join(self.plotdir,imagebase+"-se-pan-flux.png")
        os.system('cp '+fluxcomp+' '+self.pancomp_png)
        
        pass
    def get_zp_magcomp_secondpass(self):
        ''' get the final plot of inst mag vs panstarrs mag'''
        # check that png file exists
        # display png file
        #print('imagename = ',self.imagename)
        try:
            imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
            #print('imagebase = ',imagebase)        
            pancomp = glob.glob(self.zpdir+'/'+'f'+imagebase+"*se-pan-flux.png")[0]
            self.pancomp2_png = self.plotdir+'/'+'f'+imagebase+"-se-pan-flux.png"
            #print('pancomp_png = ',self.pancomp_png)
            os.system('cp '+pancomp+' '+self.pancomp2_png)
        except:
            self.pancomp2_png = None
        pass

    def get_html_data(self):
        labels = ['Date Obs','UT Time','Filter','ZP<br>(AB mag)','Max Exptime <br> (minutes)','PSF FWHM <br> (arcsec)','SE FWHM <br> (arcsec)']
        try:
            filter = self.imheader['FILTER']
        except KeyError:
            print()
            print("WARNING: could not get filter for image ",self.imagename)
            print()
            filter = 'None'
        data = [self.dateobs,self.utobs,\
                filter,\
                "{:.1f}".format(self.zp),\
                "{:.1f}".format(self.exptime/60)]
        if self.fwhm_arcsec is not None:
            data.append("{:.2f}".format(self.fwhm_arcsec))
        if self.sefwhm_arcsec is not None:
            data.append("{:.2f}".format(self.sefwhm_arcsec))
        return labels,data
    def get_gredshift_filter_curve(self):
        redshift = vmain['vr'][self.keepflag]/3.e5
        header_filter = self.imheader['FILTER']
        #print('filter from header = ',header_filter,self.filter)
        if header_filter.find('ha4') > -1:
            filter=4
        elif header_filter.find('Ha+4nm') > -1:
            # TODO - need to check that this is in fact the same filter as on HDI
            filter=4
        elif header_filter.find('Ha4nm') > -1:
            # TODO - need to check that this is in fact the same filter as on HDI
            filter=4
        elif header_filter.find('Ha6657') > -1:
            filter='intha6657'
        elif header_filter.find('Halpha') > -1:
            filter='inthalpha'
        myfilter = ft.filter_trace(filter)
        self.gals_filter_png = os.path.join(self.plotdir,'galaxies_in_filter.png')
        corrections = myfilter.get_trans_correction(redshift,outfile=self.gals_filter_png)
        filter_keepflag = corrections < 10 # this is a crazy big cut, but we can adjust with halphagui
        self.corrections = corrections[filter_keepflag]
        print()
        print(f"number of galaxies before filter cut = {np.sum(self.keepflag)}")
        self.keepflag[self.keepflag] = filter_keepflag
        print(f"number of galaxies AFTER filter cut = {np.sum(self.keepflag)}")        
        #self.gals_filter_png = os.path.join(self.plotdir,'galaxies_in_filter.png')
        #os.rename('galaxies_in_filter.png',self.gals_filter_png)
        #pass
    
class pointing():

    def __init__(self,rimage=None,haimage=None,psfdir=None,zpdir=None,outdir=None,fratiodir=None):
        '''
        INPUT:
        * rband image
        * halpha image
        * psf directory
        * directory with ZP images
        '''

        self.rimage = rimage
        self.haimage = haimage
        self.csimage = haimage.replace('.fits','-CS.fits')
        czimage = haimage.replace('.fits','-CS-ZP.fits')
        if os.path.exists(czimage):
            self.czimage = czimage
        else:
            self.czimage = None

        self.psfdir = psfdir
        self.zpdir = zpdir
        self.fratiodir = fratiodir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir
        header = fits.getheader(self.rimage)
        try:
            dateobs = header['DATE-OBS']
        except KeyError:
            try:
                dateobs = header['DATE']
            except:
                print('error: could not find date in image header')

        self.pointing_name = os.path.basename(self.rimage).replace('.fits','').replace('-r','').replace('-R','').replace('-shifted','')
        
        self.build_psf_names()
        print()        
        print("getting rband image")        
        self.get_rband_image()
        print()        
        print("getting halpha image")        
        self.get_halpha_image()
        print()        
        print("getting CS")        
        self.get_cs_image()
        print()        
        print("getting CS-ZP image")        
        self.get_cz_image()
        print()        
        print("getting galaxy cutouts")        
        self.get_gal_cutouts()
        print()        
        print("getting galaxy table")
        self.get_params_for_gal_table()
        self.write_gal_table()
        print()
        print("getting pointing position")
        self.plot_pointing_position()
        print()        
        print("getting filter ratio plot")        
        self.get_filter_ratio_plot()

    def build_psf_names(self):
        self.rpsf_image = self.psfdir+'/'+os.path.basename(self.rimage).replace('-shifted','').split('.fits')[0]+'-psf.fits'
        if os.path.exists(self.rpsf_image):
            self.rpsf_flag=True
        else:
            self.rpsf_image = self.psfdir+'/'+os.path.basename(self.rimage).split('.fits')[0]+'-psf.fits'
            if os.path.exists(self.rpsf_image):
                self.rpsf_flag=True
            else:
                self.rpsf_flag=False
                print('could not find r psf image: ',self.rpsf_image)
            
        self.hapsf_image = self.psfdir+'/'+os.path.basename(self.haimage).split('.fits')[0]+'-psf.fits'
        if os.path.exists(self.hapsf_image):
            self.hapsf_flag=True
        else:
            self.hapsf_flag=False
            print('could not find ha psf image: ',self.hapsf_image)            
      

    def get_rband_image(self):
        ''' initiate an instance of coadd image class '''
        
        if os.path.exists(self.rimage):
            outprefix = self.outdir
            filter='r'
            if self.rpsf_flag:
                self.r = coadd_image(self.rimage,psfimage=self.rpsf_image,plotdir=outprefix,zpdir=self.zpdir,filter=filter)
            else:
                self.r = coadd_image(self.rimage,psfimage=None,plotdir=outprefix,zpdir=self.zpdir,filter=filter)
            self.rcoadd_flag=True
            self.r.generate_plots()
        else:
            self.rcoadd_flag=False
        
    
    def get_halpha_image(self):
        ''' initiate an instance of coadd image class '''
        if os.path.exists(self.haimage):
            outprefix = self.outdir
            filter='ha'
            if self.hapsf_flag:
                self.ha = coadd_image(self.haimage,psfimage=self.hapsf_image,plotdir=outprefix,zpdir=self.zpdir,filter=filter)
            else:
                self.ha = coadd_image(self.haimage,psfimage=None,plotdir=outprefix,zpdir=self.zpdir,filter=filter)
            self.hacoadd_flag=True
            self.ha.generate_plots()
        else:
            self.hacoadd_flag=False
            print('could not find ha image : ',self.haimage)
        

    def get_cs_image(self):
        ''' initiate an instance of coadd image class '''
        if os.path.exists(self.csimage):
            outprefix = self.outdir
            filter='CS'
            self.cs = coadd_image(self.csimage,psfimage=None,plotdir=outprefix,zpdir=None,filter=filter)
            print('\tmaking cs plots')
            self.cs.generate_plots()            
            self.cscoadd_flag=True
            #print()
            #print('getting galaxy cutouts')
            print('\tcs image getting gal cutouts')            
            self.get_gal_cutouts()
        else:
            self.cscoadd_flag=False
            print('could not find CS ha image : ',self.csimage)
        
    def get_cz_image(self):
        ''' initiate an instance of coadd image class '''
        if (self.czimage is not None) and os.path.exists(self.czimage):
            outprefix = self.outdir
            filter='CS'
            self.cz = coadd_image(self.csimage,psfimage=None,plotdir=outprefix,zpdir=None,filter=filter)
            self.cz.generate_plots()
            self.czcoadd_flag=True
            #print()
            #print('getting galaxy cutouts')
            #self.get_gal_cutouts()
        else:
            self.czcoadd_flag=False
            print('could not find ZP CS ha image : ',self.czimage)
        

    def get_gal_cutouts(self,size=90):
        """ get cutouts galaxies in FOV """
        self.gal_cutout_images = []


        
        # loop over galaxies in FOV
        gindex=np.arange(len(self.r.galfov_imx))
        galnames = Table(self.r.cat)['prefix'][self.cs.keepflag]
        galra = Table(self.r.cat)['RA'][self.cs.keepflag]
        galdec = Table(self.r.cat)['DEC'][self.cs.keepflag]        

        ##
        # set size to 2.5 time size in coadd images
        ##
        galsizes = Table(self.r.cat)['radius'][self.cs.keepflag]*2
        
        #galsizes = size#self.rcat['radius']/.4*2
        if 'INT' in self.rimage:
            pixscale = 0.33
        elif 'HDI' in self.rimage:
            pixscale = 0.425
        elif 'BOK' in self.rimage:
            pixscale = 0.4533
        elif 'MOS' in self.rimage:
            pixscale = 0.425
        rimdata,rimheader = fits.getdata(self.rimage,header=True)
        himdata,himheader = fits.getdata(self.haimage,header=True)
        cimdata,cimheader = fits.getdata(self.csimage,header=True)
        if self.czimage is not None:
            czimdata,czimheader = fits.getdata(self.czimage,header=True)


        # loop to display other images
        if self.czimage is not None:
            images = [rimdata,himdata,cimdata,czimdata]
            imtitles = ['r','ha','cs from filt ratio','cs from ZP ratio']
        else:
            images = [rimdata,himdata,cimdata]
            imtitles = ['r','ha','cs ha']                
            
        #sizes = (galsizes/pixscale*3.,galsizes/pixscale*3.)

        rowchange = np.arange(4,50,4)
        nrow = 1
        for i in range(len(rowchange)):
            if len(gindex) > rowchange[i]:
                nrow += 1
            else:
                break

        # columns: legacy, r, halpha, cs, CS-zp
        ncol = 5
        nrow = np.sum(self.cs.keepflag)
        # change to one row per galaxy
        figsize = (12,3*np.sum(self.cs.keepflag))            
        plt.figure(figsize=figsize)
        plt.subplots_adjust(top=.95,right=.95,left=.05,bottom=.05)        
        for j in range(len(galra)):
            #print(sizes[j][0])
            try:
                imsize = galsizes[j]/pixscale*2.
            except IndexError:
                #print('hey rose - problem accessing sizes ',sizes)
                # set the default size to 60 arcsec
                imsize = 60/pixscale
            imsize_arcsec = imsize*pixscale
            # get legacy cutout
            # TODO - finish this next line
            ax = plt.subplot(nrow,ncol,5*j+1)            
            jpeg_name = get_legacy_jpg(galra[j],galdec[j],galid=galnames[j],pixscale=1,imsize=imsize_arcsec,subfolder=self.outdir)
            #print(jpeg_name)
            # plot jpg
            t = Image.open(jpeg_name)
            plt.imshow(t,origin='upper')
            plt.title(galnames[j][:20])# cutting names to avoid ridiculously long NED names
            plt.ylabel('arcsec')
            position = (self.r.galfov_imx[j],self.r.galfov_imy[j])                
            for k in range(len(images)):
                #print("displaying cutout ",imtitles[k],imsize)
                ax = plt.subplot(nrow,ncol,5*j+k+2)            
                cutout = Cutout2D(images[k],position,imsize)
                #if k > 1 :
                #    display_image(cutout.data,csimage=True)
                #else:
                display_image(cutout.data)                    
                plt.title(imtitles[k])
                # remove tick labels on greyscale images - png is in arcsec
                #plt.xticks([],[])
                #plt.yticks([],[])                
        imname = f"{self.pointing_name}-gal-cutouts.png"
        outfile = os.path.join(self.outdir,imname)
        plt.savefig(outfile)
        plt.close()
            
        # add image name to list
        self.gal_cutout_figname = imname
        
    def get_filter_ratio_plot(self):
        ''' display the filter ratio png '''
        # DONE, but need to test - TODO - populate this!
        plotdir = os.path.join(os.path.dirname(self.haimage),'plots-filterratio/')

        # need to construct plot name from the halpha image

        # for example: *Halpha-filter-ratio.png
        plotname = os.path.basename(self.haimage.replace('.fits','-filter-ratio.png'))

        filename = os.path.join(plotdir,plotname)
        #print("looking for filter ratio plot ",filename)
        if os.path.exists(filename):
            s = f'cp {filename} {os.path.join(self.outdir,os.path.basename(filename))}'
            os.system(s)
            self.filter_ratio_plot = os.path.basename(filename)
        else:
            self.filter_ratio_plot = None
        


    def get_params_for_gal_table(self):
        ''' setup arrays for galaxy table '''
        self.groupgals = self.r.cat[self.cs.keepflag]
        self.vffil = vffil[self.cs.keepflag]

    def write_gal_table(self):
        """ write out a list of with galaxies in FOV """
        # this will be good for checking which galaxies are observed in Halpha
        # in a way that is independent of running the halpha gui
        # just need the VFIDs

        # use the root file for the output name
        # keep the coadd name
        
        #outfile =
        outfile = os.path.join(self.outdir,self.pointing_name+'-galsFOV.csv')
        gals = self.r.cat['VFID'][self.cs.keepflag]

        output = open(outfile,'w')
        for i in range(len(gals)):
            s = f"{gals[i]}, {self.pointing_name} \n"
            output.write(s)
        output.close()
        pass

    def plot_pointing_position(self):
        ''' plot position of pointing wrt vf sample '''
        plt.figure(figsize=(8,8))
        plt.scatter(vmain['RA'],vmain['DEC'],marker='.',s=6,c=vmain['vr'],vmin=500,vmax=3200)
        plt.colorbar(fraction=0.046,pad=.04)
        px = self.r.imheader['CRVAL1']
        py = self.r.imheader['CRVAL2']
        boxsize=60
        plt.plot(px,py,'ks',markersize=14,mec='r',alpha=.6)
        #rect= plt.Rectangle((px-boxsize/2.,py-boxsize/2.), boxsize, boxsize,fill=True, color='k',lw=2,alpha=.5)
        #self.position_plot = self.outdir+'/positions.png'
        plt.gca().invert_xaxis()
        plt.xlabel('RA (deg)',fontsize=16)
        plt.ylabel('DEC (deg)',fontsize=16)        
        plt.title("NB: size of box is not to scale")
        self.position_plot = self.outdir+'/positions.png'        
        plt.savefig(self.position_plot)
        
        plt.close()
        
        # zoom
        zoombox=5
        plt.figure(figsize=(8,8))
        sp = plt.scatter(vmain['RA'],vmain['DEC'],marker='.',s=100,c=vmain['vr'],vmin=500,vmax=3200,label='VF')
        coflag = vmain['COflag']
        plt.plot(vmain['RA'][coflag],vmain['DEC'][coflag],'ro',markersize=12,mfc='None',label="CO")

        px = self.r.imheader['CRVAL1']
        py = self.r.imheader['CRVAL2']
        #print(px,py)
        pscale_deg = self.r.pscale/3600
            
        boxsizex=self.r.imheader['NAXIS1']*np.abs(float(pscale_deg))
        boxsizey=self.r.imheader['NAXIS2']*np.abs(float(pscale_deg))
        
        # correct the x (RA) dimension of the box for the cosine of the declination
        # boxes should appear bigger at you approach the north celestial pole
        boxsizex = boxsizex/np.cos(np.radians(py))

        #boxsizex,boxsizey=2,2
        rect= plt.Rectangle((px-boxsizex/2.,py-boxsizex/2.), boxsizex, boxsizey,fill=True, color='k',lw=2,alpha=.2)
        plt.gca().add_artist(rect)
        plt.colorbar(sp,fraction=0.046,pad=.04)        
        plt.xlabel('RA (deg)',fontsize=16)
        plt.ylabel('DEC (deg)',fontsize=16)        
        plt.axis([px-zoombox/2,px+zoombox/2,\
                  py-zoombox/2,py+zoombox/2])
        plt.gca().invert_xaxis()
        plt.legend()
        self.position_plot_zoom = self.outdir+'/positions-zoom.png'
        plt.savefig(self.position_plot_zoom)
        plt.close()
class build_html_pointing():

    def __init__(self,pointing,outdir,next=None,previous=None):
        self.pointing = pointing
        outfile = self.pointing.pointing_name+'.html'
        outfile = os.path.join(outdir,outfile)
        self.html = open(outfile,'w')
        self.next = next
        self.previous = previous
        self.htmlhome = 'index.html'
        self.build_html()
    def build_html(self):
        self.write_header()
        self.write_navigation_links()
        try:
            self.write_gal_table()
        except:
            print("problem writing gal table for ",self.pointing.pointing_name)
            
        try:
            self.write_pointing_location()
        except:
            print("problem writing pointing location for ",self.pointing.pointing_name)

        try:
            self.write_rband_header()
            self.write_rband_table()            
        except:
            print("problem writing rband header/table for ",self.pointing.pointing_name)

        
        try:
            self.write_rband_psf()
        except:
            print("problem writing rband psf for ",self.pointing.pointing_name)

        try:
            if self.pointing.r.zp_flag:
                self.write_rband_zp()
        except:
            print("problem writing rband zp for ",self.pointing.pointing_name)

                 
        #self.write_rband_div()
        try:
            self.write_ha_header()
            self.write_ha_table()
        except:
            print("problem writing ha header/table for ",self.pointing.pointing_name)

        try:
            self.write_ha_psf()
        except:
            print("problem writing ha psf for ",self.pointing.pointing_name)

        try:
            if self.pointing.ha.zp_flag:
                self.write_ha_zp()
        except:
            print("problem writing rband zp for ",self.pointing.pointing_name)

        try:
            self.write_cs_header()
            self.write_cs_table()
        except:
            print("problem writing cs header/table for ",self.pointing.pointing_name)

        # adding galaxy table again near the cutout images for easy comparison
        try:
            self.write_gal_table()
        except:
            print("problem writing gal table for ",self.pointing.pointing_name)

        try:
            self.write_cutouts_table()
        except:
            print("problem writing cutouts table for ",self.pointing.pointing_name)
            
        # add functions to write out CS image
        # first row
        # CS image, fratio plot, zoom in on galaxies in CS image?
        
        #self.write_ha_div()

        #self.write_footer()
        self.write_navigation_links()        
        self.close_html()
    def write_header(self):
        # title
        # home link
        # previous
        # next
        self.html.write('<html><body>\n')
        self.html.write('<style type="text/css">\n')
        self.html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        self.html.write('</style>\n')

        # Top navigation menu--
        #self.html.write('<h1>{}</h1>\n'.format(self.pointing.pointing_name))

        #self.html.write('<a href="../{}">Home</a>\n'.format(htmlhome))
        #self.html.write('<br />\n')
        #self.html.write('<a href="../{}">Next ({})</a>\n'.format(nexthtmlgalaxydir1, nextgalaxy[ii]))
        #self.html.write('<br />\n')
        #self.html.write('<a href="../{}">Previous ({})</a>\n'.format(prevhtmlgalaxydir1, prevgalaxy[ii]))
    def write_navigation_links(self):
        # Top navigation menu--
        self.html.write('<h1>{}</h1>\n'.format(self.pointing.pointing_name))

        self.html.write('<a href="../{}">Home</a>\n'.format(self.htmlhome))
        self.html.write('<br />\n')
        
        if self.previous is not None:
            previous = self.previous.replace('-shifted','')
            previoushtml = "{}.html".format(previous)
            self.html.write('<a href="../{}/{}">Previous ({})</a>\n'.format(previous,previoushtml, previous))
            self.html.write('<br />\n')        
        if self.next is not None:
            thisnext = self.next.replace('-shifted','')
            nexthtml = "{}.html".format(thisnext)
            self.html.write('<a href="../{}/{}">Next ({})</a>\n'.format(thisnext,nexthtml,thisnext))
        self.html.write('<p>')
        self.html.write('<p>')
        self.html.write('<p>')        



    def write_gal_table(self):
       # Add the properties of each galaxy.
        self.html.write('<h2>Galaxies in FOV</h2>\n')
        self.html.write('<table>\n')
        self.html.write('<tr>\n')
        self.html.write('<th>VFID</th>\n')
        self.html.write('<th>Galaxy</th>\n')
        #self.html.write('<th>Morphology</th>\n')
        self.html.write('<th>RA</th>\n')
        self.html.write('<th>Dec</th>\n')
        self.html.write('<th>vr<br>(km/s)</th>\n')
        self.html.write('<th>Filter<br>Cor</th>\n')                
        self.html.write('<th>D(25)<br>(arcmin)</th>\n')
        self.html.write('<th>CO</th>\n')
        self.html.write('<th>A100</th>\n')
        self.html.write('<th>Filament<br> Member</th>\n')
        self.html.write('<th>Nearest<br> Filament</th>\n')        
        self.html.write('</tr>\n')
        for i,g in enumerate(self.pointing.groupgals):
            #if '031705' in gal['GALAXY']:
            #    print(groupgal['GALAXY'])
            self.html.write('<tr>\n')
            self.html.write('<td>{}</td>\n'.format(g['VFID']))
            self.html.write('<td>{}</td>\n'.format(g['NEDname']))
            #typ = groupgal['MORPHTYPE'].strip()
            #if typ == '' or typ == 'nan':
            #    typ = '...'
            #self.html.write('<td>{}</td>\n'.format(typ))
            self.html.write('<td>{:.6f}</td>\n'.format(g['RA']))
            self.html.write('<td>{:.6f}</td>\n'.format(g['DEC']))
            self.html.write('<td>{:.0f}</td>\n'.format(g['vr']))
            try:
                self.html.write('<td>{:.2f}</td>\n'.format(self.pointing.ha.corrections[i]))
            except IndexError:
                print("WARNING: size mismatch between r and halpha catalogs!!!")
            self.html.write('<td>{:.3f}</td>\n'.format(g['radius']*2/60.))
            self.html.write('<td>{:.0f}</td>\n'.format(g['COflag']))
            self.html.write('<td>{:.0f}</td>\n'.format(g['A100flag']))
            self.html.write('<td>{}</td>'.format(self.pointing.vffil['filament_member'][i]))
            self.html.write('<td>{}</td>'.format(self.pointing.vffil['filament'][i]))            
            
            #if np.isnan(groupgal['PA']):
            #    pa = 0.0
            #else:
            #    pa = groupgal['PA']
            #self.html.write('<td>{:.2f}</td>\n'.format(pa))
            #self.html.write('<td>{:.3f}</td>\n'.format(1-groupgal['BA']))
            self.html.write('</tr>\n')
        self.html.write('</table>\n')

        # write out a table with group galaxies
    def write_pointing_location(self):
        ''' write table with rband coadd and location in sky '''
        labels=['Location','Zoom <br>5x5deg','Filter <b> Transmission']
        images=[os.path.basename(self.pointing.position_plot),\
                os.path.basename(self.pointing.position_plot_zoom),\
                os.path.basename(self.pointing.ha.gals_filter_png)]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)

        
    def write_rband_header(self):
        #self.html.write('<hr>')
        self.html.write('<h2>r-band Coadd</h2>\n')
    def write_rband_table(self):
        labels,data = self.pointing.r.get_html_data()
        buildweb.write_text_table(self.html,labels,data)

        
    def write_rband_psf(self):
        '''  make table with rband psf results '''
        labels=['Coadd','Stars','PSF']
        images=[os.path.basename(self.pointing.r.coadd_png)]
        if self.pointing.r.psf_allstars_png is not None:
            images.append(os.path.basename(self.pointing.r.psf_allstars_png))
        if self.pointing.r.psf_png is not None:
            images.append(os.path.basename(self.pointing.r.psf_png))
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)

    def write_rband_zp(self):
        ''' make table with rband zp fit '''
        labels=['Fit','Residuals']#,'Residual<br> Surface']
        images = [os.path.basename(self.pointing.r.pancomp_png),\
                  os.path.basename(self.pointing.r.zpplot_png)]#,\
                  #os.path.basename(self.pointing.r.flux_comp)]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)
        #if self.pointing.r.pancomp2_png is not None:
        #    labels=['2nd Pass Fit','Residuals','Residual<br> Surface']        
        #    images2 = [os.path.basename(self.pointing.r.pancomp2_png),\
        #               os.path.basename(self.pointing.r.zpplot2_png),\
        #               os.path.basename(self.pointing.r.zpsurf2_png)]
        #    buildweb.write_table(self.html,labels=labels,images=images2)

        

    def write_ha_header(self):
        ''' write start of halpha section '''
        #self.html.write('<hr>')
        self.html.write('<h2>Halpha Coadd</h2>\n')

    def write_ha_table(self):
        labels,data = self.pointing.ha.get_html_data()
        buildweb.write_text_table(self.html,labels,data)

    def write_ha_image(self):
        ''' write table with rband coadd and location in sky '''
        labels=['Coadd','Location','Zoom <br>5x5deg']
        images=[os.path.basename(self.pointing.ha.coadd_png),\
                os.path.basename(self.pointing.position_plot),\
                os.path.basename(self.pointing.position_plot_zoom)]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)

        
    def write_ha_psf(self):
        '''  make table with rband psf results '''
        labels=['Coadd','Stars','PSF']
        images = [os.path.basename(self.pointing.ha.coadd_png)]
        if self.pointing.ha.psf_allstars_png is not None:
            images.append(os.path.basename(self.pointing.ha.psf_allstars_png))
        if self.pointing.ha.psf_png is not None:
            images.append(os.path.basename(self.pointing.ha.psf_png))
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)

    def write_ha_zp(self):
        ''' make table with rband zp fit '''
        labels=['ZP Fit','Residuals']#,'Residual<br> Surface']
        images = [os.path.basename(self.pointing.ha.pancomp_png),\
                  os.path.basename(self.pointing.ha.zpplot_png)]
                  #os.path.basename(self.pointing.ha.zpsurf_png)]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)

    def write_cs_header(self):
        ''' write start of halpha section '''
        #self.html.write('<hr>')
        self.html.write('<h2>Continuum-Subtracted Halpha Coadd</h2>\n')

    def write_cs_table(self):
        ''' make table with rband zp fit '''

        labels=['Cont-Sub Image w/Filter Ratio']
        images = [os.path.basename(self.pointing.cs.coadd_png)]
        if self.pointing.cz.coadd_png is not None:
            images.append(os.path.basename(self.pointing.cs.coadd_png))
            labels.append('Cont-Sub Image w/ZPs')
        images.append(os.path.basename(self.pointing.filter_ratio_plot))
        labels.append('Filter Ratio')
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)
        

        #if self.pointing.ha.pancomp2_png is not None:
        #    labels=['2nd Pass Fit','Residuals','Residual<br> Surface']        
        #    images2 = [os.path.basename(self.pointing.ha.pancomp2_png),\
        #               os.path.basename(self.pointing.ha.zpplot2_png),\
        #               os.path.basename(self.pointing.ha.zpsurf2_png)]
        #    buildweb.write_table(self.html,labels=labels,images=images2)
    def write_cutouts_table(self):
        labels=['CS Cutout Images']#,'Residual<br> Surface']
        images = [self.pointing.gal_cutout_figname]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)

    def write_ha_div(self):
        # write coadd
        # write psf
        # write table with filter, zp, fwhm
        #self.html.write('<h2>Halpha Coadd</h2>\n')
        self.html.write('<table width="90%">\n')

        #self.html.write('<tr><th>Coadd<br>{}</th> <th>PSF images</th> <th>ZP Calib Orig</th> <th>ZP Calib Final</th></tr></p>\n'.format(os.path.basename(self.pointing.haimage)))
        self.html.write('<tr><th>Coadd</th> <th>PSF images</th> <th>ZP Calib Orig</th> <th>ZP Calib Final</th></tr></p>\n')
        
        pngfile = os.path.basename(self.pointing.ha.coadd_png)
        psfpng = os.path.basename(self.pointing.ha.psf_png)
        images = [pngfile,psfpng,self.pointing.ha.zpsurf_png,self.pointing.ha.pancomp_png]        

        self.html.write('<tr>')
        for i in images:
            self.html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
        self.html.write('</tr>\n')            
        self.html.write('</table>\n')
    
    def write_ha_table(self):
        # write coadd
        # write psf
        # write table with filter, zp, fwhm
        #write_coadd_prop_table(self.html,self.pointing.ha.imheader['FILTER'],self.pointing.ha.zp,self.pointing.ha.fwhm_arcsec)
        labels,data = self.pointing.ha.get_html_data()
        buildweb.write_text_table(self.html,labels,data)


        pass

    def write_footer(self):
        self.html.write('<br /><br />\n')
        self.html.write('<a href="../../{}">Home</a>\n'.format('index.html'))
        self.html.write('<br />\n')
        self.html.write('<a href="../../{}">Next ({})</a>\n'.format('testing','testing'))
        self.html.write('<br />\n')
        #self.html.write('<a href="../../{}">Previous ({})</a>\n'.format(prevhtmlgalaxydir1, prevgalaxy[ii]))
        self.html.write('<br />\n')
    
    def close_html(self):
        self.html.close()
# wrap

        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')

    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--laptop',dest = 'laptop', help='set if working on laptop')
    
    parser.add_argument('--coaddir',dest = 'coaddir', help='set to coadd directory')
    parser.add_argument('--psfdir',dest = 'psfdir', help='set to coadd directory')

    parser.add_argument('--oneimage',dest = 'oneimage',default=None, help='give full path to the r-band image name to run on just one image')
    parser.add_argument('--bokonly',dest = 'bokonly',default=False,action='store_true', help='run to rebuild bok pages only')        
    
     
    args = parser.parse_args()

    newnames = True
    if newnames:
        coadd_dir = '/media/rfinn/hdata/coadds/all-virgo-coadds/'
        coadd_dir = '/data-pool/Halpha/coadds/all-virgo-coadds/'                        
        zpdir = coadd_dir+'/plots/'
        fratiodir = coadd_dir+'/plots-filterratio/'
        vtabledir = homedir+'/research/Virgo/tables-north/v2/'
        vmain = fits.getdata(vtabledir+'vf_v2_main.fits')
        #homedir = '/mnt/qnap_home/rfinn/'
        VFFIL_PATH = vtabledir+'/vf_v2_environment.fits'
        vffil = fits.getdata(VFFIL_PATH)
        #psfdir = homedir+'/data/reduced/psf-images/'

        outpathbase = '/media/rfinn/hdata/'
        outpathbase = '/data-pool/Halpha/'        
        psfdir = outpathbase+'/psf-images/'
        outdir = outpathbase+'/html_dev/coadds/'

        # get list of r-band coadded images
        a = glob.glob(coadd_dir+'VF*INT*-r-shifted.fits')
        b = glob.glob(coadd_dir+'VF*HDI*-r.fits')
        c = glob.glob(coadd_dir+'VF*HDI*-R.fits')
        d = glob.glob(coadd_dir+'VF*BOK*-r.fits')
        d = glob.glob(coadd_dir+'VF*MOS*-R.fits')                 
        rfiles = a + b + c + d

        # changing this b/c I now store the halpha image name in the r-band header
        #halpha_names = ['ha4','Halpha','Ha6657','Ha4']        
        #a = glob.glob(coadd_dir+'VF*INT*-Halpha.fits')
        #b = glob.glob(coadd_dir+'VF*HDI*-Ha4.fits')
        #c = glob.glob(coadd_dir+'VF*HDI*-ha4.fits')
        #d = glob.glob(coadd_dir+'VF*BOK*-Ha6657.fits')         
        #hfiles = a + b + c + d

        
    #hfiles.sort()
    rfiles.sort()

    # just use one image if the argument flag was set

    if args.oneimage is not None:
        # make sure that the image exists
        if not os.path.exists(args.oneimage):
            print(f"Could not find {args.oneimage} - please check the r-band coadd name you provided")
            sys.exit()
        # find index in rfiles that corresponds to image
        try:
            coadd_index = rfiles.index(args.oneimage)
            indices = [np.arange(len(rfiles))[coadd_index]]
            print('when selecting one image, indices = ',indices,rfiles[indices[0]])
            buildone(rfiles,coadd_index,coadd_dir,psfdir,zpdir,fratiodir)
        except ValueError:
            rfiles = [args.oneimage]
            indices = np.arange(len(rfiles))
    else:
        indices = np.arange(len(rfiles))

        # just rebuild the coadd pages for bok images
        if args.bokonly:
            indices=[]
            for i in range(len(rfiles)):
                if 'BOK' in rfiles[i]:
                    indices.append(i)
                
        image_pool = mp.Pool(mp.cpu_count())
        #image_pool = mp.pool.ThreadPool(mp.cpu_count())    
        myresults = [image_pool.apply_async(buildone,args=(rfiles,i,coadd_dir,psfdir,zpdir,fratiodir,)) for i in indices]
    
        image_pool.close()
        image_pool.join()
        image_results = [r.get() for r in myresults]

    ##
    # skipping mp because of pickling errors - will get back to that someday...
    ##
    #for i in range(len(rfiles)):
    #    buildone(rfiles,i,coadd_dir,psfdir,zpdir,fratiodir)

    
    #with ProcessPoolExecutor(max_workers=24) as exe:
    #             exe.map(buildone,indices,args=(rfiles,i,coadd_dir,psfdir,zpdir,fratiodir))

    # build the web index
    cwd = os.getcwd()
    os.chdir(outdir)
    print("\nBuilding coadd index\n")
    os.system("python ~/github/Virgo/programs/build_coadd_index.py")
    os.chdir(cwd)

    #buildone(rimages,i,psfdir=psfdir,zpdir=zpdir,fratiodir = fratiodir, outdir=poutdir)

    """
    for i,rimage in enumerate(rfiles[startindex:]):
        print()
        print('###################################')
        print(f'r-band image: {rimage} ({i}/{len(rfiles)})')
        print('###################################')        
        print()
        # find matching ha4 coadd
        #if rimage.find('shifted.fits') > -1:
        #    # not sure what this section is doing
        #    h1files = glob.glob(coadd_dir+'VF*-Halpha.fits')
        #    try:
        #        haimage = h1files[0]
        #    except IndexError:
                
        #        h2files = glob.glob(coadd_dir+'VF*-Ha6657.fits')
        #        
        #        haimage = h2files[0]
        #        print(haimage)
        #    if not os.path.exists(haimage):
        #        print('WHAT IS HAPPENING???')
        #        continue
        if i < 0:
            print("just kidding...")
        else:
            rheader = fits.getheader(rimage)
            haimage = os.path.join(coadd_dir,rheader['HAIMAGE'])
            if not os.path.exists(haimage):
                print("couldn't find the halpha image ",haimage)                
                continue

            

            csimage = haimage.replace('.fits','-CS.fits')
            if not os.path.exists(csimage):
                print("couldn't find the CS halpha image ",csimage)                
                continue
            
        print('###  Halpha image = ',haimage)
        # define previous gal for html links
        if i > 0:
            previous = os.path.basename(rfiles[i-1]).replace('-R.fits','').replace('-shifted','').replace('-r.fits','').replace('.fits','').replace('-R','').replace('-r','')
            #print('previous = ',previous)
        else:
            previous = None
        if i < len(rfiles)-1:
            next = os.path.basename(rfiles[i+1]).replace('-R.fits','').replace('-shifted','').replace('-r.fits','').replace('.fits','').replace('-R','').replace('-r','')
            #print('next = ',next)
        else:
            next = None
        # define pointing name - remove fits and filter information
        pname = os.path.basename(rimage).replace('-R.fits','').replace('-shifted','').replace('-r.fits','').replace('.fits','').replace('-r','').replace('-R','')
        # create a d
        poutdir = os.path.join(outdir,pname)
        print(poutdir)
        p = pointing(rimage=rimage,haimage=haimage,psfdir=psfdir,zpdir=zpdir,fratiodir = fratiodir, outdir=poutdir)
        try:
            h = build_html_pointing(p,outdir=poutdir,next=next,previous=previous)
        except:
            print("problem builing webpage for ",p)
            print("skipping for now")
            print()
        
        #try:
        #     p = pointing(rimage=rimage,haimage=haimage,psfdir=psfdir,zpdir=zpdir,outdir=poutdir)
        #     h = build_html_pointing(p,outdir=poutdir,next=next,previous=previous)
        #    
        #except:
        #    print("")
        #    print('WE HAVE A PROBLEM!!!',rimage)
        #    print("")            
        plt.close('all')

        # uncomment the next line for testing
        #break
    """
