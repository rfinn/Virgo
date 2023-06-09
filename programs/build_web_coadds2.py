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


import argparse

homedir = os.getenv("HOME")
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'

import build_web_common as buildweb

sys.path.append(homedir+'/github/halphagui/')
import filter_transmission as ft

OVERWRITE = True

###########################################################
####  FUNCTIONS
###########################################################
image_results = []
def collect_results(result):

    global results
    image_results.append(result)


def display_image(image,percent=99.5,lowrange=False,mask=None,sigclip=False,csimage=False):
    lowrange=False
    # use inner 80% of image
    xdim,ydim = image.shape
    xmin = int(.2*xdim)
    xmax = int(.8*xdim)    
    ymin = int(.2*ydim)
    ymax = int(.8*ydim)    
    if sigclip:
        clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]
    if csimage:
        try:
            norm = simple_norm(clipped_data, stretch='linear',min_percent=10,max_percent=90)
        except:
            norm = None
    elif lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    if norm == None:
        plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    else:
        plt.imshow(image,cmap='gray_r',origin='lower')
    


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
        print('INT plot prefix = ',self.intprefix)
    def generate_plots(self):
        self.get_image()
        self.make_coadd_png()
        if self.psf_flag:
            self.get_psf_image()
            if self.found_psf:
                self.make_psf_png()
                self.get_psf_allstars()
        try:
            self.get_zpplot_firstpass()
        except:
            print('WARNING: problem getting zp calibration images')
        self.zp_flag = True
        #self.get_zpimage_firstpass()
        if self.filter != 'CS':
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
            self.pscale = np.abs(self.imheader['CD1_1'])*3600
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
        # need to cut to keep the galaxies within the right filter
        self.keepflag = keepflag        
        if os.path.exists(self.coadd_png) and not OVERWRITE:
            print('Found {}.  not remaking this.'.format(self.coadd_png))
        else:
            plt.figure(figsize=(8,8))
            ax = plt.subplot(projection=wcs.WCS(self.imheader))
            plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
            if self.filter == 'CS':
                display_image(self.imdata,csimage=True)
            else:
                display_image(self.imdata)
            galsize=60/(abs(self.imheader['CD1_1'])*3600)
            buildweb.plot_vf_gals(imx,imy,keepflag,self.cat,ax,galsize=galsize)
            ax.set_xlabel('RA (deg)',fontsize=16)
            ax.set_ylabel('DEC (deg)',fontsize=16)        
            plt.savefig(self.coadd_png)

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

    def get_psf_allstars(self):
        ''' display psf image mosaic of 100 stars '''
        # check that png file exists
        # display png file
        imagebase = os.path.basename(self.imagename).replace('.fits','').replace('-shifted','')
        #print(imagebase)
        #print('plotdir = ',self.plotdir)
        zpsurf = os.path.join(self.psfdir,'plots',imagebase+"-allstars.png")
        #print('allstars source = ',zpsurf)
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
            print("Looking for : ",os.path.join(self.zpdir,imagebase+"*imsurfit-2*.png"))
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
        corrections = myfilter.get_trans_correction(redshift)
        self.corrections = corrections
        self.gals_filter_png = os.path.join(self.plotdir,'galaxies_in_filter.png')
        os.rename('galaxies_in_filter.png',self.gals_filter_png)
        pass
    
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

        self.pointing_name = os.path.basename(self.rimage).replace('R.fits','').replace('.fits','').replace('-r','').replace('-R','').replace('-shifted','')
        
        self.build_psf_names()
        self.get_rband_image()
        self.get_halpha_image()
        self.get_cs_image()
        self.get_gal_cutouts()
        self.get_params_for_gal_table()
        self.write_gal_table()
        self.plot_pointing_position()
        self.get_filter_ratio_plot()

    def build_psf_names(self):
        self.rpsf_image = self.psfdir+'/'+os.path.basename(self.rimage).replace('-shifted','').split('.fits')[0]+'-psf.fits'
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
        
        if os.path.exists(rimage):
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
            self.cs.generate_plots()
            self.cscoadd_flag=True
            self.get_gal_cutouts()
        else:
            self.cscoadd_flag=False
            print('could not find CS ha image : ',self.csimage)
        

    def get_gal_cutouts(self,size=90):
        """ get cutouts galaxies in FOV """
        self.gal_cutout_images = []

        # loop over galaxies in FOV
        gindex=np.arange(len(self.r.galfov_imx))
        galnames = Table(self.r.cat)['prefix'][self.r.keepflag]

        #galsizes = size#self.rcat['radius']/.4*2
        if self.rimage.find('INT'):
            pixscale = 0.33
        elif self.rimage.find('HDI'):
            pixscale = 0.425
        elif self.rimage.find('BOK'):
            pixscale = 0.45
        imdata,imheader = fits.getdata(self.csimage,header=True)
        size = (size/pixscale,size/pixscale)
        rowchange = np.arange(4,50,4)
        nrow = 1
        for i in range(len(rowchange)):
            if len(gindex) > rowchange[i]:
                nrow += 1
            else:
                break
            
        ncol = 4
        figsize = (12,3*nrow)            
        plt.figure(figsize=figsize)
        plt.subplots_adjust(top=.95,right=.95,left=.05,bottom=.05)        
        for j in gindex:
            
            #galsize = galsizes[j]
            position = (self.r.galfov_imx[j],self.r.galfov_imy[j])
            # get a cutout
            #print("size = ",size)
            #print("position = ",position)            
            cutout = Cutout2D(imdata,position,size)
            # display cutout

            ax = plt.subplot(nrow,ncol,j+1)

            display_image(cutout.data,csimage=True)
            plt.title(galnames[j][:20])# cutting names to avoid ridiculously long NED names
            #ax.set_xlabel('RA (deg)',fontsize=16)
            #ax.set_ylabel('DEC (deg)',fontsize=16)

            # save cutout image            
        imname = f"{self.pointing_name}-gal-cutouts.png"
        outfile = os.path.join(self.outdir,imname)
        plt.savefig(outfile)
            
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
        print("looking for filter ratio plot ",filename)
        if os.path.exists(filename):
            s = f'cp {filename} {os.path.join(self.outdir,os.path.basename(filename))}'
            os.system(s)
            self.filter_ratio_plot = os.path.basename(filename)
        else:
            self.filter_ratio_plot = None
        


    def get_params_for_gal_table(self):
        ''' setup arrays for galaxy table '''
        self.groupgals = self.r.cat[self.r.keepflag]
        self.vffil = vffil[self.r.keepflag]

    def write_gal_table(self):
        """ write out a list of with galaxies in FOV """
        # this will be good for checking which galaxies are observed in Halpha
        # in a way that is independent of running the halpha gui
        # just need the VFIDs

        # use the root file for the output name
        # keep the coadd name
        
        #outfile =
        outfile = os.path.join(self.outdir,self.pointing_name+'-galsFOV.csv')
        gals = self.r.cat['VFID'][self.r.keepflag]

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
        try:
            pscale_deg = self.r.imheader['PIXSCAL1']/3600
        except KeyError:
            pscale_deg =self.r.imheader['CD1_1']
        boxsizex=self.r.imheader['NAXIS1']*np.abs(float(pscale_deg))
        boxsizey=self.r.imheader['NAXIS2']*np.abs(float(pscale_deg))
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
        self.write_gal_table()
        self.write_pointing_location()        
        self.write_rband_header()
        self.write_rband_table()
        self.write_rband_psf()
        if self.pointing.r.zp_flag:
            self.write_rband_zp()

                 
        #self.write_rband_div()
        self.write_ha_header()
        self.write_ha_table()
        self.write_ha_psf()
        if self.pointing.ha.zp_flag:
            self.write_ha_zp()

        self.write_cs_header()
        self.write_cs_table()
        self.write_cutouts_table()
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
        images=[os.path.basename(self.pointing.r.coadd_png),\
                os.path.basename(self.pointing.r.psf_allstars_png),\
                os.path.basename(self.pointing.r.psf_png)]
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
        images = [os.path.basename(self.pointing.ha.coadd_png),\
                  os.path.basename(self.pointing.ha.psf_allstars_png),\
                  os.path.basename(self.pointing.ha.psf_png)]
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

        labels=['Continuum-Sub Image','Filter Ratio']#,'Residual<br> Surface']
        images = [os.path.basename(self.pointing.cs.coadd_png),\
                  os.path.basename(self.pointing.filter_ratio_plot)]
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
def buildone(rimages,i,coadd_dir,psfdir,zpdir,fratiodir):
    rimage = rimages[i]
    """ code to build webpage for one coadd """
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

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')

    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--laptop',dest = 'laptop', help='set if working on laptop')
    
    parser.add_argument('--coaddir',dest = 'coaddir', help='set to coadd directory')
    parser.add_argument('--psfdir',dest = 'psfdir', help='set to coadd directory')

  
     
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
        psfdir = outpathbase+'psf-images/'
        outdir = outpathbase+'/html_dev/coadds/'

        # get list of r-band coadded images
        a = glob.glob(coadd_dir+'VF*INT*-r-shifted.fits')
        b = glob.glob(coadd_dir+'VF*HDI*-r.fits')
        c = glob.glob(coadd_dir+'VF*HDI*-R.fits')
        d = glob.glob(coadd_dir+'VF*BOK*-r.fits')         
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

    # TODO - convert this to multiprocessing!!!
    
    indices = np.arange(len(rfiles))
    image_pool = mp.Pool(mp.cpu_count())
    myresults = [image_pool.apply_async(buildone,args=(rimages,i,coadd_dir,psfdir,zpdir,fratiodir),callback=collect_results) for i in indices[0:3]]
    
    image_pool.close()
    image_pool.join()
    image_results = [r.get() for r in myresults]


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
