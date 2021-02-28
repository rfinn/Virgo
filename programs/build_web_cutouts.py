#!/usr/bin/env python

'''
GOAL:
* create web page to inspect the cutouts

USAGE:
* run from cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)
* https://docs.astropy.org/en/stable/visualization/normalization.html#:~:text=The%20astropy.visualization%20module%20provides%20a%20framework%20for%20transforming,x%20represents%20the%20values%20in%20the%20original%20image%3A

TO DO:
* need to figure out how to handle repeated observations
  - don't overwrite directory

* fix how I combine unwise images when multiple images are returned
* same for galex
'''

import os
import numpy as np
import glob
from scipy.stats import scoreatpercentile
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch, PercentileInterval
from astropy.visualization import ImageNormalize, LinearStretch
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip

homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')
import rungalfit as rg #This code has galfit defined functions 


###########################################################
####  GLOBAL VARIABLES
###########################################################


VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
haimaging_path = os.path.join(homedir,'github/HalphaImaging/python3/')
#sys.path.append(haimaging_path)

residual_stretch = LinearStretch(slope=0.5, intercept=0.5) + SinhStretch() + \
    LinearStretch(slope=2, intercept=-1)
###########################################################
####  FUNCTIONS
###########################################################

def display_image(image,percentile1=.5,percentile2=99.5,stretch='asinh',mask=None,sigclip=True):
    lowrange=False
    # use inner 80% of image
    xdim,ydim = image.shape
    xmin = int(.1*xdim)
    xmax = int(.9*xdim)    
    ymin = int(.1*ydim)
    ymax = int(.9*ydim)
    if mask is not None:
        imdata = np.ma.array(image,mask=mask)
        
    else:
        imdata = image
    v1 = scoreatpercentile(imdata,percentile1)    
    v2 = scoreatpercentile(imdata,percentile2)

    if sigclip:
        clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5)#,grow=3)
        mask = mask[xmin:xmax,ymin:ymax]
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]
        mask = mask[xmin:xmax,ymin:ymax]        
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',max_percent=percentile2)
    else:
        if mask is not None:
            norm = simple_norm(clipped_data[mask], stretch='asinh',max_percent=percentile2)
        else:
            norm = simple_norm(clipped_data, stretch='asinh',max_percent=percentile2)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower',vmin=v1,vmax=v2)
    
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

def write_table(html,images,labels):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    self.html.write('<tr>')
    for i in images:
        self.html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
    self.html.write('</tr>\n')            
    self.html.write('</table>\n')

def write_text_table(html,labels,data):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    self.html.write('<tr>')
    for d in data:
        self.html.write('<td></td>'.format(d))
    self.html.write('</tr>\n')            
    self.html.write('</table>\n')


def make_png(fitsimage,outname):
    imdata,imheader = fits.getdata(fitsimage,header=True)
    plt.figure(figsize=(6,6))
    ax = plt.subplot(projection=wcs.WCS(imheader))
    plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
    display_image(imdata,sigclip=True)
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('DEC (deg)',fontsize=16)        
    plt.savefig(outname)        


def display_galfit_model(galfile,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis',zoom=None,outdir=None):
      '''
      ARGS:
      galfile = galfit output image (with image, model, residual)
      percentile1 = min percentile for stretch of image and model
      percentile2 = max percentile for stretch of image and model
      p1residual = min percentile for stretch of residual
      p2residual = max percentile for stretch of residual
      cmap = colormap, default is viridis
      '''
      # model name


      image,h = fits.getdata(galfile,1,header=True)
      model = fits.getdata(galfile,2)
      residual = fits.getdata(galfile,3)

      if zoom is not None:
         print("who's zoomin' who?")
         # display central region of image

         # get image dimensions and center
         xmax,ymax = image.shape
         xcenter = int(xmax/2)
         ycenter = int(ymax/2)

         # calculate new size to display based on zoom factor
         new_xradius = int(xmax/2/(float(zoom)))
         new_yradius = int(ymax/2/(float(zoom)))

         # calculate pixels to keep based on zoom factor
         x1 = xcenter - new_xradius
         x2 = xcenter + new_xradius
         y1 = ycenter - new_yradius
         y2 = ycenter + new_yradius
         
         # check to make sure limits are not outsize image dimensions
         if (x1 < 1):
            x1 = 1
         if (y1 < 1):
            y1 = 1
         if (x2 > xmax):
            x2 = xmax
         if (y2 > ymax):
            y2 = ymax

         # cut images to new size
         image = image[x1:x2,y1:y2]
         model = model[x1:x2,y1:y2]
         residual = residual[x1:x2,y1:y2]         
         pass
      wcs = wcs.WCS(h)
      images = [image,model,residual]
      titles = ['image','model','residual']
      v1 = [scoreatpercentile(image,percentile1),
            scoreatpercentile(image,percentile1),
            scoreatpercentile(residual,p1residual)]
      v2 = [scoreatpercentile(image,percentile2),
            scoreatpercentile(image,percentile2),
            scoreatpercentile(residual,p2residual)]
      norms = [simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(image,'asinh',max_percent=percentile2),
               simple_norm(residual,'linear',max_percent=p2residual)]
      
      outim = ['galfit_image.png','galfit_model.png','galfit_residual.png']
      if outdir is not None:
          outim = [os.path.join(outdir,f) for f in outim]
      plt.figure(figsize=(14,6))
      plt.subplots_adjust(wspace=.0)
      for i,im in enumerate(images):
          plt.figure(figsize=(6,6))          
          plt.subplot(1,1,1,projection=wcs)
          plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
          plt.imshow(im,origin='lower',cmap=cmap,vmin=v1[i],vmax=v2[i],norm=norms[i])
          plt.xlabel('RA')
          plt.ylabel('DEC')
          plt.title(titles[i],fontsize=16)
          plt.savefig(outim[i])

###########################################################
####  CLASSES
###########################################################

    
    
class cutout_dir():

    def __init__(self,cutoutdir=None,outdir=None):
        '''
        INPUT:
        * directory containing cutouts
        * output directory for png images

        This creates the png images for different cutouts
        '''
        self.gname = os.path.basename(cutoutdir)
        

        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir

        self.cutoutdir = cutoutdir
    def runall(self):
        self.get_halpha_names()
        self.get_legacy_names()
        self.get_wise_names()
        self.get_galex_names()
        self.make_png_plots()
    def get_halpha_names(self):
        search_string = os.path.join(self.cutoutdir,self.gname+'*-R.fits')
        print(search_string)
        t = glob.glob(search_string)
        print(t)
        
        self.rimage = t[0]
        self.haimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-Ha.fits'))[0]
        self.csimage = glob.glob(os.path.join(self.cutoutdir,self.gname+'*-CS.fits'))[0]
        self.maskimage = self.rimage.replace('.fits','-mask.fits')
    def get_legacy_names(self):
        ''' get names of legacy images  '''
        legdir = os.path.join(self.cutoutdir,'legacy')
        self.legacy_g = glob.glob(legdir+'/*-g.fits')[0]
        self.legacy_r = glob.glob(legdir+'/*-r.fits')[0]
        self.legacy_z = glob.glob(legdir+'/*-z.fits')[0]        
    def get_wise_names(self):
        ''' get names of wise images  '''
        # need to fix this to check for combined unwise images
        wisedir = os.path.join(self.cutoutdir,'unwise')
        self.w1 = glob.glob(wisedir+'/*-w1-img-m.fits')[0]
        self.w2 = glob.glob(wisedir+'/*-w2-img-m.fits')[0]
        self.w3 = glob.glob(wisedir+'/*-w3-img-m.fits')[0]
        self.w4 = glob.glob(wisedir+'/*-w4-img-m.fits')[0]                        
    def get_galex_names(self):
        ''' get names of galex images  '''        
        galdir = os.path.join(self.cutoutdir,'galex')
        galexfiles = glob.glob(galdir+'/*nuv*.fits')
        if len(galexfiles) > 0:
            for f in galexfiles:
                if f.find('nuv') > -1:
                    self.nuv = f
                    self.nuv_flag = True
        else:
            self.nuv_flag = False

    def define_png_names(self):
        pass
    def make_png_plots(self):
        self.fitsimages = [self.rimage,self.haimage,self.csimage,\
                      self.legacy_g,self.legacy_r,self.legacy_z,\
                      self.w1,self.w2,self.w3,self.w4]
        if self.nuv_flag:
            self.fitsimages.append(self.nuv)
        self.pngimages=[]
        for f in self.fitsimages:
            pngfile = os.path.join(self.outdir,os.path.basename(f).replace('.fits','.png'))
            make_png(f,pngfile)
            self.pngimages.append(pngfile)          
    def make_png_cs(self):
        maskdata = fits.getdata(self.maskimage)
        boolmask = maskdata < 1
        pngfile = os.path.join(self.outdir,os.path.basename(self.csimage).replace('.fits','.png'))
        csdata = fits.getdata(self.csimage)
        display_image(csdata,mask=boolmask)
    def get_params_for_gal_table(self):
        ''' setup arrays for galaxy table '''
        self.groupgals = self.r.cat[self.r.keepflag]

    def get_galfit_model(self):
        ''' read in galfit model and make png '''
        self.galfit = self.rimage.replace('.fits','-1Comp-galfit-out.fits')

        # store fit results
        display_galfit_model(self.galfit)

        outim = ['galfit_image.png','galfit_model.png','galfit_residual.png']
        
        self.galimage = os.path.join(self.outdir,outim[0])
        self.galmodel = os.path.join(self.outdir,outim[1])
        self.galresidual = os.path.join(self.outdir,outim[2])        
        # display r-band images with ellipse

        # model

        # residual

        # store fitted parameters

        t = rg.parse_galfit_1comp(self.galfit)
        
        #header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        self.xc, self.xc_err = t[0]
        self.yc, self.yc_err = t[1]
        self.mag, self.mag_err = t[2]
        self.re, self.re_err = t[3]
        self.nsersic, self.nsersic_err = t[4]
        self.BA, self.BA_err = t[5]
        self.PA, self.PA_err = t[6]
        self.sky, self.sky_err = t[7]
        self.error = t[8]
        self.chi2nu = t[9]
        
        pass
    

    def get_phot_tables(self):
        ''' read in phot tables and make plot of flux and sb vs sma '''
        cs_galfit_phot = self.csimage.replace('.fits','-GAL-phot.fits')
        cs_gphot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','-GAL-phot.fits')
        r_gphot = fits.getdata(r_galfit_phot)

        # plot enclosed flux
        plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot]
        labels = ['r','Halphax100']
        for i,t in enumerate(tabs):
            y1 = t['flux_erg']+t['flux_erg_err']
            y2 = t['flux_erg']-t['flux_erg_err']
            if i == 1:
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend()        
        self.fluxsma_png = self.gname+'-enclosed-flux.png'
        plt.savefig(self.encflux_png)
        
        # plot mag sma
        plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot]
        labels = ['r','Halpha']
        for i,t in enumerate(tabs):
            y1 = t['mag']+t['mag_err']
            y2 = t['mag']-t['mag_err']
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('magnitude (AB)',fontsize=16)
        plt.gca().set_xscale('log')        
        self.magsma_png = self.gname+'-mag-sma.png'
        plt.legend()        
        plt.savefig(self.magsma_png)
        
        # plot sb erg vs sma
        plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot]
        labels = ['r','Halphax100']
        for i,t in enumerate(tabs):
            y1 = t['sb_erg_sqarcsec']+t['sb_erg_sqarcsec_err']
            y2 = t['sb_erg_sqarcsec']-t['sb_erg_sqarcsec_err']
            if i == 1:
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('SB (erg/s/cm^2/Hz/arcsec^2)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend()
        self.sbsma_png = self.gname+'-sb-sma.png'
        plt.savefig(self.sbsma_png)
                             
        # plot mag sma
        plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot]
        labels = ['r','Halpha']
        for i,t in enumerate(tabs):
            y1 = t['sb_mag_sqarcsec']+t['sb_mag_sqarcsec_err']
            y2 = t['sb_mag_sqarcsec']-t['sb_mag_sqarcsec_err']
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i])
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.legend()        
        self.magsma_png = self.gname+'-sbmag-sma.png'
        plt.savefig(self.magsma_png)
        
    
class build_html_cutout():

    def __init__(self,cutoutdir,outdir):
        ''' pass in instance of cutout_dir class and output directory '''
        self.cutout = cutoutdir
        outfile = self.cutout.gname+'.html'
        outfile = os.path.join(outdir,outfile)
        
        self.html = open(outfile,'w')

        # for reference, this is the order of the png images
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\
        #              self.legacy_g,self.legacy_r,self.legacy_z,\
        #              self.w1,self.w2,self.w3,self.w4]

        
        self.build_html()
    def build_html(self):
        self.write_header()
        self.write_sfr_images()
        self.write_legacy_images()
        self.write_wise_images()
        self.write_halpha_images()
        self.write_galfit_results()
        self.write_phot_enclosed()
        self.write_phot_sb()
        self.write_footer()
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
        self.html.write('<h1>{}</h1>\n'.format(self.cutoutdir.gname))

        #self.html.write('<a href="../../{}">Home</a>\n'.format(htmlhome))
        #self.html.write('<br />\n')
        #self.html.write('<a href="../../{}">Next ({})</a>\n'.format(nexthtmlgalaxydir1, nextgalaxy[ii]))
        #self.html.write('<br />\n')
        #self.html.write('<a href="../../{}">Previous ({})</a>\n'.format(prevhtmlgalaxydir1, prevgalaxy[ii]))

    def write_sfr_images(self):
        ''' legacy jpeg, galex nuv, halpha, w4 '''
        self.html.write('<h2>SF Indicators</h2>\n')
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\ # 0,1,2
        #                   self.legacy_g,self.legacy_r,self.legacy_z,\ # 3,4,5
        #                   self.w1,self.w2,self.w3,self.w4] # 6,7,8,9

        if self.cutout.nuv_flag:
            images = [self.cutout.legacy_jpg,self.cutout.pngimages[-1],\
                      self.cutout.pngimages[2],self.cutout.pngimages[9]]
            labels = ['Legacy','NUV','Halpha','W4']
        else:
            images = [self.cutout.legacy_jpg,\
                      self.cutout.pngimages[2],self.cutout.pngimages[9]]
            labels = ['Legacy','Halpha','W4']
        write_table(self.html,images,labels)
        pass

    
    def write_halpha_images(self):
        '''  r, halpha, cs, and mask images '''
        images = [self.pngimages[0:3]]
        images.append(self.cutout.mask_png)
        labels = ['R','Halpha+C','CS','Mask']
        write_table(self.html,images,labels)
        pass

    def write_legacy_images(self):
        ''' jpg, g,r,z legacy images '''
        images = [self.pngimages[3:6]]
        labels = ['g','r','z']
        write_table(self.html,images,labels)

    def write_wise_images(self):
        ''' w1 - w4 images '''
        images = [self.pngimages[6:10]]
        labels = ['W1','W2','W3','W4']
        write_table(self.html,images,labels)
    
    def write_galex_images(self):
        ''' right now just nuv '''
        pass

    def write_galfit_images(self):
        ''' display galfit model and fit parameters for r-band image '''
        images = [self.cutout.galimage,self.cutout.galmodel,self.cutout.galresidual,\
                  self.cutout.mask]
        labels = ['Images', 'Model', 'Residual','Mask']
        write_table(self.html,images,labels)
        pass
    
    def write_galfit_table(self):
        ''' display galfit model and fit parameters for r-band image '''
        labels=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        data = ['{:.1f}+/-{:.1f}'.format(self.cutout.xc, self.cutout.xc_err),
                '{:.1f}+/-{:.1f}'.format(self.cutout.yc, self.yc_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.mag, self.cutout.mag_err), 
                '{:.2f}+/-{:.2f}'.format(self.cutout.re, self.cutout.re_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.nsersic, self.cutout.nsersic_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.BA, self.cutout.BA_err),
                '{:.1f}+/-{:.1f}'.format(self.cutout.PA, self.cutout.PA_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.sky, self.cutout.sky_err),
                '{}'.format(self.cutout.error),
                '{:.2f}'.format(self.cutout.chi2nu)]
        

        write_text_table(self.html,labels,data)
        pass

    def write_phot_enclosed(self):
        ''' photometry table with galfit and photutil results '''

        images = [self.cutout.enc_flux_png,self.cutout.magsma_png]
        labels = ['Enclosed Flux','Enclosed Magnitude']
        write_table(self.html,images,labels)
    
    def write_phot_sb(self):
        ''' photometry table with galfit and photutil results '''
        images = [self.cutout.sbsma_png,self.cutout.sbmagsma_png]
        labels = ['Surface Brightness','Surface Brightness']
        write_table(self.html,images,labels)

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
    # work through coadd directory
    vmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')
    cutoutdir = '/home/rfinn/research/Virgo/gui-output-2019-june/cutouts/VFID0481-NGC6307/'

    outdir = homedir+'/research/Virgo/html-dev/cutouts/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    cutout_name = os.path.basename(cutoutdir)
    outdir = os.path.join(outdir,cutout_name)
    p = cutout_dir(cutoutdir=cutoutdir,outdir=outdir)
    #h = build_html_pointing(p,outdir=outdir)
    pass
