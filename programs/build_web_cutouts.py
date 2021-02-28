#!/usr/bin/env python

'''
GOAL:
* create web page to inspect the cutouts
* run from cutouts directory


NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)

TO DO:
* need to figure out how to handle repeated observations
  - don't overwrite directory

* fix how I combine unwise images when multiple images are returned
'''

import os
import numpy as np
import glob

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord

from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip



homedir = os.getenv("HOME")
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
haimaging_path = os.path.join(homedir,'github/HalphaImaging/python3/')
#sys.path.append(haimaging_path)

###########################################################
####  FUNCTIONS
###########################################################

def display_image(image,percent=99.5,lowrange=False,mask=None,sigclip=True):
    lowrange=False
    # use inner 80% of image
    xdim,ydim = image.shape
    xmin = int(.1*xdim)
    xmax = int(.9*xdim)    
    ymin = int(.1*ydim)
    ymax = int(.9*ydim)
    if sigclip:
        clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5)#,grow=3)
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        if mask is not None:
            norm = simple_norm(clipped_data[mask], stretch='asinh',percent=percent)
        else:
            norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    
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

def make_png(fitsimage,outname):
    imdata,imheader = fits.getdata(fitsimage,header=True)
    plt.figure(figsize=(6,6))
    ax = plt.subplot(projection=wcs.WCS(imheader))
    plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
    display_image(imdata,sigclip=True)
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('DEC (deg)',fontsize=16)        
    plt.savefig(outname)        

###########################################################
####  CLASSES
###########################################################

    
    
class cutout_dir():

    def __init__(self,cutoutdir=None,outdir=None):
        '''
        INPUT:
        * directory containing cutouts
        * output directory for png images
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
        legdir = os.path.join(self.cutoutdir,'legacy')
        self.legacy_g = glob.glob(legdir+'/*-g.fits')[0]
        self.legacy_r = glob.glob(legdir+'/*-r.fits')[0]
        self.legacy_z = glob.glob(legdir+'/*-z.fits')[0]        
    def get_wise_names(self):
        # need to fix this to check for combined unwise images
        wisedir = os.path.join(self.cutoutdir,'unwise')
        self.w1 = glob.glob(wisedir+'/*-w1-img-m.fits')[0]
        self.w2 = glob.glob(wisedir+'/*-w2-img-m.fits')[0]
        self.w3 = glob.glob(wisedir+'/*-w3-img-m.fits')[0]
        self.w4 = glob.glob(wisedir+'/*-w4-img-m.fits')[0]                        
    def get_galex_names(self):
        galdir = os.path.join(self.cutoutdir,'galex')
        galexfiles = glob.glob(galdir+'/*nuv*.fits')
        if len(galexfiles) > 0:
            for f in galexfiles:
                if f.find('nuv') > -1:
                    self.nuv = f
                    self.nuv_flag = True
        else:
            self.nuv_flag = False
            
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
        display_image(csdata,boolmask)
    def get_params_for_gal_table(self):
        ''' setup arrays for galaxy table '''
        self.groupgals = self.r.cat[self.r.keepflag]
        
class build_html_pointing():

    def __init__(self,pointing,outdir):
        self.pointing = pointing
        outfile = self.pointing.pointing_name+'.html'
        outfile = os.path.join(outdir,outfile)
        self.html = open(outfile,'w')
        self.build_html()
    def build_html(self):
        self.write_header()
        self.write_gal_table()
        self.write_rband_div()
        self.write_rband_table()
        self.write_ha_div()
        self.write_ha_table()        
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
        self.html.write('<h1>{}</h1>\n'.format(self.pointing.pointing_name))

        #self.html.write('<a href="../../{}">Home</a>\n'.format(htmlhome))
        #self.html.write('<br />\n')
        #self.html.write('<a href="../../{}">Next ({})</a>\n'.format(nexthtmlgalaxydir1, nextgalaxy[ii]))
        #self.html.write('<br />\n')
        #self.html.write('<a href="../../{}">Previous ({})</a>\n'.format(prevhtmlgalaxydir1, prevgalaxy[ii]))

    def write_gal_table(self):
       # Add the properties of each galaxy.
        self.html.write('<h3>Galaxies in FOV</h3>\n')
        self.html.write('<table>\n')
        self.html.write('<tr>\n')
        self.html.write('<th>VFID</th>\n')
        self.html.write('<th>Galaxy</th>\n')
        #self.html.write('<th>Morphology</th>\n')
        self.html.write('<th>RA</th>\n')
        self.html.write('<th>Dec</th>\n')
        self.html.write('<th>D(25)<br />(arcmin)</th>\n')
        self.html.write('<th>CO</th>\n')
        self.html.write('<th>A100</th>\n')
        self.html.write('</tr>\n')
        for g in self.pointing.groupgals:
            #if '031705' in gal['GALAXY']:
            #    print(groupgal['GALAXY'])
            self.html.write('<tr>\n')
            self.html.write('<td>{}</td>\n'.format(g['VFID']))
            self.html.write('<td>{}</td>\n'.format(g['NEDname']))
            #typ = groupgal['MORPHTYPE'].strip()
            #if typ == '' or typ == 'nan':
            #    typ = '...'
            #self.html.write('<td>{}</td>\n'.format(typ))
            self.html.write('<td>{:.7f}</td>\n'.format(g['RA']))
            self.html.write('<td>{:.7f}</td>\n'.format(g['DEC']))
            self.html.write('<td>{:.4f}</td>\n'.format(g['radius']*2/60.))
            self.html.write('<td>{:.4f}</td>\n'.format(g['COflag']))
            self.html.write('<td>{:.4f}</td>\n'.format(g['A100flag']))                        
            #if np.isnan(groupgal['PA']):
            #    pa = 0.0
            #else:
            #    pa = groupgal['PA']
            #self.html.write('<td>{:.2f}</td>\n'.format(pa))
            #self.html.write('<td>{:.3f}</td>\n'.format(1-groupgal['BA']))
            self.html.write('</tr>\n')
        self.html.write('</table>\n')        
    def write_rband_div(self):
        # write coadd
        # write psf
        # write table with filter, zp, fwhm
        self.html.write('<h2>r-band Coadd</h2>\n')




        self.html.write('<table width="90%">\n')
        #self.html.write('<tr><th>Coadd<br>{}</th> <th>PSF images</th> <th>ZP Calib Orig</th> <th>ZP Calib Final</th></tr></p>\n'.format(os.path.basename(self.pointing.rimage)))
        self.html.write('<tr><th>Coadd</th> <th>PSF images</th> <th>ZP Calib Orig</th> <th>ZP Calib Final</th></tr></p>\n')        
        pngfile = self.pointing.r.coadd_png
        psfpng = self.pointing.r.psf_png
        images = [pngfile,psfpng,psfpng,psfpng]
        self.html.write('<tr>')
        for i in images:
            self.html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
        self.html.write('</tr>\n')            
        self.html.write('</table>\n')

    def write_rband_table(self):
        # write coadd
        # write psf
        # write table with filter, zp, fwhm
        write_coadd_prop_table(self.html,self.pointing.r.imheader['FILTER'],self.pointing.r.zp,self.pointing.r.fwhm_arcsec)
        pass
    def write_ha_div(self):
        # write coadd
        # write psf
        # write table with filter, zp, fwhm
        self.html.write('<h2>Halpha Coadd</h2>\n')
        self.html.write('<table width="90%">\n')

        #self.html.write('<tr><th>Coadd<br>{}</th> <th>PSF images</th> <th>ZP Calib Orig</th> <th>ZP Calib Final</th></tr></p>\n'.format(os.path.basename(self.pointing.haimage)))
        self.html.write('<tr><th>Coadd</th> <th>PSF images</th> <th>ZP Calib Orig</th> <th>ZP Calib Final</th></tr></p>\n')
        
        pngfile = self.pointing.ha.coadd_png
        psfpng = self.pointing.ha.psf_png
        images = [pngfile,psfpng,psfpng,psfpng]
        self.html.write('<tr>')
        for i in images:
            self.html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
        self.html.write('</tr>\n')            
        self.html.write('</table>\n')
    
    def write_ha_table(self):
        # write coadd
        # write psf
        # write table with filter, zp, fwhm
        write_coadd_prop_table(self.html,self.pointing.ha.imheader['FILTER'],self.pointing.ha.zp,self.pointing.ha.fwhm_arcsec)
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
