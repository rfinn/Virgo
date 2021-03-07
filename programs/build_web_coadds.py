#!/usr/bin/env python

'''
GOAL:
* create web page to inspect the coadds, zp calibration, and psf



NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)




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


###########################################################
####  FUNCTIONS
###########################################################

def display_image(image,percent=99.5,lowrange=False,mask=None,sigclip=False):
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
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    

def get_galaxies_fov(imagename,RA,DEC):
    '''
    get the galaxies in the FOV  

    INPUT:
    * imagename
    * RA - array of RA values, like v.main['RA']
    * DEC - array of DEC values, like v.main['DEC']

    RETURNS:
    * imx, imy - RA and DEC, translated to image pixel coordinates
    * keepflag - boolean array, same length as RA and DEC
    * true for galaxies that fall within FOV of image

    '''
    # get image header
    header = fits.getheader(imagename)
    # define instance of wcs
    imwcs = wcs.WCS(header)

    # convert RA and DEC from catalog to pixels coordinates
    imx,imy = imwcs.wcs_world2pix(RA,DEC,1)

    # get image dimensions
    xmax = header['NAXIS1']
    ymax = header['NAXIS2']    

    # keep galaxies that fall within image boundary
    keepflag = (imx > 1) & (imx < xmax) & (imy > 1) & (imy < ymax)

    return imx, imy, keepflag

def plot_vf_gals(imx,imy,keepflag,cat,ax,galsize=120):
    ''' plot galaxies   '''
    gindex=np.arange(len(imx))[keepflag]
    

    for j in gindex:
        rect= plt.Rectangle((imx[j]-galsize/2.,imy[j]-galsize/2.), galsize, galsize,fill=False, color='b',lw=2)
        ax.add_artist(rect)
        s='{}\n vr={:.0f}'.format(cat['VFID'][j],cat['vr'][j])
        plt.text(imx[j],imy[j]+galsize/2.,s,fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='bottom',bbox=dict(facecolor='c',alpha=.3))
        plt.text(imx[j],imy[j]-galsize/2.,cat['NEDname'][j],fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='top',bbox=dict(facecolor='c',alpha=.3))
        #if cat['COflag'][j]:
        #    size=galsize+10
        #    rect= plt.Rectangle((imx[j]-size/2.,imy[j]-size/2.), size, size,fill=False, color='g',lw=1.5)
        #    ax.add_artist(rect)
        if cat['COflag'][j]:
            size=galsize
            #rect= plt.Circle((ran-size/2.,decn-size/2.), size,fill=False, color='g')
            rect= plt.Circle((imx[j],imy[j]), size/2,fill=False, color='g',lw=1.5)
            ax.add_artist(rect)
        if cat['A100flag'][j]:
            size=galsize+20
            rect= plt.Rectangle((imx[j]-size/2.,imy[j]-size/2.), size, size,fill=False, color='c',lw=1)
            ax.add_artist(rect)
        #if cat['HAobsflag'][j]:
        #    size=galsize+2*.005
        #    rect= plt.Rectangle((imx[i]-size/2.,imy[j]-size/2.), size, size,fill=False, color='r')
        #    ax.add_artist(rect)
    pass

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

    def __init__(self,imagename,psfimage=None,plotdir=None,cat=None,zpdir=None):
        self.imagename = imagename
        if psfimage is not None:
            self.psf_flag = True
            self.psf_image = psfimage
        else:
            self.psf_flag = False

        if plotdir is None:
            self.plotdir = os.getcwd()
        else:
            self.plotdir = plotdir
        if cat is None:
            self.cat = fits.getdata(VFMAIN_PATH)

        self.zpdir = zpdir
    def generate_plots(self):
        self.get_image()
        self.make_coadd_png()
        if self.psf_flag:
            self.get_psf_image()
            if self.found_psf:
                self.make_psf_png()
        self.get_zpimage_firstpass()
        self.get_zp_magcomparison()        
    def get_image(self):
        '''  read in image, save data and header '''
        self.imdata,self.imheader = fits.getdata(self.imagename,header=True)
        self.wcs = wcs.WCS(self.imheader)
        self.xdim,self.ydim = self.imdata.shape
        self.racenter,self.deccenter = self.wcs.wcs_pix2world(self.xdim/2,self.ydim/2,1)
        self.zp = self.imheader['PHOTZP']
        self.pscale = np.abs(self.imheader['CD1_1'])*3600
    
    def make_coadd_png(self):
        ''' display image, and mark position of galaxies '''
        self.coadd_png = self.plotdir+'coadd.png'
        imx,imy,keepflag = get_galaxies_fov(self.imagename,self.cat['RA'],self.cat['DEC'])
        self.keepflag = keepflag        
        if os.path.exists(self.coadd_png):
            print('Found {}.  not remaking this.'.format(self.coadd_png))
        else:
            plt.figure(figsize=(6,6))
            ax = plt.subplot(projection=wcs.WCS(self.imheader))
            plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
            display_image(self.imdata)
            galsize=60/(abs(self.imheader['CD1_1'])*3600)
            plot_vf_gals(imx,imy,keepflag,self.cat,ax,galsize=galsize)
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
        pass
    
    def make_psf_png(self):
        ''' read in psf image, save data and header '''
        # check that png file exists
        # display png file
        norm = simple_norm(self.psfdata, 'log', percent=99.)
        plt.figure(figsize=(6,6))
        plt.subplots_adjust(right=.9,top=.95,left=.1,bottom=.05)
        plt.imshow(self.psfdata, norm=norm, origin='lower', cmap='viridis')
        plt.colorbar(fraction=.046,pad=.04)
        #plt.show()
        self.psf_png = self.plotdir+'psf.png'
        plt.savefig(self.psf_png)        

    def display_psf_allstars(self):
        ''' display psf image mosaic of 100 stars '''
        # check that png file exists
        # display png file
        
        pass

    def get_zpimage_firstpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        imagebase = self.imagename.replace('-noback-coadd.fits','')
        zpsurf = os.path.join(self.zpdir,imagebase+"-getzp-xyresidual-fitted.png")
        self.zpsurf_png = os.path.join(self.plotdir,imagebase+"-getzp-xyresidual-fitted.png")
        os.copy(zpsurf,self.zpsurf_png)
        pass
    
    def get_zpimage_finalpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file        
        pass

    def get_zp_magcomparison(self):
        ''' get the final plot of inst mag vs panstarrs mag'''
        # check that png file exists
        # display png file
        imagebase = self.imagename.replace('-noback-coadd.fits','')        
        pancomp = os.path.join(self.zpdir,imagebase+"-se-panflux.png")
        self.pancomp_png = os.path.join(self.plotdir,imagebase+"-se-panflux.png")
        os.copy(pancomp,self.pancomp_png)
        pass

    
class pointing():

    def __init__(self,rimage=None,haimage=None,psfdir=None,zpdir=None,outdir=None):
        '''
        INPUT:
        * rband image
        * halpha image
        * psf directory
        * directory with ZP images
        '''

        self.rimage = rimage
        self.haimage = haimage
        self.psfdir = psfdir
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
        self.get_rband_image()
        self.get_halpha_image()
        self.get_params_for_gal_table()
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
            outprefix = self.outdir+'/r-'
            if self.rpsf_flag:
                self.r = coadd_image(self.rimage,psfimage=self.rpsf_image,plotdir=outprefix)
            else:
                self.r = coadd_image(self.rimage,psfimage=None,plotdir=outprefix)
            self.rcoadd_flag=True
            self.r.generate_plots()
        else:
            self.rcoadd_flag=False
        
    
    def get_halpha_image(self):
        ''' initiate an instance of coadd image class '''
        if os.path.exists(self.haimage):
            outprefix = self.outdir+'/ha-'
            if self.hapsf_flag:
                self.ha = coadd_image(self.haimage,psfimage=self.rpsf_image,plotdir=outprefix)
            else:
                self.ha = coadd_image(self.haimage,psfimage=None,plotdir=outprefix)
            self.hacoadd_flag=True
            self.ha.generate_plots()
        else:
            self.hacoadd_flag=False
            print('could not find ha image : ',self.haimage)
        
    

    def get_filter_ratio_plot(self):
        ''' display the filter ratio png '''
        pass


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
        pngfile = os.path.basename(self.pointing.r.coadd_png)
        psfpng = os.path.basename(self.pointing.r.psf_png)
        images = [pngfile,psfpng,self.pointing.r.zpsurf_png,self.pointing.r.pancomp_png]
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
    virgovms=True
    intonlaptop=False
    if virgovms:
        vmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')
        homedir = '/mnt/qnap_home/rfinn/'
        
        coadd_dir = homedir+'/Halpha/reduced/virgo-coadds-HDI/'
        zpdir = homedir+'/Halpha/reduced/virgo-coadds-HDI/plots/'        
        psfdir = homedir+'/Halpha/reduced/psf-images/'
        outdir = homedir+'/research/Virgo/html-dev/coadds/'
        
    #rimage = coadd_dir+'VF-219.9485+5.3942-INT-20190530-p019-r-shifted.fits'    
    #haimage = coadd_dir+'VF-219.9523+5.4009-INT-20190530-p019-Halpha.fits'
    #hacsimage = coadd_dir+'VF-219.9485+5.3942-INT-20190530-p019-CS.fits'
    if intonlaptop:
        vmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')      
        coadd_dir = '/home/rfinn/data/reduced/virgo-coadds-HDI/'
        zpdir = coadd_dir+'/plots/'                
        psfdir = homedir+'/data/reduced/psf-images/'
        outdir = homedir+'/research/Virgo/html-dev/coadds/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # for HDI files
    rfiles = glob.glob('*-r-*coadd.fits')
    Rfiles = glob.glob('*-R-*coadd.fits')
    hfiles = glob.glob('*-ha4-*coadd.fits')        
    hfiles.sort()

    # combine r and R files
    allrfiles = rfiles+Rfiles
    allrfiles.sort()

    # loop through r filenames
    for rimage in allrfiles:
        print(rimage)
        # find matching ha4 coadd
        haimage = rimage.replace('-r-','-ha4-').replace('-R-','-ha4-')
        print('looking for ',haimage)
        if not os.path.exists(haimage):
            print("couldn't find it")
            # haimage could have a different date
            search_string = rimage.split('HDI')[1].replace('-r-','-ha4-').replace('-R-','-ha4')
            print('looking for ',search_string)
            try:
                haimage = glob.glob('*'+search_string)[0]
            except:
                if not os.path.exists(haimage):
                    print('WARNING: could not find halpha image for ',rimage,' Skipping for now.')
                    print('\t Looking for ',haimage)
                    break
                    continue
        pname = os.path.basename(rimage).replace('-shifted','').replace('.fits','').replace('-r','').replace('-R','')
        outdir = os.path.join(outdir,pname)
        p = pointing(rimage=rimage,haimage=haimage,psfdir=psfdir,zpdir=zpdir,outdir=outdir)
        h = build_html_pointing(p,outdir=outdir)

        break
