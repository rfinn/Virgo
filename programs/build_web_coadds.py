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
import sys

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord

from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip
from astropy.time import Time

homedir = os.getenv("HOME")
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'

import build_web_common as buildweb

sys.path.append(homedir+'/github/halphagui/')
import filter_transmission as ft

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
        #temp = self.imagename.replace('-shifted','').replace('.fits','')
        #pointing = temp.split('-')[-2].replace('p','pointing')
        #self.intprefix = "{}*_{}".format(pointing,temp[-1])
        #print('INT plot prefix = ',self.intprefix)
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
            self.get_zpimage_firstpass()
            self.get_zp_magcomp_firstpass()        
            self.get_zpplot_secondpass()                
            self.get_zpimage_secondpass()
            self.get_zp_magcomp_secondpass()
            self.zp_flag = True
        except:
            print('WARNING: problem getting zp calibration images')
            self.zp_flag = False

        if self.filter == 'ha':
            self.get_gredshift_filter_curve()
    def get_image(self):
        '''  read in image, save data and header '''
        self.imdata,self.imheader = fits.getdata(self.imagename,header=True)
        self.wcs = wcs.WCS(self.imheader)
        self.xdim,self.ydim = self.imdata.shape
        self.racenter,self.deccenter = self.wcs.wcs_pix2world(self.xdim/2,self.ydim/2,1)
        self.zp = self.imheader['PHOTZP']
        try:
            self.pscale = np.abs(self.imheader['PIXSCAL1'])
        except KeyError:
            self.pscale = np.abs(self.imheader['CD1_1'])*3600
        try:
            self.exptime = self.imheader['ORIGEXPT']
        except KeyError:
            self.exptime = self.imheader['EXPTIME']    
        self.sefwhm_arcsec = self.imheader['SEFWHM']
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

        # need to cut to keep the galaxies within the right filter
        self.keepflag = keepflag        
        if os.path.exists(self.coadd_png):
            print('Found {}.  not remaking this.'.format(self.coadd_png))
        else:
            plt.figure(figsize=(8,8))
            ax = plt.subplot(projection=wcs.WCS(self.imheader))
            plt.subplots_adjust(top=.95,right=.95,left=.15,bottom=.1)
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



    def get_zpplot_firstpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        #print(self.imagename)

        #cluge b/c Becky's filenames threw a wrench in naming conventions
        
        #print(imagebase)
        #print('plotdir = ',self.plotdir)
        try:
            imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
            zpsurf = glob.glob(os.path.join(self.zpdir,imagebase+"*getzp-xyresidual-fitted.png"))[0]
        except IndexError: # the above won't work for INT data b/c my naming conventions are a mess
            # and because I ran the zp calibration from diff directory
            #print('search path = ',os.path.join(self.zpdir,self.intprefix+"*getzp-xyresidual-fitted.png"))
            zpsurf = glob.glob(os.path.join(self.zpdir,self.intprefix+"*getzp*fitted*.png"))[0]
        self.zpplot_png = os.path.join(self.plotdir,imagebase+"-getzp-xyresidual-fitted.png")
        os.system('cp '+zpsurf+' '+self.zpplot_png)
        pass
    def get_zpplot_secondpass(self):
        ''' get the zp image, first pass'''
        # check that png file exists
        # display png file
        #print(self.imagename)

        imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
        #print(imagebase)
        #print('plotdir = ',self.plotdir)
        zpplot2 = glob.glob(os.path.join(self.zpdir,"f"+imagebase+"*getzp-xyresidual-fitted.png"))[0]
        self.zpplot2_png = os.path.join(self.plotdir,"f"+imagebase+"-getzp-xyresidual-fitted.png")
        os.system('cp '+zpplot2+' '+self.zpplot2_png)
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
            
            zpsurf = glob.glob(os.path.join(self.zpdir,imagebase+"*imsurfit-2*.png"))[0]
        except IndexError:
            zpsurf = glob.glob(os.path.join(self.zpdir,self.intprefix+"*imsurfit-2*.png"))[0]
        self.zpsurf_png = os.path.join(self.plotdir,imagebase+"-imsurfit-2.png")
        os.system('cp '+zpsurf+' '+self.zpsurf_png)
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
        self.zpsurf2_png = os.path.join(self.plotdir,"f"+imagebase+"-imsurfit-2-round2.png")
        os.system('cp '+zpsurf+' '+self.zpsurf2_png)
        pass
        

    def get_zp_magcomp_firstpass(self):
        ''' get the final plot of inst mag vs panstarrs mag'''
        # check that png file exists
        # display png file
        #print('imagename = ',self.imagename)
        imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
        #print('imagebase = ',imagebase)        
        pancomp = glob.glob(self.zpdir+'/'+imagebase+"*se-pan-flux.png")[0]
        self.pancomp_png = self.plotdir+'/'+imagebase+"-se-pan-flux.png"
        #print('pancomp_png = ',self.pancomp_png)
        os.system('cp '+pancomp+' '+self.pancomp_png)
        pass
    def get_zp_magcomp_secondpass(self):
        ''' get the final plot of inst mag vs panstarrs mag'''
        # check that png file exists
        # display png file
        #print('imagename = ',self.imagename)
        imagebase = os.path.basename(self.imagename).replace('-noback-coadd.fits','').replace('.fits','-fits')
        #print('imagebase = ',imagebase)        
        pancomp = glob.glob(self.zpdir+'/'+'f'+imagebase+"*se-pan-flux.png")[0]
        self.pancomp2_png = self.plotdir+'/'+'f'+imagebase+"-se-pan-flux.png"
        #print('pancomp_png = ',self.pancomp_png)
        os.system('cp '+pancomp+' '+self.pancomp2_png)
        pass

    def get_html_data(self):
        labels = ['Date Obs','UT Time','Filter','ZP<br>(AB mag)','Max Exptime <br> (minutes)','PSF FWHM <br> (arcsec)','SE FWHM <br> (arcsec)']
        data = [self.dateobs,self.utobs,\
                self.imheader['FILTER'],\
                "{:.1f}".format(self.zp),\
                "{:.1f}".format(self.exptime/60),\
                "{:.2f}".format(self.fwhm_arcsec),\
                "{:.2f}".format(self.sefwhm_arcsec)]
        return labels,data
    def get_gredshift_filter_curve(self):
        redshift = vmain['vr'][self.keepflag]/3.e5
        header_filter = self.imheader['FILTER']
        #print('filter from header = ',header_filter,self.filter)
        if header_filter.find('ha4') > -1:
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
        self.zpdir = zpdir
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
        self.get_params_for_gal_table()
        self.plot_pointing_position()

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
        
    

    def get_filter_ratio_plot(self):
        ''' display the filter ratio png '''
        pass


    def get_params_for_gal_table(self):
        ''' setup arrays for galaxy table '''
        self.groupgals = self.r.cat[self.r.keepflag]
        self.vffil = vffil[self.r.keepflag]        

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
            self.html.write('<td>{:.2f}</td>\n'.format(self.pointing.ha.corrections[i]))                        
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
        labels=['Fit','Residuals','Residual<br> Surface']
        images = [os.path.basename(self.pointing.r.pancomp_png),\
                  os.path.basename(self.pointing.r.zpplot_png),\
                  os.path.basename(self.pointing.r.zpsurf_png)]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)
        labels=['2nd Pass Fit','Residuals','Residual<br> Surface']        
        images2 = [os.path.basename(self.pointing.r.pancomp2_png),\
                   os.path.basename(self.pointing.r.zpplot2_png),\
                   os.path.basename(self.pointing.r.zpsurf2_png)]
        buildweb.write_table(self.html,labels=labels,images=images2)

        

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
        labels=['Fit','Residuals','Residual<br> Surface']
        images = [os.path.basename(self.pointing.ha.pancomp_png),\
                  os.path.basename(self.pointing.ha.zpplot_png),\
                  os.path.basename(self.pointing.ha.zpsurf_png)]
        buildweb.write_table(self.html,labels=labels,images=images,images2=None)
        labels=['2nd Pass Fit','Residuals','Residual<br> Surface']        
        images2 = [os.path.basename(self.pointing.ha.pancomp2_png),\
                   os.path.basename(self.pointing.ha.zpplot2_png),\
                   os.path.basename(self.pointing.ha.zpsurf2_png)]
        buildweb.write_table(self.html,labels=labels,images=images2)


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
    # work through coadd directory

    # hard coding paths for now b/c this is only run in a few cases
    # could make the coadd_dir an argument that you pass in at some point

    virgovms=False

    # to rebuild all, set intonlaptop, then hdi, then ngc
    intonlaptop=True
    laptop=True
    hdi=False
    ngc=False
    if virgovms:
        vmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')
        homedir = '/mnt/qnap_home/rfinn/'
        coadd_dir = homedir+'/Halpha/reduced/virgo-coadds-HDI/'
        zpdir = homedir+'/Halpha/reduced/virgo-coadds-HDI/plots/'        
        psfdir = homedir+'/Halpha/reduced/psf-images/'
        outdir = homedir+'/research/Virgo/html-dev/coadds/'

    if laptop:
        vmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')
        #homedir = '/mnt/qnap_home/rfinn/'
        VFFIL_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main_filament_membership_allgalaxies.fits'
        vffil = fits.getdata(VFFIL_PATH)
        psfdir = homedir+'/data/reduced/psf-images/'
        outdir = homedir+'/research/Virgo/html-dev/coadds/'
        
    if hdi:
        coadd_dir = homedir+'/data/reduced/virgo-coadds-HDI/'
        zpdir = homedir+'/data/reduced/virgo-coadds-HDI/plots/'        

    if ngc:
        coadd_dir = homedir+'/data/reduced/NGC5846/'
        zpdir = homedir+'/data/reduced/NGC5846/plots/'        
        
    #rimage = coadd_dir+'VF-219.9485+5.3942-INT-20190530-p019-r-shifted.fits'    
    #haimage = coadd_dir+'VF-219.9523+5.4009-INT-20190530-p019-Halpha.fits'
    #hacsimage = coadd_dir+'VF-219.9485+5.3942-INT-20190530-p019-CS.fits'
    if intonlaptop:
        vmain = fits.getdata(homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits')    
        coadd_dir = '/home/rfinn/data/reduced/virgo-coadds-feb2019-int/'
        coadd_dir = '/home/rfinn/data/reduced/virgo-coadds-jun2019-int/'        
        zpdir = coadd_dir+'/plots/'                
        psfdir = homedir+'/data/reduced/psf-images/'
        outdir = homedir+'/research/Virgo/html-dev/coadds/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # for HDI files
    if ngc:
        # assume NGC naming convention
        allrfiles = glob.glob(coadd_dir+'nNGC*R.fits')
        print(coadd_dir+'nNGC*-R.fits')
        hfiles = glob.glob(coadd_dir+'nNGC*Ha.fits')
        print(allrfiles)
    elif hdi:
        rfiles = glob.glob(coadd_dir+'VF*-r-*coadd.fits')        
        Rfiles = glob.glob(coadd_dir+'VF*-R-*coadd.fits')
        hfiles = glob.glob(coadd_dir+'VF*-ha4-*coadd.fits')
        # combine r and R files
        allrfiles = rfiles+Rfiles
    elif intonlaptop:
        allrfiles = glob.glob(coadd_dir+'VF*-r-shifted.fits')        
        h1files = glob.glob(coadd_dir+'VF*-Halpha.fits')
        h2files = glob.glob(coadd_dir+'VF*-Ha6657.fits')        
        # combine r and R files
        hfiles = h1files+h2files

    hfiles.sort()
    allrfiles.sort()

    
    # loop through r filenames
    for i,rimage in enumerate(allrfiles):
        print(rimage)
        # find matching ha4 coadd
        if rimage.find('shifted.fits') > -1:
            h1files = glob.glob(coadd_dir+'VF*-Halpha.fits')
            try:
                haimage = h1files[0]
            except IndexError:
                
                h2files = glob.glob(coadd_dir+'VF*-Ha6657.fits')
                
                haimage = h2files[0]
                print(haimage)
            if not os.path.exists(haimage):
                print('WHAT IS HAPPENING???')
                continue
        else:
            haimage = rimage.replace('-r-','-ha4-').replace('-R-','-ha4-').replace('R.fits','Ha.fits').replace('-r-shifted.fits','-Halpha.fits')
            print('\t looking for ',haimage)
            if not os.path.exists(haimage):
                print("couldn't find it")
                # haimage could have a different date
                impath,rfile = os.path.split(rimage)
                #print(impath)
                #print(rfile)
                if intonlaptop:
                    search_string = os.path.join(impath,"VF*"+rfile.split('r')[1].replace('-r-','-ha4-').replace('-R-','-ha4-'))
                else:
                    search_string = os.path.join(impath,"VF*"+rfile.replace('-r-shifted.fits','-Ha6657.fits'))
                print('\t looking for ',search_string)
                try:
                    haimage = glob.glob(search_string)[0]
                except:
                    if not os.path.exists(haimage):
                        # assume Becky's naming convention
                        haimage = rimage.replace('-R.fits','-Ha.fits')
                        if not os.path.exists(haimage):
                            print('WARNING: could not find halpha image for ',rimage,' Skipping for now.')
                            print('\t Looking for ',haimage)
                            
                        continue
        print('###  Halpha image = ',haimage)
        # define previous gal for html links
        
        if i > 0:
            previous = os.path.basename(allrfiles[i-1]).replace('R.fits','').replace('.fits','').replace('-R','').replace('-r','')
            #print('previous = ',previous)
        else:
            previous = None
        if i < len(allrfiles)-1:
            next = os.path.basename(allrfiles[i+1]).replace('R.fits','').replace('.fits','').replace('-R','').replace('-r','')
            #print('next = ',next)
        else:
            next = None

        pname = os.path.basename(rimage).replace('R.fits','').replace('-shifted','').replace('.fits','').replace('-r','').replace('-R','')
        poutdir = os.path.join(outdir,pname)
        print(poutdir)
        try:
            p = pointing(rimage=rimage,haimage=haimage,psfdir=psfdir,zpdir=zpdir,outdir=poutdir)
            h = build_html_pointing(p,outdir=poutdir,next=next,previous=previous)
            
        except:
            print("")
            print('WE HAVE A PROBLEM!!!',rimage)
            print("")            
        plt.close('all')
        #break
