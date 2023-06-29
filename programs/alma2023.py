#!/usr/bin/env python
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from PIL import Image

# register CS halpa with legacy r image

from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp


from astropy.io import fits
from astropy.wcs import WCS

from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata.utils import Cutout2D

from astropy.coordinates import SkyCoord
from astropy.table import Table

from scipy.stats import scoreatpercentile


from astropy.cosmology import FlatLambdaCDM

# In this case we just need to define the matter density 
# and hubble parameter at z=0.

# Note the default units for the hubble parameter H0 are km/s/Mpc. 
# You can also pass an astropy `Quantity` with the units specified. 

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# TODO - add code here
def add_scale(ax,vr=1100,pscale=.331,barsize=5.,color='k'):
    '''
    ax = axis for drawing on
    vr = cosmic recession velocity
    pscale = pixel scale, in arcsec/pix
    barsize = size of marker in kpc
    '''
    z = vr/3.e5
    # get Mpc/radians
    add = cosmo.angular_diameter_distance(z)

    # convert to kpc/arcsec
    
    add_kpc_arcsec = add.value*1000*np.pi/(180*3600)
    
    # convert barsize from kpc to arcsec
    barsize_arcsec = barsize/add_kpc_arcsec
    
    # convert barsize from arcsec to pixels
    barsize_pixels = barsize_arcsec/pscale
    
    # get size of image
    x1,x2 = plt.gca().get_xlim()
    xline1 = x1 + 0.1*(x2-x1)
    y1,y2 = plt.gca().get_ylim()
    yline1 = y1 + 0.1*(y2-y1)
    # for size label
    xtext = xline1 + 0.5*barsize_pixels
    ytext = yline1 - 0.08*(y2-y1)
    
    # set up arrays for reference line
    xbar = np.array([xline1,xline1+barsize_pixels])
    ybar = np.array([yline1,yline1])  
    # plot line
    plt.plot(xbar,ybar,'r-',lw=2,color=color)
    
    # label line
    plt.text(xtext,ytext,'{:.0f} kpc'.format(barsize),horizontalalignment='center',fontsize=13,c=color)
    #plt.arrow(xline1,yline1,barsize_pixels,0)
    # draw line to show size of barsize
    
def display_gal(legacy_jpg,legacy_r,ha,gra=None,gdec=None,HImap=None,alma=False,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,shitshow=False,cmap='viridis',zoom=None,pscale=.43,ra=None,dec=None,gname=None,singlebeam=False):
    
    '''
    ARGS:
    percentile1 = min percentile for stretch of image and model
    percentile2 = max percentile for stretch of image and model
    p1residual = min percentile for stretch of residual
    p2residual = max percentile for stretch of residual
    cmap = colormap, default is viridis
    '''
    # model name
  

    # get data and header for legacy r-band image
    leg_r_data,header = fits.getdata(legacy_r,header=True)
    wcs = WCS(header) # save wcs

    # get data and header for continuum-subtracted halpha image
    image,hheader = fits.getdata(ha,header=True)
    hwcs = WCS(hheader) # save wcs

    # get the legacy color jpg image
    jpeg_data = Image.open(legacy_jpg)
    
    pngname = os.path.basename(ha).replace('.fits','.png')
    if zoom is not None:
        print("who's zoomin' who?")
        # display central region of image
        # figure out how to zoom

        ### THIS DOES NOT WORK FOR NGC5348 - PROBABLY B/C GAL IS OFF CENTER
        # get image dimensions and center
        xmax,ymax = image.shape
        if ra is not None:
            print()
            print("using ra and dec of galaxy to get center for zooming")
            print()
            xcenter,ycenter = hwcs.wcs_world2pix(ra,dec,1)
            xcenter = int(xcenter)
            ycenter = int(ycenter)
            if shitshow:
                xcenter = xmax-xcenter
                ycenter = ymax-ycenter
        else:
            xcenter = int(xmax/2)
            ycenter = int(ymax/2)

        # calculate new size to display based on zoom factor
        new_xradius = int(np.max([xcenter,xmax/2])/float(zoom))
        new_yradius = int(np.max([ycenter,ymax/2])/float(zoom))
                                                
        print("xcenter,ycenter,xmax,ymax = ",xcenter,ycenter,xmax,ymax)
        # take the max of these to account for case where original image is truncated
        new_radius = np.max([new_xradius,new_yradius])
        # calculate pixels to keep based on zoom factor
        x1 = xcenter - new_radius
        x2 = xcenter + new_radius
        y1 = ycenter - new_radius
        y2 = ycenter + new_radius
        print("new image dimensions = ",x1,x2,y1,y2)
        # check to make sure limits are not outsize image dimensions
        if (x1 < 1):
           x1 = 1
        if (y1 < 1):
           y1 = 1
        if (x2 > xmax):
            print("resetting xmax because desired value was outside image")
            x2 = xmax
        if (y2 > ymax):
            print("resetting ymax because desired value was outside image")
            y2 = ymax
        print("dimensions for halpha image = ",x1,x2,y1,y2,x2-x1,y2-y1)
        # cut images to new size
        #image = image[x1:x2,y1:y2]
        print((xcenter,ycenter),(2*new_radius,2*new_radius))
        cutout = Cutout2D(image,(xcenter,ycenter),(2*new_radius,2*new_radius),wcs=hwcs)
        image = cutout.data
        hwcs = cutout.wcs
        
        
        ## CROP JPEG IMAGE
        xmax,ymax = jpeg_data.size
        if ra is not None:
            xcenter,ycenter = wcs.wcs_world2pix(ra,dec,1)
            xcenter = int(xcenter)
            ycenter_r = int(ycenter) # maybe this fixes issue?
            # had an alignment error with one galaxy that had a cropped Halpha image
            # got rid of the following line and new I have a dec offset with the opposite sense  ...
            ycenter = int(ymax - ycenter) # maybe this fixes issue?


        else:
            xcenter = int(xmax/2)
            ycenter = int(ymax/2)
        print()
        print('for jpeg image, xcenter,ycenter = ',xcenter,ycenter,ymax/2)
        print('for fits image, xcenter,ycenter = ',xcenter,ycenter_r,ymax/2)
        print()
        # calculate new size to display based on zoom factor
        new_xradius = int(np.max([xcenter,xmax/2])/float(zoom))
        new_yradius = int(np.max([ycenter,ymax/2])/float(zoom))

        new_radius = np.max([new_xradius,new_yradius])
        # calculate pixels to keep based on zoom factor
        x1 = xcenter - new_radius
        x2 = xcenter + new_radius
        y1 = ycenter - new_radius
        y2 = ycenter + new_radius
         
        # check to make sure limits are not outsize image dimensions
        if (x1 < 1):
           x1 = 1
        if (y1 < 1):
           y1 = 1
        if (x2 > xmax):
            print("resetting xmax because desired value was outside image")
            x2 = xmax
        if (y2 > ymax):
            print("resetting ymax because desired value was outside image")
            y2 = ymax        
        
        print("dimensions for cropped jpeg image = ",x1,x2,y1,y2,x2-x1,y2-y1)
        # trying to flip x1 and x2 to see if the HI contours overlay correctly for NGC5348
        #cropped_jpg = jpeg_data.crop((x1, y1, x2, y2))
        cropped_jpg = jpeg_data.crop((x1, y1, x2, y2))
        print()
        print("shape of uncropped jpg data = ",jpeg_data.width,jpeg_data.height)
        print()
        print()
        print("shape of cropped jpg data = ",cropped_jpg.width,cropped_jpg.height)
        print()
        rcutout = Cutout2D(leg_r_data,(xcenter,ycenter_r),(x2-x1,y2-y1),wcs=wcs)
        wcscropped = rcutout.wcs
        #cropped_jpg = jpeg_data.crop((y1, x1, y2, x2))
        #cropped_jpg.save('temp.jpg')
        #jpeg_data = Image.open(cropped_jpg) 
        #return cropped_jpg

    
    v1 = [scoreatpercentile(image,percentile1)]
    #      scoreatpercentile(image,percentile1),
    #      scoreatpercentile(residual,p1residual)]
    v2 = [scoreatpercentile(image,percentile2)]
    #      scoreatpercentile(image,percentile2),
    #      scoreatpercentile(residual,p2residual)]
    norms = [simple_norm(image,'asinh',max_percent=percentile2)]
    #        simple_norm(image,'asinh',max_percent=percentile2),
    #        simple_norm(residual,'linear',max_percent=p2residual)]
               
    plt.figure(figsize=(7,4))
    plt.subplots_adjust(bottom=.01,left=.01,right=.99,top=.9,hspace=0,wspace=0)
   

    #hdu = fits.open(filename)[0]
    
 
    
    ## PLOT JPEG FROM LEGACY 

    
    if zoom is not None:
        plt.subplot(1,2,1,projection=wcscropped)
        plt.imshow(cropped_jpg)
        #plt.imshow(rcutout.data,origin="lower",norm=norms[0])
        plt.gca().set_yticks([])
        plt.gca().set_xticks([])
        ax1 = plt.gca()
        #plot_aca_centers(gname,plt.gca(),wcs)
        #plot_aca_centers(gname,plt.gca(),wcs,jpeg=True,ymax=cropped_jpg.size[1])
    else:
        print('using this one')
        plt.subplot(1,2,1,projection=wcs)
        plt.imshow(jpeg_data,origin="lower")
        ax1 = plt.gca()
        #plot_aca_centers(gname,plt.gca(),wcs)
    # plot HI contours
    if HImap is not None:
        #levels_HI = np.logspace(3,100,num=6)
        levels_HI = 2.**np.arange(10) + 1
        if gname.find('5348') > -1:
            contarr = np.array([0,1,2,3,4,5])
            levels_HI = 2.0**(contarr)*1e20
        plt.contour((HImap.data),transform=ax1.get_transform(WCS(HImap.header)),levels=levels_HI, colors='white',linestyles='-',linewidths=1)
    #plt.grid(color='white', ls='solid')
    #plt.xlabel('RA (deg)',fontsize=20)
    #plt.ylabel('Dec (deg)',fontsize=20)
    #plt.title(r'$Legacy \ grz$',fontsize=20)
    #plt.axis([x1,x2,y1,y2])
    ax1 = plt.gca()
    
    # call function to plot ACA centers
    
    
    
    if gname is not None:
        plt.text(.05,.9,gname,horizontalalignment='left',transform=ax1.transAxes,c='w',fontsize=20)
    
    ## PLOT CS HALPHA DATA
    plt.subplot(1,2,2)#,projection=wcs)
    plt.imshow(image,origin='lower',cmap=cmap,vmin=v1[0],vmax=v2[0],norm=norms[0])
    #plt.xlabel('RA (deg)',fontsize=20)
    #plt.ylabel('Dec (deg)',fontsize=20)
    #plt.title(r'$H\alpha$',fontsize=20)
    plt.gca().set_yticks([])
    plt.gca().set_xticks([])
    # call function to plot ACA centers
    # trying to plot centers on the Halpha image b/c not working on jpg
    if alma:
        #try:
        print('plotting alma centers')
        plot_alma_centers(gname,plt.gca(),hwcs,gra=gra,gdec=gdec)
        #except:
        #    print("problem plotting alma pointings")
    else:
        try:
            plot_aca_centers(gname,plt.gca(),hwcs)
        except:
            print("problem plotting ACA pointings")

    x1,x2=plt.xlim()
    y1,y2=plt.ylim()
    size=2.5/2/pscale # beam diam is 2.5 arcsec
    circ= plt.Circle((.9*x2,.1*y2),size,fill=True, color='r')
    plt.gca().add_artist(circ)
    if singlebeam:    
        iram_rad = 23/2
        kpno_rad = 55./2
        japan_rad = 7.4
        alma_rad = 50.8/2
        #beamrad_dict = {'NGC 3504':alma_rad,'NGC 4314':alma_rad,'NGC 5348':alma_rad,'NGC 5560':alma_rad,'NGC 5577':alma_rad}
        beamrad_dict = {'NGC 3504':japan_rad,'NGC 4314':iram_rad,\
                    'NGC 5348':iram_rad,'NGC 5560':kpno_rad,\
                    'NGC 5577':kpno_rad,'NGC 5566':iram_rad,\
                   'NGC 5356':kpno_rad,'UGC 09661':iram_rad,\
                   'NGC 5470':iram_rad}
        try:
            beamrad_arcsec = beamrad_dict[gname]
            size=beamrad_arcsec/pscale

            if ra is not None:
                print()
                print("using ra and dec of galaxy to get center for zooming")
                print()
                xcenter,ycenter = hwcs.wcs_world2pix(ra,dec,1)
                xcenter = int(xcenter)
                ycenter = int(ycenter)
                circ= plt.Circle((xcenter,ycenter),size,fill=False, color='w',lw=3)
                circ= plt.Circle((xcenter,ycenter),size,fill=False, color='b',lw=2)
            else:
        
                circ= plt.Circle((.5*x2,.5*y2),size,fill=False, color='w',lw=3)
                circ= plt.Circle((.5*x2,.5*y2),size,fill=False, color='b',lw=2)
            plt.gca().add_artist(circ)
        except KeyError:
            print("problem plotting beam size")
    
    plt.axis('equal')
    #plt.savefig(pngname)
    ax2 = plt.gca()

  
    
    return ax1, ax2


def plot_aca_centers(gname, ax, wcs,jpeg=False,ymax=None):
    ''' 
    GOAL:
    * read in files from Gialuca, with aca centers
    * plot positions of ACA centers with correct beam size
    
    PARAMS:
    * gname = galaxy name
    * ax = axis to plot beam circles on
    
    NOTES:
    beam sizes (diameter) from gianluca's email on 10/4/2021
    NGC 3504  beam = 13.19 ''
    NGC 4314 beam = 13.32''
    NGC 5348 beam = 11.74 ''
    NGC 5560 beam = 11.71''
    NGC 5577 beam = 11.68''
    '''
    
    beamdiam_dict = {'NGC 3504':13.19,'NGC 4314':13.32,'NGC 5348':11.74,'NGC 5560':11.71,'NGC 5577':11.68}
    beamdiam_arcsec = beamdiam_dict[gname]
    beamdiam_arcsec = 45.5
    
    # read in file containing centers
    infile = '{}_7m.pointings'.format(gname.replace('GC ',''))
    tab = Table.read(infile,format='csv',data_start=2)
    c = SkyCoord(ra=tab['RA'],dec=tab['Dec'],unit=(u.hr,u.deg))
    
    # convert centers to pixels (to plot on jpeg image)
    x,y = wcs.world_to_pixel(c)
    
    # get pixel scale
    #print(wcs)
    
    pscale = wcs.pixel_scale_matrix[1][1]
    print('pixel scale = {} deg/pix'.format(pscale))
    beamdiam_pix = beamdiam_arcsec/3600/pscale
    print('beam diam in pixels = {}'.format(beamdiam_pix))
    # plot circle at each center
    
    # flip ycoord for jpeg image
    if jpeg and (ymax is not None):
        y = ymax-y
    for xx,yy in zip(x,y):
        circ = Circle((xx,yy),beamdiam_pix,color='r',fc='None',zorder=10)
        #print('circle info = ',xx,yy,beamdiam_pix)
        ax.add_artist(circ)
        #ax.add_patch(circ)
    #plt.show()

def plot_alma_centers(gname, ax, wcs,gra=None,gdec=None,jpeg=False,ymax=None):
    ''' 
    GOAL:
    * read in files from Gialuca, with ALMA centers
    * plot positions of ALMA centers with correct beam size
    
    PARAMS:
    * gname = galaxy name
    * ax = axis to plot beam circles on
    
    NOTES:
    beam sizes (diameter) from gianluca's email on 10/4/2021
    NGC 3504  beam = 13.19 ''
    NGC 4314 beam = 13.32''
    NGC 5348 beam = 11.74 ''
    NGC 5560 beam = 11.71''
    NGC 5577 beam = 11.68''
    '''
    
    #beamdiam_dict = {'NGC 3504':13.19,'NGC 4314':13.32,'NGC 5348':11.74,'NGC 5560':11.71,'NGC 5577':11.68}
    beamdiam_dict = {'NGC 5348':50.758,\
                     'NGC 5560':50.810,\
                     'NGC 5577':50.767,\
                     'NGC 5356':50.747,\
                     'NGC 5470':50.688,\
                     'NGC 5566':50.769,\
                     'UGC 09661':50.724}
    beamdiam_arcsec = beamdiam_dict[gname]
    #beamdiam_arcsec = 50.8 # from Gianluca's README file
    
    # read in file containing centers
    try:
        infile = '{}_mosaic-12m.txt'.format(gname.replace(' ',''))
        print(infile)
        tab = Table.read(infile,format='csv',data_start=2)
        c = SkyCoord(ra=tab['RA'],dec=tab['Dec'],unit=(u.deg,u.deg))
    except:
        if gname.startswith('U'): #special case for UGC
            infile = 'mosaic_ALMA_pointings/{}_mosaic-12m.txt'.format(gname.replace(' ','').replace('0',''))
        else:
            infile = 'mosaic_ALMA_pointings/{}_moisaic-12m.txt'.format(gname.replace(' ','').replace('GC',''))
        #print(infile)
        tab = Table.read(infile,format='csv',data_start=2)
        #print(tab)
        c = SkyCoord(ra=tab['RA'],dec=tab['Dec'],unit=(u.deg,u.deg))
        
    #print(c)
    
    # convert centers to pixels (to plot on jpeg image)
    x,y = wcs.world_to_pixel(c)
    #print(x,y)
    # get pixel scale
    #print(wcs)
    
    pscale = wcs.pixel_scale_matrix[1][1]
    print('pixel scale = {} deg/pix'.format(pscale))
    beamdiam_pix = beamdiam_arcsec/3600/pscale
    print('beam diam in pixels = {}'.format(beamdiam_pix))
    # plot circle at each center
    
    # flip ycoord for jpeg image
    if jpeg and (ymax is not None):
        y = ymax-y
    for xx,yy in zip(x,y):
        #print(xx,yy)
        circ = Circle((xx,yy),beamdiam_pix/2,color='r',fc='None',zorder=10,alpha=.5)
        #print('circle info = ',xx,yy,beamdiam_pix)
        ax.add_artist(circ)
        #ax.add_patch(circ)
    #plt.show()
    
        
    
    
    
