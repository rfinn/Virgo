#!/usr/bin/env python


'''

Written by Rose Finn

PURPOSE:
- to help plan Halpha observations for Virgo Filament galaxies
- to make airmass plots
- to make finding charts

USEAGE:

source activate py36

start ipython in github/Virgo/

%run programs/kpno-halpha.p3.py


INT Feb 2019:

airmass_plots(kittpeak=False,ING=True)





To plot NE filament
make_plot()

To plot one galaxy:
zoomin(1)

To plot all galaxies that we need to observe:
plotall(delta=1)

USEFUL SITES FOR SKY CHART

http://astroweb.case.edu/jakub/TA/Query_databases.html

https://astroquery.readthedocs.io/en/latest/skyview/skyview.html

AIRMASS PLOTS

http://www.astropy.org/astropy-tutorials/Coordinates.html


'''

########################################
###### IMPORT MODULES  ########
########################################

from astropy.io import fits
from virgoCommon import *
import pylab as plt
import numpy as np
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
from astroquery.skyview import SkyView

from astropy.coordinates import EarthLocation
from astropy.time import Time

from astropy.coordinates import AltAz

from astroplan import Observer
from astroplan.plots import plot_airmass

########################################
###### RUN-SPECIFIC PARAMETERS  ########
########################################
# output file prefix
outfile_prefix = 'observing/2018March-'
max_pointing = 57

#2019
outfile_prefix = 'observing/2019Feb-INT-'
max_pointing = None

########################################
###### OTHER PARAMETERS  ########
########################################
# mass cuts for plotting NSA galaxies
minmass = 8.0
maxmass = 11.2


########################################
######  READ IN DATA TABLES  #######
########################################
tablepath = gitpath+'Virgo/tables/'
cofile = 'nsa_CO-Gianluca.virgo.fits'
co = fits.getdata(tablepath+cofile)

nsa = fits.getdata(nsa_file)
jmass = fits.getdata(mass_file)
wise = fits.getdata(wise_file)
halpha = fits.getdata(halpha_file)

CJcat = fits.getdata(tablepath+'All-virgo-20feb18_env_H070-FITS.fits')
# find CO targets that are not in NSA?

CJngcflag = CJcat.filament_name == b'1' # == 'Filament1-Group5354'


########################################
######  DEFINE SAMPLES  #######
########################################
#ngcflag =  co.filament == 'Filament1-Group5354'
ngcflag =  co.filament_name == b'1'
leoflag = co.filament_name == b'L2B'
# CO sample includes all galaxies that have been observed in CO
# not all are detection
COsample = co.CO != b''
COflag = (co.COdetected == b'1')
noCOflag = (co.COdetected == b'0')

# the following are galaxies that we still need to observe
# detected in CO but not yet observed in Halpha
ha_obs = (halpha.date_obs != b'')
ha_detect = (halpha.date_obs != b'') & (halpha.halpha == 1)
need_obs = (halpha.date_obs == b'') & (co.CO != b'') #((co.COdetected == '1') | (co.COdetected == '0'))



HIflag = (co.HI == b'1')
# NGC5353/4 Filament
radec = (nsa.RA > 192.) & (nsa.RA < 209) & (nsa.DEC > 0.) & (nsa.DEC < 50.) 
radec_flag = radec & (nsa.DEC >(2*(nsa.RA - 205.) + 20) ) & (nsa.DEC < (2*(nsa.RA - 205.) + 55))
filament = (co.filament_name !='') & (nsa.Z*3.e5 >2000.) & (nsa.Z*3.e5 < 3238.)
nsa_flag = (nsa.Z*3.e5 >1234.) & (nsa.Z*3.e5 < 3976.)
mass_flag = (jmass.MSTAR_50 > 8.3) & (jmass.MSTAR_50 < 10.2)

gas_flag = COflag | HIflag
NGCfilament = filament

# Halpha selection
# set RA and DEC as galaxies with
# CO
# no Halpha
# stellar mass between 8.5 < log(M*/Msun) < 10.  according to NSF proposal
obs_mass_flag = COsample & ~ha_obs & (jmass.MSTAR_50 > 8.5) & (jmass.MSTAR_50 < 10.) #& (nsa.SERSIC_BA > 0.2)

# resetting to COsample for 2019 observing season
filter_flag = (nsa.ZDIST*3.e5 > 2490.) & (nsa.ZDIST*3.e5 < 6000.)
obs_mass_flag = COsample & ~ha_obs #& filter_flag


########################################
###### POINTING INFO FROM MAY 2017 #####
########################################


# I sorted the pointings by RA in ipython
# sindex = argsort(pointing_dec)
# pointing_ra[sindex]
## pointing_ra = np.array([ 196.95 ,#1
##                          198.1  ,
##                          198.075,
##                          200.63 ,
##                          203.395,#5
##                          204.25 ,
##                          206.13 ,
##                          207.   ,
##                          210.348,
##                          210.048   ,#10
##                          207.183,
##                          207.89,
##                          208.3  ,
##                          207.452,
##                          207.886,#15
##                          208.48 ,
##                          208.861,
##                          208.64 ,
##                          207.39 ,
##                          206.3  ,#20
##                          207.94 ,
##                          207.1  ,
##                          202.25  ,
##                          202.34,
##                          205.124 ,#25
##                          209.272,
##                          200.752,
##                          200.73,
##                          204.434])
## pointing_dec = np.array([ 21.14  ,#1
##                           21.59  ,
##                           22.95  ,
##                           28.48  ,
##                           31.826 ,#5
##                           31.95  ,
##                           35.19  ,
##                           36.16  ,
##                           36.7909,
##                           38.352  ,#10
##                           39.387 ,
##                           39.575 ,
##                           39.6   ,
##                           39.887 ,
##                           40.2799,#15
##                           40.34  ,
##                           40.37  ,
##                           41.308 ,
##                           41.55  ,
##                           41.59  ,#20
##                           43.345 ,
##                           43.56  ,
##                           46.47  ,
##                           46.78 ,
##                           42.998,#25
##                           41.8254,
##                           23.2995,
##                           26.977,
##                           33.0055 ])

pointing_ra = nsa.RA[obs_mass_flag]
pointing_dec = nsa.DEC[obs_mass_flag]
pointing_id = nsa.NSAID[obs_mass_flag]
sorted_indices = np.argsort(pointing_ra)
pointing_dec = pointing_dec[sorted_indices]
pointing_ra = pointing_ra[sorted_indices]
pointing_id = pointing_id[sorted_indices]
nsadict = dict((a,b) for a,b in zip(pointing_id,np.arange(len(pointing_id))))

########################################
### OFFSETS TO GET MULTIPLE GALAXIES ###
########################################

####  OFFSETS  #####
# add offsets to try to get multiple galaxies in pointings
pointing_offsets_ra = np.zeros(len(pointing_ra))
pointing_offsets_dec = np.zeros(len(pointing_ra))

try:
    pointing_offsets_dec[nsadict[50209]] = -0.15
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_dec[nsadict[157256]] = 0.2
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_ra[nsadict[85510]] = -0.1
    pointing_offsets_dec[nsadict[85510]] = -0.07
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_ra[nsadict[90957]] = 0.1
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_dec[nsadict[160613]] = -0.1
except KeyError:
    print('nsa id not found in list of pointings')
    
try:
    pointing_offsets_dec[nsadict[164223]] = 0.15
except KeyError:
    print('nsa id not found in list of pointings')

pointing_ra += pointing_offsets_ra
pointing_dec += pointing_offsets_dec



def find_CO_noNSA():
    virgocat = coords.SkyCoord(nsa.RA*u.degree,nsa.DEC*u.degree,frame='icrs')
    jcat = coords.SkyCoord(CJcat.RA*u.degree,CJcat.DEC*u.degree,frame='icrs')

    index,dist2d,dist3d = jcat.match_to_catalog_sky(virgocat)

    # only keep matches with matched RA and Dec w/in 1 arcsec
    matchflag = dist2d.degree < 5./3600

    print('number in CO catalog = ',len(CJcat.RA))
    print('number with NSA matches = ',sum(matchflag))
    return matchflag
def add_detections():
    pointing_ra = nsa.RA[COflag]
    pointing_dec = nsa.DEC[COflag]
    for i in range(len(pointing_ra)):
        rect= plt.Rectangle((pointing_ra[i]-.25,pointing_dec[i]-.25), .5, .5,fill=False, color='k')
        plt.gca().add_artist(rect)


def add_CJpoints():
    flag = CJngcflag  &  (CJcat.COdetected == 1) #& CJnoNSA
    plt.plot(CJcat.RA[flag],CJcat.DEC[flag],'gs',mfc='None',mec='g',markersize=20)
    flag = CJngcflag & (CJcat.HI == 1)  #& CJnoNSA  
    plt.plot(CJcat.RA[flag],CJcat.DEC[flag],'bs',mfc='None',mec='b',markersize=25)
    

def add_pointings():
    for i in range(len(pointing_ra)):
        rect= plt.Rectangle((pointing_ra[i]-.25,pointing_dec[i]-.25), .5, .5,fill=False, color='g',lw=2)
        plt.gca().add_artist(rect)
        plt.text(pointing_ra[i]-.26,pointing_dec[i],i+1,fontsize=16,clip_on=True)
 

# match with stellar mass, NSA, WISE

def plot_positions(plotsingle=True, flag = need_obs,plotha=True):
    #if plotsingle:
    #    plt.figure()
    #flag = need_obs
    plt.scatter(nsa.RA[flag],nsa.DEC[flag],s=30,c=jmass.MSTAR_50[flag],vmin=minmass,vmax=maxmass,cmap='jet',marker='o',label='CO')
    plt.scatter(nsa.RA[noCOflag],nsa.DEC[noCOflag],s=80,c=jmass.MSTAR_50[noCOflag],vmin=minmass,vmax=maxmass,cmap='jet',marker='x',label='no CO')
    plt.scatter(nsa.RA[HIflag],nsa.DEC[HIflag],s=100,c=jmass.MSTAR_50[HIflag],vmin=minmass,vmax=maxmass,cmap='jet',marker='+',label='HI')
    #plt.scatter(nsa.RA[gas_flag],nsa.DEC[gas_flag],s=50,c=jmass.MSTAR_50[gas_flag],vmin=8,vmax=11)
    if plotha:
        plt.plot(nsa.RA[ha_obs],nsa.DEC[ha_obs],'cs',mfc='None',markersize=8)
    if plotsingle:
        plt.colorbar(label='$ \log_{10}(M_*/M_\odot) $',fraction=0.08)
        plt.legend()
        #plt.gca().invert_xaxis()
def add_nsa():
    plt.scatter(nsa.RA[nsa_flag],nsa.DEC[nsa_flag],s=10,c=jmass.MSTAR_50[nsa_flag],vmin=minmass,vmax=maxmass,alpha=.5,cmap='jet',marker='v',label='NSA')

def make_plot_2018(plotsingle=True):
    if plotsingle:
        plt.figure(figsize=(9,4))
    add_nsa()
    plot_positions(plotsingle=plotsingle,flag=obs_mass_flag)
    #add_detections()
    add_pointings()
    add_CJpoints()
    #plt.axis([190,212,15,50])
    plt.gca().invert_xaxis()
    if plotsingle:
        plt.xlabel('$RA \ (deg) $')
        plt.ylabel('$DEC \ (deg) $')
    
def make_plot(plotsingle=True):
    if plotsingle:
        plt.figure()
    add_nsa()
    plot_positions(plotsingle=plotsingle)
    #add_detections()
    add_pointings()
    add_CJpoints()
    plt.axis([190,212,15,50])
    plt.gca().invert_xaxis()
    if plotsingle:
        plt.xlabel('$RA \ (deg) $')
        plt.ylabel('$DEC \ (deg) $')

# FOV of HDI is 0.5 x 0.5 deg
# draw squares with size = FOV

def plotzoom1():
    make_plot()
    plt.axis('equal')
    plt.axis([201.55,202.55,46,47.5])
    plt.gca().invert_xaxis()

def zoomin(ngal,delta=1.5,plotsingle=True):
    make_plot(plotsingle=plotsingle)
    plt.axis([pointing_ra[ngal-1]-delta, pointing_ra[ngal-1]+delta, pointing_dec[ngal-1]-delta, pointing_dec[ngal-1]+delta])
    plt.gca().invert_xaxis()
    if plotsingle:
        plt.savefig(outfile_prefix+'%i-zoomin.png'%(ngal))
        
def plotall(delta=1.5,nrow=3,ncol=3):
    sorted_index = np.argsort(pointing_ra)
    count = 0
    nfig = 0
    
    if max_pointing != None:
        pointing_range = list(range(max_pointing))
    else:
        pointing_range = list(range(len(pointing_ra)))
        
    for i in pointing_range:
        if count == 0:
            plt.figure(figsize=(10,8))
            plt.subplots_adjust(wspace=.2,hspace=.4)
            allax = []
        plt.subplot(nrow,ncol,count+1)
        zoomin(sorted_index[i]+1,plotsingle=False,delta=delta)
        plt.title('NSAID '+str(pointing_id[i]))
        allax.append(plt.gca())
        count += 1
        if count == nrow*ncol:
            plt.colorbar(ax = allax,label='$ \log_{10}(M_*/M_\odot) $',fraction=0.08)
            plt.legend()
            plt.savefig(outfile_prefix+'plotall-%i.png'%(nfig))
            nfig += 1
            count = 0

    if count != 0:
        plt.colorbar(ax = allax,label='$ \log_{10}(M_*/M_\odot) $',fraction=0.08)
        plt.legend()
        plt.savefig(outfile_prefix+'plotall-%i.png'%(nfig)) 
def fix_project(ra,dec,b):
    a = np.cos(np.radians(dec))*np.cos(np.radians(ra-b.header['CRVAL1']))
    f = 1./b.header['CDELT2']*(180./np.pi)/(np.sin(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))+a*np.cos(np.radians(b.header['CRVAL2'])))
    decn = -f*(np.cos(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))-a*np.sin(np.radians(b.header['CRVAL2'])))
    ran = -f*np.cos(np.radians(dec))*np.sin(np.radians(ra-b.header['CRVAL1']))

    ran = b.header['CRVAL1']+(ran)*b.header['CDELT1']
    decn = b.header['CRVAL2']-(decn)*b.header['CDELT2']    
    return ran,decn
        
def finding_chart_all():
    if max_pointing != None:
        pointing_range = range(1,max_pointing+1)
    else:
        pointing_range = range(1,len(pointing_ra)+1)
                               
    for i in pointing_range:
        plt.close('all')
        finding_chart(i)

def finding_chart(npointing,delta_image = .25,offset_ra=0.,offset_dec=0.,plotsingle=True,ING_flag=False):
    i = npointing-1
    galsize=0.033
    center_ra = pointing_ra[i]+offset_ra
    center_dec = pointing_dec[i] + offset_dec
    print('center ra, dec = ',center_ra,center_dec)
    if plotsingle:
        if delta_image > .3:
            plt.subplots_adjust(right=.9,top=.9,bottom=.1,left=.1)
            plt.figure(figsize=(10,9))
        else:
            plt.figure(figsize=(6,6))
    ax=plt.gca()
    pos = coords.SkyCoord(center_ra,center_dec,frame='icrs',unit='degree')
    delta_imagex=2.*delta_image
    delta_imagey=2.*delta_image
    if ING_flag:
        pos = coords.SkyCoord(center_ra,center_dec,frame='icrs',unit='degree')
        delta_imagex=.8 #image width in deg
        delta_imagey=.8 # image width in deg
        
        xout = SkyView.get_images(pos,survey=['DSS'],height=delta_imagex*u.degree,width=delta_imagey*u.degree)
    else:
        xout = SkyView.get_images(pos,survey=['DSS'],height=2*delta_image*u.degree,width=2.*delta_image*u.degree)
    b=xout[0][0]
    ax.imshow(xout[0][0].data,interpolation='none',aspect='equal',cmap='gray_r',extent=[b.header['CRVAL1']-(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL1']+(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL2']+(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2'],
                                                           b.header['CRVAL2']-(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2']])

    if ING_flag:
        # add footprint of WFC chips
        plot_INT_footprint(center_ra,center_dec)
    else:
        rect= plt.Rectangle((center_ra-delta_image,center_dec-delta_image), 0.5, 0.5,fill=False, color='k')
        plt.gca().add_artist(rect)

    #add_cameras()
    # find galaxies on FOV
    gals = (nsa.RA > (pos.ra.value-delta_imagex/2.)) & (nsa.RA < (pos.ra.value+delta_imagex/2)) & (nsa.DEC > (pos.dec.value-delta_imagey/2.)) & (nsa.DEC < (pos.dec.value+delta_imagey/2.))
    print('pointing ',i+1,' ngal = ',np.sum(gals))
    gindex=np.arange(len(nsa.RA))[gals]
    print('Pointing %02d Galaxies'%(npointing),': ',nsa.NSAID[gals])
    for j in gindex:
        ran,decn=fix_project(nsa.RA[j],nsa.DEC[j],b)
        rect= plt.Rectangle((ran-galsize/2.,decn-galsize/2.), galsize, galsize,fill=False, color='c')
        plt.gca().add_artist(rect)
        s='%i\n vr=%i'%(nsa.NSAID[j],nsa.ZDIST[j]*3.e5)
        plt.text(ran,decn+galsize/2.,s,fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='bottom')
        plt.text(ran,decn-galsize/2.,co.NEDname[j].decode("utf-8"),fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='top')
        if COflag[j]:
            size=galsize-.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='g')
            plt.gca().add_artist(rect)
        if HIflag[j]:
            size=galsize+.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='b')
            plt.gca().add_artist(rect)
        if ha_obs[j]:
            size=galsize+.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='r')
            plt.gca().add_artist(rect)
        #plt.legend(['filament','CO','Halpha'])
    plt.title('Pointing %02d: %s'%((i+1),pos.to_string('hmsdms')))
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.gca().invert_yaxis()
    if plotsingle:
        plt.savefig(outfile_prefix+'Pointing%02d.png'%(i+1))

def plot_INT_footprint(center_ra,center_dec):
    #using full detector sizes for now because 
    detector_dra = 2048.*0.33/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    detector_ddec = 4100.*0.33/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    # draw footprint of chip 4
    rect= plt.Rectangle((center_ra-detector_dra/2.,center_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)
    # draw footpring of chip 3
    # assuming chip 3 is WEST of chip 4
    offset_ra = -1.*detector_dra-17./3600. # 17 arcsec gap in RA between 
    offset_dec = 9.5/3600. # 9.5 arcsec offset toward N
    rect= plt.Rectangle((center_ra+offset_ra-detector_dra/2.,center_dec+offset_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)

    # draw footpring of chip 1
    # assuming chip 1 is EAST of chip 4
    offset_ra = detector_dra+22.7/3600. # 17 arcsec gap in RA between 
    offset_dec = -3.18/3600. # 9.5 arcsec offset toward N
    rect= plt.Rectangle((center_ra+offset_ra-detector_dra/2.,center_dec+offset_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)

    # draw footpring of chip 2
    # assuming chip 2 is NORTH of chip 4
    offset_ra = detector_dra/2.-19.2/3600. # hard to explain
    offset_dec =  detector_ddec/2.+23.8/3600.# hard to explain
    # this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
    rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -1.*detector_ddec, detector_dra,fill=False, color='k')
    plt.gca().add_artist(rect)

    # adding guide camera
    offset_ra = -1.5*detector_dra-(22.7+98.1)/3600. # hard to explain
    offset_dec =  (3.18+649.9)/3600.# hard to explain
    # this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
    rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -7./60., 7./60,fill=False, color='k')
    plt.gca().add_artist(rect)
    
        
def add_cameras():
    # add guidestar cameras
    # North camera
    delta_ra = 0
    delta_dec = 2610./3600
    width_ra = 3.3/60.
    width_dec = 2.2/60.
    rect= plt.Rectangle((center_ra+delta_ra - 0.5*width_ra,center_dec+delta_dec - 0.5*width_dec), width_ra, width_dec,fill=False, color='k')
    plt.gca().add_artist(rect)
    # South camera
    delta_ra = -25./3600
    delta_dec = -2410./3600
    width_ra = 3.3/60.
    width_dec = 2.2/60.
    rect= plt.Rectangle((center_ra+delta_ra - 0.5*width_ra,center_dec+delta_dec - 0.5*width_dec), width_ra, width_dec,fill=False, color='k')
    plt.gca().add_artist(rect)

    
def finding_chart_with_guide_stars(npointing,offset_ra=0.,offset_dec=0.):
    finding_chart(npointing, delta_image=.75,offset_ra=offset_ra, offset_dec=offset_dec,plotsingle=False)
    plt.savefig(outfile_prefix+'Pointing%02d-guiding.png'%(npointing))

def make_all_platinum(ING_flag=False):
    if max_pointing != None:
        pointing_range = range(max_pointing)
    else:
        pointing_range = range(len(pointing_ra))
    for i in pointing_range:
        plt.close('all')
        if ING_flag:
            platinum_finding_chart_ING(i+1)
        else:
            platinum_finding_chart(i+1,ING_flag=ING_flag)

def platinum_finding_chart(npointing,offset_ra=0.,offset_dec=0.,ING_flag=False):
    fig = plt.figure(figsize = (10,6))
    grid = plt.GridSpec(2,3,hspace=.4,wspace=.2,left=.05)
    hdi = fig.add_subplot(grid[:,:-1])
    finding_chart(npointing,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False,ING_flag = ING_flag)
    south_camera = fig.add_subplot(grid[1,2])
    show_guide_camera(npointing,south_camera=True,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
    north_camera = fig.add_subplot(grid[0,2])
    show_guide_camera(npointing,south_camera=False,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
    plt.savefig(outfile_prefix+'Pointing%03d-platinum.png'%(npointing))

def platinum_finding_chart_ING(npointing,offset_ra=0.,offset_dec=0.):
    fig = plt.figure(figsize = (8,8.))
    #grid = plt.GridSpec(2,3,hspace=.4,wspace=.2,left=.05)
    #hdi = fig.add_subplot(grid[:,:-1])
    finding_chart(npointing,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False,ING_flag = True)
    plt.savefig(outfile_prefix+'Pointing%03d-platinum.png'%(npointing))

def guide_cameras(npointing,offset_ra=0,offset_dec=0):
    i = npointing - 1

    plt.figure(figsize=(10,5))
    plt.subplots_adjust(wspace=.3)
    plt.subplot(1,2,1)
    show_guide_camera(npointing,south_camera=False,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
    plt.subplot(1,2,2)
    show_guide_camera(npointing,south_camera=True,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)

    mytitle = "Guide Cameras: Pointing %02d (NSAID %i)"%((i+1),pointing_id[i])
    ax = plt.gca()
    plt.text(-.1,1.15,mytitle,transform=ax.transAxes,horizontalalignment='center',fontsize=18)

def show_guide_camera(npointing,south_camera=True,offset_ra=0,offset_dec=0,plotsingle=True):
    '''
    offset_ra in arcsec
    offset_dec in arcsec

    camera = 0 for south
    camera = 1 for north
    '''
    i = npointing-1
    if south_camera:
        # South camera
        delta_ra = -25./3600
        delta_dec = -2410./3600
        width_ra = 3.3/60.
        width_dec = 2.2/60.
        outputfile=outfile_prefix+'Pointing%02d-guiding-south.png'%(npointing)
        mytitle = "Pointing %02d South Camera (3.3'x2.2'): NSAID %i"%((i+1),pointing_id[i])
    else:
        # offsets and dimensions of 
        # North camera
        delta_ra = 0
        delta_dec = 2610./3600
        width_ra = 3.3/60.
        width_dec = 2.2/60.
        outputfile=outfile_prefix+'Pointing%02d-guiding-north.png'%(npointing)
        mytitle = "Pointing %02d North Camera (3.3'x2.2'): NSAID %i"%((i+1),pointing_id[i])
    center_ra = pointing_ra[i]+offset_ra/3600. + delta_ra
    center_dec = pointing_dec[i] + offset_dec/3600. + delta_dec

    delta_image = width_ra
    if plotsingle:
        plt.figure(figsize=(6,6))
    ax = plt.gca()
    pos = coords.SkyCoord(center_ra,center_dec,frame='icrs',unit='degree')
    xout = SkyView.get_images(pos,survey=['DSS'],height=2*delta_image*u.degree,width=2.*delta_image*u.degree)
    b=xout[0][0]
    ax.imshow(xout[0][0].data,interpolation='none',aspect='equal',cmap='gray_r',extent=[b.header['CRVAL1']-(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL1']+(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL2']+(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2'],
                                                           b.header['CRVAL2']-(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2']])

    rect= plt.Rectangle((center_ra - 0.5*width_ra,center_dec - 0.5*width_dec), width_ra, width_dec,fill=False, color='k')
    plt.gca().add_artist(rect)
    if plotsingle:
        plt.title(mytitle)
    else:
        if south_camera:
            plt.title("South Camera (3.3'x2.2')")
        else:
            plt.title("North Camera (3.3'x2.2')")
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.gca().invert_yaxis()
    if plotsingle:
        plt.savefig(outputfile)



        
def airmass_plots(kittpeak=True,ING=False):

    observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
    if kittpeak:
        print('plotting airmass curves for Kitt Peak')
        observing_location = EarthLocation.of_site('Kitt Peak')
        observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
    elif ING:
        print('plotting airmass curves for INT')
        observing_location = EarthLocation.of_site(u'Roque de los Muchachos')
        observer_site = Observer.at_site("Roque de los Muchachos", timezone="GMT")
        
    #observing_time = Time('2017-05-19 07:00')  # 1am UTC=6pm AZ mountain time
    #observing_time = Time('2018-03-12 07:00')  # 1am UTC=6pm AZ mountain time
    #aa = AltAz(location=observing_location, obstime=observing_time)


    #for i in range(len(pointing_ra)):

    start_time = Time('2017-03-12 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
    end_time = Time('2017-03-12 14:00:00')

    # for run starting 2019-Feb-04 at INT
    start_time = Time('2019-02-04 19:00:00') # INT is on UTC
    end_time = Time('2019-02-05 07:00:00')

    delta_t = end_time - start_time
    observing_time = start_time + delta_t*np.linspace(0, 1, 75)
    nplots = int(sum(obs_mass_flag)/8.)
    print(nplots)
    for j in range(nplots):
        plt.figure()
        legend_list = []
        for i in range(8):
            pointing_center = coords.SkyCoord(pointing_ra[8*j+i]*u.deg, pointing_dec[8*j+i]*u.deg, frame='icrs')
            if i == 3:
                plot_airmass(pointing_center,observer_site,observing_time,brightness_shading=True)
            else:
                plot_airmass(pointing_center,observer_site,observing_time)
            legend_list.append('Pointing %02d'%(8*j+i+1))
    
        plt.legend(legend_list)
        #plt.ylim(0.9,2.5)
        plt.gca().invert_yaxis()
        plt.subplots_adjust(bottom=.15)
        #plt.axvline(x=7*u.hour,ls='--',color='k')
        plt.axhline(y=2,ls='--',color='k')
        plt.savefig(outfile_prefix+'airmass-%02d.png'%(j+1))
        

    ##     delta_hours = np.linspace(0, 12, 100)*u.hour
    ##     full_night_times = observing_time + delta_hours
    ##     full_night_aa_frames = AltAz(location=observing_location, obstime=full_night_times)
    ##     full_night_aa_coos = pointing_center.transform_to(full_night_aa_frames)

    ##     plt.plot(delta_hours, full_night_aa_coos.secz,label='Pointing %02d'%(i+1))
    ##     plot_airmass(pointing_center, observing_location, observing_time)
    ##     plt.xlabel('Hours from 6 pm AZ time')
    ##     plt.ylabel('Airmass [Sec(z)]')
    ##     plt.ylim(0.9,2.5)
    ##     plt.tight_layout()

def check_CO():
    '''
    look for galaxies that are in Gianluca's CO file but not in the NSA catalog
    '''
    matchflag = find_CO_noNSA()
    CJnoNSA = ~matchflag
    for i in range(len(CJnoNSA)):
        if CJnoNSA[i]:
            print('%18s %18s %.8f %.8f'%(CJcat.source_name[i], CJcat.NEDname[i],CJcat.RA[i],CJcat.DEC[i]))
        
    '''
    results in
    number in CO catalog =  227
    number with NSA matches =  226
           UGC8656          UGC 08656 205.12995833 42.99388889
    '''
