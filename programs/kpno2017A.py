#!/usr/bin/env python


'''

Written by Rose Finn

PURPOSE:
- to help plan Halpha observations for Virgo Filament galaxies
- to make airmass plots
- to make finding charts

USEAGE:


USEFUL SITES FOR SKY CHART

http://astroweb.case.edu/jakub/TA/Query_databases.html

https://astroquery.readthedocs.io/en/latest/skyview/skyview.html

AIRMASS PLOTS

http://www.astropy.org/astropy-tutorials/Coordinates.html


'''


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

tablepath = gitpath+'Virgo/tables/'
cofile = 'CO-HI_virgo.fits'
co = fits.getdata(tablepath+cofile)

nsa = fits.getdata(nsa_file)
jmass = fits.getdata(mass_file)
wise = fits.getdata(wise_file)

ngcflag =  co.filament == 'Filament1-Group5354'
COflag = ngcflag & (co.CO_DETECT == 1)
noCOflag = ngcflag & (co.CO_DETECT == 0)

HIflag = ngcflag & (co.HI == 1)
# NGC5353/4 Filament
radec = (nsa.RA > 192.) & (nsa.RA < 209) & (nsa.DEC > 0.) & (nsa.DEC < 50.) 
radec_flag = radec & (nsa.DEC >(2*(nsa.RA - 205.) + 20) ) & (nsa.DEC < (2*(nsa.RA - 205.) + 55))
filament = radec_flag & (nsa.Z*3.e5 >2000.) & (nsa.Z*3.e5 < 3238.)
mass_flag = (jmass.MSTAR_50 > 8.3) & (jmass.MSTAR_50 < 10.2)

gas_flag = COflag | HIflag
NGCfilament = filament

minmass=8
maxmass=11.5

CJcat = fits.getdata(tablepath+'CO-MasterFile-2017May15.fits')
# find CO targets that are not in NSA?

CJngcflag = CJcat.filament == 'Filament1-Group5354'


# I sorted the pointings by RA in ipython
# sindex = argsort(pointing_dec)
# pointing_ra[sindex]
pointing_ra = np.array([ 196.95 ,#1
                         198.1  ,
                         198.075,
                         200.63 ,
                         203.395,#5
                         204.25 ,
                         206.13 ,
                         207.   ,
                         210.348,
                         210.048   ,#10
                         207.183,
                         207.89,
                         208.3  ,
                         207.452,
                         207.886,#15
                         208.48 ,
                         208.861,
                         208.64 ,
                         207.39 ,
                         206.3  ,#20
                         207.94 ,
                         207.1  ,
                         202.25  ,
                         202.34,
                         205.124 ,#25
                         209.272,
                         200.752,
                         200.73,
                         204.434])
pointing_dec = np.array([ 21.14  ,#1
                          21.59  ,
                          22.95  ,
                          28.48  ,
                          31.826 ,#5
                          31.95  ,
                          35.19  ,
                          36.16  ,
                          36.7909,
                          38.352  ,#10
                          39.387 ,
                          39.575 ,
                          39.6   ,
                          39.887 ,
                          40.2799,#15
                          40.34  ,
                          40.37  ,
                          41.308 ,
                          41.55  ,
                          41.59  ,#20
                          43.345 ,
                          43.56  ,
                          46.47  ,
                          46.78 ,
                          42.998,#25
                          41.8254,
                          23.2995,
                          26.977,
                          33.0055 ])

def find_CO_noNSA():


    virgocat = coords.SkyCoord(nsa.RA*u.degree,nsa.DEC*u.degree,frame='icrs')
    jcat = coords.SkyCoord(CJcat.RA*u.degree,CJcat.DEC*u.degree,frame='icrs')

    index,dist2d,dist3d = jcat.match_to_catalog_sky(virgocat)

    # only keep matches with matched RA and Dec w/in 1 arcsec
    matchflag = dist2d.degree < 5./3600

    print 'number in CO catalog = ',len(CJcat.RA)
    print 'number with NSA matches = ',sum(matchflag)
    return matchflag
def add_detections():
    pointing_ra = nsa.RA[COflag]
    pointing_dec = nsa.DEC[COflag]
    for i in range(len(pointing_ra)):
        rect= plt.Rectangle((pointing_ra[i]-.25,pointing_dec[i]-.25), .5, .5,fill=False, color='k')
        plt.gca().add_artist(rect)


def add_CJpoints():
    flag = CJngcflag  &  (CJcat.CO_DETECT == 1) #& CJnoNSA
    plt.plot(CJcat.RA[flag],CJcat.DEC[flag],'gs',mfc='None',mec='g',markersize=20)
    flag = CJngcflag & (CJcat.HI == 1)  #& CJnoNSA  
    plt.plot(CJcat.RA[flag],CJcat.DEC[flag],'bs',mfc='None',mec='b',markersize=25)
    

def add_pointings():
    for i in range(len(pointing_ra)):
        rect= plt.Rectangle((pointing_ra[i]-.25,pointing_dec[i]-.25), .5, .5,fill=False, color='g',lw=2)
        plt.gca().add_artist(rect)
        plt.text(pointing_ra[i]-.26,pointing_dec[i],i+1,fontsize=16,clip_on=True)
 

# match with stellar mass, NSA, WISE

def plot_positions(plotsingle=True):
    #plt.figure()
    plt.scatter(nsa.RA[COflag],nsa.DEC[COflag],s=30,c=jmass.MSTAR_50[COflag],vmin=minmass,vmax=maxmass,cmap='jet',marker='o',label='CO')
    plt.scatter(nsa.RA[noCOflag],nsa.DEC[noCOflag],s=80,c=jmass.MSTAR_50[noCOflag],vmin=minmass,vmax=maxmass,cmap='jet',marker='x',label='no CO')
    plt.scatter(nsa.RA[HIflag],nsa.DEC[HIflag],s=100,c=jmass.MSTAR_50[HIflag],vmin=minmass,vmax=maxmass,cmap='jet',marker='+',label='HI')
    #plt.scatter(nsa.RA[gasflag],nsa.DEC[gasflag],s=50,c=jmass.MSTAR_50[COflag],vmin=8,vmax=11)
    if plotsingle:
        plt.colorbar(label='$ \log_{10}(M_*/M_\odot) $')
        plt.legend()

def add_nsa():
    plt.scatter(nsa.RA[NGCfilament],nsa.DEC[NGCfilament],s=10,c=jmass.MSTAR_50[NGCfilament],vmin=minmass,vmax=maxmass,alpha=.5,cmap='jet',marker='v',label='NSA')

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

def plotall():
    plt.figure(figsize=(12,10))
    sorted_index = np.argsort(pointing_ra)
    for i in range(len(pointing_ra)):
        plt.subplot(5,5,i+1)
        zoomin(sorted_index[i]+1,plotsingle=False,delta=1)

def fix_project(ra,dec,b):
    a = np.cos(np.radians(dec))*np.cos(np.radians(ra-b.header['CRVAL1']))
    f = 1./b.header['CDELT2']*(180./np.pi)/(np.sin(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))+a*np.cos(np.radians(b.header['CRVAL2'])))
    decn = -f*(np.cos(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))-a*np.sin(np.radians(b.header['CRVAL2'])))
    ran = -f*np.cos(np.radians(dec))*np.sin(np.radians(ra-b.header['CRVAL1']))

    ran = b.header['CRVAL1']+(ran)*b.header['CDELT1']
    decn = b.header['CRVAL2']-(decn)*b.header['CDELT2']    
    return ran,decn
        
def finding_chart_all():
    for i in range(len(pointing_ra)):
        finding_chart(i)

def finding_chart(npointing):
    i = npointing-1
    galsize=0.033
    plt.figure(figsize=(6,6))
    ax=plt.gca()
    pos = coords.SkyCoord(pointing_ra[i],pointing_dec[i],frame='icrs',unit='degree')
    xout = SkyView.get_images(pos,survey=['DSS'],height=0.5*u.degree,width=0.5*u.degree)
    b=xout[0][0]
    ax.imshow(xout[0][0].data,interpolation='none',aspect='equal',cmap='gray_r',extent=[b.header['CRVAL1']-(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL1']+(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL2']+(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2'],
                                                           b.header['CRVAL2']-(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2']])
    # find galaxies on FOV

    gals = (filament | ngcflag) & (nsa.RA > (pointing_ra[i]-0.5)) & (nsa.RA < (pointing_ra[i]+0.5)) & (nsa.DEC > (pointing_dec[i]-0.5)) & (nsa.DEC < (pointing_dec[i]+0.5))
    print 'pointing ',i+1,' ngal = ',np.sum(gals)
    gindex=np.arange(len(nsa.RA))[gals]
    print 'Pointing %02d Galaxies'%(npointing),': ',nsa.NSAID[gals]
    for j in gindex:
        ran,decn=fix_project(nsa.RA[j],nsa.DEC[j],b)
        rect= plt.Rectangle((ran-galsize/2.,decn-galsize/2.), galsize, galsize,fill=False, color='c')
        plt.gca().add_artist(rect)
        plt.text(ran,decn+galsize/2.,nsa.NSAID[j],fontsize=10,clip_on=True,horizontalalignment='center')
        plt.text(ran,decn-galsize/2.,co.NED_name[j],fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='top')
        if COflag[j]:
            size=galsize-.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='g')
            plt.gca().add_artist(rect)
        if HIflag[j]:
            size=galsize+.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='b')
            plt.gca().add_artist(rect)
        #plt.legend(['filament','CO','Halpha'])
    plt.title('Pointing %02d: %s'%((i+1),pos.to_string('hmsdms')))
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.gca().invert_yaxis()
    plt.savefig('observing/2017May-Pointing%02d.png'%(i+1))

def airmass_plots():

    kpno = Observer.at_site("Kitt Peak", timezone="US/Mountain")
    observing_location = EarthLocation.of_site('Kitt Peak')
    observing_time = Time('2017-05-19 07:00')  # 1am UTC=6pm AZ mountain time
    #aa = AltAz(location=observing_location, obstime=observing_time)


    #for i in range(len(pointing_ra)):

    start_time = Time('2015-05-19 01:00:00')
    end_time = Time('2015-05-19 14:00:00')
    delta_t = end_time - start_time
    observing_time = start_time + delta_t*np.linspace(0, 1, 75)
    for j in range(3):
        plt.figure()
        legend_list = []
        for i in range(8):
            pointing_center = coords.SkyCoord(pointing_ra[8*j+i]*u.deg, pointing_dec[8*j+i]*u.deg, frame='icrs')
            if i == 3:
                plot_airmass(pointing_center,kpno,observing_time,brightness_shading=True)
            else:
                plot_airmass(pointing_center,kpno,observing_time)
            legend_list.append('Pointing %02d'%(8*j+i+1))
    
        plt.legend(legend_list)
        #plt.ylim(0.9,2.5)
        plt.gca().invert_yaxis()
        plt.subplots_adjust(bottom=.15)
        #plt.axvline(x=7*u.hour,ls='--',color='k')
        plt.axhline(y=2,ls='--',color='k')
        plt.savefig('observing/2017May-airmass-%02d.png'%(j+1))
        

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
matchflag = find_CO_noNSA()
CJnoNSA = ~matchflag
