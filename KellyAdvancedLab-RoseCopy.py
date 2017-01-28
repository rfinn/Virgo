
# coding: utf-8

# ## Kelly's Advanced Lab Project ##
# 
# 
# ### Project Goals ###
# 
# * Review the literature on Virgo filaments - what has been done already?
# * identify galaxies in the filaments surrounding the Virgo Cluster
#     * this will require an analysis of their redshifts
# * make a catalog that characterizes their
#     * stellar mass
#     * star formation rate from 22-micron WISE flux and GALEX UV flux
#     * bulge-to-total ratio
#     * disk scale length
# * make a poster showing SDSS color images of filament galaxies
# * compile other information on these galaxies that you will get from the literature
#     * existing Halpha observations
#     * radio observations
# * compare galaxy properties to matched field and cluster samples
#     * create a matched field sample
#     * use Becky's Virgo cluster sample
# * make a plot of SFR vs stellar mass for these galaxies in comparison to field and cluster samples.
# * prepare proposal to get Halpha imaging of Virgo filament galaxies
# 
# ### To do 9/26 ###
# 
# * get python 2 installed
# * read in NSA catalog
# * make a plot showing the region around Virgo, showing galaxies with:
#     * raflag = (nsa.RA > 150.) & (nsa.RA < 220.) 
#     * decflag= (nsa.DEC > -10.) & (nsa.DEC < 50.) 
#     * velflag = (nsa.ZDIST*3.e5 > 1000.) & (nsa.ZDIST*3.e5 < 3000.)
#     * vflag = raflag & decflag & velflag
# * color-code points according to recession velocity
# 
# 
# ### To do 10/10 ###
# * after we get Simard files:
#      * Make plot of cluster, color coded by bulge to total ratio
#      * make color-color plots of filaments, color coded by bulge to total ratio
#      
# ### To do 10/24 ###
# * make NUV-r vs r-22 plot for Virgo and all filaments combined
#     * color code by stellar mass
#     
# ### To do 10/31 ###
# * make NUV-r vs r-22 plot for Virgo and all filaments combined
#     * color code by stellar mass
#     * with error bars- excluding bad points
# * propagating the 22 micron error- see camera roll
# 
'''
GOAL:
  - to convert RA,DEC,redshift to supergalactic coordinates
  - to use these coordinates to select filament galaxies around Virgo
  
PROCEDURE:
  - 
USAGE:

NOTES:

'''

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import sys
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from mpl_toolkits.mplot3d import Axes3D
from astropy.cosmology import WMAP9 as cosmo
from astropy.constants import c

#Get current path so program can tell if this is being run on Kelly's or Rose's computer
mypath=os.getcwd()
if mypath.find('rfinn') > -1:
    print "Running on Rose's computer"
    #agcfile='/Users/rfinn/idl/programs/idl_alfa/agctotal.sav'
    gitpath='/Users/rfinn/github/'
elif mypath.find('kelly') > -1:
    print "Running on Kellys's computer"
    gitpath='/Users/kellywhalen/Github/'


#### READ IN DATA #####

#NSA Table
NSA_file = gitpath+'Virgo/nsa_v0_1_2_with_AGC.Virgo.fits'
NSA_file='/Users/rfinn/research/nsa/nsa_v0_1_2.fits'
nsa = fits.getdata(NSA_file)

# read in WISE catalog
#wisefile=gitpath+'Virgo/nsa_v0_1_2_wise.Virgo.fits'
#wise=fits.getdata(wisefile)

#read in John's stellar masses
#massfile=gitpath+'Virgo/nsa_v1_2_fsps_v2.4_miles_chab_charlot_sfhgrid01.Virgo.fits'
#jmass=fits.getdata(massfile)

## #Read in Simard Tables
## # table 1 = B/T using n=4 for bulge
## simard1file = gitpath+'Virgo/Simard1ToNSA.fits'
## simard1 = fits.getdata(simard1file)

## # table 2 - let n be a free parameter
## simard2file = gitpath+'Virgo/Simard2ToNSA.fits'
## simard2 = fits.getdata(simard1file)


## # pure sersic profile
## simard3file = gitpath+'Virgo/Simard3ToNSA.fits'
## simard3 = fits.getdata(simard3file)



#set flags
raflag = (nsa.RA > 115.) & (nsa.RA < 240.)
decflag= (nsa.DEC > -30.) & (nsa.DEC < 60.)
velflag =  (nsa.Z*3.e5 < 3000.) # & (nsa.Z*3.e5 > 1000.) 
vflag = raflag & decflag & velflag 

nsa = nsa[vflag]


# NGC5353/4 Filament - selection I used in NSF proposal
radec = (nsa.RA > 192.) & (nsa.RA < 209) & (nsa.DEC > 0.) & (nsa.DEC < 50.) 
radec_flag = radec & (nsa.DEC >(2*(nsa.RA - 205.) + 20) ) & (nsa.DEC < (2*(nsa.RA - 205.) + 55))
filament = radec_flag & (nsa.Z*3.e5 >2000.) & (nsa.Z*3.e5 < 3238.)
plt.scatter(nsa.RA[filament],nsa.DEC[filament],c=nsa.Z[filament]*3.e5,zorder=20,s=20,vmin=1000,vmax=3000,lw=0.5)
xl = np.linspace(196,230,100)
yl = (2*(xl - 205.) + 20)
#plt.plot(xl,yl,'r-')
NGCfilament = filament


# wise flag
#wflag = vflag & (wise.W4MPRO > 0.1) & (wise.W4SNR > 2.)

# define nsa RA and Dec as SkyCoord
nsa_sc = SkyCoord(nsa.RA*u.degree, nsa.DEC*u.degree)

# convert RA and DEC to galactic
# l = galcoords.l.degree
# b = galcoords.b.degree
################################################
#####  FLOW MODEL  ###############
################################################
# using mould+00 method to correct redshifts
# described in detail in appendix A
#
# http://adsabs.harvard.edu/abs/2000ApJ...529..786M

V_H = nsa.Z*c.to('km/s')


# 1. correction of observed heliocentric velocity to centroid of local group
# VLG = VH - 79 cos l cos b + 296 sin l cos b - 36 sin b

V_LG = - 79.*u.km/u.second*np.cos(nsa_sc.galactic.l.radian)*np.cos(nsa_sc.galactic.b.radian) + 296.*u.km/u.second*np.sin(nsa_sc.galactic.l.radian)*np.cos(nsa_sc.galactic.b.radian) - 36.*u.km/u.second*np.sin(nsa_sc.galactic.b.radian)


# 2. Correction for Virgo infall
# from Mould+2000, ApJ, 529, 786
# not sure if we are using this equation exactly right

V_fid = 200.*u.km/u.second # infall of LG into Virgo?
# Virgo coordinates given by Mould+2000
Virgo = SkyCoord('12h28m19s', '+12d40m00s', frame='fk5',equinox='J1950.') # epoch = 1950
Virgo = Virgo.transform_to(FK5(equinox='J2000'))
### need to fix this to use spherical distance
#theta = np.sqrt((nsa_sc.ra.radian - Virgo.ra.radian)**2 + (nsa_sc.dec.radian - Virgo.dec.radian)**2)
theta = Virgo.separation(nsa_sc).radian 
# cluster radius in deg
gamma = 2.

V_a = 1035.*u.km/u.second # recession vel of Virgo from Mould+2000
V_o = nsa.Z*c.to('km/s') + V_LG # recession velocities of the galaxies
r_oa = np.sqrt(V_o**2 + V_a**2 - 2.*V_o*V_a*np.cos(theta))
V_infall = V_fid*(np.cos(theta) + (V_o - V_a*np.cos(theta))/r_oa*(r_oa/V_a)**(1-gamma))



# 3. Correction for GA infall

V_fid = 400.*u.km/u.second # infall of LG into Virgo?
# GA coordinates given by Mould+2000
GA = SkyCoord('13h20m00s', '+44d00m00s', frame='fk5',equinox='J1950.') # epoch = 1950
GA = GA.transform_to(FK5(equinox='J2000'))
theta = GA.separation(nsa_sc).radian
#theta = np.sqrt((nsa_sc.ra.radian - GA.ra.radian)**2 + (nsa_sc.dec.radian - GA.dec.radian)**2)

gamma = 2.
V_a = 4380.*u.km/u.second # recession vel of Virgo from Mould+2000
V_o = nsa.Z*c.to('km/s') + V_LG # recession velocities of the galaxies
r_oa = np.sqrt(V_o**2 + V_a**2 - 2.*V_o*V_a*np.cos(theta))
V_GA = V_fid*(np.cos(theta) + (V_o - V_a*np.cos(theta))/r_oa*(r_oa/V_a)**(1-gamma))

# 4. Correction for Shapley supercluster infall.

# Final, correction cosmic velocity is
#
#  Vcosmic = VH + Vc,LG + Vin,Virgo + Vin,GA + Vin,Shap + ...

# ## Tranforming to Supergalactic Coordinates##
# 
# Looking to match the plots shown in Kim+2016
# 
# https://arxiv.org/abs/1611.00437

#Plot of Virgo Cluster in galactic coordinates
V_cosmic = V_H + V_LG + V_infall #+V_GA
V = V_cosmic
#V = nsa.Z*c.to('km/s')
#V = nsa.ZDIST*c.to('km/s')

SGX = V/cosmo.H(0)*np.cos(nsa_sc.supergalactic.sgl.radian)*np.cos(nsa_sc.supergalactic.sgb.radian)
SGY = V/cosmo.H(0)*np.sin(nsa_sc.supergalactic.sgl.radian)*np.cos(nsa_sc.supergalactic.sgb.radian)
SGZ = V/cosmo.H(0)*np.sin(nsa_sc.supergalactic.sgb.radian)

H0 = 74.*u.km/u.second/u.Mpc
SGX = V/H0*np.cos(nsa_sc.supergalactic.sgl.radian)*np.cos(nsa_sc.supergalactic.sgb.radian)
SGY = V/H0*np.sin(nsa_sc.supergalactic.sgl.radian)*np.cos(nsa_sc.supergalactic.sgb.radian)
SGZ = V/H0*np.sin(nsa_sc.supergalactic.sgb.radian)

# latitude and longitude of Virgo based on NED
#l = 102.93 # longitude, deg
#b = -2.73 # latitude, deg

distance_to_virgo = 16.5*u.Mpc #Mpc
SGX_Virgo = distance_to_virgo*np.cos(Virgo.supergalactic.sgl.radian)*np.cos(Virgo.supergalactic.sgb.radian)
SGY_Virgo = distance_to_virgo*np.sin(Virgo.supergalactic.sgl.radian)*np.cos(Virgo.supergalactic.sgb.radian)
SGZ_Virgo = distance_to_virgo*np.sin(Virgo.supergalactic.sgb.radian)

# difference 
DSGX = SGX - SGX_Virgo
DSGY = SGY - SGY_Virgo
DSGZ = SGZ - SGZ_Virgo



# don't use objects w/in 3.6 Mpc of Virgo (Kim+2016)

dist = np.sqrt(DSGX**2 + DSGY**2 + DSGZ**2) 
dist_flag = dist.value > 3.6


# In[7]:

kim_vflag = nsa.Z*3.e5 < 3300.  
kim_raflag = (nsa.RA > 115.) & (nsa.RA < 240.)
kim_decflag = (nsa.DEC > -35.) & (nsa.DEC < 60.)
kim_flag = kim_vflag & kim_raflag & kim_decflag & dist_flag 

d_split = 16.
#In front of cluster
SGYfront_flag = (SGY > 4. * u.Mpc) & (SGY < d_split * u.Mpc)
SGYback_flag = (SGY > d_split * u.Mpc) & (SGY < 40. * u.Mpc)
#SGYback_flag = (SGY > 21. * u.Mpc) & (SGY < 27. * u.Mpc)

SGYvirgo =  (SGY > 4.*u.Mpc) & (SGY < d_split*u.Mpc)



# ** First plot 2-D projections in SGX-SGZ plane, and SGX-SGZ **

# In[8]:
def compare_vel():
    plt.figure()
    plt.plot(nsa.Z*c.to('km/s'),V_cosmic,'b.')
    xl = np.linspace(500,3500,20)
    plt.plot(xl,xl,'k-')
    plt.xlabel('nsa.Z*c')
    plt.ylabel('V_cosmic')
def plotxzplane():
    plt.figure(figsize=(6,4))
    plt.scatter(SGX[kim_flag & SGYvirgo],SGZ[kim_flag & SGYvirgo],alpha=.5,c=SGY[kim_flag & SGYvirgo])
    plt.colorbar(fraction=.08,label='SGY (Mpc)')
    plt.subplots_adjust(bottom=.15)
    #plt.axis([-15,10,-10,20])
    plt.xlabel('$\Delta SGX \ (Mpc)$')
    plt.ylabel('$\Delta SGZ \ (Mpc)$')
    plt.title('Region around Virgo Cluster (Compare to Kim Fig 2)')
    plt.axis([-15,18,-17,12])
    plt.axvline(x=0,ls='--',color='k')
    plt.axhline(y=0,ls='--',color='k')

    ## plt.figure(figsize=(6,4))
    ## plt.scatter(SGX[kim_flag & SGYfront_flag],SGZ[kim_flag & SGYfront_flag],alpha=.5,c=SGY[kim_flag & SGYfront_flag])
    ## plt.colorbar(fraction=.08)
    ## plt.axis([-15,10,-10,20])
    ## plt.xlabel('SGX (Mpc)')
    ## plt.ylabel('SGZ (Mpc)')
    ## plt.title('Region In Front of Virgo Cluster')
    ## plt.axvline(x=0,ls='--',color='k')
    ## plt.axhline(y=0,ls='--',color='k')

    plt.figure(figsize=(6,4))
    plt.subplots_adjust(bottom=.15)
    plt.scatter(DSGX[kim_flag & SGYback_flag],DSGZ[kim_flag & SGYback_flag],alpha=.5,c=SGY[kim_flag & SGYback_flag])
    #plt.scatter(SGX[NGCfilament],SGZ[NGCfilament],alpha=1,c=SGY[NGCfilament],s=60)
    plt.colorbar(fraction=.08,label='SGY (Mpc)')
    plt.axis([-15,10,-10,20])
    plt.xlabel('$\Delta SGX \ (Mpc)$')
    plt.ylabel('$\Delta SGZ \ (Mpc)$')
    plt.title('Region Behind Virgo Cluster (Compare to Kim Fig 3)')
    plt.axvline(x=0,ls='--',color='k')
    plt.axhline(y=0,ls='--',color='k')
    plt.axis([-18,7,-10,14])


def plotxyplane():
    plt.figure(figsize=(6,4))
    plt.scatter(DSGX[kim_flag & SGYvirgo],DSGY[kim_flag & SGYvirgo],alpha=.5,c=DSGZ[kim_flag & SGYvirgo])
    plt.colorbar(fraction=.08,label='SGZ (Mpc)')
    plt.subplots_adjust(bottom=.15)
    #plt.axis([-15,10,-10,20])
    plt.xlabel('$\Delta SGX \ (Mpc)$')
    plt.ylabel('$\Delta SGY \ (Mpc)$')
    plt.title('Region around Virgo Cluster')
    plt.axis([-15,18,-17,12])
    plt.axvline(x=0,ls='--',color='k')
    plt.axhline(y=0,ls='--',color='k')



def plotv3d():
    fig = plt.figure(figsize = (10,8))
    ax = fig.add_subplot(111, projection='3d')
    flag = kim_flag & (SGZ < 12.*u.Mpc) & (SGZ > -17.*u.Mpc) & (SGX > -15.*u.Mpc) & (SGX < 17.*u.Mpc)
    ax.scatter(DSGX[flag], DSGY[flag], DSGZ[flag])
    plt.title('Region Around Virgo Cluster')
    plt.xlabel('SGX (Mpc)')
    plt.ylabel('SGY (Mpc)')
    #plt.zlabel('SGZ (Mpc)')

def plotpositions():
    # plot region covered by these fits files
    plt.figure(figsize=(12,4))
    plt.subplot(1,2,1)
    flag = SGY > 4.*u.Mpc
    plt.scatter(nsa.RA[flag],nsa.DEC[flag],c=nsa.Z[flag]*c.to('km/s'),s=15,edgecolors='None')
    plt.colorbar(fraction=.08)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.subplot(1,2,2)
    plt.hist(nsa.Z*c.to('km/s')/1000,bins=50,color='0.5')
    plt.xlabel('V_H/1000')




