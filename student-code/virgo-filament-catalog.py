#Load in all packages
import numpy as np
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

def VirgoCatalog(NSA_file):


    #read in .fits file
    nsa = fits.getdata(NSA_file)
    nsadict=dict((a,b) for a,b in zip(nsa.NSAID,np.arange(len(nsa.NSAID))))
    
    #Pick out galaxies in the general region of virgo
    raflag = (nsa.RA > 115.) & (nsa.RA < 240.)
    decflag= (nsa.DEC > -30.) & (nsa.DEC < 60.)
    velflag =  (nsa.Z*3.e5 < 3000.) 
    vflag = raflag & decflag & velflag 

    # define nsa RA and Dec as SkyCoord
    nsa_sc = SkyCoord(nsa.RA*u.degree, nsa.DEC*u.degree)
    
    # convert helio-centric velocity to units of km/s
    V_H = nsa.Z*c.to('km/s')

    # 1. correction of observed heliocentric velocity to centroid of local group
    # VLG = VH - 79 cos l cos b + 296 sin l cos b - 36 sin b

    V_LG = - 79.*u.km/u.second*np.cos(nsa_sc.galactic.l.radian)*np.cos(nsa_sc.galactic.b.radian) + 296.*u.km/u.second*np.sin(nsa_sc.galactic.l.radian)*np.cos(nsa_sc.galactic.b.radian) - 36.*u.km/u.second*np.sin(nsa_sc.galactic.b.radian)

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
    #V_a = 1221.*u.km/u.second # recession vel of Virgo as calculated by Gianluca
    V_o = nsa.Z*c.to('km/s') + V_LG # recession velocities of the galaxies
    r_oa = np.sqrt(V_o**2 + V_a**2 - 2.*V_o*V_a*np.cos(theta))
    V_infall = V_fid*(np.cos(theta) + (V_o - V_a*np.cos(theta))/r_oa*(r_oa/V_a)**(1-gamma))
    
    V_fid = 400.*u.km/u.second # infall of LG into Virgo?
    # GA coordinates given by Mould+2000
    GA = SkyCoord('13h20m00s', '+44d00m00s', frame='fk5',equinox='J1950.') # epoch = 1950
    GA = GA.transform_to(FK5(equinox='J2000'))
    theta = GA.separation(nsa_sc).radian
    #theta = np.sqrt((nsa_sc.ra.radian - GA.ra.radian)**2 + (nsa_sc.dec.radian - GA.dec.radian)**2)

    gamma = 2.
    V_a = 4380.*u.km/u.second # recession vel of Great Attractor from Mould+2000
    V_o = nsa.Z*c.to('km/s') + V_LG # recession velocities of the galaxies
    r_oa = np.sqrt(V_o**2 + V_a**2 - 2.*V_o*V_a*np.cos(theta))
    V_GA = V_fid*(np.cos(theta) + (V_o - V_a*np.cos(theta))/r_oa*(r_oa/V_a)**(1-gamma))
    
    #Add up all corrections
    V_cosmic = V_H + V_LG + V_infall +V_GA
    V = V_cosmic
    
    #Convert to supergalactic coordinates
    # using H0 = 74 to match Kim+2016 paper
    #H0 = 74.*u.km/u.second/u.Mpc
    H0 = 100.*u.km/u.second/u.Mpc #Kim's email
    SGX = V/H0*np.cos(nsa_sc.supergalactic.sgl.radian)*np.cos(nsa_sc.supergalactic.sgb.radian)
    SGY = V/H0*np.sin(nsa_sc.supergalactic.sgl.radian)*np.cos(nsa_sc.supergalactic.sgb.radian)
    SGZ = V/H0*np.sin(nsa_sc.supergalactic.sgb.radian)
    
    #Getting distance from virgo
    distance_to_virgo = 16.5*u.Mpc #Mpc
    SGX_Virgo = distance_to_virgo*np.cos(Virgo.supergalactic.sgl.radian)*np.cos(Virgo.supergalactic.sgb.radian)
    SGY_Virgo = distance_to_virgo*np.sin(Virgo.supergalactic.sgl.radian)*np.cos(Virgo.supergalactic.sgb.radian)
    SGZ_Virgo = distance_to_virgo*np.sin(Virgo.supergalactic.sgb.radian)

    # difference 
    DSGX = SGX - SGX_Virgo
    DSGY = SGY - SGY_Virgo
    DSGZ = SGZ - SGZ_Virgo

    dist = np.sqrt(DSGX**2 + DSGY**2 + DSGZ**2) 
    dist = dist / u.Mpc
    
    distarray = np.array(dist)
    vdistflag = distarray < 50 #What is the radius of virgo??

    VirgoTestSet = nsa[vdistflag]

    fits.writeto('VirgoCatalog.fits', VirgoTestSet, clobber=True)
