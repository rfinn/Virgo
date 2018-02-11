#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
from astropy import units as u

from virgoCommon import *

hyperfile = 'tables/hyperleda/HyperLedaVirgo.fits'
hleda = fits.getdata(hyperfile)

nsa = fits.getdata(nsa_file)

nsacat = SkyCoord(nsa.RA*u.degree,nsa.DEC*u.degree,frame='icrs')
hledacat = SkyCoord(hleda.al2000*u.hourangle,hleda.de2000*u.degree,frame='icrs')

# match Simard+2011 Table 1 to NSA
index,dist2d,dist3d = nsacat.match_to_catalog_sky(hledacat)

matchflag = dist2d.degree < 1./3600

matchedarray1=np.zeros(len(nsa),dtype=hleda.dtype)
matchedarray1[matchflag] = hleda[index[matchflag]]


diff = matchedarray1['vopt'] - nsa.Z*3.e5
plt.figure()
#plt.plot(nsa.Z[matchflag]*3.e5,diff,'k.')
plt.plot(matchedarray1['vopt'][matchflag],diff[matchflag],'k.')
plt.xlabel('NSA c*z (km/s)')
plt.ylabel('Hyperleda vopt - NSA c*z (km/s)')

cutoff = 200. # 300 km/s
flag = (abs(diff) > cutoff) & matchflag
jindex = np.arange(len(nsa.RA))
jindex = jindex[flag]
for i in jindex:
    print nsa.NSAID[i],nsa.RA[i],nsa.DEC[i],nsa.ISDSS[i],nsa.Z[i]*3.e5,hleda['vopt'][index[i]]
