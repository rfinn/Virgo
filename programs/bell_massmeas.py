#/Users/grudnick/anaconda/envs/python3env/bin/python

from astropy.table import Table
from matplotlib import pyplot as plt
import numpy as np
import os

tabledir = '/Users/grudnick/Work/Virgo_outskirts/Catalogs/v0-03Jul2020/'
plotdir = '/Users/grudnick/Work/Virgo_outskirts/Plots/'

maintab = Table.read(tabledir+'vf_north_v0_main.fits')
nsav1tab = Table.read(tabledir+'vf_north_v0_nsa.fits')
nsav0tab = Table.read(tabledir+'vf_north_v0_nsa_v0.fits')

#the indices for bands
fuvind = 0
nuvind = 1
uind = 2
gind = 3
rind = 4
iind = 5
zind =6

#calculate g-r color for v0 catalog
g_rv0 = nsav0tab['ABSMAG'][:,gind] - nsav0tab['ABSMAG'][:,rind]
g_rv1 = nsav1tab['SERSIC_ABSMAG'][:,gind] - nsav1tab['SERSIC_ABSMAG'][:,rind]

#correct to H0=70
#catalog absmag are listed as M - 5 log10 h
h=0.7
r_absmagv0 = nsav0tab['ABSMAG'][:,rind] + 5 * np.log10(h)
r_absmagv1 = nsav1tab['SERSIC_ABSMAG'][:,rind] + 5 * np.log10(h)

#calculate M/L using Bell et al. 2003
#v0
logMLrv0 = -0.306 + 1.097 * g_rv0
#convert to Kroupa (2002) IMF
logMLrv0 = logMLrv0 - 0.15

#v1
logMLrv1 = -0.306 + 1.097 * g_rv1
#convert to Kroupa (2002) IMF
logMLrv1 = logMLrv1 - 0.15

#compute L in r-band
Msunr = 4.67 # from Bell+2003, different from http://mips.as.arizona.edu/~cnaw/sun.html value of 4.65
#compute luminosity of sun in r-band
#flux of sun in r-band at 10pc in erg/s/cm^2/Hz
fsunr = 10**(-0.4 * (4.50 + 48.60))
#luminosity of sun in r-band in erg/s/Hz
pc_cgs = 3.086e18
Lsunr =  fsunr * 4. * np.pi * (10.0 * pc_cgs)**2

#measure luminosity of galaxy in units of solar luminosity
#flux of galaxy at 10pc 
fgalr_v0 =  10**(-0.4 * (r_absmagv0 + 48.60))
Lgalr_v0 =  fgalr_v0 * 4. * np.pi * (10.0 * pc_cgs)**2
Lgalr_v0_Lsun = Lgalr_v0 / Lsunr

fgalr_v1 =  10**(-0.4 * (r_absmagv1 + 48.60))
Lgalr_v1 =  fgalr_v1 * 4. * np.pi * (10.0 * pc_cgs)**2
Lgalr_v1_Lsun = Lgalr_v1 / Lsunr

#compute Mstar
Mstar_v0 = 10**logMLrv0 * Lgalr_v0_Lsun
Mstar_v1 = 10**logMLrv1 * Lgalr_v1_Lsun
#print(Mstar_v0)
    
plt.figure()
plt.scatter(np.log10(Lgalr_v0_Lsun),np.log10(Mstar_v0),c=g_rv0,vmin=-0.5, vmax=1.5,s=5)
cb = plt.colorbar()
cb.set_label('g-r')
ax = plt.gca()
plt.xlim(6.5,11.5)
plt.ylim(7.0,11.5)
plt.xlabel('log$_{10}~L_r$')
plt.ylabel('log$_{10}~M_\star$')
plt.title('NSA v0')
#plt.show()
plt.savefig(plotdir + 'l_mstar.v0.png')

plt.figure()
plt.scatter(np.log10(Lgalr_v1_Lsun),np.log10(Mstar_v1),c=g_rv1,vmin=-0.5, vmax=1.5,s=5)
cb = plt.colorbar()
cb.set_label('g-r')
ax = plt.gca()
plt.xlim(6.5,11.5)
plt.ylim(7.0,11.5)
plt.xlabel('log$_{10}~L_r$')
plt.ylabel('log$_{10}~M_\star$')
plt.title('NSA v1')
plt.savefig(plotdir + 'l_mstar.v1.png')
