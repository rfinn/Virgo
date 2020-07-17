#/Users/grudnick/anaconda/envs/python3env/bin/python

from astropy.table import Table
from matplotlib import pyplot as plt
import numpy as np
import os
import matplotlib
matplotlib.use('MACOSX')

tabledir = '/Users/grudnick/Work/Google Drive/VirgoFilamentCollaboration/vf-tables/north-only/v0-03Jul2020/'
plotdir = '/Users/grudnick/Work/Virgo_outskirts/Plots/'
tabdir = '/Users/grudnick/Work/Virgo_outskirts/Rfinn_github/Virgo/tables/'

maintab = Table.read(tabledir+'vf_north_v0_main.fits')
nsav1tab = Table.read(tabledir+'vf_north_v0_nsa.fits')
nsav0tab = Table.read(tabledir+'vf_north_v0_nsa_v0.fits')
envtab = Table.read(tabledir+'vf_north_v0_main_env_prop_H0_74_0.fits')

#the indices for bands
fuvind = 0
nuvind = 1
uind = 2
gind = 3
rind = 4
iind = 5
zind =6

#use the vcosmic velocity to derive the absolute magnitude from the observed flux

#make observed magnitudes from the catalog nmgy fluxes
nsav0tab.add_column(22.5 - 2.5 * np.log10(nsav0tab['NMGY'][:,gind]),name = 'mg')
nsav0tab.add_column(22.5 - 2.5 * np.log10(nsav0tab['NMGY'][:,rind]),name = 'mr')
nsav1tab.add_column(22.5 - 2.5 * np.log10(nsav1tab['SERSIC_NMGY'][:,gind]),name = 'mg')
nsav1tab.add_column(22.5 - 2.5 * np.log10(nsav1tab['SERSIC_NMGY'][:,rind]),name = 'mr')

#calculate g-r color for v0 catalog
#g_rv0 = nsav0tab['ABSMAG'][:,gind] - nsav0tab['ABSMAG'][:,rind]
#g_rv1 = nsav1tab['SERSIC_ABSMAG'][:,gind] - nsav1tab['SERSIC_ABSMAG'][:,rind]
g_rv0 = nsav0tab['mg'] - nsav0tab['mr'] 
g_rv1 = nsav1tab['mg'] - nsav1tab['mr'] 

#compute absmag in R-band using Vcosmic and H0=74
h=0.74
dgal = envtab['Vcosmic'] / (100.0 * h) * 1.e6    #distance in pc
r_absmagv0 = nsav0tab['mr'] - 5 * np.log10(dgal / 10.0)
r_absmagv1 = nsav1tab['mr'] - 5 * np.log10(dgal / 10.0)

#correct to H0=70
#catalog absmag are listed as M - 5 log10 h
#h=0.7
#r_absmagv0 = nsav0tab['ABSMAG'][:,rind] + 5 * np.log10(h)
#r_absmagv1 = nsav1tab['SERSIC_ABSMAG'][:,rind] + 5 * np.log10(h)

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
#first compute the luminosity of the sun in the r-band
Msunr = 4.67 # from Bell+2003, different from http://mips.as.arizona.edu/~cnaw/sun.html value of 4.65
#compute luminosity of sun in r-band
#flux of sun in r-band at 10pc in erg/s/cm^2/Hz
fsunr = 10**(-0.4 * (4.50 + 48.60))
#luminosity of sun in r-band in erg/s/Hz
pc_cgs = 3.086e18
Lsunr =  fsunr * 4. * np.pi * (10.0 * pc_cgs)**2

#measure luminosity of galaxy in units of solar luminosity
#flux of galaxy at 10pc in cgs units
fgalr_v0 =  10**(-0.4 * (r_absmagv0 + 48.60))
#luminosity in erg/s/Hz
Lgalr_v0 =  fgalr_v0 * 4. * np.pi * (10.0 * pc_cgs)**2
#luminosity in units of Lsol
Lgalr_v0_Lsun = Lgalr_v0 / Lsunr
lgLgalr_v0_Lsun = np.log10(Lgalr_v0_Lsun)

#flux of galaxy at 10pc in cgs units
fgalr_v1 =  10**(-0.4 * (r_absmagv1 + 48.60))
#luminosity in erg/s/Hz
Lgalr_v1 =  fgalr_v1 * 4. * np.pi * (10.0 * pc_cgs)**2
#luminosity in units of Lsol
Lgalr_v1_Lsun = Lgalr_v1 / Lsunr
lgLgalr_v1_Lsun = np.log10(Lgalr_v1_Lsun)

#compute Mstar
Mstar_v0 = 10**logMLrv0 * Lgalr_v0_Lsun
Mstar_v1 = 10**logMLrv1 * Lgalr_v1_Lsun
lgMstar_v0 = np.log10(10**logMLrv0 * Lgalr_v0_Lsun)
lgMstar_v1 = np.log10(10**logMLrv1 * Lgalr_v1_Lsun)

#reset bad absolute magnitude values to 0
ibadmagv0 = np.where((nsav0tab['ABSMAG'][:,rind] > -1) | (nsav0tab['ABSMAG'][:,gind] > -1))
ibadmagv1 = np.where((nsav1tab['SERSIC_ABSMAG'][:,rind] > -1) | (nsav1tab['SERSIC_ABSMAG'][:,gind] > -1 ))
g_rv0[ibadmagv0] = -99
lgMstar_v0[ibadmagv0] = -99
logMLrv0[ibadmagv0] = -99
lgLgalr_v0_Lsun[ibadmagv0] = -99

g_rv1[ibadmagv1] = -99
lgMstar_v1[ibadmagv1] = -99
logMLrv1[ibadmagv1] = -99
lgLgalr_v1_Lsun[ibadmagv1] = -99

#make new tables
nsav0tab.add_column(g_rv0,name = 'g-r')
nsav0tab.add_column(logMLrv0,name = 'logM_L_r')
nsav0tab.add_column( lgLgalr_v0_Lsun,name = 'logL_r')
nsav0tab.add_column( lgMstar_v0,name = 'logMstar')
newtab_v0 = Table([maintab['VFID'], nsav0tab['g-r'],nsav0tab['logM_L_r'],nsav0tab['logL_r'],nsav0tab['logMstar']])
#print(newtab_v0['logMstar'])
newtab_v0.write(tabledir + 'vf_north_v0_nsa_v0_bellmasses.fits', overwrite=True)
                      
nsav1tab.add_column(g_rv1,name = 'g-r')
nsav1tab.add_column(logMLrv1,name = 'logM_L_r')
nsav1tab.add_column( lgLgalr_v1_Lsun,name = 'logL_r')
nsav1tab.add_column( lgMstar_v1,name = 'logMstar')
newtab_v1 = Table([maintab['VFID'],nsav1tab['g-r'],nsav1tab['logM_L_r'],nsav1tab['logL_r'],nsav1tab['logMstar']])
#print(newtab_v1['logMstar'])
newtab_v1.write(tabledir + 'vf_north_v0_nsa_bellmasses.fits', overwrite=True)
                     

    
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

plt.figure()
plt.scatter(lgMstar_v0,lgMstar_v1,c=g_rv1,vmin=-0.5, vmax=1.5,s=5)
cb = plt.colorbar()
cb.set_label('g-r')
ax = plt.gca()
plt.xlim(6.5,11.5)
plt.ylim(7.0,11.5)
plt.xlabel('log$_{10}~M_\star v0r$')
plt.ylabel('log$_{10}~M_\star v1$')
plt.title('NSA v1')
plt.savefig(plotdir + 'mstar_comp.png')

