#!/usr/bin/env python

from astropy.table import Table
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord

from matplotlib import pyplot as plt
import numpy as np
v8file = '/home/rfinn/research/Virgo/V8k-catalog/EDDtable05Aug2021150224.xml'
v8file = '/home/rfinn/research/Virgo/V8k-catalog/EDDtable05Aug2021.csv'
v8file = '/home/rfinn/research/Virgo/V8k-catalog/EDDtable05Aug2021-2.fits'
colnames = ['pgc','Dist','Nd','DM','eDM',\
            'C','T','L','M','S','N','H','I','F',\
            'DM2','eD2','SNIa','Ns','DMsn','DMsp','eDsp',\
            'DM6d','eD6d','Mt',\
            'RAJ','DeJ','Glon','Glat','SGL','SGB',\
            'Ty','Asf','Btot','Ks',\
            'Vhel','Vgsr','Vls','Vcmb','Vmod',\
            'Name','Nest','Ndgp',\
            'DMgp','eDgp','Dgp','Abell','GroupName','NV',\
            '1PGC','Glongp','Glatgp','SGLgp','SGBgp','lgLgp',\
            'cf','sigp','R2t',\
            'Vhgp','Vggp','Vlsgp','Vcgp','Vmgp','Vrms',\
            'bwMass12','L_Mass12','LDC','HDC','2M++','MKgp','Icnt',\
            'l','b','sgL2','sgB2',\
            'absB','Diam','BAratio','PA','V_GSR',\
            'sgX','sgY','sgZ',\
            'objname','GrpID','Filament',\
            'al1950','de1950','test1','test2']
#v8 = ascii.read(v8file,delimiter=',',data_start=0,comment='\s*#',names=colnames)
#v8 = Table.read(v8file,delimiter=',',data_start=0,comment='\s*#',names=colnames)
v8 = Table.read(v8file,format='fits')#,delimiter=',')
flag = (v8['Glon'] != '        ') & (v8['Dist'] != '      ')
v8 = v8[flag]
c1 = SkyCoord(l=v8['Glon'],b=v8['Glat'],frame="galactic",unit='deg')

c1_radec = c1.icrs

# read in VF catalog

vmain = Table.read('/home/rfinn/research/Virgo/tables-north/v1/vf_north_v1_main.fits',format='fits')
steer = Table.read('/home/rfinn/research/Virgo/tables-north/v1/vf_north_v1_steer17.fits',format='fits')
c2 = SkyCoord(vmain['RA'],vmain['DEC'],frame='icrs',unit='deg')
# match by RA and DEC


idx, d2d, d3d = c2.match_to_catalog_sky(c1_radec)
matchflag = d2d < 15./3600*u.deg


# keep only close matched

vmain = vmain[matchflag]
steer = steer[matchflag]
matched_v8k = v8[idx][matchflag]

v8k_coords_matched2_vf = c1_radec[idx][matchflag]
# compare distances

plt.figure()
flag2 = (matched_v8k['Dist'] != '      ') & (steer['Steerflag'])
plt.plot(steer['Dmedian'][flag2],np.array(matched_v8k['Dist'][flag2],'f'),'b.')
xl = np.linspace(0,70,100)
plt.plot(xl,xl,'k--')
plt.xlabel('Median Distance from Steer+2017 (Mpc)')
plt.ylabel('Distance from Cosmicflows3 (Mpc)')
plt.savefig('distance-comparison-v8k-vf.png')
# calculate dispersion

diff = steer['Dmedian'][flag2] - np.array(matched_v8k['Dist'][flag2],'f')
std = np.std(diff)
print('difference in distances = {:.2f} +/- {:.2f} Mpc'.format(np.mean(diff),std))

# in log scale
logdiff = np.log10(steer['Dmedian'][flag2]) - np.log10(np.array(matched_v8k['Dist'][flag2],'f'))
std = np.std(logdiff)
print('log diff in distances = {:.3f} +/- {:.3f} (logMpc) amd'.format(np.mean(logdiff),std))





