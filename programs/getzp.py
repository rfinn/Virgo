#!/usr/bin/env python

'''
potentially useful references

SDSS
https://astroquery.readthedocs.io/en/latest/api/astroquery.sdss.SDSSClass.html#astroquery.sdss.SDSSClass.query_photoobj

http://skyserver.sdss.org/dr8/en/help/docs/realquery.asp#cleanStars

SELECT TOP 10 u,g,r,i,z,ra,dec, flags_r 
FROM Star 
WHERE 
ra BETWEEN 180 and 181 AND dec BETWEEN -0.5 and 0.5
AND ((flags_r & 0x10000000) != 0)
-- detected in BINNED1
AND ((flags_r & 0x8100000c00a4) = 0)
-- not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
-- SATURATED, or BAD_COUNTS_ERROR
AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2))
-- not DEBLEND_NOPEAK or small PSF error
-- (substitute psfmagerr in other band as appropriate)
AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)
-- not INTERP_CENTER or not COSMIC_RAY

https://www.sdss.org/dr12/algorithms/photo_flags_recommend/


GAIA
https://gea.esac.esa.int/archive/documentation/GDR1/Data_processing/chap_cu5phot/sec_phot_calibr.html

https://www.cosmos.esa.int/web/gaia/dr2-known-issues


Pan-STARRS
https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/

https://panstarrs.stsci.edu/




from astroquery.sdss import SDSS


query = 'SELECT TOP 10 ra, dec, u,g,r,i,z, flags_r FROM Star WHERE (clean = 1) AND ra BETWEEN 180 and 181 AND dec BETWEEN -0.5 and 0.5 AND ((flags_r & 0x10000000) != 0) AND ((flags_r & 0x8100000c00a4) = 0) AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2)) AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)'

t = SDSS.query_sql(query, data_release=14)
'''

import argparse
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.io import fits

from astroquery.vizier import Vizier


def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=20,
                    maxsources=10000):
    """
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag',
                             'ymag', 'e_ymag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag)},
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/349/ps1")[0]




parser = argparse.ArgumentParser(description ='Run sextractor, get Pan-STARRS catalog, and then computer photometric ZP')
parser.add_argument('--image', dest = 'image', default = None, help = 'Image for ZP calibration')
parser.add_argument('--instrument', dest = 'instrument', default = None, help = 'HDI = h, KPNO mosaic = m, INT = i')
parser.add_argument('--nsigma', dest = 'nsigma', default = 2.0, help = 'number of std to use in iterative rejection of ZP fitting')
parser.add_argument('--d',dest = 'd', default ='~/github/HalphaImaging/astromatic', help = 'Locates path of default config files.  Default is ~/github/HalphaImaging/astromatic')
args = parser.parse_args()

os.system('cp ' +args.d + '/default.* .')




# Example query
# print(panstarrs_query(12.345, 67.89, 0.1))


#####
# Run Source Extractor on image to measure magnitudes
####
t = args.image.split('.fits')
froot = t[0]
if args.instrument == 'h':
    os.system('sex ' + args.image + ' -c default.sex.HDI -CATALOG_NAME ' + froot + '.cat')

elif args.instrument == 'i':
    os.system('sex ' + args.image + ' -c default.sex.INT -CATALOG_NAME ' + froot + '.cat')

elif args.instrument == 'm':
    os.system('sex ' + args.image + ' -c default.sex.HDI -CATALOG_NAME ' + froot + '.cat')



#####
# Read in Source Extractor catalog
####

secat_filename = froot+'.cat'
secat = fits.getdata(secat_filename,2)

#####
# get max/min RA and DEC for the image
####

minRA = min(secat['ALPHA_J2000'])
maxRA = max(secat['ALPHA_J2000'])
minDEC = min(secat['DELTA_J2000'])
maxDEC = max(secat['DELTA_J2000'])
centerRA = 0.5*(minRA + maxRA)
centerDEC = 0.5*(minDEC + maxDEC)
radius = np.sqrt((maxRA- centerRA)**2 + (maxDEC - centerDEC)**2)
print(radius)

radius = max((maxRA - minRA), (maxDEC - minDEC))


#####
# get Pan-STARRS catalog over the same region
####

pan = panstarrs_query(centerRA, centerDEC, radius)

#####
# match Pan-STARRS1 data to Source Extractor sources
# remove any objects that are saturated or non-linear in our r-band image
####

pancoords = SkyCoord(pan['RAJ2000'],pan['DEJ2000'],frame='icrs')
secoords = SkyCoord(secat['ALPHA_J2000']*u.degree,secat['DELTA_J2000']*u.degree,frame='icrs')

index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)

# only keep matches with matched RA and Dec w/in 5 arcsec
matchflag = dist2d.degree < 5./3600


matchedarray1=np.zeros(len(pancoords),dtype=secat.dtype)
matchedarray1[matchflag] = secat[index[matchflag]]

fitflag = matchflag & (matchedarray1['FLAGS'] == 0)  & (pan['rmag'] > 14.) & (pan['rmag'] < 17.) & (matchedarray1['CLASS_STAR'] > 0.95)

#####
# Calculate Johnson R
####
R = pan['rmag'] + (-0.153)*(pan['rmag']-pan['imag']) - 0.117


#####
# Show location of residuals
#####



#####
# Solve for the zeropoint
####

# plot Pan-STARRS r mag on x axis, observed R-mag on y axis


def fitzp():
    plt.figure(figsize=(8,6))
    flag = fitflag
    plt.plot(pan['rmag'][flag],matchedarray1['MAG_AUTO'][flag],'bo')
    plt.errorbar(pan['rmag'][flag],matchedarray1['MAG_AUTO'][flag],xerr= pan['e_rmag'][fitflag],yerr=matchedarray1['MAGERR_AUTO'][fitflag],fmt='none')
    plt.plot(pan['rmag'][flag],matchedarray1['MAG_BEST'][flag],'ro',label='MAG_BEST')
    plt.plot(pan['rmag'][flag],matchedarray1['MAG_PETRO'][flag],'go',label='MAG_PETRO')
    plt.plot(pan['rmag'][flag],matchedarray1['MAG_ISO'][flag],'ko',label='MAG_ISO')
    plt.xlabel('Pan-STARRS r',fontsize=16)
    plt.ylabel('SE R-band MAG_AUTO',fontsize=16)
    c = np.polyfit(pan['rmag'][flag],matchedarray1['MAG_AUTO'][flag],1)
    xl = np.linspace(14,17,10)
    yl = np.polyval(c,xl)
    plt.plot(xl,yl,'k--')
    #plt.plot(xl,1.2*yl,'k:')
    print(c)
    
    yfit = np.polyval(c,pan['rmag'])
    residual = np.zeros(len(fitflag))
    residual[flag] = (yfit[flag] - matchedarray1['MAG_AUTO'][flag])/yfit[flag]
    bestc = np.array([0,0],'f')
    plt.figure()
    plt.scatter(matchedarray1['X_IMAGE'][flag],matchedarray1['Y_IMAGE'][flag],c = abs(residual[flag]))
    plt.colorbar()
    delta = 100.
    
    x = pan['rmag'][fitflag]
    x = R[fitflag]
    y = matchedarray1['MAG_AUTO'][fitflag]
    while delta > 1.e-3:
    
        plt.figure(figsize=(8,8))
        plt.subplot(2,1,1)
    
        print('number of points retained = ',sum(flag))
        plt.plot(x,y,'bo',label='MAG_AUTO')
        plt.xlabel('Pan-STARRS r',fontsize=16)
        plt.ylabel('SE R-band MAG',fontsize=16)
        c = np.polyfit(x,y,1)
        xl = np.linspace(14,17,10)
        yl = np.polyval(c,xl)
        s = 'fit: y = %.2f PAN + %.2f'%(c[0],c[1])
        plt.plot(xl,yl,'k--',label=s)
        plt.legend()
        print(c)

        yfit = np.polyval(c,x)
        residual = (yfit - y)/yfit 
        plt.subplot(2,1,2)
        s = 'std = %.4f'%(np.std(residual))
        plt.plot(x,residual, 'ko',label=s)
        plt.xlabel('Pan-STARRS r',fontsize=16)
        plt.ylabel('YFIT - SE R-band MAG_AUTO',fontsize=16)
        plt.legend()
        plt.axhline(y=0,color='r')
    
        # check for convergence
        print(bestc[1],c[1])
        delta = abs(bestc[1] - c[1])
        bestc = c
        flag =  (abs(residual) < 2.0*np.std(residual))
        x = x[flag]
        y = y[flag]
    return x,y


