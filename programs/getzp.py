#!/usr/bin/env python

'''
USAGE:

from within ipython:

%run ~/github/Virgo/programs/getzp.py --image pointing031-r.coadd.fits --instrument i --filter r

then:
fitzp()

The y intercept is -1*ZP


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
    FOUND THIS ON THE WEB 
    https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/

    
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


class getzp():
    def __init__(self, image, instrument='h', filter='r', astromatic_dir = '~/github/HalphaImaging/astromatic/'):

        self.image = image
        self.astrodir = astromatic_dir
        self.instrument = instrument
        self.filter = filter
    def getzp(self):
        self.runse()
        self.get_panstarrs()
        self.match_coords()
        self.fitzp()
    def runse(self):
        #####
        # Run Source Extractor on image to measure magnitudes
        ####

        os.system('cp ' +args.d + '/default.* .')
        t = self.image.split('.fits')
        froot = t[0]
        if self.instrument == 'h':
            os.system('sex ' + args.image + ' -c default.sex.HDI -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0')

        elif self.instrument == 'i':
            os.system('sex ' + args.image + ' -c default.sex.INT -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0')

        elif self.instrument == 'm':
            os.system('sex ' + args.image + ' -c default.sex.HDI -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0')

        # clean up SE files
        # skipping for now in case the following command accidentally deletes user files
        # os.system('rm default.* .')


        ###################################
        # Read in Source Extractor catalog
        ###################################

        secat_filename = froot+'.cat'
        self.secat = fits.getdata(secat_filename,2)

        ###################################
        # get max/min RA and DEC for the image
        ###################################

        minRA = min(self.secat['ALPHA_J2000'])
        maxRA = max(self.secat['ALPHA_J2000'])
        minDEC = min(self.secat['DELTA_J2000'])
        maxDEC = max(self.secat['DELTA_J2000'])
        self.centerRA = 0.5*(minRA + maxRA)
        self.centerDEC = 0.5*(minDEC + maxDEC)
        #radius = np.sqrt((maxRA- centerRA)**2 + (maxDEC - centerDEC)**2)
        #print(radius)

        # NOTE - radius is with width of the rectangular
        # search area.  This is not a circular search, as I orginally thought.
        self.radius = max((maxRA - minRA), (maxDEC - minDEC))
    def get_panstarrs(self):

        ###################################
        # get Pan-STARRS catalog over the same region
        ###################################

        self.pan = panstarrs_query(self.centerRA, self.centerDEC, self.radius)
    def match_coords(self):
        ###################################
        # match Pan-STARRS1 data to Source Extractor sources
        # remove any objects that are saturated or non-linear in our r-band image
        ###################################

        pancoords = SkyCoord(self.pan['RAJ2000'],self.pan['DEJ2000'],frame='icrs')
        secoords = SkyCoord(self.secat['ALPHA_J2000']*u.degree,self.secat['DELTA_J2000']*u.degree,frame='icrs')

        index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)

        # only keep matches with matched RA and Dec w/in 5 arcsec
        matchflag = dist2d.degree < 5./3600


        self.matchedarray1=np.zeros(len(pancoords),dtype=self.secat.dtype)
        self.matchedarray1[matchflag] = self.secat[index[matchflag]]

        ###################################
        # remove any objects that are saturated, have FLAGS set, galaxies, must have 14 < r < 17 according to Pan-STARRS
        ###################################


        fitflag = matchflag & (self.matchedarray1['FLAGS'] == 0)  & (self.pan['rmag'] > 14.) & (self.pan['rmag'] < 17.) & (self.matchedarray1['CLASS_STAR'] > 0.95)

        if self.filter == 'R':
            ###################################
            # Calculate Johnson R
            ###################################
            self.R = self.pan['rmag'] + (-0.153)*(self.pan['rmag']-self.pan['imag']) - 0.117
        else:
            self.R = self.pan['rmag']
        ###################################
        # Show location of residuals
        ###################################
    def plot_fitresults(self):#, polyfit_results = self.bestc):
        # plot best-fit results
        yfit = np.polyval(polyfit_results,self.x)
        residual = (yfit - self.y)/yfit 
        plt.figure(figsize=(8,8))

        plt.subplot(2,1,1)
        plt.plot(x,y,'bo',label='MAG_AUTO')
        plt.xlabel('Pan-STARRS r',fontsize=16)
        plt.ylabel('SE R-band MAG',fontsize=16)
        xl = np.linspace(14,17,10)
        yl = np.polyval(self.polyfit_results,xl)
        s = 'fit: y = %.2f PAN + %.2f'%(polyfit_results[0],polyfit_results[1])
        plt.plot(xl,yl,'k--',label=s)
        plt.legend()
        
        plt.subplot(2,1,2)
        s = 'std = %.4f'%(np.std(residual))
        plt.plot(self.x,residual, 'ko',label=s)
        plt.xlabel('Pan-STARRS r',fontsize=16)
        plt.ylabel('YFIT - SE R-band MAG_AUTO',fontsize=16)
        plt.legend()
        plt.axhline(y=0,color='r')

    def fitzp(plotall=False):
        ###################################
        # Solve for the zeropoint
        ###################################

        # plot Pan-STARRS r mag on x axis, observed R-mag on y axis
        flag = self.fitflag
        c = np.polyfit(self.pan['rmag'][flag],self.matchedarray1['MAG_AUTO'][flag],1)
        if plotall:
            plt.figure(figsize=(8,6))

            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_AUTO'][flag],'bo')
            plt.errorbar(self.pan['rmag'][flag],self.matchedarray1['MAG_AUTO'][flag],xerr= self.pan['e_rmag'][flag],yerr=self.matchedarray1['MAGERR_AUTO'][flag],fmt='none')
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_BEST'][flag],'ro',label='MAG_BEST')
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_PETRO'][flag],'go',label='MAG_PETRO')
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_ISO'][flag],'ko',label='MAG_ISO')
            plt.xlabel('Pan-STARRS r',fontsize=16)
            plt.ylabel('SE R-band MAG_AUTO',fontsize=16)

            xl = np.linspace(14,17,10)
            yl = np.polyval(c,xl)
            plt.plot(xl,yl,'k--')
            #plt.plot(xl,1.2*yl,'k:')
            #print(c)
    
        yfit = np.polyval(c,self.pan['rmag'])
        residual = np.zeros(len(flag))
        residual[flag] = (yfit[flag] - self.matchedarray1['MAG_AUTO'][flag])/yfit[flag]
        self.bestc = np.array([0,0],'f')
        
        plt.figure()
        plt.scatter(self.matchedarray1['X_IMAGE'][flag],self.matchedarray1['Y_IMAGE'][flag],c = abs(residual[flag]))
        plt.colorbar()
        delta = 100.
        
        #self.x = pan['rmag'][fitflag]
        self.x = self.R[flag]
        self.y = self.matchedarray1['MAG_AUTO'][flag]
        while delta > 1.e-3:
            c = np.polyfit(self.x,self.y,1)
            print('number of points retained = ',sum(flag))
            yfit = np.polyval(c,x)
            self.residual = (yfit - y)/yfit 
            if plotall:
                self.plot_fitresults(polyfit_results = c)
    
            # check for convergence
            print(bestc[1],c[1])
            delta = abs(bestc[1] - c[1])
            self.bestc = c
            flag =  (abs(residual) < 2.0*np.std(residual))
            self.x = x[flag]
            self.y = y[flag]
        self.plot_fitresults()
        
        
    def update_header(self):
        print('working on this')
        # add best-fit ZP to image header



if __name__ == '__main__':


    parser = argparse.ArgumentParser(description ='Run sextractor, get Pan-STARRS catalog, and then computer photometric ZP\n \n from within ipython: \n %run ~/github/Virgo/programs/getzp.py --image pointing031-r.coadd.fits --instrument i \n\n then:\n x,y = fitzp() \n \n The y intercept is -1*ZP. \n \n x and y data are returned in case you want to make additional plots.', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--image', dest = 'image', default = 'test.coadd.fits', help = 'Image for ZP calibration')
    parser.add_argument('--instrument', dest = 'instrument', default = None, help = 'HDI = h, KPNO mosaic = m, INT = i')
    parser.add_argument('--filter', dest = 'filter', default = 'R', help = 'filter (R or r; use r for Halpha)')
    parser.add_argument('--nsigma', dest = 'nsigma', default = 2.0, help = 'number of std to use in iterative rejection of ZP fitting')
    parser.add_argument('--d',dest = 'd', default ='~/github/HalphaImaging/astromatic', help = 'Locates path of default config files.  Default is ~/github/HalphaImaging/astromatic')
    args = parser.parse_args()

    zp = getzp(args.image, instrument=args.instrument, filter=args.filter, astromatic_dir = args.d)
    zp.getzp()



