#!/usr/bin/env python

from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
input = '/home/rfinn/research/Virgo/supersample/clean_sample.fits'
input = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink.fits'


class catalog():
    def __init__(self,input_table):
        self.cat = fits.getdata(input_table)
        print('got here')
        self.cat = Table(self.cat)
        self.cat.rename_column('RA-COMBINED','RA')
        self.cat.rename_column('DEC-COMBINED','DEC')        
        self.cat.rename_column('HL-AGC-NSA-VEL','VEL')
        self.define_ngcfilament()
    def define_ngcfilament(self):
        
        radec = (self.cat['RA'] > 192.) & (self.cat['RA'] < 209) & (self.cat['DEC'] > 0.) & (self.cat['DEC'] < 50.) 
        radec_flag = radec & (self.cat['DEC'] >(2*(self.cat['RA'] - 205.) + 20) ) & (self.cat['DEC'] < (2*(self.cat['RA'] - 205.) + 55))
        filament = radec_flag & (self.cat['VEL'] >2000.) & (self.cat['VEL'] < 3238.)
        nsa_flag = (self.cat['VEL'] >1234.) & (self.cat['VEL'] < 3976.)
        #mass_flag = (jmass.MSTAR_50 > 8.3) & (jmass.MSTAR_50 < 10.2)

        # used ~(max vel of INT 197 - 250) for low-z end of filter gap
        # used ~(min vel of INT 227 + 250) for high-z end of filter gap 
        #INTvflag =   (self.cat['VEL'] < 2100.) #| (self.cat['VEL'] > 2700.) #


        #gas_flag = COflag | HIflag
        self.NGCfilament = filament

    def plot_ngcfilament(self):
        plt.figure(figsize=(4,8))

        plt.scatter(self.cat['RA'][self.NGCfilament],self.cat['DEC'][self.NGCfilament],c = self.cat['VEL'][self.NGCfilament],s=2)
        #plt.axis([204,210,35,45])        
        plt.gca().invert_xaxis()
        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')

        sitelle_fov = 11#arcmin
        rect = Rectangle((208,40),sitelle_fov/60,sitelle_fov/60,fill=False, color='k')
        plt.gca().add_artist(rect)
        rect = Rectangle((206.9,36.1),sitelle_fov/60,sitelle_fov/60,fill=False, color='k')
        plt.gca().add_artist(rect)
        plt.colorbar(label='vr (km/s)')

if __name__ == '__main__':
    c = catalog(input)
