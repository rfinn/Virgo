#!/usr/bin/env python


"""
GOAL: 

* create a SFR-Mstar relation, with cutouts of halpha morphologies


APPROACH:

- could break SFR-Mstar plane into subplots, and plot an image in each subplot


"""


from matplotlib import pyplot as plt


import os
import sys
from astropy.table import Table,hstack
from astropy.io import fits
import numpy as np

from scipy.stats import scoreatpercentile
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch, PercentileInterval
from astropy.visualization import ImageNormalize
from astropy.visualization import LinearStretch,SinhStretch
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip



homedir = os.getenv("HOME")
sys.path.append(os.path.join(homedir,'github/Virgo/programs/'))

from readtablesv2 import vtables
from build_web_cutouts2 import display_image
H0 = 70.




cutoutdir = '/data-pool/Halpha/halphagui-output-20230713/cutouts/'

#cutoutdir = os.path.join(homedir,'/research/Virgo/plots/')

def construct_dirname(prefix,tel,dateobs,pointing):
    return f"{prefix}-{tel}-{dateobs}-{pointing}"
    pass
def get_imagenames(prefix,tel,dateobs,pointing):
    """get CS halpha cutout name """
    dirname = construct_dirname(prefix,tel,dateobs,pointing)
    cs_imname = os.path.join(dirname,dirname+'-CS.fits')
    
    return cs_imname
        
class vplots(vtables):
    def plot_sfr_mstar(self):
        """plot SFR vs mstar   """

        plt.figure(figsize=(8,6))
        plt.scatter(self.magphys['logMstar'],self.magphys['logSFR'])
        plt.xlabel('$\log_{10}(M_\star/M_\odot)$',fontsize=22)
        plt.xlabel('$\log_{10}(SFR/M_\odot/yr)$',fontsize=22)        
    def get_hadetect(self,snr=3):
        """ define halpha detection flag """

        # require SNR > 3
        self.hadetect = (self.halpha['HF_TOT']/self.halpha['HF_TOT_ERR']) > snr
        pass
    def plot_sfr_mstar_hacutout(self,xmin=8,xmax=10,ymin=-2,ymax=1,nbins=5):
        """plot halpha  cutouts   """
        
        
        dx = (xmax-xmin)/nbins
        xedge = np.linspace(xmin,xmax,nbins)

        dy = (ymax-ymin)/nbins
        yedge = np.linspace(ymin,ymax,nbins)
        # step over x
        plt.figure(figsize=(12,10))
        nplot = 1
        for i in range(nbins-1):
            print('nplot = ',nplot)
            # step over y
            for j in range(nbins-1):
                plt.subplot(nbins,nbins,nplot)
                xlo = xedge[j]
                xhi = xedge[j+1]
                ylo = yedge[i]
                yhi = yedge[i+1]
                print(xlo,xhi,ylo,yhi)

                # select galaxies in dx and dy range
                flag = (self.magphys['logMstar'] > xlo) & \
                  (self.magphys['logMstar'] < xhi) & \
                  (self.magphys['logSFR'] > ylo) & \
                  (self.magphys['logSFR'] > yhi)
                # select halpha detections in this bin
                flag = flag  & self.hadetect
                if np.sum(flag) > 1:
                    vfids = np.arange(len(self.magphys))[flag]

                    # now select one vfid at random
                    vfindices = np.random.choice(vfids,size=len(vfids))
                    # this will check all galaxies in the mstar-sfr bin
                    # in a random order
                    # once one image is found, break out of loop and
                    # move on to the next mstar-sfr bin
                    for k in vfindices:
                        imname = get_imagenames(v.main['prefix'][k],\
                                            v.halpha['TEL'][k],\
                                            v.halpha['DATE-OBS'][k],\
                                            v.halpha['POINTING'][k].split('-')[-1])
                        maskname = imname.replace('-CS.fits','-R-mask.fits')
                        imname = os.path.join(cutoutdir,imname)
                        maskname = os.path.join(cutoutdir,maskname)                        
                        print(imname)
                        print(maskname)
                        # check if image exist
                        if os.path.exists(imname):
                            print("FOUND AN IMAGE!!!!")
                            print()
                            imdata = fits.getdata(imname)
                            maskdata = fits.getdata(maskname)
                            display_image(imdata,mask=maskname)
                            # remove ticks on axes
                            plt.xticks([],[])
                            plt.yticks([],[])
                            nplot += 1
                            break
                    # increment plot in the subplot
                    print("did not find an image")
                    nplot += 1
                else:
                    print("no galaxies in this bin")
                    nplot += 1

        plt.savefig('sfr-mstar-hamorph.png')
        plt.savefig('sfr-mstar-hamorph.pdf')            

        


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v2/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_v2_', help = 'prefix for tables; default is vf_v2')                               
    args = parser.parse_args()

    if args.tabledir.startswith('/home/rfinn/'):
        homedir = os.getenv("HOME")
        args.tabledir = args.tabledir.replace('/home/rfinn',homedir)
    v = vplots(args.tabledir,args.tableprefix)
    v.read_all()
    v.get_hadetect(snr=3)
    v.plot_sfr_mstar_hacutout()
    
