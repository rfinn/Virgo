#!/usr/bin/env python

"""
Testing display parameters for images 
- all images are coming out grey in web pages.
- need to increase contrast
- copied test image from draco
"""

import os
import numpy as np
import glob
import sys

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table, Column
from astropy import wcs
from astropy.coordinates import SkyCoord

from astropy.visualization import simple_norm
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip
from astropy.time import Time

def display_image(image,percent=99,lowrange=False,mask=None,sigclip=False,csimage=False):
    lowrange=False
    # use inner 80% of image
    xdim,ydim = image.shape
    xmin = int(.2*xdim)
    xmax = int(.8*xdim)    
    ymin = int(.2*ydim)
    ymax = int(.8*ydim)    
    if sigclip:
        clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]
    if csimage:
        try:
            norm = simple_norm(clipped_data, stretch='asinh',min_percent=10,max_percent=90)
        except:
            print("error getting norm")
            norm = None
    elif lowrange:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    if norm == None:
        plt.imshow(image,cmap='gray_r',origin='lower')
    else:
        plt.imshow(image,cmap='gray_r',origin='lower', norm=norm)
    
