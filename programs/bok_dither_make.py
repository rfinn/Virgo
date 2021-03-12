from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle  
from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt

def bok_dither_make(exptime,objroot,band,ra,dec,outfile,dithsize = 50,nexp=1):

    '''Written by Gregory Rudnick 10 March 2021
    
    PURPOSE: 

    Generate a list of dithers for Bok telescope observations.  This
    will use a 5 point dither pattern

    INPUT PARAMETERS:

    exptime: the exposure time in seconds

    objroot: the root of the object name.  This will have a dither identifier appended

    band:  'r' or 'Ha+4nm'

    ra,dec: The start of the pointing in decimal units

    dithsize: the dimensions of the side of the box containing the dithers (arcsec).

    nexp: the number of exposures at each dither position

    outfile: The name of the script written by the code

    OUTPUT

    Will write a file with the dithers

    Make a plot of dither positions

    EXAMPLE CALLING SEQUENCE

    import bok_dither_make as bkd
    bkd.bok_dither_make(240.0,'test','r',23.456,9.001,'test.out')

    '''

    #dither offsets relative to *initial* pointing
    radith = np.array([0.0,
                           dithsize/2.0,
                           -dithsize/2.0,
                           -dithsize/2.0,
                           dithsize/2.0])

    decdith = np.array([0.0,
                            dithsize/2.0,
                            dithsize/2.0,
                           -dithsize/2.0,
                           -dithsize/2.0])

    #convert to degrees
    radith /= 3600.0
    decdith /= 3600.0
    
    #apply dither
    radithdeg = ra + radith
    decdithdeg = dec + decdith

    #convert to sexigesimal
    RAdithhms = Angle(radithdeg,unit='deg').to_string(unit=u.hour,sep='')
    DECdithhms = Angle(decdithdeg,unit='deg').to_string(unit=u.degree,sep='')

    #open output file
    fo = open(outfile, "w")

    for idith,val in enumerate(radith):
        dithname = objroot + '_dith' + str(idith)

        #make strings for positions
        rastr = "{:0>9.2f}".format(float(RAdithhms[idith]))
        #DECdithhms[idith]
        if float(DECdithhms[idith]) < 0:
            decstr = "-{:0>8.1f}".format(abs(float(DECdithhms[idith])))
        else:
            decstr = "+{:0>8.1f}".format(float(DECdithhms[idith]))
            
        print(rastr,decstr)
        #print(idith,RAdithhms[idith],DECdithhms[idith])
        fo.write('obs {} object {} {} {} {} {} 2000.0\n'.format(exptime,dithname,nexp,band,rastr,decstr))

    fo.close()


    
