from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle  
from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt

def bok_dither_make(objroot,ra,dec,dithsize = 50,nexp=1, nloop = 1):

    '''Written by Gregory Rudnick 10 March 2021
    
    PURPOSE: 

    Generate a list of dithers for Bok telescope observations.  This
    will use a 5 point dither pattern

    Current version hardcodes output for two filters

    INPUT PARAMETERS:

    objroot: the root of the object name.  This will have a dither identifier appended

    ra,dec: The start of the pointing in decimal units

    dithsize: the dimensions of the side of the box containing the dithers (arcsec).

    nexp: the number of exposures at each dither position

    nloop: the number of repeats of the dither pattern, with a dithsize/4 offset between each set

    OUTPUT

    Will write a file with the dithers

    Make a plot of dither positions

    EXAMPLE CALLING SEQUENCE

    import bok_dither_make as bkd
    bkd.bok_dither_make(240.0,'test',23.456,9.001)

    TODO:

    build possibility to repeat dither pattern with offsets in center between dither patterns

    '''

    #dither offsets relative to *initial* pointing
    #include small additions to some offset points to make sure that we aren't
    #offsetting directly along rows and column
    dithpad = 5.0

    #what is the offset between adjacent loops in arcsec
    loopoffset = 20.0
    
    radith = np.array([])
    decdith = np.array([])
    
    radithbase = np.array([0.0,
                               dithsize/2.0,
                               -dithsize/2.0,
                               -dithsize/2.0 - dithpad,
                               dithsize/2.0 + dithpad])

    decdithbase = np.array([0.0,
                                dithsize/2.0,
                                dithsize/2.0 + dithpad,
                                -dithsize/2.0,
                                -dithsize/2.0 - dithpad])
    
    for iloop in range (nloop):
        
        radith = np.append(radith, radithbase + iloop * loopoffset)

        decdith = np.append(decdith, decdithbase + iloop * loopoffset)

    #convert to degrees
    radith /= 3600.0
    decdith /= 3600.0

    #convert dec to radians
    decrad = dec * np.pi / 180.
    
    #apply dither
    radithdeg = ra + radith /  np.cos(decrad)   #account for cos(dec) term in ra offset
    decdithdeg = dec + decdith

    #plot dithers
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)

    ax.plot(radithdeg,decdithdeg,'ro')
    for i in range(len(radithdeg)):
        ax.text(radithdeg[i],decdithdeg[i],i)
    plt.savefig('test.png')
    
    #convert to sexigesimal
    RAdithhms = Angle(radithdeg,unit='deg').to_string(unit=u.hour,sep='')
    DECdithhms = Angle(decdithdeg,unit='deg').to_string(unit=u.degree,sep='')

    #exposure time for each
    exptime = {
        "r" : 180.0,
        "Ha+4nm" : 300.0
        }
    filts = ["r", "Ha+4nm"]
    #output files
    for filt in filts:
        #output prefix
        outfile_pre = "dither_" + objroot + "_t" + str(exptime[filt]) + "_nexp" + str(nexp) + "_dsize" + str(dithsize)

        outfile = outfile_pre + "_filt_" + filt + ".txt"
        fo = open(outfile, "w")

        #compute dithers, assumed to be the same for each filter
        for idith,val in enumerate(radith):
            filtname = filt.replace('+4nm','4')
            dithname = objroot + '_' + filtname

            #make strings for positions
            rastr = "{:0>9.2f}".format(float(RAdithhms[idith]))
            #DECdithhms[idith]
            if float(DECdithhms[idith]) < 0:
                decstr = "-{:0>8.1f}".format(abs(float(DECdithhms[idith])))
            else:
                decstr = "+{:0>8.1f}".format(float(DECdithhms[idith]))
            
            print(rastr,decstr)
            #print(idith,RAdithhms[idith],DECdithhms[idith])
            #write to each filter
            fo.write('obs {} object {} {} {} {} {} 2000.0\n'.format(exptime[filt],dithname,nexp,filt,rastr,decstr))
       
        fo.close()


    
