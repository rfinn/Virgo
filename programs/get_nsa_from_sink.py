#!/usr/bin/env python

from astropy.io import fits
from astropy.table import Table

infile1 = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink.fits'
infile2 = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink_v2.fits'
outfile1 = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink_NSAonly.fits'
outfile2 = '/home/rfinn/research/Virgo/supersample/smart_kitchen_sink_v2_NSAonly.fits'
# read in kitchen_sink.fits
infiles = [infile1, infile2]
outfiles = [outfile1, outfile2]

for i in range(len(infiles)):
    
    cat = fits.getdata(infiles[i])
    
    # keep rows 9302 onward (NSA only)
    #cat = cat[9302:]

    # keep NSA only
    HLflag = cat['objname'] != ""
    AGCflag = cat['AGCnr'] > 0
    NSAflag = cat['NSAID'] > 0
    flag = ~HLflag & ~AGCflag & NSAflag

    fits.writeto(outfiles[i],cat[flag],overwrite=True)

# get list of NSA galaxies in v2 but not in v1
# match on NSAID


# need to inspect these images

