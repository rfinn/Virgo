#!/usr/bin/env python

'''
GOAL:
- read in csv version of agc
- write out fits

'''

from astropy.io import ascii, fits
import os
from astropy.table import Table

homedir = os.getenv('HOME')
agc_csv = homedir+'/research/AGC/agcnorthminus1.csv'
agc_fits = homedir+'/research/AGC/agcnorthminus1.fits'

ac = ascii.read(agc_csv, delimiter=',')
ac = Table(ac)
ac.write(agc_fits,format='fits')

