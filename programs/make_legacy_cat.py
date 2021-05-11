#!/usr/bin/env python

'''
GOAL:
* merge north and south legacy catalogs into one

* north for dec >= 32.375
* south for dec < 32.375


'''

from astropy.table import Table
import numpy as np
import os
homedir = os.getenv("HOME")
tabledir = os.path.join(homedir,'research','Virgo','tables-moustakas',"")
vfdir = os.path.join(homedir,'research','Virgo','tables-north',"v1","")
northcat = Table.read(tabledir+'vf_north_v1_main_dr9north.fits')
southcat = Table.read(tabledir+'vf_north_v1_main_dr9south.fits')


legcat = Table(np.empty(len(northcat),dtype=northcat.dtype))

vmain = Table.read(vfdir+"vf_north_v1_main.fits")
nflag = vmain['DEC'] >= 32.375
legcat[nflag] = northcat[nflag]
legcat[~nflag] = southcat[~nflag]
legcat.write(vfdir+'vf_north_v1_legacy_dr9.fits',overwrite=True)
