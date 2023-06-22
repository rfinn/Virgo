#!/usr/bin/env python

"""
GOAL:
* gather up the content from the csv files that are created when making webpage for coadds
* csv files contain list of galaxies in FOV
* this will allow us to get estimate of halpha sample and cross-match with main table


USEAGE:
* run on draco from /data-pool/html_dev/coadds directory

TESTING


"""

import os
from astropy.table import Table, Column
import numpy as np

os.chdir('/home/rfinn/research/Virgo-dev/html-dev/coadds/')
t = os.listdir()
csvfiles = []
for d in t:
    csvfiles.append(os.path.join(d,d+'-galsFOV.csv'))

# merge all csv files
vfid = []
coadd = []

for f in csvfiles:
    if f.startswith('VF-'):
        input = open(f,'r')
        for line in input:
            s = line.rstrip().split(',')
            vfid.append(s[0])
            coadd.append(s[1])        

unique_ids = set(vfid)

# read in vf main

vmain = Table.read('/home/rfinn/research/Virgo/tables-north/v2/vf_v2_main.fits')

# create a new table with full list of vfids

newtab = Table([vmain['VFID']],names=['VFID'])

newcol = Column(np.zeros(len(vmain),'bool'),name='HaObsFlag')

newtab.add_column(newcol)
# for creaate
coadd1 = ["" for x in range(len(vmain))]
coadd2 = ["" for x in range(len(vmain))]
coadd3 = ["" for x in range(len(vmain))]

# make a dictionary with vfid

vfid_dict = dict((a,b) for a,b in zip(vmain['VFID'],np.arange(len(vmain))))
ncoadds = np.zeros(len(vmain),'i')
for v,c in zip(vfid,coadd):
    i = vfid_dict[v]
    
    newtab['HaObsFlag'][i] = True
    
    if coadd1[i] == "":
        coadd1[i] = c
        ncoadds[i] = 1
    elif coadd2[i] == "":
        coadd2[i] = c
        ncoadds[i] = 2
    elif coadd3[i] == "":
        coadd3[i] = c
        ncoadds[i] = 3        
col1 = Column(coadd1,name='Coadd1')
col2 = Column(coadd2,name='Coadd2')
col3 = Column(coadd3,name='Coadd3')
col4 = Column(ncoadds,name='ncoadd')
newtab.add_columns([col1,col2,col3,col4])


outfile = '/home/rfinn/research/Virgo/tables-north/v2/vf_v2_halpha_obs.fits'
newtab.write(outfile,format='fits',overwrite=True)
