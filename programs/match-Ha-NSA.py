#!/usr/bin/env python

'''
read in Halpha file
- downloaded Halpha file from https://docs.google.com/spreadsheets/d/1vmY5RrzM_LU2rPkebBq3chwj85vPlUcTKR1US2OvchM/edit?usp=sharing
- removed the last line
- will read it in and write it out as fits file

- copy downloaded file to "latest"
cp Observing-Summary-Halpha-good-2019May28.csv Observing-Summary-Halpha-good-latest.csv

%run programs/match-Ha-NSA.py --write-fits



'''

from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import argparse
import os

parser = argparse.ArgumentParser(description ='Match the Halpha observations with Virgo NSA catalog')

parser.add_argument('--table-path', dest = 'tablepath', default = 'github/Virgo/tables/', help = 'path to github/Virgo/tables, relative to home directory')
parser.add_argument('--write-fits',dest = 'writefits', action='store_true',help='write out fits version of observing summary file?')
parser.add_argument('--input',dest = 'input', default='Observing-Summary-Halpha-good-latest.csv',help='write out fits version of observing summary file?')
        
args = parser.parse_args()

homedir = os.getenv("HOME")

if args.writefits:
    ### Read in csv and write out fits
    infile = homedir+'/'+args.tablepath+args.input
    outfile = infile.replace('csv','fits')
    hadat = np.recfromcsv(infile)
    fits.writeto(outfile,hadat,overwrite=True)


vdat = fits.getdata(homedir+'/'+args.tablepath +'nsa.virgo.fits')
hdat = fits.getdata(homedir+'/'+args.tablepath + 'Observing-Summary-Halpha-good-latest.fits')

nsadict=dict((a,b) for a,b in zip(vdat.NSAID,np.arange(len(vdat.NSAID))))
# match by NSAID
index = np.zeros(len(hdat.nsa_id),'i')
matchflag = np.zeros(len(hdat.nsa_id),'bool')
for i in range(len(hdat.nsa_id)):
    try:
        index[i] = nsadict[int(hdat.nsa_id[i])]
        matchflag[i] = 1
    except:
        print('could not match NSAID = ',hdat.nsa_id[i])
        #print 'recession velocity = ',vdat.Z[nsadict[hdat.nsa_id[i]]]*3.e5
    





# write out line-matched catalog
outfile= args.tablepath + 'nsa_Halpha.virgo.fits'
matchedarray1=np.zeros(len(vdat),dtype=hdat.dtype)
matchedarray1[index[matchflag]] = hdat[matchflag]

withcoords  = np.lib.recfunctions.append_fields(matchedarray1,['RA','DEC','NSAID'],[vdat.RA,vdat.DEC,vdat.NSAID],dtypes=[vdat.RA.dtype,vdat.DEC.dtype,vdat.NSAID.dtype],asrecarray=True)

fits.writeto(outfile,withcoords,overwrite=True)
