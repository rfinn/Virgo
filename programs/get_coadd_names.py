#!/usr/bin/env python
"""
run from directory that contains the coadds
"""
import glob
import argparse

# get list of r-band coadded images
parser = argparse.ArgumentParser(description ='run the halpha gui in automatic mode.  Run this from the directory where you want the output data stored.  For example: /home/rfinn/research/Virgo/gui-output-NGC5846/')
parser.add_argument('--bokonly',dest = 'bokonly', default = False, action='store_true', help = 'set this to create file with bok filenames only.  needed this for rerunning gui after alignment.')
args = parser.parse_args()



a = glob.glob('VF*INT*-r-shifted.fits')
b = glob.glob('VF*HDI*-r.fits')
c = glob.glob('VF*HDI*-R.fits')

#d = glob.glob('VF*BOK*-r.fits')
#######################################################
# updating to use the shifted r-band images
#######################################################    
d = glob.glob('VF*BOK*-r-shifted.fits')

e = glob.glob('VF*MOS*-R.fits')         
rfiles = a + b + c + d + e

if args.bokonly:
    rfiles = d
rfiles.sort()
print(f"number of targets = {len(rfiles)}")

# write out as a csv file
outfile = open('virgo-coadds.csv','w')
outfile2 = open('virgo-coadds-fullpath.txt','w')
outfile3 = open('virgo-coadds-fullpath-test.txt','w')
if args.bokonly:
    outfile = open('virgo-bok-coadds.csv','w')
    outfile2 = open('virgo-bok-coadds-fullpath.txt','w')
    outfile3 = open('virgo-bok-coadds-fullpath-test.txt','w')
    
coadd_dir = '/data-pool/Halpha/coadds/all-virgo-coadds'
for i in range(len(rfiles)):
    #basname = rfiles[i].replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    outfile.write(f"{rfiles[i]}\n")
    outfile2.write(f"{coadd_dir}/{rfiles[i]}\n")
    if i < 2:
        outfile3.write(f"{coadd_dir}/{rfiles[i]}\n")
outfile.close()
outfile2.close()
outfile3.close()
