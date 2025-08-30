#!/usr/bin/env python
"""
run from directory that contains the coadds

this will create the list of coadded images that I can then use 
as input for other files


"""
import glob
import argparse

# get list of r-band coadded images
parser = argparse.ArgumentParser(description ='Get the list of coadded images.')
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
outfile4 = open('virgo-coadds-test.csv','w')
outfile2 = open('virgo-coadds-fullpath.txt','w')
outfile3 = open('virgo-coadds-fullpath-test.txt','w')
if args.bokonly:
    outfile = open('virgo-bok-coadds.csv','w')
    outfile2 = open('virgo-bok-coadds-fullpath.txt','w')
    outfile3 = open('virgo-bok-coadds-fullpath-test.txt','w')
    outfile4 = open('virgo-bok-coadds-test.csv','w')    
coadd_dir = '/data-pool/Halpha/coadds/all-virgo-coadds'
for i in range(len(rfiles)):
    #basname = rfiles[i].replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    outfile.write(f"{rfiles[i]}\n")
    outfile2.write(f"{coadd_dir}/{rfiles[i]}\n")
    if i < 2:
        outfile3.write(f"{coadd_dir}/{rfiles[i]}\n")
        outfile4.write(f"{rfiles[i]}\n")
outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()
