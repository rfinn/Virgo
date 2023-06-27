#!/usr/bin/env python
"""
run from directory that contains the coadds
"""
import glob
# get list of r-band coadded images

a = glob.glob('VF*INT*-r-shifted.fits')
b = glob.glob('VF*HDI*-r.fits')
c = glob.glob('VF*HDI*-R.fits')
d = glob.glob('VF*BOK*-r.fits')         
rfiles = a + b + c + d

rfiles.sort()
print(f"number of targets = {len(rfiles)}")

# write out as a csv file
outfile = open('virgo-coadds.csv','w')
for i in range(len(rfiles)):
    basname = rfiles[i].replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    outfile.write(f"{basname} \n")
outfile.close()
