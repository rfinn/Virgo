#!/usr/bin/env python
# get list of r-band coadded images
a = glob.glob(coadd_dir+'VF*INT*-r-shifted.fits')
b = glob.glob(coadd_dir+'VF*HDI*-r.fits')
c = glob.glob(coadd_dir+'VF*HDI*-R.fits')
d = glob.glob(coadd_dir+'VF*BOK*-r.fits')         
rfiles = a + b + c + d

rfiles.sort()
print(f"number of targets = {len(rfiles)}")

# write out as a csv file
outfile = open('virgo-coadds.csv','w')
for i in range(len(rfiles)):
    outfile.write(f"{rfiles[i]} \n")
outfile.close()
