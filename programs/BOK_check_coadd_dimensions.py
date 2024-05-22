#!/usr/bin/env python

from astropy.io import fits
import glob


rfiles = glob.glob("*r-shifted.fits")

nbad=0
for r in rfiles:
    rheader = fits.getheader(r)
    rdata = fits.getdata(r)
    himage = rheader['HAIMAGE']
    rnaxis1,rnaxis2 = rdata.shape


    # get dimensions of halpha image
    hdata = fits.getdata(himage)

    hnaxis1,hnaxis2 = hdata.shape    

    if (rnaxis1 == hnaxis1) & (rnaxis2 == hnaxis2):
        continue
    else:
        print("HOLD UP: ",r)
        print(f"\tr-band dimensions = {rnaxis1},{rnaxis2}")
        print(f"\thalpha dimensions = {hnaxis1},{hnaxis2}")             
        nbad += 1

        # realign images
        im1 = himage
        im2 = r.replace("-shifted.fits",".fits")
        weight2 = im2.replace(".fits",".weight.fits")
        os.system(r"python ~/github/HalphaImaging/python3/INT_align_images.py --image1 {im1} --image2 {im2} --weight2 {weight2}"

print
print(f"number of problems = {nbad} out of {len(rfiles)}")
if nbad > 0:
    print(":(")
else:
    print(":)")
