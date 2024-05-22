#!/usr/bin/env python

from astropy.io import fits
import glob


rfiles = glob.glob("*r-shifted.fits")
badfiles = open('redo_alignment_bok.txt','w')
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
        # add halpha image name to file of images to realign
        badfiles.write(f"{himage}\n")

badfiles.close()
print(f"number of problems = {nbad} out of {len(rfiles)}")
if nbad > 0:
    print(":(")
else:
    print(":)")

print("to fix the problems, type:\n\n parallel --eta python ~/github/HalphaImaging/python3/BOK_align_images_wrapper.py  :::: redo_alignment_bok.txt"
