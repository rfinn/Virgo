#!/usr/bin/env python

'''
Run this from a directory from a given night

it will remove original BIAS and SKYFLATs

then remove ORIGINAL in target-r, target-Halpha directories

you can do this after the data is calibrated

goal is to save some disk space by removing copies of intermediate stage images

'''

import os
import glob
directories = ['BIAS','SKYFLAT-r','SKYFLAT-Halpha','SKYFLAT-Ha6657','target-r','target-Halpha','target-Ha6657']
subdir = ['ORIGINALS','SPLIT_IMAGES','MASK_IMAGES']
for d in directories:
    for s in subdir:
        path = '/'.join([d,s,'r*.fit*'])
        print(path)
        rmfiles = glob.glob(path)
        for f in rmfiles:
            os.remove(f)

# remove input files for coadd
files = os.listdir()
for f in files:
    if os.path.isdir(f):
        
        coadd_dirs = ['coadd_r','coadd_Halpha']
        for c in coadd_dirs:
            subdir = '/'.join([f,c])
            if os.path.exists(subdir):
                rmfiles = glob.glob('/'.join([subdir,'r*.fits']))
                for f in rmfiles:
                    os.remove(f)
        

