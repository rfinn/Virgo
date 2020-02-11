#!/usr/bin/env python

'''
GOAL:
- sort INT WFC data according to target and filter
-

RELEVANT HEADER KEYWORDS:

    OBSTYPE: TARGET, SKY (skyflat), BIAS, FOCUS

    IMAGETYP: object, sky, focus, zero

    OBJECT: name that we assigned, i.e. pointing-021
    - need to strip spaces

    CATNAME:  similar to OBJECT, but capitalized? not all are the same.  use OBJECT

    WFFBAND - filter; need to strip spaces

    WFFPOS - filter position (can use this to double check)

PROCEDURE:

get list of filter, imagetype, and object names

for flats, create a directory for each filter called FLAT-HALPHA.  Move flat to appropriate directory

for bias, create BIAS, and move files there

for focus, create FOCUS directory, move files there

for science frames, sort them by filter into, e.g. SCIENCE-r or SCIENCE-Halpha
- these will get sorted by object once bias subtraction and flatfielding is done.

After basic calibration, we can sort science objects, make directory for each unique set of OBJECT-WFFBAND.  move files to appropriate directory.

'''
#import ccdproc
from ccdproc import ImageFileCollection
import os
import sys

##### SET PATHS  ####
homedir = os.getenv("HOME")
working_dir = os.getcwd()
# path to theli scripts
theli_path = homedir+'/theli/scripts/Linux_64/'


# get list of directory names
t=[f.path for f in os.scandir(os.getcwd()) if f.is_dir() ] 
dirnames = []
for d in t:
        if d.find('junk') > -1:
                continue
        dirnames.append(os.path.basename(d))

# start reduction script

# split files
for d in dirnames:
    if os.path.exists(d):
        command_string = theli_path+'process_split_WFC@INT.sh '+working_dir+' '+d
        os.system(command_string)


# process BIAS frames
for d in dirnames:
    if os.path.exists(d):
        if d == 'BIAS':
            # process bias frames
            command_string = theli_path+'parallel_manager.sh process_bias_para.sh '+working_dir+' BIAS'
            os.system(command_string)
            
# process Flat frames
for d in dirnames:
    if os.path.exists(d):
        if d.find('FLAT'):
            # process flats
            command_string = theli_path+'parallel_manager.sh process_flat_para.sh '+working_dir+' BIAS '+d
            command_string = theli_path+'create_flat_ratio.sh '+working_dir+' '+d
            command_string = theli_path+'parallel_manager.sh create_norm_para.sh '+working_dir+' '+d

# calibrate images
# need to apply the correct flat to each image

flat_images = []

for d in dirnames:
    if d.find('FLAT'):
        flat_images.append(d)

for f in flat_images:
    ffilter = f.split('-')[1]
    for d in dirnames:
        #skip bias and flat images
        if d.find('BIAS') > -1:
            continue
        if d.find('FLAT') > -1:
            continue
        if os.path.exists(d):
            ofilter = d.split('-')[1]
            if ofilter == ffilter:
                command_string =  theli_path+'/parallel_manager.sh process_science_para.sh '+working_dir+' BIAS '+f+' '+d
                os.system(command_string)

# create a preview image
# calibrate weights
for d in dirnames:
        #skip bias and flat images
        if d.find('BIAS') > -1:
            continue
        if d.find('FLAT') > -1:
            continue
        if os.path.exists(d):
            ofilter = d.split('-')[1]
            command_string = theli_path+'make_album_WFC@INT.sh '+working_dir+' '+d+' OFC'
            os.system(command_string)
            command_string = theli_path+'create_tiff.sh '+working_dir+' '+d+' OFC'
            os.system(command_string)
            command_string = theli_path+'parallel_manager.sh create_global_weights_para.sh '+working_dir+' SKYFLAT-'+str(ofilter)+'_norm '+d 
            os.system(command_string)

# create weights
for d in dirnames:
        #skip bias and flat images
        if d.find('BIAS') > -1:
            continue
        if d.find('FLAT') > -1:
            continue
        if os.path.exists(d):
            command_string = theli_path+'transform_ds9_reg.sh '+working_dir+' '+d
            os.system(command_string)
            command_string = theli_path+'parallel_manager.sh create_weights_para.sh '+working_dir+' '+d+' OFC'
            os.system(command_string)


