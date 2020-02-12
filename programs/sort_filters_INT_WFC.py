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

ic = ImageFileCollection(os.getcwd(),keywords='*',glob_include='r*.fit*')

# get list of filter, imagetype, and object names

filters = ic.values('wffband')
object_names = ic.values('object')

# create a unique set of object-filter
filefilterlist = ["{}-{}".format(str(i).replace(" ","").replace("-","").replace("_",""),str(j).replace(" ","")) for i,j in zip(object_names,filters)]

t = ['target']*len(filefilterlist)   
targetfilterlist = ["{}-{}".format(str(i).replace(" ","").replace("-","").replace("_",""),str(j).replace(" ","")) for i,j in zip(t,filters)]

# find bias and change the name to BIAS (without filter)
# replace Skybackground with SKYFLAT
for i,d in enumerate(filefilterlist):
    if d.find('Bias') > -1:
        filefilterlist[i] = 'BIAS'
    if d.find('Skybackground') > -1:
        print('found it')
        filefilterlist[i] = filefilterlist[i].replace("Skybackground", "SKYFLAT")

# make directory for BIAS
# make directory for each flat/filter combination

for d in dirnames:
    if d.find('--') > -1:
        continue
    if os.path.exists(d):
        continue
    else:
        if d.find('BIAS') > -1:
            os.mkdir(d)
        elif d.find('FLAT') > -1:
            os.mkdir(d)

# move bias frames to subdirectory
# move flats there
# move files to appropriate subdirectory
for f,d in zip(ic.files,filefilterlist):
    if (filefilterlist.find('FLAT') > -1) | (filefilterlist.find('BIAS') > -1):
        try:
            os.rename(f,d+'/'+f)
        except:
            print("Error moving file {} to directory {}".format(f,d))
            

# make target directories
t = set(filter)
unique_filters = [f.replace(' ','') for f in t]
for f in unique_filters:
    fdir = 'target'+f
    if os.path.exists(fdir):
        continue
    else:
        os.mkdir(fdir)
        
# group remaining targets into directories according to filter
for f,d in zip(ic.files,targetfilterlist):
    try:
        os.rename(f,d+'/'+f)
    except:
        print("Error moving file {} to directory {}".format(f,d))






        
