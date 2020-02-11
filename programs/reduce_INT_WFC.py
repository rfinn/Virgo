#!/usr/bin/env python

'''
GOAL:
- Call theli scripts from the command line

USAGE:
- run from top directory for a particular night, after you have run sort_INT_WFC.py
-


'''
#import ccdproc
from ccdproc import ImageFileCollection
import os
import stat
import sys
import subprocess


###  SET ENVIRONMENT VARIABLES


os.environ["INSTRUMENT"] = "WFC@INT"


##### SET PATHS  ####
homedir = os.getenv("HOME")
data_dir = os.getcwd()
# path to theli scripts
if data_dir.find('Users/rfinn') > -1: # running on macbook
    theli_path = homedir+'/software/theli/gui-2.10.5/scripts/'
    # make sure python 2 is called instead of python 3
    p = os.getenv("PATH")
    os.environ["PATH"] = '/usr/bin/:'+p

else: # assume this is Virgo machine
    theli_path = homedir+'/gui-2.10.5/scripts/'
    # make sure python 2 is called instead of python 3
    #p = os.getenv("PATH")
    #os.environ["PATH"] = '/usr/local/anaconda2/bin/:'+p


# append path to scripts so it knows where they are
p = os.getenv("PATH")
os.environ["PATH"] = p+':'+theli_path

os.environ["INSTRUMENT"] = "WFC@INT"


# get list of directory names
t=[f.path for f in os.scandir(os.getcwd()) if f.is_dir() ] 
dirnames = []
for d in t:
        if d.find('junk') > -1:
                continue
        dirnames.append(os.path.basename(d))

# start reduction script

def split_files():
    # change to theli dir because all scripts call other programs in a relative way

    # split files
    i=0
    for d in dirnames:
        if os.path.exists(d):
            print('splitting files in ',d)
            os.chdir(theli_path)
            #os.chmod('process_split_WFC@INT.sh', stat.S_IEXEC)
            command_string = './process_split_WFC@INT.sh '+data_dir+' '+d
            rc = subprocess.call(command_string, shell=True)
            #os.system(command_string)
            os.chdir(data_dir)
            i += 1
            if i > 0:
                return
def process_bias():
    # process BIAS frames
    for d in dirnames:
        if os.path.exists(d):
            if d == 'BIAS':
                # process bias frames
                print('processing bias frames in', d)
                command_string = theli_path+'parallel_manager.sh process_bias_para.sh '+data_dir+' BIAS'
                os.system(command_string)
def process_flats():
    # process Flat frames
    for d in dirnames:
        if os.path.exists(d):
            if d.find('FLAT'):
                # process flats
                print('processing flatfield images in ',d)
                command_string = theli_path+'parallel_manager.sh process_flat_para.sh '+data_dir+' BIAS '+d
                command_string = theli_path+'create_flat_ratio.sh '+data_dir+' '+d
                command_string = theli_path+'parallel_manager.sh create_norm_para.sh '+data_dir+' '+d

def calibrate_images():
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
                    print('processing images in ',d)
                    command_string =  theli_path+'/parallel_manager.sh process_science_para.sh '+data_dir+' BIAS '+f+' '+d
                    os.system(command_string)

def preview_weights():
    # create a preview image
    # calibrate weights
    for d in dirnames:#skip bias and flat images
        if d.find('BIAS') > -1:
            continue
        if d.find('FLAT') > -1:
            continue
        if os.path.exists(d):
            print('creating preview and calibrating weights for ',d)
            ofilter = d.split('-')[1]
            command_string = theli_path+'make_album_WFC@INT.sh '+data_dir+' '+d+' OFC'
            os.system(command_string)
            command_string = theli_path+'create_tiff.sh '+data_dir+' '+d+' OFC'
            os.system(command_string)
            command_string = theli_path+'parallel_manager.sh create_global_weights_para.sh '+data_dir+' SKYFLAT-'+str(ofilter)+'_norm '+d 
            os.system(command_string)

def create_weights():
    
    # create weights
    for d in dirnames:
        #skip bias and flat images
        if d.find('BIAS') > -1:
            continue
        if d.find('FLAT') > -1:
            continue
        if os.path.exists(d):
            print('creating weights for ',d)
            command_string = theli_path+'transform_ds9_reg.sh '+data_dir+' '+d
            os.system(command_string)
            command_string = theli_path+'parallel_manager.sh create_weights_para.sh '+data_dir+' '+d+' OFC'
            os.system(command_string)

# create astrometric reference catalog
'''
This calls scripts/create_astrorefcat_fromWEB.sh

which then calls sqlcl.py - so need to make sure that python is python2



'''
def create_astroref_cat():
    pass


if __name__ == '__main__':
    print('good luck!')
    
    split_files()
    #process_bias()
    #process_flats()
    #calibrate_images()
    #preview_weights()
    #create_weights()
