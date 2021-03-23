#!/usr/bin/env python

'''
GOAL:

gather zp plots from individual pointing directories in 
/home/rfinn/data/reduced/scratch-int-feb2019/attempt2

cp subdir/plots/*.png to 
/home/rfinn/data/reduced/virgo-coadds-feb2019-int

making it general b/c I might need to do the same for the June INT data

'''

import os
import argparse


parser = argparse.ArgumentParser(description ='gather plots from ZP calibration of INT images')
parser.add_argument('--path2subdirs', dest = 'path2subdirs', default = '/home/rfinn/data/reduced/scratch-int-feb2019/attempt2/', help = 'path to folder that has subdirectory for each pointing.  default is /home/rfinn/data/reduced/scratch-int-feb2019/attempt2.')
parser.add_argument('--outdir', dest = 'outdir', default = '/home/rfinn/data/reduced/virgo-coadds-feb2019-int/plots/', help = 'output directory for ZP plots.  the default is /home/rfinn/data/reduced/virgo-coadds-feb2019-int')

args = parser.parse_args()


# create outdir if it doesn't already exist
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

dirlist = os.listdir(args.path2subdirs)
print(dirlist)
dirlist.sort()

for d in dirlist:
    if os.path.isdir(os.path.join(args.path2subdirs,d)):
        print('######################################')
        print('###### WORKING ON ',d)        
        print('######################################')
        source_dir = os.path.join(args.path2subdirs,d,'plots')
        command_string = "cp {}/*.png {}/.".format(source_dir,args.outdir)
        print(command_string)
        os.system(command_string)
