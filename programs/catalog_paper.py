#! /usr/bin/env python

'''
GOAL: Plots for catalog paper
* 3D plots
* density vs distance
  * for each filament
  * calculate density within 0.5, 1 and 1.5 Mpc


'''
import os
import numpy as np
from matplotlib import pyplot as plt
from readtables import vtables


######################################################
##### SET PATHS
######################################################

homedir = os.getenv('HOME')
tablepath = homedir+'/research/APPSS/tables/'
tablepath = homedir+'/research/Virgo/tables-north/v1/'
spinepath = '/home/rfinn/research/Virgo/tables-north/spines/'

######################################################
##### FUNCTIONS
######################################################
def plotspines(ymin,ymax):
    sfiles = glob.glob(spinepath+"*.fits")
    for f in sfiles:
        spine  = Table.read(f)
        flag = (spine['SGY'] > ymin) & (spine['SGY'] < ymax)
        fname = os.path.basename(f).replace('filament_spine_','').replace('.fits','').replace('_Filament','')
        plt.plot(spine['SGX'][flag],spine['SGZ'][flag],lw=3,zorder=15,label=fname)



######################################################
##### CLASSES
######################################################

class vfcatalogs(vtables):
    ''' inherits vtables from read_tables, so will read in catalogs when initiated   '''
    def get_filaments(self):
        pass

class filament():
    ''' class for each filament.  pass in filament name and the vf env table   '''
    def __init__(self,name=None,envcat=None):
        # set name
        # read in spine
        pass

    def calc_density(self):
        ''' calculate dens within 0.5,1,1.5 Mpc '''
        pass

    def get_spine(self):
        ''' get spine for this filament '''

        # check if spine file is available

        # then open, read spine

        pass
    
        
