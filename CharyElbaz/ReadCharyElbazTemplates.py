#!/usr/bin/env python
'''
    Written by Rose A. Finn, November 20, 2012

    PURPOSE
      program reads in the chary_elbaz.sav file using idlsave

 
    CALLING SEQUENCE
      import ReadCharyElbazTemplates
      ce=ReadCharyElbazTemplates.charyelbaz(savefile)

      where savefile is the location of chary_elbaz.save
      e.g.
      savefile='/home/rfinn/research/SEDs/templates/CharyElbaz01/chary_elbaz_codes/chary_elbaz.save'

      to plot SEDs
      plot(ce.lamb,ce.nulnuinlsun[:,1])

      or to plot all SEDs
      ce.plotseds()
    INPUT PARAMETERS
      location of chary_elbaz.save file

    OUTPUT PARAMETERS
      class charyelbaz with the following associated vectors:
        self.lamb
        self.nulnu_iras25
        self.nulnu_iras100
        self.nulnu_iras12
        self.nulnu_iras60
        self.nulnuinlsun
        self.lir_sanders
        self.lir
        self.nulnu_lw3
        self.nulnu_lw2

    REQUIRED PYTHON MODULES
      idlsave
      pylab
      os
      scipy
NOTES:
2018-06-11 - updated to use scipy.io.readsav instead of idlsave

'''

#import idlsave
from scipy.io import readsav
from pylab import *
import os

class charyelbaz:
    def __init__(self,infile):
        mypath=os.getcwd()
        if mypath.find('/Users/rfinn') > -1:
            print "Running on Rose's mac pro or laptop"
            homedir='/Users/rfinn/'
        elif mypath.find('Users/kellywhalen') > -1:
            print "Running on Kelly's Laptop"
            homedir='/Users/kellywhalen/Github/Virgo/'

        cefile=readsav(infile)
        self.nulnu_iras25=cefile['nulnu_iras25']
        self.nulnu_iras100=cefile['nulnu_iras100']
        self.nulnu_iras12=cefile['nulnu_iras12']
        self.nulnu_iras60=cefile['nulnu_iras60']
        self.nulnuinlsun=cefile['nulnuinlsun']
        self.lir_sanders=cefile['lir_sanders']
        self.lir=cefile['lir']
        self.nulnu_lw3=cefile['nulnu_lw3']
        self.nulnu_lw2=cefile['nulnu_lw2']
        self.lamb=cefile['lambda']
        #
        # convert all to double-precision arrays
        #
        self.lamb=array(self.lamb,'d')
        self.nulnu_iras25=array(self.nulnu_iras25,'d')
        self.nulnu_iras100=array(self.nulnu_iras100,'d')
        self.nulnu_iras12=array(self.nulnu_iras12,'d')
        self.nulnu_iras60=array(self.nulnu_iras60,'d')
        self.nulnuinlsun=array(self.nulnuinlsun,'d')
        self.lir_sanders=array(self.lir_sanders,'d')
        self.lir=array(self.lir,'d')
        self.nulnu_lw3=array(self.nulnu_lw3,'d')
        self.nulnu_lw2=array(self.nulnu_lw3,'d')

    def plotseds(self):
        figure(figsize=(10,10))
        for i in range(len(self.nulnuinlsun[1])):
            plot(self.lamb,self.nulnuinlsun[:,i])
        ax=gca()
        ax.set_yscale('log')
        ax.set_xscale('log')
        title('Chary & Elbaz SEDs')
        xlabel(r'$\mathrm{Wavelength \ (\mu m)}$',fontsize=16)
        ylabel(r'$\mathrm{\nu L_\nu/L_\odot}$',fontsize=16)
        
        
                
