#!/usr/bin/env python

'''
 written by Rose A. Finn on 2/3/2013
 updated on 2/22/2016

 GOAL:
    - this program is only used to calculate local density and write out local density files

 METHOD:
    - read in NSA catalog
    - calculate various local density estimates
    - write output

 USAGE:

   from w/in ipython 

   %run ~/Dropbox/pythonCode/LCSCalcLocalDensity.py

   from command line

   python ~/Dropbox/pythonCode/LCSCalcLocalDensity.py
 

 UPDATES:

   
'''


from pylab import *
import os
from astropy.io import fits
from LCScommon import *

#import ReadAGCsav
from pyraf import iraf
import mystuff as my
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-v',"--velocitycut",help='velocity cut to use when calculating local density.  Galaxies within +/- vcut will be used.  Default value is 1500 km/s.',default=1500.,type=float)
parser.add_argument('-m',"--magnitudecut",help='r-band absolute magnitude cut.  Galaxies with mag less than Mr will be used for calculating density.',default=-18.3,type=float)

args = parser.parse_args()


class cluster:
    def __init__(self,clustername):
        self.prefix=clustername

        # read NSA table for each cluster
        infile=homedir+'research/NSA/'+self.prefix+'_NSA.Fits'
        self.ndat=fits.getdata(infile)
        self.nsadir=homedir+'research/NSA/'

        self.cra=clusterRA[self.prefix]
        self.cdec=clusterDec[self.prefix]
        self.cz=clusterz[self.prefix] 
        self.biweightvel=clusterbiweightcenter[self.prefix]
        self.biweightscale=clusterbiweightscale[self.prefix]
        self.r200=2.02*(self.biweightscale)/1000./sqrt(OmegaL+OmegaM*(1.+self.cz)**3)*H0/70. # in Mpc
        self.r200deg=self.r200*1000./my.DA(self.cz,h)/3600.

        self.cdMpc=self.biweightvel/H0
        self.cdcm=self.cdMpc*3.e24
        self.csigma=self.biweightscale
        self.mcl=my.clusterMass(self.csigma,self.cz,h)
        self.AngDistance=my.DA(self.cz,h)
        self.sdssflag=(self.ndat.ISDSS > -1) & (self.ndat.ABSMAG[:,4] < args.magnitudecut)


    def localdensity5(self):
        # get local density by summing mass w/in 300 kpc
        # NOTE:
        #  - this measure will be unreliable for galaxies near the edge of the 3 deg area
        #  - this will be fine for galaxies on 24um image

        DA=self.AngDistance
        x=self.ndat.RA
        y=self.ndat.DEC
        xref=self.ndat.RA[self.sdssflag]
        yref=self.ndat.DEC[self.sdssflag]
        n1=6
        sigma5=zeros(len(x),'d')
        sigma10=zeros(len(x),'d')
        for i in range(len(x)):
            deltav=abs(self.ndat.ZDIST[self.sdssflag]-self.ndat.ZDIST[i])*3.e5
            d=sqrt((x[i]-xref)**2+(y[i]-yref)**2)*3600./1000.*DA#d in Mpc
            flag = (deltav < args.velocitycut) 

            #print i,sum(vflag)
            d=d[flag]#apply a velocity cut of +/- 1500 km/s
            d.sort()#sort in ascending order, zeroth element is distance from galaxy to itself
            sig=0
            if len(d) < n1:
                print 'not enough points to calculate local density'
                print 'only have ',len(d),' galaxies w/in 1500 km/s'
                print 'skipping to next entry'
                continue
            else:
                sigma5[i]=1./(4.*pi)*(1.*5)/(d[5])**2
            try:
                sigma10[i]=1./(4.*pi)*(1.*10)/(d[10])**2
            except IndexError:
                continue
        self.sigma_5=sigma5
        self.sigma_10=sigma10
        

    def localdensity(self):
        # get local density by summing mass w/in 300 kpc
        # NOTE:
        #  - this measure will be unreliable for galaxies near the edge of the 3 deg area
        #  - this will be fine for galaxies on 24um image

        DA=self.AngDistance
        x=self.ndat.RA
        y=self.ndat.DEC
        xref=self.ndat.RA[self.sdssflag]
        yref=self.ndat.DEC[self.sdssflag]

        n1=3
        n2=6
        sigma=zeros(len(x),'d')
        for i in range(len(x)):
            deltav=abs(self.ndat.ZDIST[self.sdssflag]-self.ndat.ZDIST[i])*3.e5
            d=sqrt((x[i]-xref)**2+(y[i]-yref)**2)*3600./1000.*DA#d in Mpc
            flag = (deltav < args.velocitycut) 

            d=d[flag]#apply a velocity cut of +/- 1500 km/s
            d.sort()#sort in ascending order, zeroth element is distance from galaxy to itself
            sig=0
            if len(d) < n1:
                print 'not enough points to calculate local density'
                print 'only have ',len(d),' galaxies w/in 1500 km/s'
                print 'skipping to next entry'
                continue
            else:
                sigma[i]=1./(4.*pi)*d[1]
                #sigma5[i]=1./(4.*pi)*(1.*5)/(d[5])**2

            

        self.sigma_nn=sigma
        
        return
    def localdensitybymass(self):#find local density, using nearest neighbors n1 through n2
        DA=self.AngDistance #kpc/arcsec
        angDist_300kpc=300./DA/3600. # angular dist corresponding to 300 kpc, in degrees
        x=self.ndat.RA
        y=self.ndat.DEC
        # should probably use SDSS as the reference sample, but for now, using NSA

        xref=self.ndat.RA[self.sdssflag]
        yref=self.ndat.DEC[self.sdssflag]
        self.rhomass=zeros(len(x),'d')
        for i in range(len(x)):
            deltav=abs(self.ndat.ZDIST[self.sdssflag]-self.ndat.ZDIST[i])*3.e5
            d=sqrt((x[i]-xref)**2+(y[i]-yref)**2)*3600./1000.*DA#d in Mpc
            flag = (deltav < args.velocitycut) 
            d=d[flag]#apply a velocity cut of +/- 1500 km/s
            neighbor_flag=d < .3

            self.rhomass[i]=sum(self.ndat.MASS[neighbor_flag])
            #print i, self.rhomass[i],sum(neighbor_flag)
    def matchnsa2sdssphot(self):
        return
    def writeoutput(self):
        # append flag to NSA table 
        # create new table
        #ldat=atpy.Table()
        #ldat.add_column('On24ImageFlag',self.On24ImageFlag)

        # will add more columns here as time permits
        c1 = fits.Column(name='NSAID',array=self.ndat.NSAID,unit='',format='J')#,description='NSA ID')
        c2 = fits.Column(name='RA',array=self.ndat.RA,unit='deg',format='E')#,description='RA')
        c3 = fits.Column(name='DEC',array=self.ndat.DEC,unit='deg',format='E')#,description='DEC')
        c4 = fits.Column(name='Z',array=self.ndat.Z,unit='',format='E')#,description='redshift')
        c5 = fits.Column(name='ZDIST',array=self.ndat.ZDIST,unit='',format='E')#,description='NSA ZDIST')
        c6 = fits.Column(name='RHOMASS',array=self.rhomass,unit='Msun',format='D')#,description='Mass of galaxies w/in 300 kpc')
        c7 = fits.Column(name='SIGMA_NN',array=self.sigma_nn,unit='Ngal/Mpc^2',format='E')#, description='3rd to 6th nearest neighbor')
        c8 = fits.Column(name='SIGMA_5',array=self.sigma_5,unit='Ngal/Mpc^2',format='E')#, description='to 5th nearest neighbor')
        c9 = fits.Column(name='SIGMA_10',array=self.sigma_10,unit='Ngal/Mpc^2',format='E')#, description='to 10th nearest neighbor')
        coldefs=fits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9])
        ldat=fits.new_table(coldefs)
        # write out table
        outfile=homedir+'research/LocalClusters/NSAmastertables/LocalDensityTables/'+self.prefix+'_localdensity.fits'
        ldat.writeto(outfile,clobber=True)


if __name__ == '__main__':


    for cname in clusternames:
        cl=cluster(cname)
        print  '\n ',cl.prefix,'\n'
        cl.localdensity()
        cl.localdensity5()
        cl.localdensitybymass()
        cl.writeoutput()

