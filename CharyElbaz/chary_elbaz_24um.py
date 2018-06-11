#!/usr/bin/env python

'''
    Written by Rose Finn, November 27, 2012

    PURPOSE: 
      This module is a python translations of the chary & elbaz 2001 code for
      converting the observed 24um flux to the total IR luminosity (and SFR)
      
    CALLING SEQUENCE

      This program assumes that github/Virgo is off the home directory.
      
      ceLir,ceSFR=chary.chary_elbaz_24um(redshift[flag],self.mipsflux[flag])

    INPUT PARAMETERS
      redshift - either a number or an array
      mipsflux - observed mips flux in micro-Jy
      
    OUTPUT PARAMETERS
      Lir - total IR luminosity (8-1000 um) in solar luminosities
      SFRir - SFR corresponding to IR luminosity

    EXAMPLES
      import chary_elbaz_24um as chary
      ceLir,ceSFR=chary.chary_elbaz_24um(redshift,mipsflux)

      or in ipython, in the directory where you unpacked the tar file:
        > ipython -pylab
        > %run chary_elbaz_24um.py
        > redshift=arange(.1,.6,.1)
        > mipsflux = 5000.*ones(len(redshift))
        > Lir, SFR = chary_elbaz_24um(redshift,mipsflux)
        > figure()
        > plot(redshift, Lir,'ro')
        > figure()
        > plot(redshift, SFR,'ro')

        # to plot Chary & Elbaz 2001 SEDs
        > ce.plotseds()

    PROCEDURE

    REQUIRED PYTHON MODULES
        scipy
        pylab
        atpy
        idlsave (required by ReadCharyElbazTemplates.py)

    ADDITIONAL REQUIRED MODULES
        ReadCharyElbazTemplates.py
        astrofuncs.py

    NOTES
      Required Files
        z_response directory, which contains filter traces
        
      You must change the chary_path *in this file * to the directory that contains
      the z_response directory (probably where you unpacked the tar files) e.g.
      chary_path='/Users/rfinn/idl/programs/chary_elbaz_codes/'




      UPDATES: 
      - created readzresponse(file) after I upgraded to maverickcs and atpy and atropy.io.ascii would no longer read the z_response files (12/16/13)

'''

import astrofuncs
# read in chary_elbaz.save file
import ReadCharyElbazTemplates

from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.integrate import trapz
from pylab import *
#import atpy
import numpy as np
import os



# python magic to get environment variable
homedir = os.environ['HOME']

chary_path=homedir+'/github/Virgo/CharyElbaz/CharyElbaz01/chary_elbaz_codes/'

savefile=chary_path+'chary_elbaz.save'
ce=ReadCharyElbazTemplates.charyelbaz(savefile)
#read and fit iras filter responses

def readzresponse(file):
    input=open(file,'r')
    x=[]
    y=[]
    for line in input:
        if line.find('#') > -1:
            continue
        t=line.split()
        x.append(t[0])
        y.append(t[1])
    input.close()
    x=np.array(x,'f')
    y=np.array(y,'f')
    return x,y


infile=chary_path+'z_response/iras12.txt'
#mdata=atpy.Table(infile,type='ascii')
#mdata=ascii.read(infile)#,type='ascii')
#wave12=mdata['col1'] # wavelength
#resp12=mdata['col2'] # filter response
wave12,resp12=readzresponse(infile)

infile=chary_path+'z_response/iras25.txt'
#mdata=atpy.Table(infile,type='ascii')
#wave25=mdata['col1'] # wavelength
#resp25=mdata['col2'] # filter response
wave25,resp25=readzresponse(infile)

infile=chary_path+'z_response/iras60.txt'
#mdata=atpy.Table(infile,type='ascii')
#wave60=mdata['col1'] # wavelength
#resp60=mdata['col2'] # filter response
wave60,resp60=readzresponse(infile)

infile=chary_path+'z_response/iras100.txt'
#mdata=atpy.Table(infile,type='ascii')
#wave100=mdata['col1'] # wavelength
#resp100=mdata['col2'] # filter response
wave100,resp100=readzresponse(infile)

def chary_elbaz_24um(redshift,SmicroJy, H0=70., o_m=0.3, o_l=0.7):
    '''
    PURPOSE
    CALLING SEQUENCE
    INPUT PARAMETERS
    OUTPUT PARAMETERS
    EXAMPLES
    PROCEDURE
    NOTES

    '''
    try:
        t=len(redshift)
        redshift=array(redshift,dtype='d')  # allows redshift to be a single element or an array of values
        SmicroJy=array(SmicroJy, dtype='d')
    except:
        redshift=array([redshift],dtype='d')  # allows redshift to be a single element or an array of values
        SmicroJy=array([SmicroJy], dtype='d')

    #print redshift, type(redshift), SmicroJy
    try:
        ngals=len(redshift) 
        if ngals != len(SmicroJy):
            print len(redshift),len(SmicroJy)
            print'ERROR in dimensions'
            return
    except TypeError:
        print 'looks like I got a single number'
    
    #print redshift
    distance=astrofuncs.dL(redshift, H0/100.)

    lambda_0 = 23.675               # microns
    c_light_um = 2.99793e14         # microns/s
    
    Slambda_seds_obs=zeros([105,1366],'d')
    Snu_seds_obs=zeros([105,1366],'d')

    #readcol, 'z_response/mips24lg.dat', l_r24, r_r24, format='D,D', /silent
    #
    # replaced with the following 3 lines
    infile=chary_path+'z_response/mips24lg.dat'
    #mdata=atpy.Table(chary_path+'z_response/mips24lg.dat',type='ascii')
    #l_r24=mdata['col1'] # wavelength
    #r_r24=mdata['col2'] # filter response
    l_r24,r_r24=readzresponse(infile)
    #print l_r24
    NN_r24 = len(l_r24)
    Slambda_seds_24um = zeros([105, NN_r24],'f')

    T_corps_noir = 10000            # K
    corps_noir = astrofuncs.planck(l_r24*1e4, T_corps_noir)
    #print 'corps_noir[0:10] = ',corps_noir[0:10]
    num_seds = zeros(105,'d')
    #print l_r24
    func4=interp1d(l_r24*1e4,corps_noir) # save for later use w/in loops
    
    #    denom_seds = int_tabulated(l_r24, r_r24*corps_noir)
    # the following command replaces the int_tabulated() function in idl
    denom_seds=trapz(r_r24*corps_noir, l_r24)

    Slambda_seds_mips24 = zeros(105,'f')
    Snu = zeros([ngals,1366],'f')
    luminosity =  zeros([ngals,1366],'f')
    LIR_Sanders=zeros(ngals,'f')
    SFR_Sanders=zeros(ngals,'f')
    S24um=zeros(ngals,'f')

    for igal in range(ngals):

        #if (keyword_set(verbose) and igal mod 10 eq 0):
        #    print,igal,'/',ngals

        # SmicroJy in 10^-29 erg/s/cm^2/Hz
        Slambda = 1e-29*1e-4*(c_light_um/lambda_0**2)*SmicroJy[igal] # in erg/s/cm^2/Angstrom

        lambda_obs   = ce.lamb*(1e0+redshift[igal])

        for i in range(105):
            Snu_seds_obs[i,:] = (1.e0+redshift[igal])*ce.nulnuinlsun[:,i]/(1.e-32*4*pi*(distance[igal]*3.0856e22)**2*(3.e14/ce.lamb)/3.826e26) # in microJy
            c0=1./(4*pi*(distance[igal]*3.0856e24)**2)
            c2=c0*(1e4*ce.lamb)/(1.e0+redshift[igal]) # erg/s/cm**2/Angstrom
            Slambda_seds_obs[i,:]= 3.826e33*ce.nulnuinlsun[:,i]/(4.*pi*(distance[igal]*3.0856e24)**2*(1.e4*ce.lamb))/(1.e0+redshift[igal]) # erg/s/cm**2/Angstrom
            ## spectrum of  each template interpolated at the 24um response curve wavlelengths
            # Slambda_seds_24um[i,:] = interpol(Slambda_seds_obs[i, *], lambda_obs, l_r24)
            func3=interp1d(lambda_obs,Slambda_seds_obs[i,:])
            Slambda_seds_24um[i,:]=func3(l_r24)

            ## integration over the 24 response curve
            # num_seds[i] = int_tabulated(l_r24, r_r24*Slambda_seds_24um(i, *))
            num_seds[i]=trapz(r_r24*Slambda_seds_24um[i,:], l_r24)

        ## to scale to a black body
        norm_seds = num_seds/denom_seds 

        ## model flux density at 24um (over the filter)
        #for i=0,104 do Slambda_seds_mips24[i] = norm_seds[i]*interpol(corps_noir, l_r24*1e4, lambda_0*1e4) # in erg/s/cm**2/Angstrom
        #print lambda_0
        Slambda_seds_mips24 = norm_seds*func4(lambda_0*1.e4)

        # ind_sed=interpol(arange(105,dtype='f'),Slambda_seds_mips24,Slambda) ## get index for this observed flu
        if Slambda < min(Slambda_seds_mips24):
            print 'Warning: attempt to extrapolate the 24um luminosity outside the model limits! (flux too low)'
            ind_sed=0
            scalefactor=1.*Slambda/min(Slambda_seds_mips24)
        elif Slambda > max(Slambda_seds_mips24):
            print 'Warning: attempt to extrapolate the 24um luminosity outside the model limits! (flux too high)'
            ind_sed=104
            scalefactor=1.*Slambda/max(Slambda_seds_mips24)
        else:
            func5=interp1d(Slambda_seds_mips24,arange(105))
            ind_sed=func5(Slambda)
            scalefactor=1
        #  Snu[igal,:]=reform(interpolate(Snu_seds_obs,ind_sed,findgen(1366),/grid))
        # this command interpolates between seds to get an sed appropriate for ind_sed (of length 1366)
        # then it 'reform's the array to make it into a 1-d array
        
        # I am going to rewrite it to pull the nearest sed
        nearest_sed=int(round(ind_sed))
        Snu[igal,:]=Snu_seds_obs[nearest_sed]


        f12_rf = iras(ce.lamb, Snu[igal,:]/(1e0+redshift[igal]),band=12)*scalefactor
        f25_rf = iras(ce.lamb, Snu[igal,:]/(1e0+redshift[igal]),band=25)*scalefactor
        f60_rf = iras(ce.lamb, Snu[igal,:]/(1e0+redshift[igal]),band=60)*scalefactor
        f100_rf = iras(ce.lamb, Snu[igal,:]/(1e0+redshift[igal]),band=100)*scalefactor


        LIR_Sanders[igal] = 1.8e-14*1e-6*(13.48*f12_rf + 5.16*f25_rf + 2.58*f60_rf + f100_rf)*4*pi*(distance[igal]*3.0856e22)**2/3.826e26
        SFR_Sanders[igal] = 1.7217e-10*LIR_Sanders[igal]

        S24um[igal] = SmicroJy[igal]

    return LIR_Sanders,SFR_Sanders

def iras(mylambda, flux, band=12): # band is 12, 25, 60, 100
    '''
    PURPOSE
      Integrated observed flux over a given IRAS band
    CALLING SEQUENCE
      iras(mylambda, flux, band)

      where band is the iras band and can be 12, 25, 60 or 100

    INPUT PARAMETERS
      Wavelength - in microns
      flux - in Snu

    OUTPUT PARAMETERS
    EXAMPLES
    PROCEDURE
    NOTES

    '''
    
    if band == 12:
        lambda0 = 12.e0
        wave=wave12
        resp=resp12
    elif band == 25:
        lambda0 = 25.e0
        wave=wave25
        resp=resp25
    elif band == 60:
        lambda0 = 60.e0
        wave=wave60
        resp=resp60 
    elif band == 100:
        lambda0 = 100.e0
        wave=wave100
        resp=resp100
    else:
        print 'Error reading band.  should be 12, 25, 60 or 100'
        return

    # flux_wave = interpol(flux, mylambda, wave)
    func1=interp1d(mylambda,flux)
    flux_wave=func1(wave)

    fct_numerateur = flux_wave*resp/wave**2
    # int_numerateur = int_tabulated(wave, fct_numerateur)
    int_numerateur = trapz(fct_numerateur, wave)
    
    fct_denumerateur = resp/(wave*lambda0)
    # int_denumerateur = int_tabulated(wave, fct_denumerateur)

    int_denumerateur = trapz(fct_denumerateur, wave)

    result = int_numerateur/int_denumerateur

    return result

