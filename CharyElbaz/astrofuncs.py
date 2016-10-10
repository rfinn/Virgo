#!/usr/bin/env python
from pylab import *
from scipy.integrate import romberg
from scipy.integrate import quad
cl=3.e5
OmegaL = 0.7
OmegaM = 0.3

def planck(wave,temp):
    ''' based on idlutil planck function
    PURPOSE: 
      To calculate the Planck function in units of ergs/cm2/s/A  

    CALLING SEQUENCE: 
      bbflux = PLANCK( wave, temp) 

    INPUT PARAMETERS: 
        WAVE   Scalar or vector giving the wavelength(s) in **Angstroms**
              at which the Planck function is to be evaluated.
        TEMP   Scalar giving the temperature of the planck function in degree K

    OUTPUT PARAMETERS:
        BBFLUX - Scalar or vector giving the blackbody flux (i.e. !pi*Intensity)
              in erg/cm^2/s/A in at the specified wavelength points.

    EXAMPLES:
      To calculate the blackbody flux at 30,000 K every 100 Angstroms between
      2000A and 2900 A

      import idlutils
      bbflux = idlutils.planck(wave,30000)

      If a star with a blackbody spectrum has a radius R, and distance,d, then
      the flux at Earth in erg/cm^2/s/A will be bbflux*R^2/d^2

    PROCEDURE:
      The wavelength data are converted to cm, and the Planck function
      is calculated for each wavelength point. See Allen (1973), Astrophysical
      Quantities, section 44 for more information.

    NOTES:
      See the procedure planck_radiance.pro in 
      ftp://origin.ssec.wisc.edu/pub/paulv/idl/Radiance/planck_radiance.pro
      for computation of Planck radiance given wavenumber in cm-1 or  
      wavelength in microns 
    ''' 
    w = wave / 1.e8                   # Angstroms to cm    ;constants appropriate to cgs units.
    c1 =  3.7417749e-5                # =2*!DPI*h*c*c       
    c2 =  1.4387687                  # =h*c/k
    val =  c2/w/temp
    bbflux=c1/(w**5*(exp(val)-1))
    return bbflux*1.e-8              # Convert to ergs/cm2/s/A

def dL2(z,h):
    c=3.e5
    H0=100.*h
    if type(z) == ndarray:
	n=len(z)
	DL=zeros(n,'d')
	for i in range(n):
	    s=romberg(func,0.,z[i])#,tol=1.e-6)
	    DL[i]=c/H0*(1+z[i])*s #(Mpc/h)
    else:
	s=romberg(func,0.,z)#,tol=1.e-6)
	DL=c/H0*(1+z)*s #(Mpc/h)
    return DL
def dL(z,h):
    c=3.e5
    H0=100.*h

    try:#multiple objects
	n=len(z)
	DL=zeros(n,'d')
	for i in range(n):
	    s=romberg(func,0.,z[i])#,tol=1.e-6)
	    DL[i]=c/H0*(1+z[i])*s #(Mpc/h)
    except TypeError:
	s=romberg(func,0.,z)#,tol=1.e-6)
	DL=c/H0*(1+z)*s #(Mpc/h)
    return DL


def dLcm(z,h):
    DL=dL(z,h)*1.e6*3.08568025e18#convert to cm
    return DL

def lookbackt(z,h):
    c=3.e5
    H0=100.*h
    tH=9.78*h
    try:
        n=len(z)
        t=zeros(n,'d')
        for i in range(n):
            t[i]=quad(funct,0.,z[i])
    except TypeError:
        t=quad(funct,0.,z)
        print 't = ',t[0]
    t=t*tH #time in Gyr
    return t,2*t

def E(z):
    E = N.sqrt(OmegaL + OmegaM*(1.+z)**3)
    return E

def func(x):
    return 1.0/(sqrt(OmegaM*(1+x)**3 + OmegaL))

def funct(x):#for lookback time
    return 1.0/(sqrt(OmegaM*(1+x)**3 + OmegaL)*(1+x))

def DA(z,h):#uses numerical integration technique gaussian quadrature
    c=3.e5
    H0=100.*h
    s=romberg(func,0.,z) 
    DA=c/H0/(1+z)*s     #(Mpc/radian)
    DA=DA*1000./206264  #(kpc/arcsec)    
    return DA
