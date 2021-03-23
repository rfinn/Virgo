#!/usr/bin/env python

'''
GOAL
* measure density profiles of filaments in cylindrical volumes along length
* fit exponential profiles


'''

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import os
homedir = os.getenv("HOME")
sys.path.append(homedir+'/github/Virgo/programs/')

import virgoCommon


#########################################################
####  READ IN VIRGO TABLES
#########################################################
import readtables
TABLEDIR = '/home/rfinn/research/Virgo/tables-north/v1/'
TABLEPREFIX = 'vf_north_v1_'
v = readtables.vtables(TABLEDIR,TABLEPREFIX) 
v.read_filaments() # sets up v.fil table
v.read_env() # sets up v.env table

#########################################################
####  FUNCTIONS
#########################################################

def profile_fun(x,x0,a,b):
    """
    PURPOSE:
    analytical function for the density profile
    
    INPUT:
    x : radial distance
    x0 : characteristic exponential scale
    a  : coeffiecient of the exponential function, a+b is the density value at x<<x0
    b : density value at x>>x0
        
    OUTPUT:
    density value: rho = a*exp(-x/x0)+b
    """
    result = a*np.exp(-x/x0)+b
    return result


def density_profile_fit_proc(rho,distance,N=None,p0=[1.,1.,0.]):
    """
    PURPOSE:
    macro that performs the density profile fit
    
    INPUT:
    rho: number density (counts / volume)
    distance : radial distance from the filament spine
    N : number of sources counted within 'distance'
    p0 : guess for the best fit parameters : p0 = [x0,a,b] as in 'profile_fun'


    OUTPUT:
    best fit parameter array and associated error array:  
                                                          
    par_best = [x0,a,b]                                   
    err_par_best = [err_x0,err_a,err_b]                   
                                                          
    adopted function: density profile as in 'profile_fun'
    """
    p0 = np.array(p0)
    rho      = np.array(rho)
    distance = np.array(distance)
    if N is not None: N = np.array(N)
    ####
    if N is not None:
        err_rho = np.sqrt(N)/N*rho
        cond = (err_rho == 0)
        err_rho[cond] = 1./N[cond]*rho[cond]
    else:
        err_rho = np.array([1.]*len(rho))
    ######
    par_best, par_cov = curve_fit(profile_fun,xdata=distance, ydata=rho, p0=p0, sigma=err_rho) #best fit parameters and covariance matrix
    err_par_best = np.diag(par_cov)**0.5 #best fit parameter errors
    return par_best,err_par_best


def test_best_exp_fit():
    """
    PURPOSE:
    test macro that perform a test with an exponential fit for the profile
    """
    #####
    #x and y values for the test
    x_values = np.arange(0,6,6./30)
    y_values = profile_fun(x_values,x0=1.5,a=1.,b=0.5)+(np.random.uniform(size=len(x_values))-0.5)/2.
    ######
    print('I find the best fit curve parameters')
    print('insert N in density_profile_fit_proc, below, for the real case')
    par_best,err_par_best  = density_profile_fit_proc(rho=y_values,distance=x_values,N=None,p0=[1.,1.,0.]) #I find the best fit curve parameters
    ####
    print('x0 = '+str(round(par_best[0],3))+' +/- '+str(round(err_par_best[0],3)))
    print('a = '+str(round(par_best[1],3))+' +/- '+str(round(err_par_best[1],3)))
    print('b = '+str(round(par_best[2],3))+' +/- '+str(round(err_par_best[2],3)))
    ######
    #best fit curve
    xcurve = np.arange(-1,8,0.001)
    ycurve = profile_fun(xcurve,x0=par_best[0],a=par_best[1],b=par_best[2])
    ######
    plt.close('all')
    plt.scatter(x_values,y_values)
    plt.plot(xcurve,ycurve,color='red')
    plt.xlabel('distance')
    plt.ylabel('density')
    plt.show()
    return


def plot_profiles():
    os.chdir(plotdir) 
    plt.figure(figsize=(10,8))
    plt.subplots_adjust
    for nfil,f in enumerate(virgoCommon.filaments):
        print('##############################')
        print('###  ',f)
        print('##############################')
        if nfil > len(virgoCommon.mycolors)-1:
            ls = '--'
            marker='^'
            alpha=.7
        else:
            ls = '-'
            marker='o'
            alpha=1
        plt.plot(dist,dens,label=f,marker=marker,ls=ls,alpha=alpha)
    plt.gca().set_yscale('log')
    plt.ylim(.06,15)
    plt.legend(fontsize=10,loc='upper right')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(r'$Distance \ (h_{74}^{-1} \ Mpc)$',fontsize=20)
    plt.ylabel(r'$ \rho _{gal}\ (h_{74}^{-1} \ Mpc)^{-3}$',fontsize=20)


    plt.savefig('density-profiles.png')
    plt.savefig('density-profiles.pdf')    

#########################################################
####  CLASSES
#########################################################

class filament():
    ''' class for profiles, to measure density profile and fit exp to each '''
    def __init__(self,filament_name,marker=None,ls=None,alpha=None):
        ''' take filament name as input '''
        self.filname = filament_name
        if marker is not None:
            self.marker = marker
        if ls is not None:
            self.ls = ls
        if alpha is not None:
            self.alpha = alpha
        self.name = filament_name
        pass
    def measure_profile(self,dmin=0.25,dmax=3,stepsize=.2):
        '''
        PURPOSE: to measures the density profile 
        
        INPUT:
        * dmin - first bin for computing density
        * dmax
        * stepsize

        OUTPUT:
        * self.ngal - number of galaxies in each bin
        * self.density - ngal/volume
        * self.distance - distance bins, ymax of bin

        '''
        dist = np.arange(dmin,dmax+stepsize,stepsize)
        dens = np.zeros(len(dist),'f')
        length = virgoCommon.fil_lengths[self.filname]
        self.length = length
        colname = "dist_3D_{}".format(self.filname)
        self.ngal = np.zeros(len(dist),'d')
        self.density = np.zeros(len(dist),'d')
        self.density_err = np.zeros(len(dist),'d')
        self.volume = np.zeros(len(dist),'d')                
        try:
            t = v.fil[colname]
        except KeyError:
            print('\t WARNING: no column named ',colname)
            return
        for i,d in enumerate(dist):
            flag =  (v.fil[colname] < d) & (v.env['flag_clus'] == 0)
            self.ngal[i] = sum(flag)
            self.volume[i] = (np.pi*d**2*length)
            self.density[i] = sum(flag)/self.volume[i]
            self.density_err[i] = np.sqrt(self.ngal[i])/self.volume[i]
            
        self.distance = dist
    def calc_local_density(self):

        self.local_ngal = (self.ngal[1:]-self.ngal[0:-1])
        self.local_volume = (self.volume[1:]-self.volume[0:-1])
        self.local_density = self.local_ngal/self.local_volume
        self.local_distance = (self.distance[1:]+self.distance[0:-1])/2
        self.local_density_err = np.sqrt(self.local_ngal)/self.local_volume        
    def fit_profile(self):
        '''
        PURPOSE: use Gianluca's functions to fit an exponential profile to each
        '''
        self.par_best,self.err_par_best =  density_profile_fit_proc(self.density,self.distance,N=self.ngal,p0=[1.,1.,0.])
    def fit_local_profile(self):
        '''
        PURPOSE: use Gianluca's functions to fit an exponential profile to each
        '''
        self.par_best_local,self.err_par_best_local =  density_profile_fit_proc(self.local_density,self.local_distance,N=self.local_ngal,p0=[1.,1.,0.])

    def plot_density(self):
        '''
        plot best-fit profile
        '''
        plt.plot(self.distance,self.density,marker='o')
        plt.errorbar(self.distance,self.density,yerr=self.density_err)
    def plot_local_density(self):
        '''
        plot best-fit profile
        '''
        plt.plot(self.local_distance,self.local_density,marker='o')
        plt.errorbar(self.local_distance,self.local_density,yerr=self.local_density_err)
        
    def plot_exp_fit(self):
        xl = np.linspace(min(self.distance),max(self.distance),100)
        x0 = self.par_best[0]
        a = self.par_best[1]
        b = self.par_best[2]              

        yl = a*np.exp(-1*xl/x0)+b
        plabel = '{:.1f}exp(-x/{:.1f})+{:.1f}'.format(a,x0,b)
        plt.plot(xl,yl,'k--',label=plabel)
                         
        pass
    def plot_local_exp_fit(self):
        xl = np.linspace(min(self.local_distance),max(self.local_distance),100)
        x0,a,b = self.par_best_local

        yl = a*np.exp(-1*xl/x0)+b
        plabel = '{:.1f}exp(-x/{:.1f})+{:.1f}'.format(a,x0,b)
        plt.plot(xl,yl,'k--',label=plabel)
                         
        pass
    def measure_concentration(self):
        self.conc = self.density[4]/self.density[-1]
    def plot_single(self):
        plt.figure(figsize=(10,8))
        self.plot_density()
        self.plot_exp_fit()
        ftitle = "{} length={:.1f}".format(self.name,self.length)
        plt.title(ftitle,fontsize=20)
        plt.xlabel(r'$Distance \ (Mpc/h) $',fontsize=20)
        plt.ylabel(r'$\rho(<r) \ (Mpc/h)^{-3} $',fontsize=20)        
        plt.legend(fontsize=18)
    def plot_localdens(self,plotsingle=True):
        if plotsingle:
            plt.figure(figsize=(10,8))
            fsize=20
        else:
            fsize=12
        self.plot_local_density()
        self.plot_local_exp_fit()
        ftitle = "{} (length={:.1f} Mpc)".format(self.name,self.length)
        plt.title(ftitle,fontsize=fsize)
        plt.xlabel(r'$Distance \ (Mpc/h) $',fontsize=fsize)
        plt.ylabel(r'$\rho(r) \ (Mpc/h)^{-3} $',fontsize=fsize)        
        plt.legend(fontsize=fsize)
        
#########################################################
####  MAIN PROGRAM
#########################################################

if __name__ == '__main__':
    #initiate filament classes
    testing = np.arange(len(virgoCommon.filaments))
    allfilaments = []
    alllengths = []
    allscalelengths = []
    allerr = []
    allconc = []
    plt.figure(figsize=(10,14))
    plt.subplots_adjust(hspace=.5,wspace=.3)
    for i in testing:
        f = filament(virgoCommon.filaments[i])
        f.measure_profile(dmax=6)
        f.calc_local_density()
        #plt.figure()
        #f.plot_density()
        
        f.fit_profile()
        f.fit_local_profile()
        #f.plot_exp_fit()
        #plt.title(f.name)
        #f.plot_single()
        plt.subplot(7,2,i+1)
        f.plot_localdens(plotsingle=False)
        if i < len(testing)-2:
            #plt.xticks([])
            plt.xlabel('')

        #plt.axis([0,6.1,.01,9])

        plt.xlim(0.,4)
        plt.gca().set_yscale('log')
        f.measure_concentration()
        alllengths.append(f.length)
        allscalelengths.append(f.par_best[0])
        allerr.append(f.err_par_best[0])
        allconc.append(f.conc)
        allfilaments.append(f)
    # measure density profile and fit exp for each filament
    
    # plot profiles
    pass
