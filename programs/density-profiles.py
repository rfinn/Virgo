#!/usr/bin/env python

'''
GOAL
* measure density profiles of filaments in cylindrical volumes along length
* fit exponential profiles


'''

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
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
v.read_all()
#v.read_filaments() # sets up v.fil table
#v.read_env() # sets up v.env table

#########################################################
####  FUNCTIONS
#########################################################

def run_spearman(x,y):
    t = spearmanr(x,y)
    print('spearman rank test: rho={:.2f}, p={:.4f}'.format(t[0],t[1]))
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
        self.name = self.filname.replace('_Filament','').replace('Virgo_','').replace('_',' ')        
        if marker is not None:
            self.marker = marker
        if ls is not None:
            self.ls = ls
        if alpha is not None:
            self.alpha = alpha

        pass
    def measure_profile(self,dmin=0.25,dmax=3,stepsize=.2,magcut=False,samemag=False):
        '''
        PURPOSE: to measures the density profile 
        
        INPUT:
        * dmin - first bin for computing density
        * dmax
        * stepsize
        * magcut - boolean, apply a mag cut when calculating density.
        * samemag - boolean; True=use the same Mr=-15.7 cut for all; False=use Mr cut tuned to each filament

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
            flag_mag = v.rphot['M_r']<v.rphot['M_lim_values_fil']
            if magcut:
                if samemag:
                    flag = flag & (v.rphot['M_r'] < -15.7)
                else:
                    flag = flag & flag_mag
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
        self.conc_err = self.density_err[4]/self.density[-1]        
    def plot_single(self):
        plt.figure(figsize=(10,8))
        self.plot_density()
        self.plot_exp_fit()
        ftitle = "{} length={:.1f}".format(self.name,self.length)
        # updating to remove length, remove _Filament, replace_ with " "
        #ftitle = "{} length={:.1f}".format(self.name,self.length)
        ftitle = "{} length={:.1f}".format(self.name.replace("_Filament","").replace("Virgo_","").replace("_"," "))        
        plt.title(ftitle,fontsize=20)
        plt.xlabel(r'$r \ (h^{-1}~Mpc) $',fontsize=20)
        plt.ylabel(r'$\rho(<r) \ (h^{-1}~Mpc)^{-3} $',fontsize=20)        
        #plt.legend(fontsize=18)
    def plot_localdens(self,plotsingle=True,showx=False,showy=False):
        if plotsingle:
            plt.figure(figsize=(10,8))
            fsize=20
        else:
            fsize=12
        self.plot_local_density()
        self.plot_local_exp_fit()
        #ftitle = "{} (length={:.1f} Mpc)".format(self.name,self.length)
        ftitle = "{}".format(self.name.replace("_Filament","").replace("Virgo_","").replace("_"," "))                
        #plt.title(ftitle,fontsize=fsize)
        plt.text(0.95,.85,ftitle,transform=plt.gca().transAxes,horizontalalignment='right',fontsize=12)
        if showx:
            plt.xlabel(r'$r \ (h^{-1}~Mpc) $',fontsize=fsize)
        else:
            #plt.xticks([])
            plt.gca().tick_params(labelbottom=False)
        if showy:
            plt.ylabel(r'$\rho \ (h^{-1}~Mpc)^{-3} $',fontsize=fsize)
        else:
            plt.yticks([],[])
        #plt.legend(fontsize=fsize)

class virgofilaments:
    def __init__(self,filament_list):
        '''  pass in list of filament instances from main program '''
        if args.magcut & args.samemag:
            self.suffix = 'fixedMr'
        elif args.magcut & ~args.samemag:
            self.suffix = 'filament_Mr'
        else:
            self.suffix = 'no_Mr_cut'
        self.allfilaments = filament_list
        
    def plot_scale_length_vs_length(self,ymax=None):
        ''' plot scale length vs length '''
        plt.figure(figsize=(8,6))

        # update to use BV colors
        for i,f in enumerate(allfilaments):
            plt.plot(f.length,f.par_best_local[0],'bo',color=virgoCommon.mycolors[i],markersize=14,label=f.name)
            plt.errorbar(f.length,f.par_best_local[0],yerr=f.err_par_best_local[0],color=virgoCommon.mycolors[i])
        plt.xlabel("$r \ (h^{-1}~Mpc)$",fontsize=20)
        plt.ylabel("$r_0 \ (h^{-1}~Mpc)$",fontsize=20)
        plt.legend()
        if ymax is not None:
            plt.ylim(0,ymax)
        plt.savefig('scale_length_vs_L_'+self.suffix+'.pdf')
        plt.savefig('scale_length_vs_L_'+self.suffix+'.png')
        
    def plot_central_density_vs_length(self,ymax=None):
        ''' plot central density vs length '''
        plt.figure(figsize=(8,6))

        # update to use BV colors
        for i,f in enumerate(allfilaments):
            plt.plot(f.length,f.par_best_local[1],'bo',color=virgoCommon.mycolors[i],markersize=14,label=f.name)
            plt.errorbar(f.length,f.par_best_local[1],yerr=f.err_par_best_local[1],color=virgoCommon.mycolors[i])
        plt.xlabel("$r \  (h^{-1}~Mpc)$",fontsize=20)
        plt.ylabel("$a \ (h~Mpc^{-1})^3$",fontsize=20)
        plt.legend()
        if ymax is not None:
            plt.ylim(0,ymax)
        plt.savefig('central_density_vs_L_'+self.suffix+'.pdf')
        plt.savefig('central_density_vs_L_'+self.suffix+'.png')

    def plot_density_contrast_vs_length(self,ymax=None):
        ''' plot central density vs length '''
        plt.figure(figsize=(8,6))

        # update to use BV colors
        for i,f in enumerate(allfilaments):
            y = f.par_best_local[1]/f.par_best_local[2]
            yerr = f.err_par_best_local[1]/f.par_best_local[2]            
            plt.plot(f.length,y,'bo',color=virgoCommon.mycolors[i],markersize=14,label=f.name)
            plt.errorbar(f.length,y,yerr=yerr,color=virgoCommon.mycolors[i])
        plt.xlabel("Length",fontsize=20)
        plt.ylabel("Density Contrast",fontsize=20)
        plt.legend()
        if ymax is not None:
            plt.ylim(0,ymax)
        plt.savefig('density_contrast_vs_L_'+self.suffix+'.pdf')
        plt.savefig('density_contrast_vs_L_'+self.suffix+'.png')

    def plot_central_density_contrast_vs_length(self,ymax=None):
        ''' plot central density vs length '''
        plt.figure(figsize=(8,6))

        # update to use BV colors
        for i,f in enumerate(allfilaments):
            y = f.density[4]/f.par_best_local[2]
            yerr = f.density_err[4]/f.par_best_local[2]            
            plt.plot(f.length,y,'bo',color=virgoCommon.mycolors[i],markersize=14,label=f.name)
            plt.errorbar(f.length,y,yerr=yerr,color=virgoCommon.mycolors[i])
        plt.xlabel("$r \ (h^{-1}~Mpc)$",fontsize=20)
        plt.ylabel(r"$\rho (r < 1~h^{-1}~Mpc)/b$",fontsize=20)
        plt.legend()
        if ymax is not None:
            plt.ylim(0,ymax)
        plt.savefig('central_density_contrast_vs_L_'+self.suffix+'.pdf')
        plt.savefig('central_density_contrast_vs_L_'+self.suffix+'.png')
    def plot_concentration_vs_length(self,ymax=None):
        ''' plot central density vs length '''
        plt.figure(figsize=(8,6))

        # update to use BV colors
        for i,f in enumerate(allfilaments):
            y = f.conc
            yerr = f.conc_err        
            plt.plot(f.length,y,'bo',color=virgoCommon.mycolors[i],markersize=14,label=f.name)
            plt.errorbar(f.length,y,yerr=yerr,color=virgoCommon.mycolors[i])
        plt.xlabel("$Length \ (h^{-1}~Mpc)$",fontsize=20)
        plt.ylabel("Concentration",fontsize=20)
        plt.legend()
        if ymax is not None:
            plt.ylim(0,ymax)
        plt.savefig('concentration_vs_L_'+self.suffix+'.pdf')
        plt.savefig('concentration_vs_L_'+self.suffix+'.png')

#########################################################
####  MAIN PROGRAM
#########################################################

if __name__ == '__main__':
    ###########################
    ##### SET UP ARGPARSE
    ###########################
    import argparse
    parser = argparse.ArgumentParser(description ='Program to plot radial density profiles of filaments')
    parser.add_argument('--magcut', dest = 'magcut', default = False,action='store_true', help = 'apply Mr cut when counting galaxies.  default is False')
    parser.add_argument('--samemag', dest = 'samemag', default = False, action='store_true', help = 'Set this to use the same Mr cut for each galaxy (-15.7).  otherwise use Mr cut that is tuned to each galaxy.')
    parser.add_argument('--stepsize', dest = 'stepsize', default = 0.2, help = 'Set step size for binning in radial direction.  the default is 0.2')    

    args = parser.parse_args()
    
    #initiate filament classes
    testing = np.arange(len(virgoCommon.filaments))
    allfilaments = []
    alllengths = []
    allscalelengths = []
    alldensitycontrasts = []
    allerr = []
    allconc = []
    plt.figure(figsize=(10,12))
    plt.subplots_adjust(hspace=.1,wspace=.15,bottom=.07,left=.1,top=.95,right=.95)
    for i in range(len(virgoCommon.filaments)):
        f = filament(virgoCommon.filaments[i])
        f.measure_profile(dmax=6,magcut=args.magcut,samemag=args.samemag,stepsize=float(args.stepsize))
        f.calc_local_density()
        #plt.figure()
        #f.plot_density()
        
        f.fit_profile()
        f.fit_local_profile()
        #f.plot_exp_fit()
        #plt.title(f.name)
        #f.plot_single()
        plt.subplot(7,2,i+1)
        plotx=True
        ploty=True
        if i < (len(testing)-2):
            plotx=False
        if i%2 == 1:
            ploty=False
        f.plot_localdens(plotsingle=False,showy=ploty,showx=plotx)
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
        alldensitycontrasts.append(f.par_best_local[1]/f.par_best_local[2])
    # measure density profile and fit exp for each filament
    
    # plot profiles
    pass
