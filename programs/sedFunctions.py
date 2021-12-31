#!/usr/bin/env python

'''
GOAL:
* This program contains tools to visualize:
  * the fluxes fed into magphys
  * the resulting SED and pdfs from magphys
* In the magphys_sed class, I have translated much of plot_sed.pro, the idl program that comes in the magphys tarball

USEAGE:

From the directory containing the *.sed and *.fit files:

effective_wavelengths = np.array([ 0.1516,0.2267,0.48623,0.64606,0.91993,3.40025,4.65201,12.81034,22.37528],'d')
%run ~/github/Virgo/programs/sedFunctions.py
s = magphys_sed('2681',all_effective_wavelengths)
s.plot_sed()
s.plot_histograms()


TO DO:
* allow user to pass the filter file to magphys_sed, and read the effective wavelengths from there


WRITTEN BY:
Rose Finn, August 2021
'''


from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from astropy.table import Table

def plotsed(wave,flux,flux_err,title=None):
    ''' 
    plot sed of one galaxy 
    
    INPUT:
    * wave - effective wavelength of photometric points in microns
    * flux - flux of photometric points in Jy
    * flux_err - error in flux, in Jy

    OUTPUT:
    * sed plot
    '''
    plt.figure()
    frequency = 3.e8/(wave*1.e-4)
    plt.plot(wave,frequency*flux,'bo',markersize=6)
    plt.errorbar(wave,frequency*flux,yerr=frequency*flux_err,fmt='.',color='b')
    plt.xlabel('Wavelength (um)',fontsize=18)
    plt.ylabel('nu Flux (Jy-Hz) ',fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)    
    if title is not None:
        plt.title(title,fontsize=20)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.axhline(y=0,ls='--',color='k')
def plotsed_ABmag(wave,flux,flux_err,title=None):
    ''' plot sed of one galaxy in AB mag'''
    plt.figure()
    frequency = 3.e8/(wave*1.e-4)
    flux = 22.5 - 2.5*np.log10(flux/3.631e-6) # convert from Jy to nanomaggy, then to AB mag
    plt.plot(wave,flux,'bo',markersize=6)
    #plt.errorbar(wave,flux,yerr=flux_err,fmt='.',color='b')
    plt.xlabel('Wavelength (um)',fontsize=18)
    plt.ylabel('AB Mag',fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)    
    if title is not None:
        plt.title(title,fontsize=20)
    plt.gca().set_xscale('log')
    #plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    #plt.axhline(y=0,ls='--',color='k')
def plotallsed(wave,flux,flux_err,titles=None,ABmag=False):
    ''' 
    use this for a table with flux values for multiple galaxies 

    INPUT: 
    * wave - effective wavelengths of the L photometric points in microns
    * flux - (NxL) array for flux values in Jy, where N is the number of gals and L is number of phtometric points
    * flux_err - error in the flux values in Jy

    OPTIONAL INPUTS:
    * titles - a length N array with titles for each SED plot; this is usually the gal id
    * ABmag - plot SED in terms of AB mag instead of lambda F_lambda

    OUTPUT:
    * produces an SED plot for each galaxy

    '''
    wave = np.array(wave,'f')
    for i,f in enumerate(flux):
        if ABmag:
            plotsed_ABmag(wave,flux[i],flux_err[i],title='VFID'+str(titles[i]))
        else:
            plotsed(wave,flux[i],flux_err[i],title='VFID'+str(titles[i]))
            
def parse_pdf(lines):
    ''' special function to parse the pdf output from magphys galid.fit file  '''
    a = np.zeros(len(lines),'d')
    b = np.zeros(len(lines),'d')
    for i,l in enumerate(lines):
        t = l.split()
        a[i] = float(t[0])
        b[i] = float(t[1])
    return a,b
        

class magphys_sed():
    
    '''
    class to visualize output from magphys
    '''
    
    def __init__(self,galid,wavelengths):
        self.sed_file = '{}.sed'.format(galid)
        self.fit_file = '{}.fit'.format(galid)
        self.lambda_eff = np.array(wavelengths,'d')
        self.galid = galid
    def plot_sed(self,plot_unattenuated=True):
        fig = plt.figure(figsize=(8,6))
        plt.subplots_adjust(left=.2)
        gs = fig.add_gridspec(2,1, height_ratios=(7,2),left=.1,right=.9,bottom=.2,top=.9,hspace=.1)
        ax = fig.add_subplot(gs[0,0])
        resid_ax = fig.add_subplot(gs[1,0],sharex=ax)
        
        self.read_sed_file(self.sed_file,ax1=ax,plot_unattenuated=plot_unattenuated)
        self.read_fit_file(self.fit_file,ax1=ax,ax2=resid_ax)

        ax.set_xlim(.07,100)
        # these are the limits from plot_sed.pro
        # but my y values are offset considerably
        #plt.ylim(7.1,12)
        ax.set_ylim(7,10.5)                
        ax.set_xscale('log')
        ax.tick_params(axis='x',labelbottom=False)
        resid_ax.set_xscale('log')
        resid_ax.set_ylabel('(obs-model)/obs')
        resid_ax.axhline(y=0,ls='--',color='k')
        plt.sca(ax)
        resid_ax.set_xlabel(r'$Wavelength \ (\mu m) \ [observed \ frame]$',fontsize=14)
        plt.ylabel(r'$log(\lambda L_\lambda/L_\odot) $',fontsize=14)        
        plt.legend()#loc='lower left')
        s = 'VFID{}: logMstar = {:.2f}, logSFR = {:.2f}'.format(self.galid,np.log10(self.mstar),np.log10(self.sfr))
        plt.title(s,fontsize=14)
        outfile = 'VFID{}-magphys-sed.png'.format(self.galid)
        
        plt.savefig(outfile)
    def plot_sed_noresidual(self,plot_unattenuated=True):
        fig = plt.figure(figsize=(8,6))
        plt.subplots_adjust(left=.15)
        ax = plt.gca()
        #gs = fig.add_gridspec(2,1, height_ratios=(7,2),left=.1,right=.9,bottom=.2,top=.9,hspace=.1)
        #ax = fig.add_subplot(gs[0,0])
        #resid_ax = fig.add_subplot(gs[1,0],sharex=ax)
        
        self.read_sed_file(self.sed_file,ax1=ax,plot_unattenuated=plot_unattenuated)
        self.read_fit_file(self.fit_file,ax1=ax,plotresid=False)

        ax.set_xlim(.07,100)
        # these are the limits from plot_sed.pro
        # but my y values are offset considerably
        #plt.ylim(7.1,12)
        ax.set_ylim(8.5,10.25)                
        ax.set_xscale('log')
        #ax.tick_params(axis='x',labelbottom=False)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        ax.set_xlabel(r'$Wavelength \ (\mu m) \ [observed \ frame]$',fontsize=20)
        plt.ylabel(r'$log(\lambda L_\lambda/L_\odot) $',fontsize=20)        
        plt.legend(['model','NGC3938','Predicted'],fontsize=16)#loc='lower left')
        s = 'VFID{}: logMstar = {:.2f}, logSFR = {:.2f}'.format(self.galid,np.log10(self.mstar),np.log10(self.sfr))
        #plt.title(s,fontsize=14)
        outfile = 'VFID{}-magphys-sed-noresidual.png'.format(self.galid)
        
        plt.savefig(outfile)
    def read_sed_file(self,sed_file,ax1=None,plot_unattenuated=True):
        self.sed = Table.read(sed_file,data_start=2,format='ascii')
        self.wave = np.array(self.sed['col1'],'d')

        self.lum_at = np.array(self.sed['col2'],'d')
        self.lum_un = np.array(self.sed['col3'],'d')

        self.sed_wave = 10.**self.wave#lambda in AA
        #L_at=np.log10(sed_wave*10**lum_at) #total attenuated SED
        # doing the conversion in log space instead of the eqn above, which is from plot_sed
        self.L_at=self.wave+self.lum_at #total attenuated SED        
        #L_un=np.log10(sed_wave*10**lum_un) #unattenuated SED
        # doing the conversion in log space instead of the eqn above, which is from plot_sed        
        self.L_un=self.wave+self.lum_un #unattenuated SED
        self.sed_wave_um = self.sed_wave/1.e4 # convert to microns    

        ax1.plot(self.sed_wave_um,self.L_at,label='attenuated')
        if plot_unattenuated:        
            ax1.plot(self.sed_wave_um,self.L_un,label='unattenuated')
        

    def read_fit_file(self,fit_file,ax1=None,ax2=None,plotresid=True):
        in1 = open(fit_file,'r')
        fit_lines = in1.readlines()
        # flux is given as L_nu in units of L_sun/Hz
        # I am keeping syntax/variable names similar to plot_sed
        # but this is really a luminosity, or luminosity density?
        flux = np.array(fit_lines[2].split(),'d') # observed L
        e_flux = np.array(fit_lines[3].split(),'d') # error
        # not sure exactly what these are yet, except for z of course...
        i_sfh,i_ir,chi2,z = np.array(fit_lines[8].split(),'d')
        fmu_sfh,fmu_ir,mu,tauv,ssfr,mstar,Ldust,T_W,T_C_ISM,\
            xi_Ctot,xi_PAHtot,xi_MIRtot,xi_Wtot,\
            tvism,Mdust,SFR = np.array(fit_lines[10].split(),'d')
        # model flux is given as L_nu in units of L_sun/Hz
        # again, keeping naming convention like plot_sed.pro, but this is a L not f
        p_flux = np.array(fit_lines[12].split(),'d')
        in1.close()

        # convert to luminosities using eqns from plot_sed.pro
        frequency = 3.e8/(self.lambda_eff*1.e-6) # freq in Hz; lambda_eff is in microns
        L_flux=np.log10((1.+z)*flux*frequency) # log of nu L_nu
        # haven't checked all this yet, but assuming errors are correct
        L_eflux_lo=np.log10((1.+z)*flux*frequency)-np.log10((1.+z)*flux*frequency-e_flux*(1.+z)*frequency)
        L_eflux_hi=-np.log10((1.+z)*flux*frequency)+np.log10((1.+z)*flux*frequency+e_flux*(1.+z)*frequency)

        # plot_sed.pro has the difference between the two luminosities to show the residuals
        # I am just using the predicted luminosity
        #L_pflux=np.log10((1.+z)*flux*3e+14/self.lambda_eff)-np.log10((1.+z)*p_flux*3e+14/self.lambda_eff) # from plot_sed.pro
        L_pflux=np.log10((1.+z)*p_flux*frequency) # log of nu L_nu, predicted
        # show residuals in linear scale, as fraction of observed flux
        resid = (10.**L_flux - 10.**L_pflux)/10.**L_flux
        # approximate error...
        resid_err = e_flux*(1+z)*frequency/10.**L_flux
        ax1.plot(self.lambda_eff,L_flux,'bs',label='VFID '+self.galid)
        yerr = np.array((L_eflux_lo,L_eflux_hi))
        #print('yerr shape = ',yerr.shape)
        #print(yerr)
        ax1.errorbar(self.lambda_eff,L_flux,yerr=yerr,fmt='None',color='b')        
        ax1.plot(self.lambda_eff,L_pflux,'k^',label='Predicted')
        if plotresid:
            ax2.plot(self.lambda_eff,resid,'ko',label='residuals')
            ax2.errorbar(self.lambda_eff,resid,yerr=resid_err,fmt='none',color='k')                
        self.mstar = mstar
        self.sfr = SFR

    def plot_mstar_hist(self,ymin=0,ymax=1,xmin=8,xmax=12):
        in1 = open(self.fit_file,'r')
        fit_lines = in1.readlines()        
        mstar = parse_pdf(fit_lines[209:269])
        in1.close()

        plt.figure(figsize=(8,6))
        plt.subplots_adjust(bottom=.2)
        plt.subplots_adjust(wspace=.1,hspace=.8)
        plt.fill_between(mstar[0],mstar[1])
        plt.axis([xmin,xmax,ymin,ymax])
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)        
        plt.xlabel(r'$\rm log_{10}(M_\star/M_\odot)$',fontsize=18)
        plt.ylabel(r'$\rm PDF$',fontsize=18)        

        outfile = 'VFID{}-mstar-pdfs.png'.format(self.galid)
        plt.savefig(outfile)
    def plot_mstar_sfr_hist(self,ymin=0,ymax=1,xmin=8,xmax=12):
        in1 = open(self.fit_file,'r')
        fit_lines = in1.readlines()        
        mstar = parse_pdf(fit_lines[209:269])
        SFR = parse_pdf(fit_lines[619:679])        
        in1.close()

        plt.figure(figsize=(8,6))
        plt.subplots_adjust(wspace=.01,hspace=.8,bottom=.2)
        plt.subplot(1,2,1)
        plt.fill_between(mstar[0],mstar[1])
        plt.axis([xmin,xmax,ymin,ymax])
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)        
        plt.xlabel(r'$\rm log_{10}(M_\star/M_\odot)$',fontsize=18)
        plt.ylabel(r'$\rm PDF$',fontsize=18)        
        
        plt.subplot(1,2,2)
        plt.fill_between(SFR[0],SFR[1])
        #plt.axis([xmin,xmax,ymin,ymax])
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(-.65,.05)
        plt.ylim(ymin,ymax)
        plt.yticks([],[])
        plt.xlabel(r'$\rm log_{10}(SFR/(M_\odot/yr))$',fontsize=18)
        #plt.ylabel(r'$\rm PDF$',fontsize=18)        

        outfile = 'VFID{}-mstar-sfr-pdfs.png'.format(self.galid)
        plt.savefig(outfile)
        
    def plot_histograms(self):
        ''' plot liklihood histograms of fitted parameters   '''

        # following plot_sed.pro
        # make a 2x6 grid of histograms
        in1 = open(self.fit_file,'r')
        fit_lines = in1.readlines()
        fmu_SFR = parse_pdf(fit_lines[16:36])
        fmu_IR = parse_pdf(fit_lines[39:59])
        mu = parse_pdf(fit_lines[62:82])
        tau_V = parse_pdf(fit_lines[85:133])
        sSFR = parse_pdf(fit_lines[136:206])
        Mstar = parse_pdf(fit_lines[209:269])
        Ld_tot = parse_pdf(fit_lines[272:332])
        Tc_ISM = parse_pdf(fit_lines[335:345])
        Tw_BC = parse_pdf(fit_lines[348:378])
        xi_C_tot = parse_pdf(fit_lines[381:401])
        xi_W_tot = parse_pdf(fit_lines[450:470])
        tau_V_ISM = parse_pdf(fit_lines[473:553])
        Mdust = parse_pdf(fit_lines[556:616])                
        SFR = parse_pdf(fit_lines[619:679])
        in1.close()

        allhist = [fmu_SFR,fmu_IR,mu,tau_V,sSFR,Mstar,\
                   Ld_tot,Tc_ISM,Tw_BC,xi_C_tot,xi_W_tot,\
                   tau_V_ISM,Mdust,SFR]
        allhist_names = ['fmu_SFR','fmu_IR','mu','tau_V','sSFR','Mstar',\
                   'Ld_tot','Tc_ISM','Tw_BC','xi_C_tot','xi_W_tot',\
                   'tau_V_ISM','Mdust','SFR']
        plt.figure(figsize=(8,6))
        plt.subplots_adjust(wspace=.1,hspace=.8)
        for i,h in enumerate(allhist):
            plt.subplot(3,5,i+1)
            t = plt.fill_between(h[0],h[1])
            plt.xlabel(allhist_names[i])
            plt.ylim(0,1)
            if (i == 0) | (i == 5) | (i == 10) :
                plt.yticks(np.linspace(0,1,3))
            else:
                plt.yticks([])

    
        outfile = 'VFID{}-magphys-pdfs.png'.format(self.galid)
        plt.savefig(outfile)
