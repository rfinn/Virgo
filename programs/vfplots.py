#!/usr/bin/env python
'''
GOAL:
- make plots!




'''

from matplotlib import pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import os
homedir = os.getenv("HOME")
os.sys.path.append(homedir+'/github/Virgo/programs/')
from readtables import *
plotdir = homedir+'/research/Virgo/plots/'

### FROM NED
virgo_ra,virgo_dec = 187.697083, 12.336944
virgo_vr = 1079.25

h=.74

# using colors from matplotlib default color cycle
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

class vfplots(readfulltables):
    def plotpositions(self):
        '''plot DEC vs RA, color coded by recession velocity.'''
        plt.figure(figsize=(12,8))
        plt.scatter(self.main['RA'],self.main['DEC'],c = self.main['vr'],s=10)
        plt.colorbar()
    def plothasize(self):
        '''plot ha R90/ rband R90 vs R-band C30'''
        plt.figure()
        plt.plot(self.ha['C30'],self.ha['HR_F75']/self.ha['R_F75'],'b.')
    def compare_nsa(self):
        '''compare r-band sersic fit from galfit with nsa'''
        plt.figure(figsize=(10,6))
        flag = self.ha['HAobsflag']
        yvars = ['GAL_RE','GAL_N','GAL_BA','GAL_PA']
        xvars = ['SERSIC_TH50','SERSIC_N','SERSIC_BA','SERSIC_PHI']
        plt.subplots_adjust(hspace=.3,wspace=.3)
        for i,x in enumerate(xvars):
            plt.subplot(2,2,i+1)
            if i == 0:
                plt.plot(self.nsa0[x][flag],self.ha[yvars[i]][flag]*.43,'b.')
            else:
                plt.plot(self.nsa0[x][flag],self.ha[yvars[i]][flag],'b.')
            x1,x2 = plt.xlim()
            xl = np.linspace(x1,x2,100)
            plt.plot(xl,xl,'k--')
            plt.xlabel('NSA '+xvars[i])
            plt.ylabel('GALFIT '+yvars[i])
            if i == 0:
                plt.ylim(plt.xlim())
    def HIdef_env(self):
        '''plot HI def for filament, field, cluster'''
        rmax = 3.6
        #virgo_flag = (np.abs(self.env['distSGY_Virgo']) < rmax) & (np.abs(self.env['distSGZ_Virgo']) < rmax) & (np.abs(self.env['distSGX_Virgo']) < rmax)
        virgo_flag = (self.env['distSGY_Virgo']**2 + self.env['distSGZ_Virgo']**2+ self.env['distSGX_Virgo']**2) < rmax**2
        print('number of cluster members = ',sum(virgo_flag))
        plt.figure(figsize=(12,8))
        plt.subplots_adjust(left=.05,right=.95,hspace=.35)
        envflag = [(self.filmemb['filament']=='---'),(self.filmemb['filament']!='---')]
        filaments = set(self.filmemb['filament'])
        filaments.remove('Virgo_Northern_Filament')
        filaments = list(filaments)
        allax=[]
        for i in range(len(filaments)+1):
            plt.subplot(3,4,i+1)            
            plt.plot(self.main['RA'],self.main['DEC'],'k.',alpha=.01)
            if i < len(filaments):
                flag = (self.filmemb['filament'] == filaments[i]) & (self.a100['HIdef_flag']) & ~virgo_flag
            else:
                flag =  (self.a100['HIdef_flag']) & virgo_flag

            plt.scatter(self.main['RA'][flag],self.main['DEC'][flag],c=self.a100['HIdef_bos'][flag],s=3,vmin=-.5,vmax=.5)
            if i < len(filaments):
                if filaments[i] == '---':
                    mytitle = 'Field'
                    plt.title(mytitle)
                else:
                    mytitle=filaments[i]
                    plt.title(mytitle)
            else:
                mytitle = 'Virgo'
                plt.title(mytitle)
            plt.axis([100,280,-5,40])
            ax = plt.gca()
            allax.append(ax)
            ax.invert_xaxis()

            # print fraction of HI def galaxies
            a = np.sum(self.a100['HIdef_bos'][flag] > 0.3)
            b = np.sum(flag)
            s = mytitle+': %.2f (%i/%i)'%(a/b,a,b)
            print(s)
            plt.title(s,fontsize=10)
            #if i%4 != 0:
            #    plt.yticks(())
            #if i < 8:
            #    plt.xticks(())
        
        cb=plt.colorbar(ax=allax,fraction=.08)
        cb.set_label('HI def')
        plt.savefig(plotdir+'HIdef-filaments.pdf')
        plt.savefig(plotdir+'HIdef-filaments.png')
    def virgo_phasespace(self):
        '''show phasespace diagram relative to center of virgo'''
        dr = np.sqrt((virgo_ra - self.main['RA'])**2 \
                     + (virgo_dec - self.main['DEC'])**2)

        dv = self.main['vr'] - virgo_vr

        plt.figure()
        plt.plot(dr,dv,'k.',alpha=.1)
        rmax = 3.6*h
        #virgo_flag = (np.abs(self.env['distSGY_Virgo']) < rmax) & (np.abs(self.env['distSGZ_Virgo']) < rmax) & (np.abs(self.env['distSGX_Virgo']) < rmax)
        virgo_flag = (self.env['distSGY_Virgo']**2 + self.env['distSGZ_Virgo']**2+ self.env['distSGX_Virgo']**2) < rmax**2
        plt.axhline(y=0)        
        plt.plot(dr[virgo_flag],dv[virgo_flag],'bo',alpha=.3)


    def HIdef_env(self):
        # histogram of Boselli HIdef
        rmax = 3.6
        virgo_flag = (self.env['distSGY_Virgo']**2 + self.env['distSGZ_Virgo']**2+ self.env['distSGX_Virgo']**2) < rmax**2
        print('number of cluster members = ',sum(virgo_flag))
        filament_flag = (self.filmemb['filament']!='---') & (~virgo_flag)
        spiral_flag = (self.hl['t'] > 0)

        plt.figure(figsize=(10,6))
        nbin=np.linspace(-2,2,20)
        #t = plt.hist(self.a100['HIdef_bos'][self.a100['HIdef_flag']],bins=nbin,histtype='step',label='all',normed=True,lw=2)

        # virgo
        flag = self.a100['HIdef_flag'] & (virgo_flag) & spiral_flag
        t = plt.hist(self.a100['HIdef_bos'][flag],bins=nbin,histtype='step',label='Virgo',normed=True,lw=2)
        plt.axvline(x=np.median(self.a100['HIdef_bos'][flag]),c=mycolors[0],ls='--')        
        flag = self.a100['HIdef_flag'] & (filament_flag) & spiral_flag
        t = plt.hist(self.a100['HIdef_bos'][flag],bins=nbin,histtype='step',label='Filaments',normed=True,lw=2)
        plt.axvline(x=np.median(self.a100['HIdef_bos'][flag]),c=mycolors[1],ls='--')                

        flag = self.a100['HIdef_flag'] & (~virgo_flag & ~filament_flag) & spiral_flag
        t = plt.hist(self.a100['HIdef_bos'][flag],bins=nbin,histtype='step', label='Field',normed=True,lw=2)
        plt.axvline(x=np.median(self.a100['HIdef_bos'][flag]),c=mycolors[2],ls='--')        

        plt.xlabel("HI Deficiency (Boselli & Gavazzi)",fontsize=14)
        plt.ylabel('Normalized Histogram',fontsize=14)
        plt.axvline(x=0.3,ls='--',c='k',label='Deficient')
        plt.axvline(x=0.,ls='-',c='k')
        plt.legend()
                

    def compare_velocities(self):
        '''compare recession velocity with flow-corrected velocities'''
        plt.figure()
        plt.plot(self.main['vr'], self.env['Vcosmic'],'k.',alpha=.1,label='Vcosmic')
        
        plt.plot(self.main['vr'][self.main['A100flag']], self.a100['Dist'][self.main['A100flag']]*74,'b.',alpha=.1,label='ALFALFA')        
        xl = np.linspace(500,3500)
        plt.plot(xl,xl,'r--',label='1:1')
        plt.xlabel('Our vr (km/s)',fontsize=14)
        plt.ylabel('Flow-corrected vr (km/s)',fontsize=14)        
        plt.ylim(-500,4000)
        leg = plt.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(.8)
    def compare_flow_velocities(self):
        '''compare recession velocity with flow-corrected velocities'''
        plt.figure()
        plt.plot(self.env['Vcosmic'][self.main['A100flag']],self.a100['Dist'][self.main['A100flag']]*74,'k.',alpha=.1,label='Vcosmic')
        
        #plt.plot(self.main['vr'][self.main['A100flag']], ,'b.',alpha=.1,label='ALFALFA')        
        xl = np.linspace(500,3500)
        plt.plot(xl,xl,'r--',label='1:1')
        plt.xlabel('Vcosmic (km/s)',fontsize=14)
        plt.ylabel('ALFALFA flow-corrected vr (km/s)',fontsize=14)        
        plt.ylim(-500,4000)
        leg = plt.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(.8)

    def compare_a100_distance(self):
        '''plot a100 distance versus distance from Vcosmic'''
        plt.figure()
        flag = self.main['A100flag']
        plt.plot(self.a100['Dist'][flag], self.env['Vcosmic'][flag]/70,'k.',alpha=.1,label='full sample')
        x1,x2 = plt.xlim()
        x1,x2 = 0,60
        xl = np.linspace(x1,x2,100)
        plt.plot(xl,xl,'r--',label='1:1')
        plt.xlabel('ALFALFA Distance (Mpc)',fontsize=14)
        plt.ylabel('Vcosmic Distance (Mpc)',fontsize=14)        
        plt.ylim(x1,x2)
        plt.xlim(x1,x2)        
        plt.legend()

    def sitelle_sample(self):
        '''
        select galaxies with NUV and W3 emission; 
        identify galaxies with multiple neighbors within r = 10arcmin
        focus on Virgo 3 and NGC filament
        '''

        flag = (self.nsa['SERSIC_ABSMAG'][:,1] < 0) &\
            (self.unwise['W3SNR'] > 5)

    def HIdef_morphology(self):
        flag = self.main['A100flag'] & (self.mstar['logMstar'] > -99)
        y = self.a100['HIdef_bos']
        x = self.hl['t']
        mycolors = [self.env['n5th_2D'],self.mstar['logMstar']]
        ylabels = ['5th_2D','logMstar']
        v1 = [.01,8.]
        v2 = [40,10.5]
        plt.figure(figsize=(12,4))
        plt.subplots_adjust(wspace=0.3)
        for i,c in enumerate(mycolors):
            plt.subplot(1,2,i+1)
            plt.scatter(x[flag],y[flag],c=mycolors[i][flag],alpha=.7,vmin=v1[i],vmax=v2[i],s=10)
            cb = plt.colorbar()
            cb.set_label(ylabels[i])
            plt.ylim(-1.5,2)            
            plt.ylabel('HI Deficiency')
            plt.xlabel('T type')
        plt.savefig(plotdir+'/HIdef_morphology_fullsample.png')
        pass
        
class coplots(readCOtables):
    def plotpositions(self):
        plt.figure(figsize=(12,8))
        plt.scatter(self.main['RA'],self.main['DEC'],c = self.main['vr'],s=10)
        plt.colorbar()

if __name__ == '__main__':
    v = vfplots()
    #v.plotpositions()
