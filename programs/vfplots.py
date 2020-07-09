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
        rmax = 3.6
        #virgo_flag = (np.abs(self.env['distSGY_Virgo']) < rmax) & (np.abs(self.env['distSGZ_Virgo']) < rmax) & (np.abs(self.env['distSGX_Virgo']) < rmax)
        virgo_flag = (self.env['distSGY_Virgo']**2 + self.env['distSGZ_Virgo']**2+ self.env['distSGX_Virgo']**2) < rmax**2
        plt.axhline(y=0)        
        plt.plot(dr[virgo_flag],dv[virgo_flag],'bo',alpha=.3)




        
class coplots(readCOtables):
    def plotpositions(self):
        plt.figure(figsize=(12,8))
        plt.scatter(self.main['RA'],self.main['DEC'],c = self.main['vr'],s=10)
        plt.colorbar()

if __name__ == '__main__':
    v = vfplots()
    #v.plotpositions()
