#!/usr/bin/env python
'''
GOAL:
- make plots!




'''

from matplotlib import pyplot as plt
import numpy as np
import os
homedir = os.getenv("HOME")
os.sys.path.append(homedir+'/github/Virgo/programs/')
from readtables import *

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
        
        
class coplots(readCOtables):
    def plotpositions(self):
        plt.figure(figsize=(12,8))
        plt.scatter(self.main['RA'],self.main['DEC'],c = self.main['vr'],s=10)
        plt.colorbar()

if __name__ == '__main__':
    v = vfplots()
    #v.plotpositions()
