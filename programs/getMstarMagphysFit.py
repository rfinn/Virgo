#!/usr/bin/env python

import glob

from astropy.table import Table
import numpy as np

flist = glob.glob('*.fit')

vfid = []
mstar = []
sfr = []
for f in flist:
    t = f.split('.')
    vfid.append(t[0])
    in1 = open(f,'r')
    fit_lines = in1.readlines()        
    bestfit_model = fit_lines[10]
    print(bestfit_model)
    mstar.append(np.log10(bestfit_model[5]))
    sfr.append(bestfit_model[2])    
    in1.close()
