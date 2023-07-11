#!/usr/bin/env python

'''
GOAL:
* combine tables that are created from the halpha gui
* astropy.table vstack is really slow, so I am going to try to do it more manually

DATA:
* I moved all the most recent files from HDI, INT feb 2019 and INT jun 2019 to ~/research/Virgo/halpha-tables-20210311


PROCEDURE:
* count total lines in input files
* create a new empty array with the right number of rows and dataype of the individual tables
* write individual tables into the combined table


OUTPUT:
* combined table
'''

from datetime import date
import glob
import numpy as np

from astropy.io import fits
from astropy.table import Table

# get input files
flist1 = glob.glob("v*.fits")
flist2 = glob.glob("VF*rfinn*.fits")
flist3 = glob.glob("nNGC*.fits")
flist = flist1 + flist2 + flist3
flist = flist2 = glob.glob("VF*data-rfinn*.fits")
flist.sort()
print(f"Found {len(flist)} files to stack")
print()
# get total number of lines
print('getting total number of lines')
nlines=0
for f in flist:
    t = fits.getdata(f)
    nlines += len(t)

# create output table

outtab = np.zeros(nlines,dtype=t.dtype)

# loop through input tables and write the rows into the output table
print('filling in output table')
startindex=0
for f in flist:
    t = fits.getdata(f)
    n = len(t)
    outtab[startindex:startindex+n] = t
    startindex += n

# write the output table
print('writing output table')
today = date.today()
str_date_today = today.strftime('%Y-%b-%d')
outtab = Table(outtab)
outtab_name = 'halphagui-output-combined-{}.fits'.format(str_date_today)

outtab.write(outtab_name,overwrite=True)
