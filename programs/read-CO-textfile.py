#!/usr/bin/env python

'''
PURPOSE
- read in table from Francoise and Pascale
  - latest is ~/github/tables/CO-masterfile-2017May09.dat
- match with observation data
  - latest is in ~/proposals/NSF2016/all-co10.tab
  - latest is in ~/proposals/NSF2016/all-co21.tab
- write out a fits table with combined data from both tables
  - output is ~/github/tables/CO-MasterFile-2017MonDD.fits



USEAGE:
- this is set up to run from ~/github/Virgo directory

in ipython -pylab

%run read-CO-textfie.py



'''

import numpy as np
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import time

gitpath = '/Users/rfinn/git/Virgo/'
# I removed some extra lines from Pascale's file
infile = 'tables/Jablonka-MasterFile-Finn.list'

infile = gitpath+'tables/Master_file_fm.dat'
infile = gitpath+'tables/Combes-2017May09-masterfile.dat'
infile = gitpath+'tables/CO-masterfile-2017May09.dat'

'''
UGCA201 - changed 1? to 2 for CO 
NGC3595 - changed HI detect from xxxx to xxx
'''


filaments = ['Filament1','Filament1-complement','N5353group','Filament3','LeoII-B','LeoII-A','LeoMinorFilament4','VirgoIIIcenter','VirgoIIIWest']

inf1 = open(infile,'r')
NEDname=[]
SDSS_ID = []
RA = []
DEC = []
vr = []
CO = []
CO_DETECT = []
HI =[]
filament= []
HImass = []
H2mass = []
upperlimit=[]
i=0
for line in inf1:
    i+=1
    if i < 14: # skipp first 15 header lines - this may change with future versions - check!!
        continue
    if len(line) < 1: # skip blank lines
        continue
    if line.startswith('#'):
        t = line.split('#')
        j,current_filament  = t
        current_filament = "".join(current_filament.replace("&","-").split())
        i += 1
        continue
    t=line.split()
    if len(t) == 6:
        NEDname.append(t[0])
        SDSS_ID.append(t[1])
        RA.append(t[2])
        DEC.append(t[3])
        CO.append(t[4])
        CO_DETECT.append(t[5])
        filament.append(current_filament)
        HI.append('-1')
    elif len(t) == 7:
        NEDname.append(t[0])
        SDSS_ID.append(t[1])
        RA.append(t[2])
        DEC.append(t[3])
        CO.append(t[8])
        CO_DETECT.append(t[9])
        HI.append(t[6])
        filament.append(current_filament)
    elif len(t) == 14:
        NEDname.append(t[0])
        SDSS_ID.append(t[1])
        RA.append(t[2]+':'+t[3]+':'+t[4])
        DEC.append(t[5]+':'+t[6]+':'+t[7])
        vr.append(t[8])
        CO.append(t[9])
        CO_DETECT.append(t[10].replace('xxx','-1'))
        HI.append(t[11].replace('xxx','-1').replace('N','0'))
        HImass.append(t[12].replace('xxx','-1'))
        h2 = t[13].replace('xxx','-1')
        if h2.find('<') > -1:
            upperlimit.append('1')
            H2mass.append(h2.replace('<',''))  # needed to edit catalog- some have 4 x's
        else:
            upperlimit.append('0')
            H2mass.append(h2)
        filament.append(current_filament)

coords = SkyCoord(RA,DEC, unit= (u.hour,u.deg))

# read in observation files
co10 = np.recfromcsv('/Users/rfinn/proposals/NSF2016/all-co10.tab',delimiter='')

matchflag_co10 = np.zeros(len(RA),'bool')
matchindex_co10 = np.zeros(len(RA),'i')

for j in range(len(co10['source'])):    
    for i in range(len(NEDname)):
        if co10['source'][j].lower().find(NEDname[i].lower()) > -1:
            matchflag_co10[i] = 1
            matchindex_co10[i] = j
            break
    if not(matchflag_co10[i]):
        print('no match for ',co10['source'][j])

area_co10 = np.zeros(len(RA),'f')
area_err_co10 = np.zeros(len(RA),'f')
V_co10 = np.zeros(len(RA),'f')
V_err_co10 = np.zeros(len(RA),'f')
DV_co10 = np.zeros(len(RA),'f')
DV_err_co10 = np.zeros(len(RA),'f')
Tpeak_co10 = np.zeros(len(RA),'f')
Tpeak_err_co10 = np.zeros(len(RA),'f')
rms_co10 = np.zeros(len(RA),'f')

area_co10[matchflag_co10] = co10['area'][matchindex_co10[matchflag_co10]] 
area_err_co10[matchflag_co10] = co10['err'][matchindex_co10[matchflag_co10]] 
V_co10[matchflag_co10] = co10['v'][matchindex_co10[matchflag_co10]] 
V_err_co10[matchflag_co10] = co10['err_1'][matchindex_co10[matchflag_co10]] 
DV_co10[matchflag_co10] = co10['dv'][matchindex_co10[matchflag_co10]] 
DV_err_co10[matchflag_co10] = co10['err_2'][matchindex_co10[matchflag_co10]] 
Tpeak_co10[matchflag_co10] = co10['tpeak'][matchindex_co10[matchflag_co10]] 
Tpeak_err_co10[matchflag_co10] = co10['err_3'][matchindex_co10[matchflag_co10]] 
rms_co10[matchflag_co10] = co10['rms'][matchindex_co10[matchflag_co10]] 

col12 = fits.Column(name='CO10_area',format='E',array=area_co10)
col13 = fits.Column(name='CO10_area_err',format='E',array=area_err_co10)
col14 = fits.Column(name='CO10_V',format='E',array=V_co10)
col15 = fits.Column(name='CO10_V_err',format='E',array=V_err_co10)
col16 = fits.Column(name='CO10_DV',format='E',array=DV_co10)
col17 = fits.Column(name='CO10_DV_err',format='E',array=DV_err_co10)
col18 = fits.Column(name='CO10_Tpeak',format='E',array=Tpeak_co10)
col19 = fits.Column(name='CO10_Tpeak_err',format='E',array=Tpeak_err_co10)
col20 = fits.Column(name='CO10_rms',format='E',array=rms_co10)

# repeat for CO 2-1 file
co21 = np.recfromcsv('/Users/rfinn/proposals/NSF2016/all-co21.tab',delimiter='')

matchflag_co21 = np.zeros(len(RA),'bool')
matchindex_co21 = np.zeros(len(RA),'i')

for j in range(len(co21['source'])):    
    for i in range(len(NEDname)):
        if co21['source'][j].lower().find(NEDname[i].lower()) > -1:
            matchflag_co21[i] = 1
            matchindex_co21[i] = j
            break
    if not(matchflag_co21[i]):
        print('no match for ',co21['source'][j])

area_co21 = np.zeros(len(RA),'f')
area_err_co21 = np.zeros(len(RA),'f')
V_co21 = np.zeros(len(RA),'f')
V_err_co21 = np.zeros(len(RA),'f')
DV_co21 = np.zeros(len(RA),'f')
DV_err_co21 = np.zeros(len(RA),'f')
Tpeak_co21 = np.zeros(len(RA),'f')
Tpeak_err_co21 = np.zeros(len(RA),'f')
rms_co21 = np.zeros(len(RA),'f')

area_co21[matchflag_co21] = co21['area'][matchindex_co21[matchflag_co21]] 
area_err_co21[matchflag_co21] = co21['err'][matchindex_co21[matchflag_co21]] 
V_co21[matchflag_co21] = co21['v'][matchindex_co21[matchflag_co21]] 
V_err_co21[matchflag_co21] = co21['err_1'][matchindex_co21[matchflag_co21]] 
DV_co21[matchflag_co21] = co21['dv'][matchindex_co21[matchflag_co21]] 
DV_err_co21[matchflag_co21] = co21['err_2'][matchindex_co21[matchflag_co21]] 
Tpeak_co21[matchflag_co21] = co21['tpeak'][matchindex_co21[matchflag_co21]] 
Tpeak_err_co21[matchflag_co21] = co21['err_3'][matchindex_co21[matchflag_co21]] 
rms_co21[matchflag_co21] = co21['rms'][matchindex_co21[matchflag_co21]] 

col21 = fits.Column(name='CO21_area',format='E',array=area_co21)
col22 = fits.Column(name='CO21_area_err',format='E',array=area_err_co21)
col23 = fits.Column(name='CO21_V',format='E',array=V_co21)
col24 = fits.Column(name='CO21_V_err',format='E',array=V_err_co21)
col25 = fits.Column(name='CO21_DV',format='E',array=DV_co21)
col26 = fits.Column(name='CO21_DV_err',format='E',array=DV_err_co21)
col27 = fits.Column(name='CO21_Tpeak',format='E',array=Tpeak_co21)
col28 = fits.Column(name='CO21_Tpeak_err',format='E',array=Tpeak_err_co21)
col29 = fits.Column(name='CO21_rms',format='E',array=rms_co21)



# write out a fits file


col1 = fits.Column(name='NED_name',format='20A',array=np.array(NEDname))
col2 = fits.Column(name='SDSS_ID',format='20A',array=np.array(SDSS_ID))
col3 = fits.Column(name='RA',format='E',array=coords.ra.deg)
col4 = fits.Column(name='DEC',format='E',array=np.array(coords.dec.deg))
col5 = fits.Column(name='CO',format='2A',array=np.array(CO))
col6 = fits.Column(name='CO_DETECT',format='J',array=np.array(CO_DETECT))
col7 = fits.Column(name='HI',format='J',array=np.array(HI))
col8 = fits.Column(name='HImass',format='E',array=np.array(HImass))
col9 = fits.Column(name='H2mass',format='F',array=np.array(H2mass))
col10 = fits.Column(name='H2mass_upperlimit',format='I',array=np.array(upperlimit))
col11 = fits.Column(name='filament',format='20A',array=np.array(filament))

cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29])

tbhdu = fits.BinTableHDU.from_columns(cols)

date_stamp = time.strftime("%Y%b%d")
tbhdu.writeto('tables/CO-MasterFile-'+date_stamp+'.fits',overwrite=True)
#tbhdu.writeto('tables/Jablonka-MasterFile-'+date_stamp+'.fits',clobber=True)
#tbhdu.writeto('tables/Jablonka-MasterFile.fits',clobber=True)

