#!/usr/bin/env python


'''

Written by Rose Finn

PURPOSE:
- to help plan Halpha observations for Virgo Filament galaxies
- to make airmass plots
- to make finding charts

USEAGE:

source activate py36

start ipython in github/Virgo/

%run programs/kpno-halpha.p3.py


To make airmass plots:

   airmass_plots(kittpeak=False,ING=True)  # INT Feb 2019:

   or

   airmass_plots(kittpeak=False,MLO=True)  # MLO April 2019


To make finding charts:

   platinum_finding_chart(1,MLO=True,ING=False) # makes 1

   make_all_platinum(MLO=True, KPNO=False,ING=False)

   make_all_platinum(KPNO=True,ING=False,MLO=False,startnumber=None) # makes all, starting at number XX
   
To plot NE filament
make_plot()

To plot one galaxy:
zoomin(1)

To plot all galaxies that we need to observe:
plotall(delta=1)

USEFUL SITES FOR SKY CHART

http://astroweb.case.edu/jakub/TA/Query_databases.html

https://astroquery.readthedocs.io/en/latest/skyview/skyview.html

AIRMASS PLOTS

http://www.astropy.org/astropy-tutorials/Coordinates.html


910,910 is where object is placed in CCD4

'''

########################################
###### IMPORT MODULES  ########
########################################

from astropy.io import fits
from virgoCommon import *
import pylab as plt
import numpy as np
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
from astroquery.skyview import SkyView

from astropy.coordinates import EarthLocation
from astropy.time import Time

from astropy.coordinates import AltAz

from astroplan import Observer
from astroplan.plots import plot_airmass

########################################
###### RUN-SPECIFIC PARAMETERS  ########
########################################
telescope_run = 'MLO-2019Apr-'

# output file prefix
outfile_prefix = 'observing/'+telescope_run
max_pointing = 57

#2019
outfile_prefix = 'observing/2019Feb-INT-'
max_pointing = None

# Mt Laguna Instrument
outfile_prefix = '/Users/rfinn/Dropbox/Research/Virgo/finding-charts/'+telescope_run


########################################
###### OTHER PARAMETERS  ########
########################################
# mass cuts for plotting NSA galaxies
minmass = 8.0
maxmass = 11.2

###############################################
##### Set moretargets to be true to
##### look for lower mass sources at early RA
###############################################
moretargets = False

########################################
######  READ IN DATA TABLES  #######
########################################
tablepath = gitpath+'Virgo/tables/'
cofile = 'nsa_CO-Gianluca.virgo.fits'
co = fits.getdata(tablepath+cofile)

nsa = fits.getdata(nsa_file)
jmass = fits.getdata(mass_file)
wise = fits.getdata(wise_file)
if moretargets:
    halpha = fits.getdata(halpha_file)
else:
    halpha = fits.getdata('/Users/rfinn/github/Virgo/tables/nsa_Halpha.virgo.2019Feb04.fits')
CJcat = fits.getdata(tablepath+'All-virgo-20feb18_env_H070-FITS.fits')
# find CO targets that are not in NSA?

CJngcflag = CJcat.filament_name == b'1' # == 'Filament1-Group5354'


########################################
######  DEFINE SAMPLES  #######
########################################
#ngcflag =  co.filament == 'Filament1-Group5354'
ngcflag =  co.filament_name == b'1'
leoflag = co.filament_name == b'L2B'
# CO sample includes all galaxies that have been observed in CO
# not all are detection
COsample = co.CO != b''
COflag = (co.COdetected == b'1')
noCOflag = (co.COdetected == b'0')

# the following are galaxies that we still need to observe
# detected in CO but not yet observed in Halpha
ha_obs = (halpha.date_obs != b'')
# galaxies that we have observed AND detected in Halpha
ha_detect = (halpha.date_obs != b'') & (halpha.halpha == 1)

# galaxies that have NOT been observed in Halpha, but have been observed in CO
need_obs = (halpha.date_obs == b'') & (co.CO != b'') #((co.COdetected == '1') | (co.COdetected == '0'))

HIflag = (co.HI == b'1')

########################################
######  DEFINE FILAMENTS  #######
########################################
# NGC5353/4 Filament
radec = (nsa.RA > 192.) & (nsa.RA < 209) & (nsa.DEC > 0.) & (nsa.DEC < 50.) 
radec_flag = radec & (nsa.DEC >(2*(nsa.RA - 205.) + 20) ) & (nsa.DEC < (2*(nsa.RA - 205.) + 55))
filament = (co.filament_name !='') & (nsa.Z*3.e5 >2000.) & (nsa.Z*3.e5 < 3238.)
nsa_flag = (nsa.Z*3.e5 >1234.) & (nsa.Z*3.e5 < 3976.)
mass_flag = (jmass.MSTAR_50 > 8.3) & (jmass.MSTAR_50 < 10.2)

gas_flag = COflag | HIflag
NGCfilament = filament

# Halpha selection
# set RA and DEC as galaxies with
# CO
# no Halpha
# stellar mass between 8.5 < log(M*/Msun) < 10.  according to NSF proposal
obs_mass_flag = COsample & ~ha_obs #& (jmass.MSTAR_50 > 8.5) #& (jmass.MSTAR_50 < 10.) #& (nsa.SERSIC_BA > 0.2)

# resetting to COsample for 2019 observing season
filter_flag = (nsa.Z*3.e5 > 2490.) & (nsa.Z*3.e5 < 6000.)
obs_mass_flag = COsample & ~ha_obs #& filter_flag


# SELECTING LOWER MASS TARGETS FOR INT RUN IN FEB 2019
more_targets_flag = (nsa.Z*3.e5 < 2300.) & (nsa.RA > 115.) & (nsa.RA < 140.) & ~COsample & ~ha_obs & (nsa.DEC > 20.) & (nsa.DEC < 40.) & (jmass.MSTAR_50 > 8.5)

if moretargets == True:
    print('setting sample to more_targets_flag')
    obs_mass_flag = more_targets_flag

########################################
###### COLOR CODE FOR SCATTER PLOTS
########################################
mycolor=jmass.MSTAR_50
v1=minmass
v2=maxmass
mylabel='$ \log_{10}(M_*/M_\odot) $'


mycolor=nsa.Z*3.e5
v1=1000
v2=3000
mylabel='$ v_r \ (km/s) $'


########################################
###### POINTING INFO FROM MAY 2017 #####
########################################


# keep only objects that have been observed in CO but NOT observed in Halpha
pointing_ra = nsa.RA[obs_mass_flag]
pointing_dec = nsa.DEC[obs_mass_flag]
pointing_id = nsa.NSAID[obs_mass_flag]

# sort by RA
sorted_indices = np.argsort(pointing_ra)
pointing_dec = pointing_dec[sorted_indices]
pointing_ra = pointing_ra[sorted_indices]
pointing_id = pointing_id[sorted_indices]
nsadict = dict((a,b) for a,b in zip(pointing_id,np.arange(len(pointing_id))))

########################################
### OFFSETS TO GET MULTIPLE GALAXIES ###
########################################

####  OFFSETS  #####
# add offsets to try to get multiple galaxies in pointings
pointing_offsets_ra = np.zeros(len(pointing_ra))
pointing_offsets_dec = np.zeros(len(pointing_ra))

#### KPNO OFFSETS
'''
try:
    pointing_offsets_dec[nsadict[50209]] = -0.15
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_dec[nsadict[157256]] = 0.2
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_ra[nsadict[85510]] = -0.1
    pointing_offsets_dec[nsadict[85510]] = -0.07
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_ra[nsadict[90957]] = 0.1
except KeyError:
    print('nsa id not found in list of pointings')

try:
    pointing_offsets_dec[nsadict[160613]] = -0.1
except KeyError:
    print('nsa id not found in list of pointings')
    
try:
    pointing_offsets_dec[nsadict[164223]] = 0.15
except KeyError:
    print('nsa id not found in list of pointings')
'''

'''
##################################################
############ INT WFC OFFSETS
##################################################
try:
    pointing_offsets_ra[nsadict[50207]] = -8./60
    pointing_offsets_dec[nsadict[50207]] = 12./60

    pointing_offsets_ra[nsadict[64280]] = -8.5/60
    pointing_offsets_dec[nsadict[64280]] = -3./60

    pointing_offsets_ra[nsadict[64353]] = 14.5/60
    pointing_offsets_dec[nsadict[64353]] = 12./60

    pointing_offsets_ra[nsadict[135051]] = 90./3600
    pointing_offsets_dec[nsadict[135051]] = 90./3600

    pointing_offsets_ra[nsadict[15877]] = 4./60
    pointing_offsets_dec[nsadict[15877]] = -7./60

    pointing_offsets_ra[nsadict[156774]] = -4./60
    pointing_offsets_dec[nsadict[156774]] = -2./60

    pointing_offsets_ra[nsadict[135296]] = -3./60
    pointing_offsets_dec[nsadict[135296]] = 3./60

    pointing_offsets_ra[nsadict[135465]] = 0./60
    pointing_offsets_dec[nsadict[135465]] = 2./60

    pointing_offsets_ra[nsadict[50379]] = 0./60
    pointing_offsets_dec[nsadict[50379]] = 3./60

    pointing_offsets_ra[nsadict[135527]] = 0./60
    pointing_offsets_dec[nsadict[135527]] = -2.5/60

    pointing_offsets_ra[nsadict[135602]] = 0./60
    pointing_offsets_dec[nsadict[135602]] = -2.5/60

    pointing_offsets_ra[nsadict[135606]] = 4./60
    pointing_offsets_dec[nsadict[135606]] = 0./60

    pointing_offsets_ra[nsadict[50569]] = -.5/60
    pointing_offsets_dec[nsadict[50569]] = 3/60

    pointing_offsets_ra[nsadict[135797]] = -6./60
    pointing_offsets_dec[nsadict[135797]] = 4./60

    pointing_offsets_ra[nsadict[47220]] = -7.5/60
    pointing_offsets_dec[nsadict[47220]] = -2/60

    pointing_offsets_ra[nsadict[135852]] = -9./60
    pointing_offsets_dec[nsadict[135852]] = -2/60

    pointing_offsets_ra[nsadict[135862]] = 5/60
    pointing_offsets_dec[nsadict[135862]] = 3/60

    pointing_offsets_ra[nsadict[157256]] = 0./60
    pointing_offsets_dec[nsadict[157256]] = 13/60

    pointing_offsets_ra[nsadict[64909]] = -2/60
    pointing_offsets_dec[nsadict[64909]] = 2.5/60

    pointing_offsets_ra[nsadict[136042]] = 0./60
    pointing_offsets_dec[nsadict[136042]] = 1/60

    pointing_offsets_ra[nsadict[85513]] = 0./60
    pointing_offsets_dec[nsadict[85513]] = 13/60

    pointing_offsets_ra[nsadict[85367]] = 0./60
    pointing_offsets_dec[nsadict[85367]] = 3./60

    pointing_offsets_ra[nsadict[157480]] = 1.5/60
    pointing_offsets_dec[nsadict[157480]] = 0./60

    pointing_offsets_ra[nsadict[157495]] = 3.5/60
    pointing_offsets_dec[nsadict[157495]] = 3./60

    pointing_offsets_ra[nsadict[107148]] = -5/60
    pointing_offsets_dec[nsadict[107148]] = -2.5/60

    pointing_offsets_ra[nsadict[85977]] = 0./60
    pointing_offsets_dec[nsadict[85977]] = -2./60

    pointing_offsets_ra[nsadict[48222]] = -4./60
    pointing_offsets_dec[nsadict[48222]] = 3.3/60

    pointing_offsets_ra[nsadict[137045]] = -5./60
    pointing_offsets_dec[nsadict[137045]] = -3/60

    pointing_offsets_ra[nsadict[107715]] = 0./60
    pointing_offsets_dec[nsadict[107715]] = 2./60

    pointing_offsets_ra[nsadict[137391]] = 10.5/60
    pointing_offsets_dec[nsadict[137391]] = -2./60

    pointing_offsets_ra[nsadict[137460]] = -8/60
    pointing_offsets_dec[nsadict[137460]] = 0./60

    pointing_offsets_ra[nsadict[107764]] = 19./60
    pointing_offsets_dec[nsadict[107764]] = -3./60

    pointing_offsets_ra[nsadict[88142]] = 5./60
    pointing_offsets_dec[nsadict[88142]] = -10./60

    pointing_offsets_ra[nsadict[137993]] = -7./60
    pointing_offsets_dec[nsadict[137993]] = 0./60

    pointing_offsets_ra[nsadict[90176]] = 0./60
    pointing_offsets_dec[nsadict[90176]] = 3./60

    pointing_offsets_ra[nsadict[138221]] = -3.6/60
    pointing_offsets_dec[nsadict[138221]] = 2.2/60

    pointing_offsets_ra[nsadict[87097]] = -3./60
    pointing_offsets_dec[nsadict[87097]] = -3./60

    pointing_offsets_ra[nsadict[87086]] = 3./60
    pointing_offsets_dec[nsadict[87086]] = 0./60

    pointing_offsets_ra[nsadict[138642]] = 4./60
    pointing_offsets_dec[nsadict[138642]] = 2./60

    pointing_offsets_ra[nsadict[159520]] = 0./60
    pointing_offsets_dec[nsadict[159520]] = -4.5/60

    pointing_offsets_ra[nsadict[159779]] = 0.75/60
    pointing_offsets_dec[nsadict[159779]] = -2./60

    pointing_offsets_ra[nsadict[101649]] = 8./60
    pointing_offsets_dec[nsadict[101649]] = 0./60

    pointing_offsets_ra[nsadict[93963]] = 11.5/60
    pointing_offsets_dec[nsadict[93963]] = -2./60

    pointing_offsets_ra[nsadict[90957]] = 8./60
    pointing_offsets_dec[nsadict[90957]] = 9./60

    pointing_offsets_ra[nsadict[92459]] = 5./60
    pointing_offsets_dec[nsadict[92459]] = -3./60

    pointing_offsets_ra[nsadict[140301]] = 0./60
    pointing_offsets_dec[nsadict[140301]] = -3.3/60

    pointing_offsets_ra[nsadict[160627]] = 4./60
    pointing_offsets_dec[nsadict[160627]] = -14/60

    pointing_offsets_ra[nsadict[117685]] = 0./60
    pointing_offsets_dec[nsadict[117685]] = -1.5/60

    pointing_offsets_ra[nsadict[118414]] = -4./60
    pointing_offsets_dec[nsadict[118414]] = 0./60

    pointing_offsets_ra[nsadict[143701]] = 0./60
    pointing_offsets_dec[nsadict[143701]] = 0./60

    pointing_offsets_ra[nsadict[163875]] = 0./60
    pointing_offsets_dec[nsadict[163875]] = 2./60

    pointing_offsets_ra[nsadict[143841]] = -2./60
    pointing_offsets_dec[nsadict[143841]] = 12./60

    pointing_offsets_ra[nsadict[144056]] = 21.5/60
    pointing_offsets_dec[nsadict[144056]] = -2./60

    pointing_offsets_ra[nsadict[67567]] = 10/60
    pointing_offsets_dec[nsadict[67567]] = -2.5/60

    pointing_offsets_ra[nsadict[17878]] = 3.5/60
    pointing_offsets_dec[nsadict[17878]] = 2./60

    pointing_offsets_ra[nsadict[164911]] = -6/60
    pointing_offsets_dec[nsadict[164911]] = 3.5/60

    pointing_offsets_ra[nsadict[165082]] = -1./60
    pointing_offsets_dec[nsadict[165082]] = -1./60

    pointing_offsets_ra[nsadict[18052]] = 3./60
    pointing_offsets_dec[nsadict[18052]] = 3./60

    pointing_offsets_ra[nsadict[165115]] = -6.5/60
    pointing_offsets_dec[nsadict[165115]] = -3./60

    id=145218
    pointing_offsets_ra[nsadict[id]] = 10./60
    pointing_offsets_dec[nsadict[id]] = 10.5/60
except KeyError:
    print('nsa id not found in list of pointings')
try: #143
    id=165200
    pointing_offsets_ra[nsadict[id]] = 9.5/60
    #pointing_offsets_dec[nsadict[id]] = 10.5/60
except KeyError:
    print('nsa id not found in list of pointings')

try: #144
    id=18153
    pointing_offsets_ra[nsadict[id]] = 4.3/60
    pointing_offsets_dec[nsadict[id]] = -12/60
except KeyError:
    print('nsa id not found in list of pointings')
try: #145
    id=145398
    pointing_offsets_ra[nsadict[id]] = 15./60
    pointing_offsets_dec[nsadict[id]] = 2./60
except KeyError:
    print('nsa id not found in list of pointings')
try: #147
    id=18301
    pointing_offsets_ra[nsadict[id]] = 7./60
    pointing_offsets_dec[nsadict[id]] = -10./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #147
    id=145554
    pointing_offsets_ra[nsadict[id]] = -2./60
    pointing_offsets_dec[nsadict[id]] = -10./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #152
    id=18363
    pointing_offsets_ra[nsadict[id]] = 6./60
    pointing_offsets_dec[nsadict[id]] = 14./60
except KeyError:
    print('nsa id not found in list of pointings')
try: #154
    id=145672
    pointing_offsets_ra[nsadict[id]] = 4./60
    pointing_offsets_dec[nsadict[id]] = -10./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #158
    id=165875
    pointing_offsets_ra[nsadict[id]] = 10./60
    pointing_offsets_dec[nsadict[id]] = 2./60
except KeyError:
    print('nsa id not found in list of pointings')
try: #164
    id=145846
    pointing_offsets_ra[nsadict[id]] = 4./60
    pointing_offsets_dec[nsadict[id]] = 11./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #166
    id=165956
    pointing_offsets_ra[nsadict[id]] = -3./60
    pointing_offsets_dec[nsadict[id]] = 0./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #167
    id=145879
    pointing_offsets_ra[nsadict[id]] = -4./60
    pointing_offsets_dec[nsadict[id]] = 12./60
except KeyError:
    print('nsa id not found in list of pointings')
try: #173
    id=166280
    pointing_offsets_ra[nsadict[id]] = -6./60
    pointing_offsets_dec[nsadict[id]] = 0./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #175
    id=121129
    pointing_offsets_ra[nsadict[id]] = 4./60
    pointing_offsets_dec[nsadict[id]] = 9./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #178
    id=166330
    #pointing_offsets_ra[nsadict[id]] = 4./60
    pointing_offsets_dec[nsadict[id]] = -3./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #189
    id=68462
    pointing_offsets_ra[nsadict[id]] = -7./60
    #pointing_offsets_dec[nsadict[id]] = -3./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #192
    id=69842
    pointing_offsets_ra[nsadict[id]] = 16./60
    pointing_offsets_dec[nsadict[id]] = 1./60
except KeyError:
    print('nsa id not found in list of pointings')

try: #195
    id=147731
    pointing_offsets_ra[nsadict[id]] = 0./60
    pointing_offsets_dec[nsadict[id]] = 3./60
except KeyError:
    print('nsa id not found in list of pointings')



##################################################
############ low mass extension of leo filaments
##################################################

try:
    id=156595
    pointing_offsets_ra[nsadict[id]] = 8.5/60
    pointing_offsets_dec[nsadict[id]] = 12./60
except KeyError:
    print('nsa id not found in list of pointings')


try:
    id=64369
    pointing_offsets_ra[nsadict[id]] = 20./60
    pointing_offsets_dec[nsadict[id]] = -10./60
except KeyError:
    print('nsa id not found in list of pointings')

try:
    id=64410
    pointing_offsets_ra[nsadict[id]] = -3./60
    pointing_offsets_dec[nsadict[id]] = 1./60
except KeyError:
    print('nsa id not found in list of pointings')

try:
    id=64403
    pointing_offsets_ra[nsadict[id]] = -7./60
    pointing_offsets_dec[nsadict[id]] = -2./60
except KeyError:
    print('nsa id not found in list of pointings')

try:
    id=50476
    pointing_offsets_ra[nsadict[id]] = 2./60
    pointing_offsets_dec[nsadict[id]] = 4./60
except KeyError:
    print('nsa id not found in list of pointings')

try:
    id=84130
    pointing_offsets_ra[nsadict[id]] = 9./60
    pointing_offsets_dec[nsadict[id]] = -1./60
except KeyError:
    print('nsa id not found in list of pointings')


try:
    id=64305
    pointing_offsets_ra[nsadict[id]] = 3./60
    pointing_offsets_dec[nsadict[id]] = 14./60
except KeyError:
    print('nsa id not found in list of pointings')

try:
    id=84294
    pointing_offsets_ra[nsadict[id]] = 0./60
    pointing_offsets_dec[nsadict[id]] = 14.5/60
except KeyError:
    print('nsa id not found in list of pointings')

##################################################
############ END OF INT WFC OFFSETS
##################################################
'''

# make a dictionary to store the offsets according to NSA ID
# format is offsets = {nsaid:[dra,ddec]}
# offsets = {135046:[5.,4.],
#            84889:[3.,2.],
#            157073:[4.,0],
#            64408:[2.,-1]
#            }

# The following are INT vallues, which I deleted from the long list above
offsets_INT = {135046:[5.,4.],
           84889:[3.,2.],
           157073:[4.,0],
           64408:[2.,-1]
           }


offsets_MLO = {135046:[5.,4.],
           84889:[3.,2.]
           }

# change this to use the offsets for the desired telescope
offsets = offsets_MLO

for key in offsets:
    pointing_offsets_ra[nsadict[key]] = offsets[key][0]/60.
    pointing_offsets_dec[nsadict[key]] = offsets[key][1]/60.

# update pointing centers to account for offsets
pointing_ra += pointing_offsets_ra
pointing_dec += pointing_offsets_dec



def find_CO_noNSA():
    virgocat = coords.SkyCoord(nsa.RA*u.degree,nsa.DEC*u.degree,frame='icrs')
    jcat = coords.SkyCoord(CJcat.RA*u.degree,CJcat.DEC*u.degree,frame='icrs')

    index,dist2d,dist3d = jcat.match_to_catalog_sky(virgocat)

    # only keep matches with matched RA and Dec w/in 1 arcsec
    matchflag = dist2d.degree < 5./3600

    print('number in CO catalog = ',len(CJcat.RA))
    print('number with NSA matches = ',sum(matchflag))
    return matchflag
def add_detections():
    pointing_ra = nsa.RA[COflag]
    pointing_dec = nsa.DEC[COflag]
    for i in range(len(pointing_ra)):
        rect= plt.Rectangle((pointing_ra[i]-.25,pointing_dec[i]-.25), .5, .5,fill=False, color='k')
        plt.gca().add_artist(rect)


def add_CJpoints():
    flag = CJngcflag  &  (CJcat.COdetected == 1) #& CJnoNSA
    plt.plot(CJcat.RA[flag],CJcat.DEC[flag],'gs',mfc='None',mec='g',markersize=20)
    flag = CJngcflag & (CJcat.HI == 1)  #& CJnoNSA  
    plt.plot(CJcat.RA[flag],CJcat.DEC[flag],'bs',mfc='None',mec='b',markersize=25)
    

def add_pointings():
    for i in range(len(pointing_ra)):
        rect= plt.Rectangle((pointing_ra[i]-.25,pointing_dec[i]-.25), .5, .5,fill=False, color='g',lw=2)
        plt.gca().add_artist(rect)
        plt.text(pointing_ra[i]-.26,pointing_dec[i],'  '+str(i+1),fontsize=8,clip_on=True)
 

# match with stellar mass, NSA, WISE

def plot_positions(plotsingle=True, flag = need_obs,plotha=True):
    #if plotsingle:
    #    plt.figure()
    #flag = need_obs
    
    plt.scatter(nsa.RA[flag],nsa.DEC[flag],s=30,c=mycolor[flag],vmin=v1,vmax=v2,cmap='jet',marker='o',label='CO')
    plt.scatter(nsa.RA[noCOflag],nsa.DEC[noCOflag],s=80,c=mycolor[noCOflag],vmin=v1,vmax=v2,cmap='jet',marker='x',label='no CO')
    plt.scatter(nsa.RA[HIflag],nsa.DEC[HIflag],s=100,c=mycolor[HIflag],vmin=v1,vmax=v2,cmap='jet',marker='+',label='HI')
    #plt.scatter(nsa.RA[gas_flag],nsa.DEC[gas_flag],s=50,c=jmass.MSTAR_50[gas_flag],vmin=8,vmax=11)
    if plotha:
        plt.plot(nsa.RA[ha_obs],nsa.DEC[ha_obs],'cs',mfc='None',markersize=8)
    if plotsingle:
        plt.colorbar(label=mylabel,fraction=0.08)
        plt.legend()
        #plt.gca().invert_xaxis()
def add_nsa():
    #plt.scatter(nsa.RA[nsa_flag],nsa.DEC[nsa_flag],s=10,c=jmass.MSTAR_50[nsa_flag],vmin=minmass,vmax=maxmass,alpha=.5,cmap='jet',marker='v',label='NSA')
    plt.scatter(nsa.RA[nsa_flag],nsa.DEC[nsa_flag],s=10,c=nsa.Z[nsa_flag]*3.e5,vmin=1500,vmax=3500,alpha=.5,cmap='jet',marker='v',label='NSA')

def make_plot_2018(plotsingle=True):
    if plotsingle:
        plt.figure(figsize=(9,4))
    add_nsa()
    plot_positions(plotsingle=plotsingle,flag=obs_mass_flag)
    #add_detections()
    add_pointings()
    add_CJpoints()
    #plt.axis([190,212,15,50])
    plt.gca().invert_xaxis()
    if plotsingle:
        plt.xlabel('$RA \ (deg) $')
        plt.ylabel('$DEC \ (deg) $')

def show_INT_feb():
    make_plot_2018()
    plt.axis([115,200,15,52])
    plt.gca().invert_xaxis()
    plt.savefig('plots/INT_feb_locations.png')
    
def make_plot(plotsingle=True):
    if plotsingle:
        plt.figure()
    add_nsa()
    plot_positions(plotsingle=plotsingle)
    #add_detections()
    add_pointings()
    add_CJpoints()
    plt.axis([190,212,15,50])
    plt.gca().invert_xaxis()
    if plotsingle:
        plt.xlabel('$RA \ (deg) $')
        plt.ylabel('$DEC \ (deg) $')

# FOV of HDI is 0.5 x 0.5 deg
# draw squares with size = FOV

def plotzoom1():
    make_plot()
    plt.axis('equal')
    plt.axis([201.55,202.55,46,47.5])
    plt.gca().invert_xaxis()

def zoomin(ngal,delta=1.5,plotsingle=True):
    make_plot(plotsingle=plotsingle)
    plt.axis([pointing_ra[ngal-1]-delta, pointing_ra[ngal-1]+delta, pointing_dec[ngal-1]-delta, pointing_dec[ngal-1]+delta])
    plt.gca().invert_xaxis()
    if plotsingle:
        plt.savefig(outfile_prefix+'%i-zoomin.png'%(ngal))
        
def plotall(delta=1.5,nrow=3,ncol=3):
    sorted_index = np.argsort(pointing_ra)
    count = 0
    nfig = 0
    
    if max_pointing != None:
        pointing_range = list(range(max_pointing))
    else:
        pointing_range = list(range(len(pointing_ra)))
        
    for i in pointing_range:
        if count == 0:
            plt.figure(figsize=(10,8))
            plt.subplots_adjust(wspace=.2,hspace=.4)
            allax = []
        plt.subplot(nrow,ncol,count+1)
        zoomin(sorted_index[i]+1,plotsingle=False,delta=delta)
        plt.title('NSAID '+str(pointing_id[i]))
        allax.append(plt.gca())
        count += 1
        if count == nrow*ncol:
            plt.colorbar(ax = allax,label='$ \log_{10}(M_*/M_\odot) $',fraction=0.08)
            plt.legend()
            plt.savefig(outfile_prefix+'plotall-%i.png'%(nfig))
            nfig += 1
            count = 0

    if count != 0:
        plt.colorbar(ax = allax,label='$ \log_{10}(M_*/M_\odot) $',fraction=0.08)
        plt.legend()
        plt.savefig(outfile_prefix+'plotall-%i.png'%(nfig)) 
def fix_project(ra,dec,b):
    a = np.cos(np.radians(dec))*np.cos(np.radians(ra-b.header['CRVAL1']))
    f = 1./b.header['CDELT2']*(180./np.pi)/(np.sin(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))+a*np.cos(np.radians(b.header['CRVAL2'])))
    decn = -f*(np.cos(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))-a*np.sin(np.radians(b.header['CRVAL2'])))
    ran = -f*np.cos(np.radians(dec))*np.sin(np.radians(ra-b.header['CRVAL1']))

    ran = b.header['CRVAL1']+(ran)*b.header['CDELT1']
    decn = b.header['CRVAL2']-(decn)*b.header['CDELT2']    
    return ran,decn
        
def finding_chart_all():
    if max_pointing != None:
        pointing_range = range(1,max_pointing+1)
    else:
        pointing_range = range(1,len(pointing_ra)+1)
                               
    for i in pointing_range:
        plt.close('all')
        finding_chart(i)

def finding_chart(npointing,delta_image = .25,offset_ra=0.,offset_dec=0.,plotsingle=True,ING=False,MLO=False,KPNO=False):
    galsize=0.033
    i = npointing-1

    center_ra = pointing_ra[i]+offset_ra
    center_dec = pointing_dec[i] + offset_dec
    print('center ra, dec = ',center_ra,center_dec)
    if plotsingle:
        if delta_image > .3:
            plt.subplots_adjust(right=.9,top=.9,bottom=.1,left=.1)
            plt.figure(figsize=(10,9))
        else:
            plt.figure(figsize=(6,6))
    ax=plt.gca()
    pos = coords.SkyCoord(center_ra,center_dec,frame='icrs',unit='degree')
    delta_imagex=2.*delta_image
    delta_imagey=2.*delta_image
    if ING:
        delta_imagex=.8 #image width in deg
        delta_imagey=.8 # image width in deg
        xout = SkyView.get_images(pos,survey=['DSS'],height=delta_imagex*u.degree,width=delta_imagey*u.degree)
    elif MLO:
        delta_imagex = 13./60. # image width in deg
        delta_imagey = 13./60 # image width in deg
    xout = SkyView.get_images(pos,survey=['DSS'],height=2*delta_image*u.degree,width=2.*delta_image*u.degree)
    b=xout[0][0]
    ax.imshow(xout[0][0].data,interpolation='none',aspect='equal',cmap='gray_r',extent=[b.header['CRVAL1']-(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL1']+(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL2']+(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2'],
                                                           b.header['CRVAL2']-(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2']])

    if ING:
        # add footprint of WFC chips
        plot_INT_footprint(center_ra,center_dec)
    elif MLO:
        delta_arcmin = 13.
        rect= plt.Rectangle((center_ra-0.5*delta_imagex,center_dec-0.5*delta_imagey), delta_arcmin/60., delta_arcmin/60.,fill=False, color='k')
        plt.gca().add_artist(rect)
         
    elif KPNO:
        rect= plt.Rectangle((center_ra-delta_image,center_dec-delta_image), 0.5, 0.5,fill=False, color='k')
        plt.gca().add_artist(rect)

    #add_cameras()
    # find galaxies on FOV
    gals = (nsa.RA > (pos.ra.value-delta_imagex/2.)) & (nsa.RA < (pos.ra.value+delta_imagex/2)) & (nsa.DEC > (pos.dec.value-delta_imagey/2.)) & (nsa.DEC < (pos.dec.value+delta_imagey/2.))
    print('pointing ',i+1,' ngal = ',np.sum(gals))
    gindex=np.arange(len(nsa.RA))[gals]
    print('Pointing %02d Galaxies'%(npointing),': ',nsa.NSAID[gals])
    print('Pointing %02d Galaxies'%(npointing),': ',nsa.Z[gals]*3.e5)
    for j in gindex:
        ran,decn=fix_project(nsa.RA[j],nsa.DEC[j],b)
        rect= plt.Rectangle((ran-galsize/2.,decn-galsize/2.), galsize, galsize,fill=False, color='c')
        plt.gca().add_artist(rect)
        s='%i\n vr=%i'%(nsa.NSAID[j],nsa.Z[j]*3.e5)
        plt.text(ran,decn+galsize/2.,s,fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='bottom')
        plt.text(ran,decn-galsize/2.,co.NEDname[j].decode("utf-8"),fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='top')
        if COflag[j]:
            size=galsize-.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='g')
            plt.gca().add_artist(rect)
        if HIflag[j]:
            size=galsize+.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='b')
            plt.gca().add_artist(rect)
        if ha_obs[j]:
            size=galsize+.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='r')
            plt.gca().add_artist(rect)
        #plt.legend(['filament','CO','Halpha'])
    if moretargets:
        plt.title('LM-Pointing %02d: %s'%((i+1),pos.to_string('hmsdms')))
    else:
        plt.title('Pointing %02d: %s'%((i+1),pos.to_string('hmsdms')))
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.gca().invert_yaxis()
    if plotsingle:
        plt.savefig(outfile_prefix+'Pointing%02d.png'%(i+1))

def plot_INT_footprint(center_ra,center_dec):
    #using full detector sizes for now because 
    detector_dra = 4100.*0.33/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    detector_ddec = 2048.*0.33/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    # draw footprint of chip 4
    rect= plt.Rectangle((center_ra-detector_dra/2.,center_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)
    # draw footprint of chip 3
    # assuming chip 3 is NORTH and a smidge WEST of chip 4
    offset_dec = detector_ddec+17./3600. # 17 arcsec gap in RA between 
    offset_ra = -9.5/3600. # 9.5 arcsec offset toward N
    rect= plt.Rectangle((center_ra+offset_ra-detector_dra/2.,center_dec+offset_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)

    # draw footprint of chip 1
    # assuming chip 1 is SOUTH and a smidge EAST of chip 4
    offset_dec = -1*detector_ddec-22.7/3600. # 17 arcsec gap in RA between 
    offset_ra = +3.18/3600. # 9.5 arcsec offset toward N
    rect= plt.Rectangle((center_ra+offset_ra-detector_dra/2.,center_dec+offset_dec-detector_ddec/2.), detector_dra, detector_ddec,fill=False, color='k')
    plt.gca().add_artist(rect)

    # draw footprint of chip 2
    # assuming chip 2 is WEST of chip 4
    offset_dec = detector_ddec/2.-detector_dra-19.2/3600. # hard to explain
    offset_ra =  -.5*detector_dra-23.8/3600.# hard to explain
    # this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
    rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -1.*detector_ddec, detector_dra,fill=False, color='k')
    plt.gca().add_artist(rect)

    # adding guide camera
    offset_dec = -2*detector_ddec-(22.7+98.1)/3600. # hard to explain
    offset_ra =  detector_dra/2-(3.18+649.9)/3600.# hard to explain
    # this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
    rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -7./60., 7./60,fill=False, color='k')
    plt.gca().add_artist(rect)
    
        
def add_cameras():
    # add guidestar cameras for HDI
    # North camera
    delta_ra = 0
    delta_dec = 2610./3600
    width_ra = 3.3/60.
    width_dec = 2.2/60.
    rect= plt.Rectangle((center_ra+delta_ra - 0.5*width_ra,center_dec+delta_dec - 0.5*width_dec), width_ra, width_dec,fill=False, color='k')
    plt.gca().add_artist(rect)
    # South camera
    delta_ra = -25./3600
    delta_dec = -2410./3600
    width_ra = 3.3/60.
    width_dec = 2.2/60.
    rect= plt.Rectangle((center_ra+delta_ra - 0.5*width_ra,center_dec+delta_dec - 0.5*width_dec), width_ra, width_dec,fill=False, color='k')
    plt.gca().add_artist(rect)

    
def finding_chart_with_guide_stars(npointing,offset_ra=0.,offset_dec=0.):
    finding_chart(npointing, delta_image=.75,offset_ra=offset_ra, offset_dec=offset_dec,plotsingle=False)
    plt.savefig(outfile_prefix+'Pointing%02d-guiding.png'%(npointing))

def make_all_platinum(KPNO=True,ING=False,MLO=False,startnumber=None):    
    if max_pointing != None:
        pointing_range = range(max_pointing)
    else:
        pointing_range = range(len(pointing_ra))
    for i in pointing_range:
        plt.close('all')
        platinum_finding_chart(i+1,KPNO=KPNO,ING=ING,MLO=MLO)

def platinum_finding_chart(npointing,offset_ra=0.,offset_dec=0.,ING=False,KPNO=False,MLO=False):
    if KPNO:
        fig = plt.figure(figsize = (10,6))
        grid = plt.GridSpec(2,3,hspace=.4,wspace=.2,left=.05)
        hdi = fig.add_subplot(grid[:,:-1])
        finding_chart(npointing,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False,ING=ING,MLO=MLO,KPNO=KPNO)
        south_camera = fig.add_subplot(grid[1,2])
        show_guide_camera(npointing,south_camera=True,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
        north_camera = fig.add_subplot(grid[0,2])
        show_guide_camera(npointing,south_camera=False,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
    else:
        fig = plt.figure(figsize = (8,8.))
        finding_chart(npointing,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False,ING=ING,MLO=MLO,KPNO=KPNO)
    if moretargets:
        plt.savefig(outfile_prefix+'Pointing%03d-lowMass-platinum.png'%(npointing))
    else:
        plt.savefig(outfile_prefix+'Pointing%03d-platinum.png'%(npointing))

def guide_cameras(npointing,offset_ra=0,offset_dec=0):
    i = npointing - 1

    plt.figure(figsize=(10,5))
    plt.subplots_adjust(wspace=.3)
    plt.subplot(1,2,1)
    show_guide_camera(npointing,south_camera=False,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
    plt.subplot(1,2,2)
    show_guide_camera(npointing,south_camera=True,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)

    mytitle = "Guide Cameras: Pointing %02d (NSAID %i)"%((i+1),pointing_id[i])
    ax = plt.gca()
    plt.text(-.1,1.15,mytitle,transform=ax.transAxes,horizontalalignment='center',fontsize=18)

def show_guide_camera(npointing,south_camera=True,offset_ra=0,offset_dec=0,plotsingle=True):
    '''
    offset_ra in arcsec
    offset_dec in arcsec

    camera = 0 for south
    camera = 1 for north
    '''
    i = npointing-1
    if south_camera:
        # South camera
        delta_ra = -25./3600
        delta_dec = -2410./3600
        width_ra = 3.3/60.
        width_dec = 2.2/60.
        outputfile=outfile_prefix+'Pointing%02d-guiding-south.png'%(npointing)
        mytitle = "Pointing %02d South Camera (3.3'x2.2'): NSAID %i"%((i+1),pointing_id[i])
    else:
        # offsets and dimensions of 
        # North camera
        delta_ra = 0
        delta_dec = 2610./3600
        width_ra = 3.3/60.
        width_dec = 2.2/60.
        outputfile=outfile_prefix+'Pointing%02d-guiding-north.png'%(npointing)
        mytitle = "Pointing %02d North Camera (3.3'x2.2'): NSAID %i"%((i+1),pointing_id[i])
    center_ra = pointing_ra[i]+offset_ra/3600. + delta_ra
    center_dec = pointing_dec[i] + offset_dec/3600. + delta_dec

    delta_image = width_ra
    if plotsingle:
        plt.figure(figsize=(6,6))
    ax = plt.gca()
    pos = coords.SkyCoord(center_ra,center_dec,frame='icrs',unit='degree')
    xout = SkyView.get_images(pos,survey=['DSS'],height=2*delta_image*u.degree,width=2.*delta_image*u.degree)
    b=xout[0][0]
    ax.imshow(xout[0][0].data,interpolation='none',aspect='equal',cmap='gray_r',extent=[b.header['CRVAL1']-(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL1']+(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],
                                                           b.header['CRVAL2']+(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2'],
                                                           b.header['CRVAL2']-(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2']])

    rect= plt.Rectangle((center_ra - 0.5*width_ra,center_dec - 0.5*width_dec), width_ra, width_dec,fill=False, color='k')
    plt.gca().add_artist(rect)
    if plotsingle:
        plt.title(mytitle)
    else:
        if south_camera:
            plt.title("South Camera (3.3'x2.2')")
        else:
            plt.title("North Camera (3.3'x2.2')")
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.gca().invert_yaxis()
    if plotsingle:
        plt.savefig(outputfile)



        
def airmass_plots(kittpeak=True,ING=False,MLO=False):

    observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
    if kittpeak:
        print('plotting airmass curves for Kitt Peak')
        observing_location = EarthLocation.of_site('Kitt Peak')
        observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
        start_time = Time('2017-03-12 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2017-03-12 14:00:00')

    elif MLO:
        print('plotting airmass curves for MLO')
        observing_location = EarthLocation.of_site(u'Palomar')
        observer_site = Observer.at_site("Palomar", timezone="US/Pacific")
        # for run starting 2019-Apr-04 at MLO
        start_time = Time('2019-04-03 01:00:00') # need to enter UTC time, MLO UTC+6?
        end_time = Time('2019-04-03 14:00:00')
        
    elif ING:
        print('plotting airmass curves for INT')
        observing_location = EarthLocation.of_site(u'Roque de los Muchachos')
        observer_site = Observer.at_site("Roque de los Muchachos", timezone="GMT")
        # for run starting 2019-Feb-04 at INT
        start_time = Time('2019-02-04 19:00:00') # INT is on UTC
        end_time = Time('2019-02-05 07:00:00')

    #observing_time = Time('2017-05-19 07:00')  # 1am UTC=6pm AZ mountain time
    #observing_time = Time('2018-03-12 07:00')  # 1am UTC=6pm AZ mountain time
    #aa = AltAz(location=observing_location, obstime=observing_time)


    #for i in range(len(pointing_ra)):



    delta_t = end_time - start_time
    observing_time = start_time + delta_t*np.linspace(0, 1, 75)
    nplots = int(sum(obs_mass_flag)/8.)
    print(nplots)
    for j in range(nplots):
        plt.figure()
        legend_list = []
        for i in range(8):
            pointing_center = coords.SkyCoord(pointing_ra[8*j+i]*u.deg, pointing_dec[8*j+i]*u.deg, frame='icrs')
            if i == 3:
                plot_airmass(pointing_center,observer_site,observing_time,brightness_shading=True)
            else:
                plot_airmass(pointing_center,observer_site,observing_time)
            legend_list.append('Pointing %02d'%(8*j+i+1))
    
        plt.legend(legend_list)
        #plt.ylim(0.9,2.5)
        plt.gca().invert_yaxis()
        plt.subplots_adjust(bottom=.15)
        #plt.axvline(x=7*u.hour,ls='--',color='k')
        plt.axhline(y=2,ls='--',color='k')
        plt.savefig(outfile_prefix+'airmass-%02d.png'%(j+1))
        

    ##     delta_hours = np.linspace(0, 12, 100)*u.hour
    ##     full_night_times = observing_time + delta_hours
    ##     full_night_aa_frames = AltAz(location=observing_location, obstime=full_night_times)
    ##     full_night_aa_coos = pointing_center.transform_to(full_night_aa_frames)

    ##     plt.plot(delta_hours, full_night_aa_coos.secz,label='Pointing %02d'%(i+1))
    ##     plot_airmass(pointing_center, observing_location, observing_time)
    ##     plt.xlabel('Hours from 6 pm AZ time')
    ##     plt.ylabel('Airmass [Sec(z)]')
    ##     plt.ylim(0.9,2.5)
    ##     plt.tight_layout()

def check_CO():
    '''
    look for galaxies that are in Gianluca's CO file but not in the NSA catalog
    '''
    matchflag = find_CO_noNSA()
    CJnoNSA = ~matchflag
    for i in range(len(CJnoNSA)):
        if CJnoNSA[i]:
            print('%18s %18s %.8f %.8f'%(CJcat.source_name[i], CJcat.NEDname[i],CJcat.RA[i],CJcat.DEC[i]))
        
    '''
    results in
    number in CO catalog =  227
    number with NSA matches =  226
           UGC8656          UGC 08656 205.12995833 42.99388889
    '''

def make_INT_catalog():
    pos=coords.SkyCoord(pointing_ra*u.degree,pointing_dec*u.degree,frame='icrs')
    outfile = open('observing/finn_virgo.cat','w')
    for i in range(len(pointing_ra)):
        ra = pos[i].ra.hms
        dec = pos[i].dec.dms
        s = 'pointing-%03d %02d %02i %02.2f %02d %02i %02.2d J2000 ! NSAID %6s \n'%(i+1,ra[0],ra[1],ra[2],dec[0],dec[1],dec[2],str(pointing_id[i]))
        outfile.write(s)
    outfile.close()

def get_more_targets():
    # selecting targets in early part of night for Feb 2019 INT run
    # need to supplement CO sample b/c we are so efficient!
    # 
    plt.figure()
    plt.scatter(nsa.RA[more_targets_flag],nsa.DEC[more_targets_flag],c=mycolor[more_targets_flag],s=10,vmin=v1,vmax=v2,cmap='jet',marker='o')
    plt.colorbar(fraction=0.08, label=mylabel)
    #plt.axis([120,240,-10,50])
    plt.gca().invert_xaxis()
    #plt.ylim(0,50)

