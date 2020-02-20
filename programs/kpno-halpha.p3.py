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

   airmass_plots(KPNO=False,ING=True)  # INT Feb 2019:

   or

   airmass_plots(KPNO=False,MLO=True)  # MLO April 2019


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
import os
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy import units as u
from astroquery.skyview import SkyView

from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time

from astropy.coordinates import AltAz

from astroplan import Observer
from astroplan.plots import plot_airmass

# prevent auto downloading of tables for airmass plots
#from astropy.utils import iers
#iers.conf.auto_download = False


########################################
###### FOR INT RUN              
###### set INGrun to True
###### otherwise, set to False
########################################

INGrun=False
HDIrun = True
########################################
###### RUN-SPECIFIC PARAMETERS  ########
########################################

if INGrun:
    telescope_run = '2019May/INT-2019May-227filter-'
    telescope_run = '2019May/INT-2019May-197filter-'
else:
    telescope_run = '2019June/MLO-2019June-'
    telescope_run = 'KPNO-2020Feb-'
    #telescope_run = '2019May/MLO-2019May-'
    run = '/2020Feb/'
    outfile_directory = outfile_directory+run

outfile_prefix = outfile_directory+telescope_run

max_pointing = None

#2019

max_pointing = None

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

cofile = 'nsa_CO-Gianluca.virgo.fits'
co = fits.getdata(tablepath+cofile)

nsa = fits.getdata(nsa_file)
jmass = fits.getdata(mass_file)
wise = fits.getdata(wise_file)
halpha = fits.getdata(halpha_file)
#if moretargets:
#    halpha = fits.getdata(halpha_file)
#else:
#    halpha = fits.getdata('/Users/rfinn/github/Virgo/tables/nsa_Halpha.virgo.2019Feb04.fits')
#    halpha = fits.getdata(tablepath+'nsa_Halpha.virgo.2019Feb04.fits')
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

# used ~(max vel of INT 197 - 250) for low-z end of filter gap
# used ~(min vel of INT 227 + 250) for high-z end of filter gap 
INTvflag =   (nsa.Z*3.e5 < 2100.) #| (nsa.Z*3.e5 > 2700.) #


gas_flag = COflag | HIflag
NGCfilament = filament

# Halpha selection
#add in targets that need to be reobserved
extra_targ_flag = (nsa.NSAID == 135136) #| (nsa.NSAID == 135129) ### Greg - why is this here? 135129 was observed on 2019-02-09 at INT

# set RA and DEC as galaxies with
# CO
# no Halpha
# stellar mass between 8.5 < log(M*/Msun) < 10.  according to NSF proposal
#obs_mass_flag = COsample & ~ha_obs #& (jmass.MSTAR_50 > 8.5) #& (jmass.MSTAR_50 < 10.) #& (nsa.SERSIC_BA > 0.2)
obs_mass_flag = (COsample & ~ha_obs) | extra_targ_flag
    
if INGrun:
    obs_mass_flag = obs_mass_flag & INTvflag

# resetting to COsample for 2019 INT observing season
#filter_flag = (nsa.Z*3.e5 > 2490.) & (nsa.Z*3.e5 < 6000.)
#obs_mass_flag = COsample & ~ha_obs #& filter_flag


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
pointing_mag = 22.5 - 2.5 * np.log10(nsa.NMGY[obs_mass_flag])


# sort by RA
sorted_indices = np.argsort(pointing_ra)
pointing_dec = pointing_dec[sorted_indices]
pointing_ra = pointing_ra[sorted_indices]
pointing_id = pointing_id[sorted_indices]
pointing_mag = pointing_mag[sorted_indices]
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

#######################
### REGULAR SAMPLE
#######################




##################################################
############ END OF INT WFC OFFSETS
##################################################
'''

# make a dictionary to store the offsets according to NSA ID
# format is offsets = {nsaid:[dra,ddec]} in arcmin
# offsets = {135046:[5.,4.],
#            84889:[3.,2.],
#            157073:[4.,0],
#            64408:[2.,-1]
#            }

# The following are INT vallues, which I deleted from the long list above
offsets_INT = {#135046:[5.,4.], # already observed
           #84889:[3.,2.],
           #157073:[4.,0],
           #64408:[2.,-1],
           147731:[0.,3],
           87097:[-3.,-3.],
           90957:[8.,9.],
           #50207:[-8.,12.],
           #64280:[-8.5,-3.],
           #64353:[14.5,12.],
           #135051:[1.5,1.5],
           15877:[4.,-7],
           #156774:[-4.,-2.],
           #135296:[-3.,3],
           135465:[0.,2.],
           #50379:[0.,3.],
           #135527:[0.,-2.5],
           #135602:[0,-2.5],
           #135606:[4.,0.],
           #50569:[-5.,3.],
           #135797:[-6.,4.],
           #47220:[-7.5,-2.],
           #135852:[-9.,-2.],
           #135862:[5.0,3.],
           #157256:[0.,13.],
           #64909:[-2.02,2.5],
           #136042:[0.,1.],
           #85513:[0.,13.],
           #85367:[0.,3.],
           #157480:[1.5,0],
           #157495:[3.5,3.],
           #107148:[-5.,-2.5],
           #85977:[0.,-2.],
           #48222:[-4.,3.3],
           #137045:[-5.,-3.],
           #107715:[0.,2.],
           137391:[10.5,-2.],
           #137460:[-8.,0.],
           #107764:[19.,-3.],
           #88142:[5.,-10.],
           #137993:[-7.,0.],
           90176:[0.,3.],
           #138221:[-3.6,2.2],
           87097:[-3.,-3.],
           #87086:[3.,0.],
           138642:[4.,2.],
           159520:[0.,-4.5],
           #159779:[0.75,-2.],
           101649:[8.,0.],
           #93963:[11.5,-2.],
           92459:[5.,-3.],
           #140301:[0.,-3.3],
           160627:[4.,-14.5],
           #117685:[0.,-1.5],
           118414:[-4.,0.],
           #143701:[0.,0.],
           #163875:[0.,2.],
           #143841:[-2.,12.],
           #144056:[21.5,-2.],
           #67567:[10.,-2.5],
           #17878:[3.5,2.],
           164911:[-6.,3.5],
           #165082:[-1.,-1.],
           18052:[3.,3.],
           165115:[-6.5,-3.],
           #145218:[10.,10.5],
           165200:[9.5,0],
           18153:[4.3,-12.],
           145398:[15.,2.],
           18301:[7.,-10.],
           #145554:[-2.,-10.],
           18363:[6.,14.],
           145672:[4.,-10.],
           165875:[10.,2.],
           145846:[18.,15.5],
           165956:[-3.,0.],
           145879:[-4.,12.],
           166280:[-6.,0.],
           #121129:[4.,9.],
           166330:[4.,-3.],
           68462:[-7.,-3.],
           69842:[16.,1.],
           147731:[0.,3.],
           135129:[0.,-3.],
           87100:[-8.,0],
           61693:[5.,0], # shift to get M51 on the chip :)
           165896:[6.7,2.5],
           15333:[7.,0],
           166297:[-5.,-2.5],
           146289:[9.,-2.5]
           
           }

# The following are INT vallues, which I deleted from the long list above
offsets_HDI = {#135046:[5.,4.], # already observed
           #84889:[3.,2.],
           #157073:[4.,0],
           #64408:[2.,-1],
           147731:[0.,3],
           #50207:[-8.,12.],
           #64280:[-8.5,-3.],
           #64353:[14.5,12.],
           #135051:[1.5,1.5],
           15333:[7.,-5],
           15877:[4.,-5],
           18052:[4.5,4.],    
           18153:[4.3,-8.],
           18301:[9.,-8.],
           18363:[4.,8.],
           19883:[0.,1.],    
           #156774:[-4.,-2.],
           #135296:[-3.,3],
           54619:[0.,1.],    
           56408:[0.,-13],# bright star very close to galaxy
           61693:[6.,-9],
           68462:[-10.,-7.],    
           80186:[1.,0.],    
           87097:[-2.,-4.],
           87100:[0.,7],
           90053:[-10.,-2.5],
           90176:[0.,3.],    
           90956:[0.,-1.],
           90957:[2.,2.],
           93977:[5.,0.],    
           101649:[0.,0.],
           102983:[-11.,7.],    
           117583:[-2,1.],# for guide star
           118414:[-6,4.],    
           125153:[1.5,0.],
           135465:[0.,0.],    
           136430:[1.2,-.5],# for guide star    
           137391:[10.5,-5.],    
           138642:[4.5,3.5],
           139741:[0.,-1.5],
           142509:[3.,0],# for guide star
           143305:[-0.5,4.0],    
           144151:[-1.,0],
           145398:[7.,4.5],
           145672:[6.,-12.5],
           145756:[2.,0],        
           145804:[7.,0],    
           #145846:[18.,15.5],
           145846:[8,8.5],
           146842:[0,-1.5],
           147100:[0,-1.5],
           147731:[13.,5.],    
           160627:[-1.,-10],
           162674:[-1.,-.5],
           163803:[-2.,-.5],        
           164911:[0,-3],
           165115:[0,0.],
           165200:[9.,-3],
           165862:[-2.,-5],    
           165875:[9.,-2.],
           165896:[6.7,-3],
           166297:[-4.,8],
           166330:[4.,-12.],# bright star to northe
           166548:[0.,-1.],# better position for guide stars
           166816:[0.,11.5],# bright star to south
           #50379:[0.,3.],
           #135527:[0.,-2.5],
           #135602:[0,-2.5],
           #135606:[4.,0.],
           #50569:[-5.,3.],
           #135797:[-6.,4.],
           #47220:[-7.5,-2.],
           #135852:[-9.,-2.],
           #135862:[5.0,3.],
           #157256:[0.,13.],
           #64909:[-2.02,2.5],
           #136042:[0.,1.],
           #85513:[0.,13.],
           #85367:[0.,3.],
           #157480:[1.5,0],
           #157495:[3.5,3.],
           #107148:[-5.,-2.5],
           #85977:[0.,-2.],
           #48222:[-4.,3.3],
           #137045:[-5.,-3.],
           #107715:[0.,2.],

           #137460:[-8.,0.],
           #107764:[19.,-3.],
           #88142:[5.,-10.],
           #137993:[-7.,0.],

           #138221:[-3.6,2.2],
           #87086:[3.,0.],

           159520:[0.,-4.5],
           #159779:[0.75,-2.],

           #93963:[11.5,-2.],
           92459:[14.5,-7.],
           #140301:[0.,-3.3],

           #117685:[0.,-1.5],

           #143701:[0.,0.],
           #163875:[0.,2.],
           #143841:[-2.,12.],
           #144056:[21.5,-2.],
           #67567:[10.,-2.5],
           #17878:[3.5,2.],

           #165082:[-1.,-1.],


           #145218:[10.,10.5],




           #145554:[-2.,-10.],




           165956:[-3.,0.],
           145879:[-4.,12.],
           166280:[-6.,0.],
           #121129:[4.,9.],


           69842:[16.,1.],

           135129:[0.,-3.],





           146289:[9.,-2.5]
           
           }

#I am commenting these out as the offsets aren't working well and I want to keep it simple by having one source per pointing.
offsets_MLO = {#87097:[3.,3.],
               # 90957:[-4.,3.5],


               # 61693:[0,-4.0],
               # 56478:[-4.5,0],
               # 54619:[0.,-3.],
               # 94217:[3.5,-1.5],
               # 165115:[0.,-3.0],
               # 165147:[0.,2.0],
               # 145398:[3.5,0.],
               # 18301:[0.,-3.0],
               # 165862:[-2.5,0.],
               # 165875:[0.,-2.5],
               # 165896:[1.0,-2.],
               # 145804:[3.0,0.],
               # 145814:[4.0,0.],
               # 15345:[-2.5,0],
               # 146289:[2.5,0.],
               # 166330:[-2.5,0.],
               # 166335:[-3.0,0.],
               # 114557:[3.0,0.],
               # 69842:[-2.5,0],
           }

# change this to use the offsets for the desired telescope
if INGrun:
    offsets = offsets_INT
elif HDIrun:
    offsets = offsets_HDI
else:
    offsets = {}
    
for key in offsets:
    try:
        pointing_offsets_ra[nsadict[key]] = offsets[key][0]/60.
        pointing_offsets_dec[nsadict[key]] = offsets[key][1]/60.
    except:
        print('problem setting offset for {} - already observed, or in INT filter gap?'.format(key))
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
    # delta_image is the half width of the image
    # delta_imagex and delta_imagey are the full width of the image in x and y directions, respectively
    delta_imagex=2.*delta_image
    delta_imagey=2.*delta_image
    if ING:
        delta_imagex=.8 #image width in deg
        delta_imagey=.8 # image width in deg
        #xout = SkyView.get_images(pos,survey=['DSS'],height=delta_imagex*u.degree,width=delta_imagey*u.degree)
    elif MLO:
        delta_imagex = 13./60. # image width in deg
        delta_imagey = 13./60 # image width in deg
    xout = SkyView.get_images(pos,survey=['DSS'],height=delta_imagex*u.degree,width=delta_imagey*u.degree)
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
    #for MLO plot galaxies outside the field of view to help with tweaking the pointing
    if MLO:
        source_pad = 8./60.
        gals = (nsa.RA > (pos.ra.value-(delta_imagex/2.+source_pad))) & (nsa.RA < (pos.ra.value+delta_imagex/2 + source_pad)) & (nsa.DEC > (pos.dec.value-(delta_imagey/2. + source_pad))) & (nsa.DEC < (pos.dec.value+delta_imagey/2. + source_pad))
    else:
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
        if noCOflag[j]:
            size=galsize-.005
            #rect= plt.Circle((ran-size/2.,decn-size/2.), size,fill=False, color='g')
            rect= plt.Circle((ran,decn), size,fill=False, color='g')
            plt.gca().add_artist(rect)
        if HIflag[j]:
            size=galsize+.005
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='b')
            plt.gca().add_artist(rect)
        if ha_obs[j]:
            size=galsize+2*.005
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
        plt.savefig(outfile_directory+telescope_run+'-Pointing%02d-NSA-%i.png'%(i+1, str(pointing_id[i])))

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

def make_all_platinum(KPNO=False,ING=False,MLO=False,startnumber=None):    
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
        plt.savefig(outfile_directory+telescope_run+'-Pointing%02d-NSA-%i-lowMass.png'%(i+1, str(pointing_id[i])))
        #plt.savefig(outfile_directory+'NSA-'+str(pointing_id[npointing-1])+'-'+telescope_run+'Pointing%03d-lowMass.png'%(npointing))            
        #plt.savefig(outfile_prefix+'Pointing%03d-lowMass-platinum.png'%(npointing))
    else:
        plt.savefig(outfile_directory+telescope_run+'-Pointing%02d-NSA-%i.png'%(i+1, str(pointing_id[i])))
        #plt.savefig(outfile_directory+'NSA-'+str(pointing_id[npointing-1])+'-'+telescope_run+'Pointing%03d.png'%(npointing))    
        #plt.savefig(outfile_prefix+'Pointing%03d-NSA-%i.png'%(npointing,pointing_id[npointing-1]))
        #plt.savefig(outfile_prefix+'Pointing%03d-platinum.png'%(npointing))

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



        
def airmass_plotsv2(KPNO=False,ING=False,MLO=False):

    observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")

    if KPNO:
        print('plotting airmass curves for Kitt Peak')
        observing_location = EarthLocation.of_site('Kitt Peak')
        observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
        start_time = Time('2017-03-12 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2017-03-12 14:00:00')
        start_time = Time('2020-02-24 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2020-02-24 14:00:00')
        
    elif MLO:
        print('plotting airmass curves for MLO')
        observing_location = EarthLocation.of_site(u'Palomar')
        observer_site = Observer.at_site("Palomar", timezone="US/Pacific")
        # for run starting 2019-Apr-04 at MLO
        #start_time = Time('2019-04-03 01:00:00') # need to enter UTC time, MLO UTC+6?
        #end_time = Time('2019-04-03 14:00:00')
        #start_time = Time('2019-05-04 01:00:00') # need to enter UTC time, MLO UTC+6?
        #end_time = Time('2019-05-04 14:00:00')
        start_time = Time('2019-06-04 01:00:00') # need to enter UTC time, MLO UTC+6?
        end_time = Time('2019-06-04 14:00:00')
        
    elif ING:
        print('plotting airmass curves for INT')
        observing_location = EarthLocation.of_site(u'Roque de los Muchachos')
        observer_site = Observer.at_site("Roque de los Muchachos", timezone="GMT")
        # for run starting 2019-Feb-04 at INT
        #start_time = Time('2019-02-04 19:00:00') # INT is on UTC
        #end_time = Time('2019-02-05 07:00:00')
        # for run starting 2019-May-29 at INT
        start_time = Time('2019-05-29 19:00:00') # INT is on UTC
        end_time = Time('2019-05-30 07:00:00')



    delta_t = end_time - start_time
    observing_time = start_time + delta_t*np.linspace(0, 1, 30)
    nplots = int(sum(obs_mass_flag)/8.)
    if (sum(obs_mass_flag)/8.) > nplots:
        remainder = sum(obs_mass_flag) - 8*nplots
        nplots += 1
        partial = True

    print(nplots)
    for j in range(nplots):
        print('nplots = ',nplots)
        plt.figure()
        legend_list = []
        if j == (nplots - 1):
            lastplot = remainder
        else:
            lastplot = 8
        for i in range(lastplot):
            print('\t galaxy number = ',i)
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
        plt.savefig(outfile_directory+'airmass-%02d.png'%(j+1))
        

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
def airmass_plots(KPNO=False,ING=False,MLO=False):

    observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")

    if KPNO:
        print('plotting airmass curves for Kitt Peak')
        observing_location = EarthLocation.of_site('Kitt Peak')
        observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
        start_time = Time('2017-03-12 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2017-03-12 14:00:00')
        start_time = Time('2020-02-24 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2020-02-24 14:00:00')
        
    elif MLO:
        print('plotting airmass curves for MLO')
        observing_location = EarthLocation.of_site(u'Palomar')
        observer_site = Observer.at_site("Palomar", timezone="US/Pacific")
        # for run starting 2019-Apr-04 at MLO
        #start_time = Time('2019-04-03 01:00:00') # need to enter UTC time, MLO UTC+6?
        #end_time = Time('2019-04-03 14:00:00')
        #start_time = Time('2019-05-04 01:00:00') # need to enter UTC time, MLO UTC+6?
        #end_time = Time('2019-05-04 14:00:00')
        start_time = Time('2019-06-04 01:00:00') # need to enter UTC time, MLO UTC+6?
        end_time = Time('2019-06-04 14:00:00')
        
    elif ING:
        print('plotting airmass curves for INT')
        observing_location = EarthLocation.of_site(u'Roque de los Muchachos')
        observer_site = Observer.at_site("Roque de los Muchachos", timezone="GMT")
        # for run starting 2019-Feb-04 at INT
        #start_time = Time('2019-02-04 19:00:00') # INT is on UTC
        #end_time = Time('2019-02-05 07:00:00')
        # for run starting 2019-May-29 at INT
        start_time = Time('2019-05-29 19:00:00') # INT is on UTC
        end_time = Time('2019-05-30 07:00:00')

    #observing_time = Time('2017-05-19 07:00')  # 1am UTC=6pm AZ mountain time
    #observing_time = Time('2018-03-12 07:00')  # 1am UTC=6pm AZ mountain time
    #aa = AltAz(location=observing_location, obstime=observing_time)


    #for i in range(len(pointing_ra)):



    delta_t = end_time - start_time
    observing_time = start_time + delta_t*np.linspace(0, 1, 30)
    nplots = int(sum(obs_mass_flag)/8.)
    if (sum(obs_mass_flag)/8.) > nplots:
        remainder = sum(obs_mass_flag) - 8*nplots
        nplots += 1
        partial = True

    print(nplots)
    for j in range(nplots):
        print('nplots = ',nplots)
        plt.figure()
        legend_list = []
        if j == (nplots - 1):
            lastplot = remainder
        else:
            lastplot = 8
        for i in range(lastplot):
            print('\t galaxy number = ',i)
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
        plt.savefig(outfile_directory+'airmass-%02d.png'%(j+1))
        

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

def make_MLO_catalog():
    #make two catalogs for an MLO run, one suitable for loading into
    #ACE and one to put in a google doc
    coord_cat = open(gitpath+'Virgo/observing/mlo_virgo.may2019.csv','w')
#    coord_cat = open(gitpath+'Virgo/observing/mlo_virgo.cat','w')
    #pointing_ra is a list of all sources that need to be observed, ordered by RA
    pos=coords.SkyCoord(pointing_ra*u.degree,pointing_dec*u.degree,frame='icrs')
    #s = '# Pointing\tNSAID\tRA\tDEC\trmag\n'
    s = '# Pointing,NSAID,RA,DEC,rmag\n'
    coord_cat.write(s)

    for i in range(len(pointing_ra)):
        ra = pos[i].ra.hms
        dec = pos[i].dec.dms
        rastr = '%02d:%02d:%02.2f'%(ra[0],ra[1],ra[2])
        decstr = '%02d:%02d:%02.2f'%(dec[0],dec[1],dec[2])
        #s = 'pointing-%03d %6s %02d %02i %02.2f %02d %02i %02.2d %02.2f\n'%(i+1,pointing_id[i],ra[0],ra[1],ra[2],dec[0],dec[1],dec[2], pointing_mag[i][4])
        s = '%03d, %7i, %11s, %11s, %02.2f\n'%(i+1,pointing_id[i], rastr, decstr, pointing_mag[i][4])
        #s = '%03d\t%7i\t%11s\t%11s\t%02.2f\n'%(i+1,pointing_id[i], rastr, decstr, pointing_mag[i][4])
        coord_cat.write(s)
    coord_cat.close()

    #make list of dither positions for each pointing
    dither_cat = open(gitpath+'Virgo/observing/dither_cat_MLO_virgo.csv','w')
    #pointing_ra is a list of all sources that need to be observed, ordered by RA
    pos=coords.SkyCoord(pointing_ra*u.degree,pointing_dec*u.degree,frame='icrs')
    s = '# Pointing, Dither, RA, DEC\n'
    dither_cat.write(s)

    #dither positions in asec relative to original position
    raoff = np.array([0, 0, 120., 120., -120., -120.])
    decoff = np.array([0, -120., 0, 120., 120., 0])
    # #convert to degress
    # raoff = raoff / 3600.
    # decoff = decoff / 3600.
    # #put in degree units
    # posoff = coords.SkyCoord(raoff * u.degree, decoff * u.degree, frame='icrs')
  
    #now write dither positions
    for i in range(len(pointing_ra)):

        #dither positions
        for jdith in range(len(raoff)):
            #apply offsets
            posdithra = pos[i].ra + raoff[jdith] / 3600. * u.degree
            posdithdec = pos[i].dec + decoff[jdith] / 3600 * u.degree
            #convert to HMS and DMS
            radith = posdithra.hms 
            decdith = posdithdec.dms
            #print out coordinates
            rastr = '%02d:%02d:%02.2f'%(radith[0],radith[1],radith[2])
            decstr = '%02d:%02d:%02.2f'%(decdith[0],decdith[1],decdith[2])
            s = '%03d, %1i, %11s, %11s\n'%(i+1,jdith, rastr, decstr)
            dither_cat.write(s)
        
    dither_cat.close()

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

