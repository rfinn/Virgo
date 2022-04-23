#!/usr/bin/env python


'''

Written by Rose Finn

PURPOSE:
- to help plan Halpha observations for Virgo Filament galaxies
- to make airmass plots
- to make finding charts

USEAGE:

start ipython in github/Virgo/

%run programs/observing-halpha.py


To make airmass plots:

   airmass_plots(KPNO=False,ING=True)  # INT Feb 2019:

   or

   airmass_plots(KPNO=False,MLO=True)  # MLO April 2019


To make finding charts:

   platinum_finding_chart(0,MLO=True,ING=False) # makes 1

   make_all_platinum(MLO=True, KPNO=False,ING=False)

   make_all_platinum(KPNO=True,ING=False,MLO=False,startnumber=None) # makes all, starting at number XX
   
To plot NE filament
make_plot()

To plot one galaxy:
zoomin(1)

OUTPUT:
* airmass plots are saved in the airmass_dir
* finding chart plots are saved in the finding_chart_dir


UPDATES
* 2021-02-26: doing major rewrite to use virgo v1 tables, in prep for Mar/Apr 2021 runs

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


from virgoCommon import *
import pylab as plt
import numpy as np
import os
import sys

from astroquery.sdss import SDSS
from astroquery.skyview import SkyView

from astropy import coordinates as coords
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
from astropy.coordinates import AltAz
from astropy.wcs import WCS

from astroplan import Observer
from astroplan.plots import plot_airmass


sys.path.append(homedir+'/github/Virgo/programs/')
import readtablesv2 as readtables


########################################
###### FOR INT RUN              
###### set INGrun to True
###### otherwise, set to False
########################################

INGrun=False
HDIrun = True

########################################
###### FOR BOK RUN              
###### set BOKrun to True
###### otherwise, set to False
########################################
BOKrun = True
HDIrun = False


########################################
###### RUN-SPECIFIC PARAMETERS  ########
########################################

## NOTE: outfile_directory is defined in virgoCommon

if INGrun:
    telescope_run = '2019May/INT-2019May-227filter-'
    telescope_run = '2019May/INT-2019May-197filter-'
elif BOKrun:
    runid = '2021Mar'
    runid = '2021Apr'
    finding_chart_dir = os.path.join(outfile_directory,runid,'finding-charts','')
    if not os.path.exists(finding_chart_dir):
        os.mkdir(finding_chart_dir)
    airmass_dir = os.path.join(outfile_directory,runid,'airmass','')
    if not os.path.exists(airmass_dir):
        os.mkdir(airmass_dir)
    telescope_run = '2021Mar/BOK-2021Mar-'
    telescope_run = '2021Apr/BOK-2021Apr-'    
    #telescope_run = '2021Apr/INT-2019May-197filter-'
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
moretargets = True
########################################
######  READ IN DATA TABLES  #######
########################################

homedir = os.getenv("HOME")
#v = readtables.vtables(homedir+'/research/Virgo/tables-north/v1/','vf_north_v1_')
v = readtables.vtables(homedir+'/research/Virgo/tables-north/v2/','vf_v2_')
v.read_all()

 

########################################
######  DEFINE SAMPLES  #######
########################################
COflag = v.main['COflag']
noCOflag = ~v.main['COflag']

# the following are galaxies that we still need to observe
# detected in CO but not yet observed in Halpha
ha_obs = v.main['HAobsflag']
# galaxies that we have observed AND detected in Halpha
ha_detect = v.main['HAobsflag'] & v.main['HAflag']

# galaxies that have NOT been observed in Halpha, but have been observed in CO
need_obs = (v.main['COflag'] & ~v.main['HAobsflag']) | (v.main['HAobsflag'] & (v.ha['FILT_COR'] > 2.))


HIflag = v.main['A100flag']

obs_mass_flag = need_obs

print(sum(obs_mass_flag))

# those observed at bok already
# the log says 2303, but it's 2302

observed_2021Mar = ['VFID2303','VFID1728','VFID2593','VFID1538','VFID2357',\
                    'VFID2821','VFID1573','VFID1572','VFID3299']#,'VFID0501'


for vf in observed_2021Mar:
    flag = v.main['VFID'] == vf
    if sum(flag) > 0:
        i = np.arange(len(v.main))[flag]
        #print(i,vf,': before setting flag, obs_mass_flag = ',obs_mass_flag[i])
        obs_mass_flag[i] = False
    else:
        print('trouble in paradise')
        print(vf)
        sys.exit()

print('after removing targets from bok mar run: ',sum(obs_mass_flag))
observed_2021Apr = ['VFID2947','VFID1901','VFID5541','VFID0638','VFID0422',\
                    'VFID1844','VFID2997']#,,'VFID2911','VFID0501'

for vf in observed_2021Apr:
    flag = v.main['VFID'] == vf
    if sum(flag) > 0:
        i = np.arange(len(v.main))[flag]
        obs_mass_flag[i] = False
    else:
        print('trouble in paradise')
        print(vf)
        sys.exit()


print('after removing targets from bok apr run: ',sum(obs_mass_flag))
# these are targets that are in the FOV of another pointing,
# so we don't need to do them separately

# first is 1728
bok_duplicates = ['VFID1683',\
                  # VFID2303
                  'VFID2315',\
                  # VFID 2593
                  'VFID2562','VFID2589','VFID2593','VFID2600','VFID2612',\
                  # VFID 1538
                  'VFID1463','VFID1473','VFID1475',\
                  # VFID 2357
                  'VFID2355','VFID2357','VFID2368',\
                  # VFID 2821
                  'VFID2771',\
                  # VFID1573
                  'VFID1548','VFID1552','VFID1561','VFID1574',\
                  'VFID1580','VFID1589','VFID1590','VFID1591',\
                  'VFID1595','VFID1596','VFID1597','VFID1606',\
                  'VFID1607','VFID1630',\
                  #VFID5541
                  'VFID5399','VFID5408','VFID5413','VFID5564',\
                  # VFID 0638
                  'VFID0568',\
                  # VFID1901
                  'VFID1834',\
                  # VFID1844
                  'VFID1822',\
                  # VFID0422
                  'VFID0377','VFID0385',\
                  # VFID3780
                  'VFID3739',\
                  # VFID0520
                  'VFID0469','VFID0474','VFID0487','VFID0496',\
                  'VFID0531','VFID0538',\
                  # VFID2911
                  'VFID2867','VFID2897','VFID2898','VFID2926',\
                  # VFID2997
                  'VFID2921','VFID2931','VFID2999',\
                  # VFID4894
                  'VFID4670','VFID4720','VFID4752','VFID4774',\
                  'VFID4827','VFID4835',\
                  # VFID6253
                  'VFID6165','VFID6195','VFID6242',\
                  # VFID3459
                  'VFID3445','VFID3454','VFID3457','VFID3469'\
                  ]
#                  ['VFID0607','VFID0611','VFID1979','VFID1955',\
#                  'VFID2315','VFID2368','VFID2766','VFID2797',\
#                  'VFID3457','VFID3459','VFID4835','VFID6132',\
#                  'VFID6379','VFID6406','VFID6438','VFID101'\
#                  ]
    
#for vf in bok_duplicates:
for vf in bok_duplicates:
    flag = v.main['VFID'] == vf
    if sum(flag) > 0:
        i = np.arange(len(v.main))[flag]
        #print(i,vf,': before setting flag, obs_mass_flag = ',obs_mass_flag[i])        
        obs_mass_flag[i] = False
    else:
        print('trouble in paradise')
        print(vf)
        sys.exit()



print('after removing galaxies that were in FOV of main targets from bok mar/apr runs: ',sum(obs_mass_flag))
'''
########################################
######  DEFINE FILAMENTS  #######
########################################
# NGC5353/4 Filament
radec = (v.main['RA'] > 192.) & (v.main['RA'] < 209) & (v.main['DEC'] > 0.) & (v.main['DEC'] < 50.) 
radec_flag = radec & (v.main['DEC'] >(2*(v.main['RA'] - 205.) + 20) ) & (v.main['DEC'] < (2*(v.main['RA'] - 205.) + 55))
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
extra_targ_flag = (v.main['VFID'] == 135136) #| (v.main['VFID'] == 135129) ### Greg - why is this here? 135129 was observed on 2019-02-09 at INT

# set RA and DEC as galaxies with
# CO
# no Halpha
# stellar mass between 8.5 < log(M*/Msun) < 10.  according to NSF proposal
#obs_mass_flag = COsample & ~ha_obs #& (jmass.MSTAR_50 > 8.5) #& (jmass.MSTAR_50 < 10.) #& (nsa.SERSIC_BA > 0.2)

obs_mass_flag = (COsample & ~ha_obs) | extra_targ_flag


if INGrun:
    obs_mass_flag = obs_mass_flag & INTvflag
'''



# SELECTING LOWER MASS TARGETS FOR INT RUN IN FEB 2019
#more_targets_flag = (nsa.Z*3.e5 < 2300.) & (v.main['RA'] > 115.) & (v.main['RA'] < 140.) & ~COsample & ~ha_obs & (v.main['DEC'] > 20.) & (v.main['DEC'] < 40.) & (jmass.MSTAR_50 > 8.5)

# FOR KPNO RUN IN FEB 2020
#more_targets_flag =  (v.main['RA'] > 115.) & (v.main['RA'] < 140.) & ~COsample & ~ha_obs & (v.main['DEC'] > 0.) & (v.main['DEC'] < 50.) & (jmass.MSTAR_50 > 8.5) &  (nsa.Z*3.e5 < 3000.)
moretargets = False
if moretargets == True:
    print('setting sample to more_targets_flag')
    obs_mass_flag = more_targets_flag

########################################
###### COLOR CODE FOR SCATTER PLOTS
########################################
mycolor=v.nsav0['MASS']
v1=minmass
v2=maxmass
mylabel='$ \log_{10}(M_*/M_\odot) $'


mycolor=v.main['vr']
v1=1000
v2=3000
mylabel='$ v_r \ (km/s) $'

########################################
###### FUNCTIONS
########################################
def convert_angle_2ra(angle,dec,dec2=None):
    ''' 
    DESCRIPTION:
    converts offset in deg to offset in ra;
    assumes angle, ra and dec are in deg;
    assumes both points are at the same declination by default;
    you can specify a second declination using dec2

    INPUT:
    * angle = angular offset in degrees
    * dec = declination of pointing in degrees
    * dec2 = optional, if angle does not correspond to a line of contant declination

    RETURNS:
    * equivalent angle in degrees of right ascension
    '''
    # a = angle in deg
    # b = 90 - dec
    # c = 90 - dec2, if different from dec
    # A = offset in longitude/RA
    # cos(a) = cos(b)*cos(c) + sin(b)*sin(c)*cos(A)
    # inverting gives
    # cos(A) = (cos(a)-cos(b)*cos(c))/[sin(b)*sin(c)]

    a = np.radians(angle)
    b = np.radians(90-dec)
    if dec2 is None:
        c = np.radians(90-dec)
    else:
        c = np.radians(90-dec2)
    delta_longitude_radians = np.arccos( (np.cos(a)-np.cos(b)*np.cos(c))/\
                                    (np.sin(b)*np.sin(c)))
    return np.degrees(delta_longitude_radians)
    
def finding_chart_all():
    if max_pointing != None:
        pointing_range = range(1,max_pointing+1)
    else:
        pointing_range = range(1,len(pointing_ra)+1)
                               
    for i in pointing_range:
        plt.close('all')
        finding_chart(i)

########################################
###### POINTING INFO FROM MAY 2017 #####
########################################


# keep only objects that have been observed in CO but NOT observed in Halpha
pointing_ra = v.main['RA'][obs_mass_flag]
pointing_dec = v.main['DEC'][obs_mass_flag]
pointing_id = v.main['VFID'][obs_mass_flag]
pointing_mag = 22.5 - 2.5 * np.log10(v.nsav0['NMGY'][:,4][obs_mass_flag])
pointing_vr = v.main['vr'][obs_mass_flag]

print('after cutting table to keep remaining targets: ',len(pointing_ra))

# sort by RA
sorted_indices = np.argsort(pointing_ra)
pointing_dec = pointing_dec[sorted_indices]
pointing_ra = pointing_ra[sorted_indices]
pointing_id = pointing_id[sorted_indices]
pointing_mag = pointing_mag[sorted_indices]
pointing_vr = pointing_vr[sorted_indices]
vfdict = dict((a,b) for a,b in zip(pointing_id,np.arange(len(pointing_id))))


print('after sorting arrays by RA: ',len(pointing_ra))

########################################
### OFFSETS TO GET MULTIPLE GALAXIES ###
########################################

####  OFFSETS  #####
# add offsets to try to get multiple galaxies in pointings
pointing_offsets_ra = np.zeros(len(pointing_ra))
pointing_offsets_dec = np.zeros(len(pointing_ra))

##################################################
############ INT WFC OFFSETS
##################################################

##################################################
############ low mass extension of leo filaments
##################################################


#######################
### REGULAR SAMPLE
#######################




##################################################
############ END OF INT WFC OFFSETS
##################################################


# make a dictionary to store the offsets according to VFID
# format is offsets = {nsaid:[dra,ddec]} in arcmin
# offsets = {135046:[5.,4.],
#            84889:[3.,2.],
#            157073:[4.,0],
#            64408:[2.,-1]
#            }

# The following are INT vallues, which I deleted from the long list above
offsets_INT = {#84889:[3.,2.],
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
# shifts in arcmin
offsets_BOK = {'VFID0377':[-10,0],\
               'VFID0422':[10,2],\
               'VFID0501':[0,-5],\
               'VFID0520':[0,-6],\
               'VFID0603':[-7,-1],\
               'VFID0638':[15,5],\
               'VFID0776':[-7,-7],\
               'VFID0783':[-5,0],\
               'VFID0788':[21,2],\
               'VFID0957':[-7,5],\
               'VFID1010':[-14,0],\
               'VFID1036':[5,12],\
               'VFID1168':[0,-4],\
               'VFID1213':[-12,-5],\
               'VFID1277':[5,-7],\
               'VFID1304':[-10,-5],\
               'VFID1376':[0,8],\
               'VFID1534':[-7,7],\
               'VFID1538':[-6,8],\
               'VFID1573':[18,-25],\
               'VFID1728':[5,0],\
               'VFID1756':[5,-5],\
               'VFID1819':[9,10],\
               'VFID1821':[10,-17],\
               'VFID1832':[0,9],\
               'VFID1901':[14,-5],\
               'VFID1944':[9,-9],\
               'VFID1956':[0,-4],\
               'VFID2068':[10,-3],\
               'VFID2098':[0,-10],\
               'VFID2162':[-3,-7],\
               'VFID2171':[-6,-7],\
               'VFID2259':[-5,-7],\
               'VFID2303':[25,-15],\
               'VFID2357':[2,-5],\
               'VFID2368':[15,2],\
               'VFID2484':[16,-6],\
               'VFID2488':[3,6],\
               'VFID2562':[5,0],\
               'VFID2593':[7,-8.5],\
               'VFID2601':[17,6],\
               'VFID2621':[2,3],\
               'VFID2661':[-5,-6],\
               'VFID2704':[12,4],\
               'VFID2762':[5,1.5],\
               'VFID2797':[8,0],\
               'VFID2821':[10,0],\
               'VFID2947':[8,-10],\
               'VFID2997':[11,4],\
               'VFID3098':[10,0],\
               'VFID3106':[7,0],\
               'VFID3119':[-3,0],\
               'VFID3272':[18,-8],\
               'VFID3299':[29,5],\
               'VFID3391':[10,-10],\
               'VFID3454':[14,-6],\
               'VFID3459':[8,-7],\
               'VFID3598':[-4,-10],\
               'VFID3714':[12,-9],\
               'VFID3780':[15,8],\
               'VFID3948':[10,0],\
               'VFID4025':[10,0],\
               'VFID4257':[13,0],\
               'VFID4279':[5,-8],\
               'VFID4796':[20,0],\
               'VFID4894':[14,4],\
               'VFID5541':[15,3],\
               'VFID5695':[9,-9],\
               #'VFID5726':[50,0],\
               'VFID5922':[10,-5],\
               'VFID5981':[20,-2],\
               'VFID6042':[15,-13],\
               'VFID6065':[2,12],\
               'VFID6115':[62,-8],\
               'VFID6127':[15,7],\
               'VFID6165':[20,-13],\
               'VFID6253':[15,4],\
               'VFID6293':[8,-5],\
               'VFID6305':[0,5],\
               'VFID6369':[54.,-6.5],\
               #'VFID6397':[40,5],\
               'VFID6406':[8,-4],\
               'VFID6425':[-1,3],\
               'VFID6426':[50,3],\
               'VFID6503':[50,12],\
               'VFID6599':[10,12],\
               'VFID6620':[0,10],\
               #'VFID2911':[11.5,-7],\
               # for alma proposal
               'VFID2911':[7,-11],\
           }

# change this to use the offsets for the desired telescope
if INGrun:
    offsets = offsets_INT
elif HDIrun:
    offsets = offsets_HDI
elif BOKrun:
    offsets = offsets_BOK
    
for key in offsets:
    try:
        pointing_offsets_ra[vfdict[key]] = offsets[key][0]/60.
        pointing_offsets_dec[vfdict[key]] = offsets[key][1]/60.
    except:
        print('problem setting offset for {} - already observed, or in INT filter gap?'.format(key))
# update pointing centers to account for offsets
pointing_ra += pointing_offsets_ra
pointing_dec += pointing_offsets_dec

if BOKrun:
    pointing_ra += -30./60 # 10' east
    pointing_dec += 20/60.

    # shift to bottom right chip
    pointing_ra += convert_angle_2ra(0.5,pointing_dec) # shift by 0.5 deg



def find_CO_noNSA():
    virgocat = coords.SkyCoord(v.main['RA']*u.degree,v.main['DEC']*u.degree,frame='icrs')
    jcat = coords.SkyCoord(CJcat.RA*u.degree,CJcat.DEC*u.degree,frame='icrs')

    index,dist2d,dist3d = jcat.match_to_catalog_sky(virgocat)

    # only keep matches with matched RA and Dec w/in 1 arcsec
    matchflag = dist2d.degree < 5./3600

    print('number in CO catalog = ',len(CJcat.RA))
    print('number with NSA matches = ',sum(matchflag))
    return matchflag
def add_detections():
    pointing_ra = v.main['RA'][COflag]
    pointing_dec = v.main['DEC'][COflag]
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
    
    plt.scatter(v.main['RA'][flag],v.main['DEC'][flag],s=30,c=mycolor[flag],vmin=v1,vmax=v2,cmap='jet',marker='o',label='CO')
    plt.scatter(v.main['RA'][noCOflag],v.main['DEC'][noCOflag],s=80,c=mycolor[noCOflag],vmin=v1,vmax=v2,cmap='jet',marker='x',label='no CO')
    plt.scatter(v.main['RA'][HIflag],v.main['DEC'][HIflag],s=100,c=mycolor[HIflag],vmin=v1,vmax=v2,cmap='jet',marker='+',label='HI')
    #plt.scatter(v.main['RA'][gas_flag],v.main['DEC'][gas_flag],s=50,c=jmass.MSTAR_50[gas_flag],vmin=8,vmax=11)
    if plotha:
        plt.plot(v.main['RA'][ha_obs],v.main['DEC'][ha_obs],'cs',mfc='None',markersize=8)
    if plotsingle:
        plt.colorbar(label=mylabel,fraction=0.08)
        plt.legend()
        #plt.gca().invert_xaxis()
def add_nsa():
    #plt.scatter(v.main['RA'][nsa_flag],v.main['DEC'][nsa_flag],s=10,c=jmass.MSTAR_50[nsa_flag],vmin=minmass,vmax=maxmass,alpha=.5,cmap='jet',marker='v',label='NSA')
    plt.scatter(v.main['RA'][nsa_flag],v.main['DEC'][nsa_flag],s=10,c=v.main['vr'],vmin=1500,vmax=3500,alpha=.5,cmap='jet',marker='v',label='NSA')

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
    #plt.axis([190,212,15,50])
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
    ''' correction for cosine(dec) term  '''

    # a = cos(dec)*cos(delta_ra)
    a = np.cos(np.radians(dec))*np.cos(np.radians(ra-b.header['CRVAL1']))
    # f 
    f = 1./b.header['CDELT2']*(180./np.pi)/\
        (np.sin(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))+a*np.cos(np.radians(b.header['CRVAL2'])))

    # -f * cos(dec)*sin(dec)-cos(dec)*cos(delta_ra)*sin(dec)
    decn = -f*(np.cos(np.radians(b.header['CRVAL2']))*np.sin(np.radians(dec))-a*np.sin(np.radians(b.header['CRVAL2'])))

    # -f * cos(dec)*sin(delta_ra)
    ran = -f*np.cos(np.radians(dec))*np.sin(np.radians(ra-b.header['CRVAL1']))

    ran = b.header['CRVAL1']+(ran)*b.header['CDELT1']
    decn = b.header['CRVAL2']-(decn)*b.header['CDELT2']    
    return ran,dec


def finding_chart(npointing,delta_image = .25,offset_ra=0.,offset_dec=0.,plotsingle=True,ING=False,MLO=False,KPNO=False,BOK=False):
    galsize=15
    i = npointing

    
    center_ra = pointing_ra[i]+offset_ra
    center_dec = pointing_dec[i] + offset_dec
    #print('center ra, dec = ',center_ra,center_dec)
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
    elif BOK: # 90prime
        delta_imagex = 1.2 # image width in deg
        delta_imagey = 1.2 # image width in deg
    #print(pos,delta_imagex)
    #print('running SkyView')
    #xout = SkyView.get_images(pos,survey=['DSS2 Red'],height=delta_imagex*u.degree,width=delta_imagey*u.degree)
    image_dimension = int(delta_imagex*60.)
    #print(image_dimension)
    #print(delta_imagex, delta_imagey)
    xout = SkyView.get_images(pos,survey=['DSS2 Red'],height=delta_imagex*u.degree,width=delta_imagey*u.degree)#,pixels=[image_dimension,image_dimension])
    b=xout[0][0]
    bwcs = WCS(b.header)
    ax = plt.subplot(projection=bwcs)
    ax.imshow(xout[0][0].data,interpolation='none',aspect='equal',cmap='gray_r',origin='lower')
    #          extent=[b.header['CRVAL1']-(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],\
    #                  b.header['CRVAL1']+(b.header['NAXIS1']-b.header['CRPIX1'])*b.header['CDELT1'],\
    #                  b.header['CRVAL2']+(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2'],\
    #                  b.header['CRVAL2']-(b.header['NAXIS2']-b.header['CRPIX2'])*b.header['CDELT2']])
    #ax.invert_yaxis()
    xra,ydec = ax.coords
    xra.set_major_formatter('d.ddd')
    ydec.set_major_formatter('d.ddd')    
    fits.writeto('latestim.fits',xout[0][0].data,header=b.header,overwrite=True)
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
    elif BOK:
        plot_BOK_footprint(center_ra,center_dec,header=b.header)
    #add_cameras()
    # find galaxies on FOV
    #for MLO plot galaxies outside the field of view to help with tweaking the pointing
    if MLO:
        source_pad = 8./60.
        gals = (v.main['RA'] > (pos.ra.value-(delta_imagex/2.+source_pad))) & (v.main['RA'] < (pos.ra.value+delta_imagex/2 + source_pad)) & (v.main['DEC'] > (pos.dec.value-(delta_imagey/2. + source_pad))) & (v.main['DEC'] < (pos.dec.value+delta_imagey/2. + source_pad))
    else:
        print(center_ra,center_dec,delta_imagex,delta_imagey)
        ramin=(center_ra - delta_imagex/2.)
        ramax=(center_ra + delta_imagex/2.)
        decmin=(center_dec - delta_imagey/2.)
        decmax=(center_dec + delta_imagey/2.)
        #print(ramin,ramax,decmin,decmax)
        ran,decn=fix_project(v.main['RA'],v.main['DEC'],b)        
        gals = (ran > ramin) & (ran < ramax) & \
            (decn > decmin) & (decn < decmax)

    print('pointing ',v.main['VFID'][npointing],' ngal = ',np.sum(gals))
    gindex=np.arange(len(v.main['RA']))[gals]
    galnames = v.main['VFID'][gindex]
    #print('Pointing %02d Galaxies'%(npointing),': ',v.main['VFID'][gals])
    #print('Pointing %02d Galaxies'%(npointing),': ',v.main['vr'][gals])
    
    for j in gindex:

        ran,decn=bwcs.wcs_world2pix(v.main['RA'][j],v.main['DEC'][j],1)        
        #print("RA=({:.6f},{:.6f}), DEC=({:.6f},{:.6f})".format(v.main['RA'][j],ran,v.main['DEC'][j],decn) )       
        rect= plt.Rectangle((ran-galsize/2.,decn-galsize/2.), galsize, galsize,fill=False, color='c')
        plt.gca().add_artist(rect)
        s='{}\n vr={:.0f}'.format(v.main['VFID'][j],v.main['vr'][j])
        plt.text(ran,decn+galsize/2.,s,fontsize=8,clip_on=True,horizontalalignment='center',verticalalignment='bottom')
        plt.text(ran,decn-galsize/2.,v.main['NEDname'][j],fontsize=8,clip_on=True,horizontalalignment='center',verticalalignment='top')
        if COflag[j]:
            size=galsize-5
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='g')
            plt.gca().add_artist(rect)
        #if noCOflag[j]:
        #    size=galsize-.005
        #    #rect= plt.Circle((ran-size/2.,decn-size/2.), size,fill=False, color='g')
        #    rect= plt.Circle((ran,decn), size,fill=False, color='g')
        #    plt.gca().add_artist(rect)
        if HIflag[j]:
            size=galsize+5
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='b')
            plt.gca().add_artist(rect)
        if ha_obs[j]:
            size=galsize+10
            rect= plt.Rectangle((ran-size/2.,decn-size/2.), size, size,fill=False, color='r')
            plt.gca().add_artist(rect)
        #plt.legend(['filament','CO','Halpha'])
    if moretargets:
        plt.title('LM-Pointing {}: {}'.format(pointing_id[npointing],pos.to_string('hmsdms')))
    else:
        gals_in_fov="VFIDs:"
        for i,g in enumerate(galnames):
            
            if i < len(galnames)-1:
                gals_in_fov += "{},".format(g.replace('VFID',""))
            else:
                gals_in_fov += "{}".format(g.replace('VFID',""))
                                            
        plt.title('Pointing {}: {}\n{}'.format(pointing_id[npointing],pos.to_string('hmsdms'),gals_in_fov))
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    #plt.gca().invert_xaxis()
    if plotsingle:
        plt.savefig(finding_chart_dir+'-Pointing-{}.png'.format(pointing_id[npointing]))
    return center_ra, center_dec
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
    

def plot_BOK_footprint(center_ra,center_dec,header=None):
    #using full detector sizes for now because
    pscale = .455
    xdim_pix = 4032
    ydim_pix = 4096
    # gap is 500", approximate 1060 pixels
    gap_dec = 60./3600
    gap_ra = convert_angle_2ra(160./3600,center_dec)
    #print('gap_ra = {:.2f} arcsec'.format(gap_ra*3600))
    #print('gap_dec = {:.2f} arcsec'.format(gap_dec*3600))    
    
    detector_dra = xdim_pix*pscale/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    detector_ddec = ydim_pix*pscale/3600. # 2154 pixels * 0.33"/pix, /3600 to get deg
    #print(detector_dra)
    # convert detector width in deg to deg of RA        
    detector_ra_width = convert_angle_2ra(detector_dra,center_dec)
    #detector_ra_width = detector_dra/np.cos(np.radians(center_dec))
    #print('detector_ra_width = {:.2f}'.format(detector_ra_width))
    offset_width_ra = detector_ra_width+gap_ra/2
    offset_width_dec = detector_ddec + gap_dec/2
    detector_offsets = [(-1*offset_width_ra,gap_dec/2),\
                        (-1*offset_width_ra,-1*offset_width_dec),\
                        (1*gap_ra/2,gap_dec/2),\
                        (1*gap_ra/2,-1*offset_width_dec)]
    # draw footprint of chip 4
    labels = ['1','2','3','4']
    for i,d in enumerate(detector_offsets):
        doffsetx,doffsety = d
        #print(doffsetx,doffsety)
        #print(center_ra+doffsetx,center_dec+doffsety)
        if header is not None:
            hwcs = WCS(header)
            # convert ra and dec of bottom left corner of ccd to x and y pixel values
            centerx,centery = hwcs.wcs_world2pix(center_ra,\
                                                 center_dec,1)
            plt.plot(centerx,centery,'rx')

            box_lower_x,box_lower_y = hwcs.wcs_world2pix(center_ra+doffsetx,\
                                                 center_dec+doffsety,1)
            #plt.plot(box_lower_x,box_lower_y,'bx')
            #plt.text(box_lower_x,box_lower_y,labels[i],fontsize=14)            
            dx = detector_dra/header['CDELT1']#/np.cos(np.radians(header['CRVAL2']))
            dy = detector_ddec/header['CDELT2']#/np.cos(np.radians(header['CRVAL2']))
            #print(centerx,centery,dx,dy)
            rect= plt.Rectangle((box_lower_x,box_lower_y), dx, dy,fill=False, color='k')
            rect2= plt.Rectangle((box_lower_x,box_lower_y), dx/2, dy/2,fill=False, color='.5')
            rect3= plt.Rectangle((box_lower_x+dx/2,box_lower_y+dy/2), \
                                 dx/2, dy/2,fill=False, color='.5')            
            
        else:
            rect= plt.Rectangle((center_ra+doffsetx,center_dec+doffsety), detector_dra, detector_ddec,fill=False, color='k')
            # try to mark amp quadrants
            rect2= plt.Rectangle((center_ra+doffsetx,center_dec+doffsety), detector_dra/2, detector_ddec/2,fill=False, color='.5')
            rect3= plt.Rectangle((center_ra+doffsetx+detector_dra/2,\
                                  center_dec+doffsety+detector_ddec/2), \
                                 detector_dra/2, detector_ddec/2,fill=False, color='.5')            

        plt.gca().add_artist(rect2)
        plt.gca().add_artist(rect3)                
        plt.gca().add_artist(rect)
    # adding guide camera
    #offset_dec = -2*detector_ddec-(22.7+98.1)/3600. # hard to explain
    #offset_ra =  detector_dra/2-(3.18+649.9)/3600.# hard to explain
    ## this chip is rotated 90 deg, so detecter_dra and detector_ddec are interchanged
    #rect= plt.Rectangle((center_ra+offset_ra,center_dec+offset_dec), -7./60., 7./60,fill=False, color='k')
    #plt.gca().add_artist(rect)
    

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

def make_all_platinum(KPNO=False,ING=False,MLO=False,BOK=False,startnumber=0):

    for i in range(startnumber,len(pointing_ra)):
        plt.close('all')
        #center_ra,center_dec = platinum_finding_chart(i,KPNO=KPNO,ING=ING,MLO=MLO,BOK=BOK)
        platinum_finding_chart(i,KPNO=KPNO,ING=ING,MLO=MLO,BOK=BOK)        

def platinum_finding_chart(npointing,offset_ra=0.,offset_dec=0.,ING=False,KPNO=False,MLO=False,BOK=False):
    if KPNO:
        fig = plt.figure(figsize = (10,6))
        grid = plt.GridSpec(2,3,hspace=.4,wspace=.2,left=.05)
        hdi = fig.add_subplot(grid[:,:-1])
        finding_chart(npointing,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False,ING=ING,MLO=MLO,KPNO=KPNO,BOK=BOK)
        south_camera = fig.add_subplot(grid[1,2])
        show_guide_camera(npointing,south_camera=True,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
        north_camera = fig.add_subplot(grid[0,2])
        show_guide_camera(npointing,south_camera=False,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False)
    else:
        fig = plt.figure(figsize = (8,8.))
        finding_chart(npointing,offset_ra=offset_ra,offset_dec=offset_dec,plotsingle=False,ING=ING,MLO=MLO,KPNO=KPNO,BOK=BOK)
    if moretargets:
        plt.savefig(finding_chart_dir+'-Pointing-{}-lowMass.png'.format(pointing_id[npointing]))
    else:
        #print('saving file ',outfile_directory+telescope_run+'Pointing-{}.png'.format(pointing_id[npointing]))
        plt.savefig(finding_chart_dir+'Pointing-{}.png'.format(pointing_id[npointing]))

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



        
def airmass_plotsv2(pointing_id,KPNO=False,ING=False,MLO=False):

    observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")

    if KPNO:
        print('plotting airmass curves for Kitt Peak')
        observing_location = EarthLocation.of_site('Kitt Peak')
        observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
        start_time = Time('2017-03-12 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2017-03-12 14:00:00')
        start_time = Time('2020-02-24 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2020-02-24 14:00:00')
        start_time = Time('2021-03-14 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2021-03-14 14:00:00')
        # April 2021 Run
        start_time = Time('2021-04-15 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2021-04-15 14:00:00')
        
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
    ngal=0
    for j in range(nplots):
        print('plot {} of {}'.format(j+1,nplots))
        plt.figure(figsize=(8,6))
        legend_list = []
        if j == (nplots - 1):
            lastplot = remainder
        else:
            lastplot = 8
        for i in range(lastplot):
            #print('\t galaxy number = ',i)
            pointing_center = coords.SkyCoord(pointing_ra[8*j+i]*u.deg, pointing_dec[8*j+i]*u.deg, frame='icrs')
            if i == 3:
                plot_airmass(pointing_center,observer_site,observing_time,brightness_shading=True)
            else:
                plot_airmass(pointing_center,observer_site,observing_time)
            legend_list.append('Pointing {}'.format(pointing_id[ngal]))
            ngal+= 1 
        
    
        plt.legend(legend_list)
        #plt.ylim(0.9,2.5)
        plt.gca().invert_yaxis()
        plt.subplots_adjust(bottom=.15)
        #plt.axvline(x=7*u.hour,ls='--',color='k')
        plt.axhline(y=2,ls='--',color='k')
        plt.savefig(airmass_dir+'airmass-%02d.png'%(j+1))
        

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
    # prevent auto downloading of tables for airmass plots
    from astropy.utils import iers
    iers.conf.auto_download = False

    observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")

    if KPNO:
        print('plotting airmass curves for Kitt Peak')
        observing_location = EarthLocation.of_site('Kitt Peak')
        observer_site = Observer.at_site("Kitt Peak", timezone="US/Mountain")
        start_time = Time('2017-03-12 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2017-03-12 14:00:00')
        
        start_time = Time('2020-02-24 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2020-02-24 14:00:00')

        # 2021 Mar run
        start_time = Time('2021-03-13 01:00:00') # UTC time, so 1:00 UTC = 6 pm AZ mountain time
        end_time = Time('2021-03-13 14:00:00')
        
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
    plt.scatter(v.main['RA'][more_targets_flag],v.main['DEC'][more_targets_flag],c=mycolor[more_targets_flag],s=10,vmin=v1,vmax=v2,cmap='jet',marker='o')
    plt.colorbar(fraction=0.08, label=mylabel)
    #plt.axis([120,240,-10,50])
    plt.gca().invert_xaxis()
    #plt.ylim(0,50)

def ngcfilament():
    ### FOCUS ON NGC 5353
    fig=plt.figure(figsize=(6,5))
    plt.plot(pointing_ra,pointing_dec,'ko',c='0.7',markersize=4,alpha=0.2)

    # NGC5353/4 Filament
    plt.scatter(pointing_ra[NGCfilament],pointing_dec[NGCfilament],c=pointing_vr[NGCfilament],zorder=20,s=20,vmin=1000,vmax=3000,lw=0.5,cmap='jet')
    
    xl = np.linspace(196,230,100)
    yl = (2*(xl - 205.) + 20)
    #plt.plot(xl,yl,'r-')
    #NGCfilament = filament
    plt.colorbar()
    plt.gca().invert_xaxis()
    #return filament

def dither_BOK(ra,dec,pointing):
    for i in range(len(ra)):
        pass
