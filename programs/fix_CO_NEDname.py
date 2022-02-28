import os
import numpy as np
import time

from astropy.io import fits
from astropy.table import Table, join, hstack, Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.ned import Ned

homedir = os.getenv("HOME")
#sys.path.append(homedir+'/github/appss/')
#from join_catalogs import make_new_cats, join_cats


#masterfile = homedir+'/research/Virgo/supersample/vf_clean_sample.fits'
#outdir = homedir+'/research/Virgo/tables/'

cofile = homedir+'/github/Virgo/tables/CO-MasterFile-2018Feb16.fits'
outfile = homedir+'/research/Virgo/tables/CO-MasterFile-2018Feb16-fixedNEDnames.fits'
nedname_field='NED_name'

cofile = homedir+'/research/Virgo/tables/All-virgo_Master_file_19Jun2019.fits'
outfile = homedir+'/research/Virgo/tables/All-virgo_Master_file_19Jun2019-fixedNEDnames.fits'
nedname_field='NEDname'

# updated catalog
cofile = homedir+'/research/Virgo/tables/galaxy_sample_prop_general_2020Oct28.fits'
outfile = homedir+'/research/Virgo/tables/galaxy_sample_prop_general_2020Oct28-fixedNEDnames.fits'

# updated catalog
cofile = homedir+'/research/Virgo/gianluca_input_tables/galaxy_sample_prop_general_mod.fits'
outfile = homedir+'/research/Virgo/gianluca_input_tables/galaxy_sample_prop_general_mod_fixedNEDnames.fits'
nedname_field='NEDname'
co = Table(fits.getdata(cofile))

realNEDname = np.zeros(len(co),dtype='|S30')

#coflag = np.zeros(len(catcoord),'bool')
          
for i in range(len(co)):
    # look up NED name
    ### FIX THE NAME FOR UGC8656
    ### THE CATALOG LIST UGC8656 BUT IT'S ACTUALLY  UGC8656 NOTES01

    try:
        time.sleep(1)
        if co[nedname_field][i] == 'UGC8656':
            t = Ned.query_object('UGC8656 NOTES01')
        elif co[nedname_field][i] == 'SHOC 206b':
            print('cheating on SHOC206b')
            t = Ned.query_object('SHOC 206a')                
        else:
            t = Ned.query_object(co[nedname_field][i])               
        realNEDname[i] = (t['Object Name'][0])

                        
    except IndexError:
        print(i,'2 NED did not like ',+str(co[nedname_field][i]))
    except:
        print(i,'2 NED did not like ',+str(co[nedname_field][i]))
# append realNEDname to co table

co.rename_column(nedname_field,'NEDnameCO')
c = Column(realNEDname)
co.add_column(c,name='NEDname')                        
co.write(outfile,format='fits',overwrite=True)
