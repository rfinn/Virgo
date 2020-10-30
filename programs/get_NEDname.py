#!/usr/bin/env python
'''
get the NED name for the cleaned sample by querying NED

INPUT: vf_clean_sample.fits

OUTPUT: vf_clean_sample_wNEDname.fits

PROCEDURE:
- first try to use the file /home/rfinn/research/Virgo/supersample/ned_names.fits
- look up each object according to its vf_clean_sample.superName 
  - check this against the ned_names.NEDinput
  - ned_names.NEDname is the corresponding ned name

'''

import os
import numpy as np
import time

from astropy.io import fits, ascii
from astropy.table import Table, join, hstack, vstack, Column, MaskedColumn 
from astropy.coordinates import SkyCoord
import astropy.units as u

from astroquery.ned import Ned

from matplotlib import pyplot as plt

from virgoCommon import *


import argparse
parser = argparse.ArgumentParser(description ='Create a crazy big catalog from HL, AGC, NSA')
parser.add_argument('--version',dest = 'version', default='v1',help='version of tables. default is v1')
parser.add_argument('--evcc',dest = 'evcc', default=False,action='store_true',help='run for evcc catalog containing galaxies not in our original table')
        
args = parser.parse_args()

if args.evcc:
    outfile_suffix = '_'+args.version+'_evcc'
else:
    outfile_suffix = '_'+args.version


class getNED:
    def __init__(self,clean_catalog):
        self.clean_a100 = Table(fits.getdata(clean_catalog))
    def get_NEDname(self):

        # skip for now until other parts are working
        # if file with NED names already exists from a previous query,
        # read the file

        self.nedfile = 'ned_names'+outfile_suffix+'.fits'
        if os.path.exists(self.nedfile):
            print('found file '+self.nedfile+'\nUsing this instead of querying NED')
            self.get_NEDname_from_file()
            self.write_clean()            
        else:
            #sys.exit() # for testing purposes
            print('querying NED names from database')
            self.get_NEDname_query()
    def get_NEDname_from_file(self):
        
        nednames = Table(fits.getdata(self.nedfile))
        self.nednames = nednames
        # do a left join of the input catalog and nednames
        # using column nednames.NEDinput and clean_a100.superName

        #self.clean_a100.rename_column('superName','NEDinput')
        self.clean_a100.add_column(self.clean_a100['superName'],name='NEDinput')

        # join returns a table that is sorted by the key columns
        # the following command gets table back into its original order
        self.newtab = myjoinleft(self.clean_a100,nednames,keys='NEDinput')

        
        # fix case for UGC09348, until we run the full query again...
        # this is where HL has the wrong coordinates/association
        flag = self.newtab['NEDinput'] == 'UGC09348'
        if (sum(flag) > 0) & self.newtab['NEDra'].mask[flag]: # hasn't been matched
            # assign ned RA and DEc from NGC 5658
            flag2 = nednames['NEDname']=='NGC 5658'
            if sum(flag2) > 0:
                self.newtab['NEDra'][flag] = nednames['NEDra'][flag2]
                self.newtab['NEDdec'][flag] = nednames['NEDdec'][flag2]
                self.newtab['NEDname'][flag] = 'UGC 09348'
            else:
                print('ruh roh')
        # fix case for PGC2586382 , until we run the full query again...
        # NED doesn't recognize 'PGC2586382', but it does know LEDA 2586382
        # it will also find the NSA ID, so I need to not just use the superName
        # I need to keep stepping through other name options if it doesn't return
        # a match to HL or AGC or NSA name
        flag = self.newtab['NEDinput']==  'PGC2586382'
        if (sum(flag) > 0) & self.newtab['NEDra'].mask[flag]: # hasn't been matched
            self.newtab['NEDra'][flag] = 226.398865 
            self.newtab['NEDdec'][flag] = 59.093766
            self.newtab['NEDname'][flag] = 'WISEA J150535.77+590537.2'

        # fix case for NGC 2793, until we run the full query again...
        # The NSA ID from version 1 doesn't get a NED match, but the
        # NSA ID from version 0 does...
        flag = self.newtab['NSAID_2']==135797
        if (sum(flag) > 0) & self.newtab['NEDra'].mask[flag]: # hasn't been matched
            self.newtab['NEDra'][flag] = 139.197125
            self.newtab['NEDdec'][flag] = 34.429806
            self.newtab['NEDname'][flag] = 'NGC 2793'
            

        # fix case for UGC 08, until we run the full query again...
        # The NSA ID from version 1 doesn't get a NED match, but the
        # NSA ID from version 0 does...
        flag = self.newtab['superName']=='UGC08656 NOTES01'
        if (sum(flag) > 0) & self.newtab['NEDra'].mask[flag]: # hasn't been matched
            self.newtab['NEDra'][flag] = 205.129878
            self.newtab['NEDdec'][flag] = 42.993819
            self.newtab['NEDname'][flag] = 'UGC 08656 NOTES01'
            
            
        # for those without a match by NEDname, look through the
        # entries with NEDinput = 'byposition'
        
        no_match_by_name = self.newtab['NEDname'].mask
        print(no_match_by_name)
        pos_match_indices = np.arange(len(self.clean_a100))[no_match_by_name]

        # shorten NED catalog to those that were matched by location
        # this should make search by location faster
        self.nednames_bypos = nednames[nednames['NEDinput'] == 'bylocation']
        #nedcoord = SkyCoord(self.nednames_bypos['NEDra'],self.nednames_bypos['NEDdec'],unit='deg',frame='icrs')            
        #self.catcoord = SkyCoord(self.clean_a100['RA'],self.clean_a100['DEC'],unit='deg',frame='icrs')            
        for i in pos_match_indices:
            d = np.sqrt((self.newtab['RA'][i]-self.nednames_bypos['NEDra'])**2 + \
                        (self.newtab['DEC'][i]-self.nednames_bypos['NEDdec'])**2)
            if (min(d) < 10./3600):
                j = np.where(d == min(d))[0][0]                
                #print(i,j,'found a match!')

                self.newtab['NEDinput'][i]='bylocation'
                self.newtab['NEDname'][i] = self.nednames_bypos['NEDname'][j]
                self.newtab['NEDra'][i] = self.nednames_bypos['NEDra'][j]
                self.newtab['NEDdec'][i] = self.nednames_bypos['NEDdec'][j]
        self.newtab.write('test.fits',format='fits',overwrite=True)
        
        # for those with no match, query NED
        
    def query_unmatched(self,startindex=0):
        ## running query for those that haven't yet been matched to an entry in the catalog
        ##

        no_match_by_name = self.clean_a100['NEDname'].mask
        pos_match_indices = np.arange(len(self.clean_a100))[no_match_by_name]
        
        
        for i in pos_match_indices:
            foundit=False
            
            time.sleep(.5)
            coord = SkyCoord(self.clean_a100['RA'][i],self.clean_a100['DEC'][i],unit=(u.deg,u.deg), frame='icrs')
            print(i,self.clean_a100['RA'][i],self.clean_a100['DEC'][i])
            
            try:
                t = Ned.query_region(coord,radius=10./3600*u.deg, equinox='J2000')
                print(t)
                #t = Ned.query_object('NSA '+str(self.clean_a100['NSAID_2'][i]))
                self.clean_a100['NEDname'][i] = t['Object Name'][0]
                self.clean_a100['NEDra'][i] = t['RA'][0]
                self.clean_a100['NEDdec'][i] = t['DEC'][0]
                self.clean_a100['NEDinput'][i] = 'bylocation'
                foundit=True
                continue
            except IndexError:
                print('IndexError')
                pass
            except:
                print(i,'NED did not like search by coordinate for this object')

            if not(foundit):
                print("oh no - could not find NED name! so sorry...",i)
                self.clean_a100['NEDname'][i] = ''
                self.clean_a100['NEDra'][i] = -999
                self.clean_a100['NEDdec'][i] = -999
                

    def write_NEDnames(self):
        nedtable = Table([self.clean_a100['NEDinput'],self.clean_a100['NEDra'],self.clean_a100['NEDdec'],self.clean_a100['NEDname']])
        nedtable.write(self.nedfile,format='fits',overwrite=True)            

            
    def get_NEDname_query(self,startindex=0):
        ## GOT BLOCKED BY NED FOR TOO MANY QUERIES
        ## TRYING ANOTHER APPROACH - TO MATCH TO CATALOG I DOWNLOADED FROM DEC
        
        # look up NED name for each galaxy
        # https://astroquery.readthedocs.io/en/latest/ned/ned.html


        # if in HL, look up HL name
        NEDid = []
        NEDra = []
        NEDdec = []
        NEDinput = []
        
        for i in range(startindex,len(self.clean_a100['objname'])):
        #for i in range(20):            
            # check if HL id exists, look up NED name
            # if yes, break
            #
            # if no NED comes back
            # check if NSA ID exists for that galaxy. if yes, look up NED name
            # if NED name is found, break
            #
            # if no NED name comes back
            # check if A100 ID exists for that galaxy. if yes, look up NED name
            # if no Ned name comes back, then no NED name!
            #
            foundit=False
            if self.clean_a100['HLflag'][i]:
                time.sleep(1)
                try:
                    t = Ned.query_object(self.clean_a100['objname'][i])
                    NEDid.append(t['Object Name'][0])
                    NEDra.append(t['RA'][0])
                    NEDdec.append(t['DEC'][0])
                    NEDinput.append(self.clean_a100['objname'][i])
                    foundit=True
                    continue
                except IndexError:
                    pass
                except:
                    print(i,'2 NED did not like ',self.clean_a100['objname'][i])
                    pass
                

            if self.clean_a100['NSAflag'][i]:
                if not(foundit):
                    time.sleep(1)
                    try:
                        t = Ned.query_object('NSA '+str(self.clean_a100['NSAID'][i]))
                        NEDid.append(t['Object Name'][0])
                        NEDra.append(t['RA'][0])
                        NEDdec.append(t['DEC'][0])
                        NEDinput.append('NSA '+str(self.clean_a100['NSAID'][i]))
                        foundit=True
                        continue
                    except IndexError:
                        pass
                    except:
                        print(i,'2 NED did not like ','NSA '+str(self.clean_a100['NSAID'][i]))
                        pass
            if self.clean_a100['A100flag'][i]:
                if not(foundit):
                    time.sleep(1)
                    try:
                        #print('AGC'+str(self.clean_a100['AGC'][i]))
                        t = Ned.query_object('AGC'+str(self.clean_a100['AGC'][i]))
                        NEDid.append(t['Object Name'][0])
                        NEDra.append(t['RA'][0])
                        NEDdec.append(t['DEC'][0])
                        NEDinput.append('AGC'+str(self.clean_a100['AGC'][i]))
                        foundit=True
                        continue
                    except IndexError:
                        pass
                    except:
                        print(i,'2 NED did not like ','AGC'+str(self.clean_a100['AGC'][i]))
                        pass
            if self.clean_a100['NSA0flag'][i]:
                if not(foundit):
                    time.sleep(1)
                    try:
                        t = Ned.query_object('NSA '+str(self.clean_a100['NSAID_2'][i]))
                        NEDid.append(t['Object Name'][0])
                        NEDra.append(t['RA'][0])
                        NEDdec.append(t['DEC'][0])
                        NEDinput.append('NSA '+str(self.clean_a100['NSAID_2'][i]))
                        foundit=True
                        continue
                    except IndexError:
                        pass
                    except:
                        print(i,'2 NED did not like ','NSA '+str(self.clean_a100['NSAID'][i]))
                        pass

            # check to make sure that the NED ra and dec
            # is not offset from object
            # we have a few cases where the NED object is wrong
            if foundit:
                distance = np.sqrt((self.clean_a100['RA'][i]-NEDra[-1])**2 + \
                                   (self.clean_a100['DEC'][i] - NEDdec[-1])**2)
                if (distance > 10./3600):
                    print('NED name is offset from source by {1:.3e} arcsec'.format(distance*3600))
                    print('resetting foundit to false')
                    foundit = False
            if not(foundit):
                # search by coordinates
                time.sleep(1)
                coord = SkyCoord(self.clean_a100['RA'][i],self.clean_a100['DEC'][i],unit=(u.deg,u.deg), frame='icrs')
                print(i,self.clean_a100['RA'][i],self.clean_a100['DEC'][i])

                try:
                    t = Ned.query_region(coord,radius=10./3600*u.deg, equinox='J2000')
                    print(t)
                    #t = Ned.query_object('NSA '+str(self.clean_a100['NSAID_2'][i]))

                    # should check if NED name already points to another galaxy
                    # and remove it if it does
                    if np.sum(np.array(NEDid) == t['Object Name'][0]) > 0:
                        # galaxy is already in the list, so position match is wrong
                        print('galaxy ', t['Object Name'][0],' is already in the list')
                    else:
                    
                        NEDid.append(t['Object Name'][0])
                        NEDra.append(t['RA'][0])
                        NEDdec.append(t['DEC'][0])
                        NEDinput.append('bylocation')

                    
                    
                        foundit=True
                        continue
                except IndexError:
                    pass
                except:
                    print(i,'NED did not like search by coordinate for this object')

            if not(foundit):
                print("oh no - could not find NED name! so sorry...",i)
                NEDid.append('')
                NEDra.append(-999)
                NEDdec.append(-999)
                NEDinput.append(-999)
            else:
                print(i,'found NED name')

        # ANOTHER CHECK TO ADD
        # find any remaining duplicates
        # if only one has NEDinput == 'bylocation', 
        # then remove the ned name for the duplicate that has 'bylocation'

        
        c1 = Column(NEDid,name='NEDname')
        c2 = Column(np.array(NEDra,'f'),name='NEDra')
        c3 = Column(np.array(NEDdec,'f'),name='NEDdec')
        c4 = Column(NEDinput,name='NEDinput')        
        self.NEDid = NEDid
        self.NEDra = NEDra
        self.NEDdec = NEDdec
        try:
            #nedtable = Table([c1,c2,c3,c4]).write('ned_names.tbl',format='ipac',overwrite=True)

            self.nedtable = Table([c1,c2,c3,c4]).write('ned_names'+outfile_suffix+'.fits',format='fits',overwrite=True)            
        except:
            print("couldn't write ned names")
            
        self.clean_a100.add_columns([c1,c2,c3,c4])
        self.clean_a100.write('vf_clean_sample_wNEDname'+outfile_suffix+'.fits',format='fits',overwrite=True)
    def get_GL_NEDname(self):
        # for galaxies with no NED match, use GL's catalog to match
        # HL name to NEDname (he did a position match for those that didn't return a NED name
        # when searching by name
        ref = fits.getdata(homedir+'/github/Virgo/tables/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')

        

    def write_clean(self):
        #self.newtab.write('vf_clean_sample_wNEDname.fits',format='fits',overwrite=True)
        # updating for v1
        self.newtab.write('vf_clean_sample_wNEDname'+outfile_suffix+'.fits',format='fits',overwrite=True)


if __name__ == '__main__':
    os.chdir('/home/rfinn/research/Virgo/supersample/')
    ###################################################################
    #### INPUT FILES
    ###################################################################
    #n = getNED('vf_clean_sample.fits') # for v0
    n = getNED('vf_clean_sample_v1.fits') # for v1    
    n.get_NEDname()
    #n.query_unmatched()
    #n.write_NEDnames()
    
    
