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
from astropy.table import Table, join, hstack, Column, MaskedColumn 
from astropy.coordinates import SkyCoord
import astropy.units as u

from astroquery.ned import Ned

from matplotlib import pyplot as plt

from virgoCommon import *

class getNED:
    def __init__(self,clean_catalog):
        self.clean_a100 = Table(fits.getdata(clean_catalog))
    def get_NEDname(self):

        # skip for now until other parts are working
        # if file with NED names already exists from a previous query,
        # read the file
        if os.path.exists('ned_names.fits'):
            print('found file ned_names.fits.\nUsing this instead of querying NED')
            self.get_NEDname_from_file()
        else:
            #sys.exit() # for testing purposes
            self.get_NEDname_query()
    def get_NEDname_from_file(self):
        
        nednames = fits.getdata('/home/rfinn/research/Virgo/supersample/ned_names_v2.fits')
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
            self.newtab['NEDra'][flag] = nednames['NEDra'][flag2]
            self.newtab['NEDdec'][flag] = nednames['NEDdec'][flag2]
            self.newtab['NEDname'][flag] = 'UGC 09348'
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
        nedtable.write('ned_names_v2.fits',format='fits',overwrite=True)            

            
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

            nedtable = Table([c1,c2,c3,c4]).write('ned_names.fits',format='fits',overwrite=True)            
        except:
            print("couldn't write ned names")
            
        self.clean_a100.add_columns([c1,c2,c3,c4])
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)
    def get_NEDname_orig(self):
        # matching to catalog I downloaded on March 25, 2020
        # read in NED file
        '''
        Search by by Parameters
        decmin = -35  ** this is different from first download
        decmax = 75
        ramax = 280.
        ramin = 100.
        vmax = 3300. (units km/s)
        vmin = 500.
        object = Galaxy
        
        RA has to be in hours
        ramin = 6.6666666667
        ramax = 18.666666667

        output = text, ascii, bar separated
        velocity lower limit = -99

        Had to run this on virgo b/c I got locked out of NED.
        can't access the website from any machine in my house!

        Dowloaded text file - there was garbage at top of file (unnecessary info) that I deleted.

        Saved as

        /Users/rfinn/github/Virgo/tables/ned-noprolog-25mar2020.txt
        '''
        
        self.read_ned()
        self.cull_ned()

        # create empty array of names
        NEDid = np.zeros(len(self.clean_a100['RA']),dtype='|S25')

        
        # set up SkyCoord for clean cat and NED
        n = SkyCoord(self.ned['RA'],self.ned['DEC'],frame='icrs',unit='deg')
        c = SkyCoord(self.clean_a100['RA'],self.clean_a100['DEC'],frame='icrs',unit='deg')        
        # find matches
        idn, d2d, d3d = c.match_to_catalog_sky(n) # matches NED to VF catalog
        #print(d2d)

        nedmatchflag = (d2d < 10.*u.arcsec) # keep matches within 15arcsec

        indices = np.arange(len(idn))
        #print(indices[0:50])
        #print(len(indices),len(NEDid),len(idn))
        #NEDid[indices[nedmatchflag]] = self.ned['Object Name'][idn[nedmatchflag]]
        # transpose matched NED names into clean_a100
        
        # for each galaxy in cleaned catalog, find ned sources within 10 arcsec
        self.NEDmultiples = np.zeros(len(self.clean_a100),'bool')

        multiples = 0
        '''
        procedure:
        - calc offset between VF galaxy and NED galaxies
        - keep those w/in 15" offset
        - if only one source, assign VF galaxy this NED id
        - if multiples
          - try to match the HL name to the matching NED sources, if HL name exists
            - set multiple flag to False if match is found
          - look for source with name starting with NGC, UGC, IC.  use this if found
          - otherwise set to first NED source in the list
        '''
        
        for i in range(len(self.clean_a100)):
            d = np.sqrt((self.clean_a100['RA'][i] - self.ned['RA'])**2 \
                        + (self.clean_a100['DEC'][i] - self.ned['DEC'])**2)
            matchflag = d < 15./3600
            if np.sum(matchflag) > 1:
                self.NEDmultiples[i] = True
                #if multiples < 30:
                #    print('number of matches = ',np.sum(matchflag),multiples)
                #    print(self.ned['Object Name', 'Redshift Flag','Redshift Points'][matchflag],'\nHLname = ',self.clean_a100['objname'][i],'\n')
                multiples += 1

                
                # if we get multiple matches, then look at names and find one with
                # closest match to objname if it exist.
                indexmultiples = np.arange(len(self.ned))[matchflag]
                if len(self.clean_a100['objname'][i]) > 1: # has a HL name
                    n1 = self.clean_a100['objname'][i] # HL name
                    foundmatch = False
                    # look for match with HL name
                    for j in indexmultiples:
                        n2 = self.ned['Object Name'][j] # NED name
                        if (n2[0:2] == n1[0:2]) & (n2[-2:] == n1[-2:]):
                            # start and end of names match
                            # so this is the correct NED name
                            NEDid[i] = n2
                            foundmatch = True
                            self.NEDmultiples[i] = False
                            print('found HL name: ',n1,n2)
                            break
                if (foundmatch == False): # continue looking for the best match
                    # if not, look for more common names like NGC, UGC, IC
                    for j in indexmultiples:
                        if (self.ned['Object Name'][j].startswith('NGC')) | \
                           (self.ned['Object Name'][j].startswith('UGC')) | \
                           (self.ned['Object Name'][j].startswith('IC')):
                            # NED name
                            # start and end of names match
                            NEDid[i] = self.ned['Object Name'][j]
                            foundmatch = True
                            print('found a nice name: ',self.clean_a100['objname'][i],NEDid[i])
                            break
                if (foundmatch == False):
                    # assign closest match
                    # set nedquery flag to true, and query these objects individually
                    closest_match = np.arange(len(self.clean_a100))[np.where(d == min(d))[0][0]]
                    NEDid[i] = self.ned['Object Name'][closest_match]
                    print('taking closest match',self.clean_a100['objname'][i],NEDid[i])
                    # could also look at number of redshifts
            elif np.sum(matchflag) == 1:
                NEDid[i] = self.ned['Object Name'][np.where(matchflag)[0][0]]
                
                        
        print('number with multiple NED matches = ',multiples)
        print('number of NED matches = {}/{}'.format(sum(nedmatchflag),len(nedmatchflag)))
        
        c1 = Column(NEDid,name='NEDname')
        c2 = Column(self.NEDmultiples,name='NEDmultiples')# check these by hand?
        self.clean_a100.add_columns([c1,c2,c3])
        self.clean_a100.write('vf_clean_sample.fits',format='fits',overwrite=True)
    def get_GL_NEDname(self):
        # for galaxies with no NED match, use GL's catalog to match
        # HL name to NEDname (he did a position match for those that didn't return a NED name
        # when searching by name
        ref = fits.getdata(homedir+'/github/Virgo/tables/nsa_HyperLeda_NED_Steer2017dist_Virgo_field_sources_extension_H0_74_0_final_Kim2016corr_inclCOsample.fits')

        

    def write_clean(self):
        #self.newtab.write('vf_clean_sample_wNEDname.fits',format='fits',overwrite=True)
        # updating for v1
        self.newtab.write('vf_clean_sample_wNEDname_v1.fits',format='fits',overwrite=True)


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
    n.write_clean()
    
