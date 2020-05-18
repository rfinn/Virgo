#!/usr/bin/env python

    def get_CO_sources(self,match_by_coords=False,match_by_name=True):
        # read in CO mastertable
        # match to main table
        cofile = homedir+'/github/Virgo/tables/CO-MasterFile-2018Feb16.fits'
        # this file has a new column with the exact NED names
        # can use this to match to mastertable NEDname column
        cofile = homedir+'/research/Virgo/tables/CO-MasterFile-2018Feb16-fixedNEDnames.fits'        
        self.co = Table(fits.getdata(cofile))

        cocoord = SkyCoord(self.co['RA'],self.co['DEC'],unit='deg',frame='icrs')
        self.catcoord = SkyCoord(self.clean_a100['RA'],self.clean_a100['DEC'],unit='deg',frame='icrs')
        if match_by_coords:
            # match co coords to mastertable 
            idx, d2d, d3d = self.catcoord.match_to_catalog_sky(cocoord)
            self.d2d = d2d
            self.idx = idx
            self.coflag = d2d < 15./3600*u.deg
            # create a new, blank table with same # of lines as mastertable
            # but with columns like the co table
            newco = Table(np.zeros(len(self.basictable),dtype=self.co.dtype))
            # add co information into the new table for the galaxies with
            # a match to the CO sample
            # NOTE: this could match multiple galaxies to the same CO source
            newco[self.coflag] = self.co[idx[self.coflag]]

            # join basic table and co table

            self.cotable = hstack([self.basictable,newco])

        if match_by_name:
            #self.co.rename_column('NED_name','NEDname')
            #np.searchsorted(names1,names2)
            
            # match the basictable and the CO table by matching
            # entries by the NEDname colums
            self.clean_a100 = join(self.clean_a100,self.co,keys='NEDname',join_type='left')

            self.coflag = ~self.clean_a100['CO'].mask

            # also check to see which CO sources were not matched
            self.testtable = join(self.co,self.clean_a100,keys='NEDname',join_type='left')
            #self.coflag = len(self.co['CO']) > 0
            self.comatchflag = ~self.testtable['VFID'].mask
            print('CO sources with no match in mastertable:')
            print(self.testtable['NEDname','NED_name'][~self.comatchflag])
        ## plot the positions of CO galaxies that weren't matched to mastertable
        plt.figure()
        plt.plot(self.testtable['RA_1'][~self.comatchflag],self.testtable['DEC_1'][~self.comatchflag],'bo')

        ## print the CO galaxies with no matches in the mastertable 
        print('number of galaxies with CO matches = ',sum(self.coflag))

        ## look for CO galaxies that were matched to multiple galaxies in the
        ## mastertable
        unique, counts = duplicates(self.cotable,'NED_name')
        print("CO sources that are matched to multiple galaxies in the mastertable:")
        print(unique[counts>1])
        
        
        #self.cotable.add_column(Column(self.coflag),name='COFlag')
        #self.cotable.write(outdir+'vf_co.fits',format='fits',overwrite=True)

        # print CO sources that are not in the table

