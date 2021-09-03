#!/usr/bin/env python

'''
GOAL:
* write out filter transmission and wavelength for legacy filters
* we will send the files to the authors of magphys so that we can run magphys on legacy photometry

REFERENCES:
https://speclite.readthedocs.io/en/latest/filters.html


NOTES:
* installed magphys on nyx


'''

import speclite.filters
from astropy.table import QTable

filter_names = ['decamDR1-g', 'decamDR1-r',  'decamDR1-z',\
                'BASS-g', 'BASS-r', 'MzLS-z']

for i,f in enumerate(filter_names):
    all_filters = speclite.filters.load_filters(f)
    t = QTable([all_filters[0].wavelength,all_filters[0].response],
               names = ('wavelength(A)','response'))
    outfile_name = f+'.ecsv'
    t.write(outfile_name,format='ascii.ecsv',overwrite=True)
