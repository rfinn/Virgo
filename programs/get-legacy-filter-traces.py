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
from astropy.table import Tables, Column

filter_names = ['decamDR1-g', 'decamDR1-r',  'decamDR1-z',\
                'BASS-g', 'BASS-r', 'MzLS-z']

for f in filter_names:
    this_filter = speclite.filters.load_filters(f)
    c1 = Column
