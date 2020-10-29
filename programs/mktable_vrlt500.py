#!/usr/bin/env python

'''
GOAL
* to create a table with the same exact columns as smart_kitchen_sink.v2.fits but for additional galaxies with vr < 500 km/s
* we set a min vr cut to our original sample, but we now want to include galaxies with vr<500 but with redshift-independent distances that put them in the vicinity of Virgo.

INPUT TABLES
* benedetta's table of extra galaxies (none of the columns from this table should be in the output table)
* Hyperleda table
* agc
* nsa v1
* nsa v0





OUTPUT

NOTES:
If I could make a version of the HL table that just has the 117 galaxies to add, then I could run mksupersample.py, using this as the HL table instead of the big one.

However, one complication is that 4 galaxies are not in HL :(.

I could input 

'''
