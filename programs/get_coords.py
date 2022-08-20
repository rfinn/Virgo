#!/usr/bin/env python

"""
Get coordinates of VFID object, given the VFID of v2 catalog

"""

def get_coords(VFID, catdir=None):
    """
    args: VFID = pass in 4-digit VFID number, e.g. VFID0000
    returns: RA, DEC in degrees
    """
    if catdir is None:
        print("Please provide a valid catalog, path to vf_v2_main.fits")
        return
    else:
        from astropy.table import Table
        cat = Table.read(catdir)
        return cat['RA'][cat['VFID_V1']==VFID][0],cat['DEC'][cat['VFID']==VFID][0]

if __name__ == '__main__':
    import os
    import sys
    homedir = os.getenv("HOME")
    # location of main VF table
    catdir = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
    vfid = sys.argv[1]
    ra,dec = get_coords(vfid,catdir)
    print('coordinates for V1 VFID')
    print(f'RA = {ra}  DEC = {dec}')

