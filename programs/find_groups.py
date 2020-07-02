#!/usr/bin/env python

'''
Run John Moustakas's group-finding code on VF


'''
from astropy.table import Table, Column, vstack, join
import sys
import numpy as np
import time
#sys.path.append('/home/rfinn/github/SGA/bin')
#import SGAbuildparent as SGAbp

from astropy.coordinates import SkyCoord

def degrees_between(ra1,dec1,ra2,dec2):
    c1 = SkyCoord(ra1, dec1,frame="icrs",unit="deg")
    c2 = SkyCoord(ra2, dec2,frame='icrs',unit="deg")

    return c1.separation(c2).deg

def build_group_catalog(cat, mfac=2.0, dmax=10.0/60.0):
    """dmax in arcmin

    Group SGA galaxies together where their circular radii would overlap.  Use
    the catalog D25 diameters (in arcmin) multiplied by a scaling factor MFAC.
    The output catalog adds the column GROUP_ID which is unique for each group.
    The column MULT_GROUP is the multiplicity of that galaxy's group.

    """
    from pydl.pydlutils.spheregroup import spheregroup

    ## RAF
    ## REPLACING ASTRONOMETRY FUNCTION WITH ASTROPY 
    #from astrometry.util.starutil_numpy import degrees_between
    #print('Starting spheregrouping.')

    nchar = np.max([len(gg) for gg in cat['GALAXY']])+6 # add six characters for "_GROUP"
    
    t0 = time.time()
    cat.add_column(Column(name='GROUP_ID', data=np.zeros(len(cat), dtype=np.int)-1))
    cat.add_column(Column(name='GROUP_NAME', length=len(cat), dtype='<U{}'.format(nchar)))
    cat.add_column(Column(name='GROUP_MULT', data=np.zeros(len(cat), dtype=np.int16)))
    cat.add_column(Column(name='GROUP_PRIMARY', data=np.zeros(len(cat), dtype=bool)))
    cat.add_column(Column(name='GROUP_RA', length=len(cat), dtype='f8')) # diameter-weighted center
    cat.add_column(Column(name='GROUP_DEC', length=len(cat), dtype='f8'))
    cat.add_column(Column(name='GROUP_DIAMETER', length=len(cat), dtype='f4'))

    #ww = np.where((parent['RA'] > 177) * (parent['RA'] < 178) * (parent['DEC'] > -1.5) * (parent['DEC'] < -0.5))[0]
    #ww = np.where((parent['RA'] > 200) * (parent['RA'] < 240) * (parent['DEC'] > 20))[0]
    #ww = np.where((parent['RA'] > 193) * (parent['RA'] < 196) * (parent['DEC'] > 26) * (parent['DEC'] < 30))[0]
    
    # Initialize a unique group number for each galaxy
    gnum = np.arange(len(cat)).astype(np.int)
    mgrp = np.ones(len(cat)).astype(np.int16)
    
    # First group galaxies within 10 arcmin, setting those to have the same
    # group number
    t0 = time.time()
    print('Spheregrouping took...', end='')
    ingroup, group_mult, firstgroup, nextgroup = spheregroup(cat['RA'], cat['DEC'], dmax)

    ngroup = np.count_nonzero(firstgroup != -1)
    for ii in np.arange(ngroup):
        #print(ii, ngroup)
        nn = group_mult[ii] # number of galaxies in this group
        if nn > 1:
            # Build INDX as the indices of all objects in this grouping
            indx = np.zeros(nn, dtype=int)
            indx[0] = firstgroup[ii]
            for jj in np.arange(nn-1):
                indx[jj+1] = nextgroup[indx[jj]]
            # Look at all pairs within this grouping to see if they should be connected.
            for jj in np.arange(nn-1):
                for kk in np.arange(jj, nn):
                    dd = degrees_between(cat['RA'][indx[jj]], cat['DEC'][indx[jj]], cat['RA'][indx[kk]], cat['DEC'][indx[kk]])
                    #dd = c1.separate(c2)
                    # If these two galaxies should be connected, make GNUM the
                    # same for them...
                    #print(dd, mfac * (cat['D25'][indx[jj]] / 60. + cat['D25'][indx[kk]] / 60.))
                    if dd < (0.5 * mfac * (cat['D25'][indx[jj]] / 60. + cat['D25'][indx[kk]] / 60.)):
                        jndx = np.where(np.logical_or(gnum[indx]==gnum[indx[jj]], gnum[indx]==gnum[indx[kk]]))[0]
                        gnum[indx[jndx]] = gnum[indx[jndx[0]]]
                        mgrp[indx[jndx]] = len(jndx)
            #print(ii, ngroup, gnum[indx], mgrp[indx])

    # Special-case the largest galaxies, looking for neighbhors
    ibig = np.where(cat['D25'] / 60. > dmax)[0]
    if len(ibig) > 0:
        for ii in np.arange(len(ibig)):
           dd = degrees_between(cat['RA'][ibig[ii]], cat['DEC'][ibig[ii]], cat['RA'], cat['DEC'])
           inear = np.where(dd < 0.5*(cat[ibig[ii]]['D25'] + cat['D25']) / 60.)[0]
           if len(inear) > 0:
               for jj in np.arange(len(inear)):
                  indx = np.where(np.logical_or(gnum==gnum[ibig[ii]], gnum==gnum[inear[jj]]))[0]
                  gnum[indx] = gnum[indx[0]]
                  mgrp[indx] = len(indx)
    print('...{:.3f} min'.format((time.time() - t0)/60))

    npergrp, _ = np.histogram(gnum, bins=len(gnum), range=(0, len(gnum)))

    print('Found {} total groups, including:'.format(len(set(gnum))))
    print('  {} groups with 1 member'.format(np.sum( (npergrp == 1) ).astype('int')))
    print('  {} groups with 2 members'.format(np.sum( (npergrp == 2) ).astype('int')))
    print('  {} group(s) with 3-5 members'.format(np.sum( (npergrp >= 3)*(npergrp <= 5) ).astype('int')))
    print('  {} group(s) with 6-10 members'.format(np.sum( (npergrp >= 6)*(npergrp <= 10) ).astype('int')))
    print('  {} group(s) with >10 members'.format(np.sum( (npergrp > 10) ).astype('int')))

    cat['GROUP_ID'] = gnum
    cat['GROUP_MULT'] = mgrp

    I = np.where(cat['GROUP_MULT'] == 1)[0]
    if len(I) > 0:
        cat['GROUP_RA'][I] = cat['RA'][I]
        cat['GROUP_DEC'][I] = cat['DEC'][I]
        cat['GROUP_DIAMETER'][I] = cat['D25'][I]
        cat['GROUP_NAME'][I] = cat['GALAXY'][I]
        cat['GROUP_PRIMARY'][I] = True

    more = np.where(cat['GROUP_MULT'] > 1)[0]
    for group in set(cat['GROUP_ID'][more]):
        I = np.where(cat['GROUP_ID'] == group)[0]
        # Compute the D25-weighted RA, Dec of the group:
        weight = cat[I]['D25']
        cat['GROUP_RA'][I] = np.sum(weight * cat[I]['RA']) / np.sum(weight)
        cat['GROUP_DEC'][I] = np.sum(weight * cat[I]['DEC']) / np.sum(weight)
        # Get the diameter of the group as the distance between the center of
        # the group and the outermost galaxy (plus the diameter of that galaxy,
        # in case it's a big one!).
        dd = degrees_between(cat['RA'][I], cat['DEC'][I], cat['GROUP_RA'][I[0]], cat['GROUP_DEC'][I[0]])
        pad = dd + cat['D25'][I] / 60.0
        cat['GROUP_DIAMETER'][I] = np.max(pad) * 60 # [arcmin]
        if cat['GROUP_DIAMETER'][I[0]] < np.max(cat['D25'][I]):
            print('Should not happen!')
            pdb.set_trace()

        # Assign the group name based on its largest member and also make this
        # galaxy "primary".
        primary = np.argmax(cat['D25'][I])
        cat['GROUP_NAME'][I] = '{}_GROUP'.format(cat['GALAXY'][I][primary])
        cat['GROUP_PRIMARY'][I[primary]] = True

    print('Building a group catalog took {:.3f} min'.format((time.time() - t0)/60))
        
    return cat





# read in vf main catalog
vfmain = Table.read('/home/rfinn/research/Virgo/tables-north/v0/vf_north_v0_main.fits')
# save a version that is RA, DEC, and radius

cat = vfmain['RA','DEC','radius','VFID']

# rename radius to D25
cat.rename_column('radius','D25')
cat.rename_column('VFID','GALAXY')
# run SGA-build-parent.build_group_catalog on catalog

gcat = build_group_catalog(cat)
