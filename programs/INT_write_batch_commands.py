#!/usr/bin/env

'''

just doing this for astrometry/coaddition

can cut and paste this into theli terminal

this doesn't work.  scamp gets stuck.  not sure what it's doing but it ran for hours and didn't even finish one pointing.  ugh.

'''

import os

DATE = '20190210'
OUTFILE = 'theli-astro-coadd-commands'
output = open(OUTFILE,'w')

cdir = os.getcwd()
DATE = os.path.basename(cdir)

t = os.listdir()
t.sort(reverse=True) # group pointings together, reverse puts r-band before Halpha

ALLCOMMANDS = False
COADD_ONLY = True
FILTER = '-Ha'
if ALLCOMMANDS:
    for f in t:
        if ((f.find('pointing') > -1) | (f.find('VF') > -1) | (f.find('NGC') > -1) | (f.find('UGC')>-1)) & (os.path.isdir(f)):

            p,filter = f.split('-')
            output.write('./create_astrorefcat_fromWEB.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC GAIA-DR1 vizier.cfa.harvard.edu \n')
            output.write('./parallel_manager.sh create_astromcats_para.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            output.write('./create_scampcats.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            output.write('./create_scamp.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            output.write('./create_stats_table.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC headers \n')
            output.write('./create_absphotom_coadd.sh /home/rfinn/data/INT/'+DATE+'/ '+f+'\n')
            output.write('./parallel_manager.sh create_skysub_para.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            output.write('./prepare_coadd_swarp.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC.sub \n')
            output.write('./parallel_manager.sh resample_coadd_swarp_para.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC.sub \n')
            output.write('./perform_coadd_swarp.sh /home/rfinn/data/INT/'+DATE+'/ '+f+'\n')
            output.write('./update_coadd_header.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')

    output.close()

if COADD_ONLY:
    for f in t:
        
        if (((f.find('pointing') > -1) | (f.find('VF') > -1) | (f.find('NGC') > -1) | (f.find('UGC')>-1)) & os.path.isdir(f)  & (f.find(FILTER) > -1)):

            p,filter = f.split('-')
            #output.write('./create_astrorefcat_fromWEB.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC GAIA-DR1 vizier.cfa.harvard.edu \n')
            
            #output.write('./parallel_manager.sh create_astromcats_para.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            #output.write('./create_scampcats.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            #output.write('./create_scamp.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            #output.write('./create_stats_table.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC headers \n')
            #output.write('./create_absphotom_coadd.sh /home/rfinn/data/INT/'+DATE+'/ '+f+'\n')
            #output.write('./parallel_manager.sh create_skysub_para.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')
            output.write('./prepare_coadd_swarp.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC.sub \n')
            output.write('./parallel_manager.sh resample_coadd_swarp_para.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC.sub \n')
            output.write('./perform_coadd_swarp.sh /home/rfinn/data/INT/'+DATE+'/ '+f+'\n')
            output.write('./update_coadd_header.sh /home/rfinn/data/INT/'+DATE+'/ '+f+' OFC \n')

    output.close()
