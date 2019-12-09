#!/usr/bin/env python

'''
GOAL:
   - to run theli commands from the command line rather than gui
   - reduce WFC camera data for Virgo filaments project

PROCEDURE:
   - execute scripts to reduce data
   - combine bias and flats
   - calibration science data
   - run astrometry
   - create mosaic

WHAT NEXT:
   - after this, should run photometric calibration and write ZP to header

COMMANDS:

# move images to subdirecotries
/usr/local/anaconda2/bin/python ~/github/Virgo/programs/sort_int.py

# split images into 4 images, one for each WFC chip
./process_split_WFC@INT.sh /home/rfinn/INT/temp/20190207 BIAS
./process_split_WFC@INT.sh /home/rfinn/INT/temp/20190207 SKYFLAT-r
./process_split_WFC@INT.sh /home/rfinn/INT/temp/20190207 pointing021-r

# process bias
./parallel_manager.sh process_bias_para.sh /home/rfinn/INT/temp/20190207/ BIAS

# process flats
./parallel_manager.sh process_flat_para.sh /home/rfinn/INT/temp/20190207 BIAS SKYFLAT-r
./create_flat_ratio.sh /home/rfinn/INT/20190207 SKYFLAT-r
./parallel_manager.sh create_norm_para.sh /home/rfinn/INT/20190207 SKYFLAT-r

# calibrate data
./parallel_manager.sh process_science_para.sh /home/rfinn/INT/20190207 BIAS SKYFLAT-r pointing021-r

# create binned preview and global weights
./make_album_WFC@INT.sh /home/rfinnn/INT/temp/20190207 pointing-r OFC
./create_tiff.sh /home/rfinnn/INT/temp/20190207 pointing-r OFC
./parallel_manager.sh create_global_weights_para.sh /home/rfinnn/INT/temp/20190207 SKYFLAT-r_norm pointing-r

# create weights
./transform_ds9_reg.sh /home/rfinnn/INT/temp/20190207 pointing-r 
./parallel_manager.sh create_weights_para.sh /home/rfinnn/INT/temp/20190207 pointing-r OFC


# get source cat


# create catalogs from images
./parallel_manager.sh create_astromcats_para.sh /home/rfinn/INT/temp/20190207 pointing021-r OFC
./create_scampcats.sh /home/rfinn/INT/temp/20190207 pointing021-r OFC

# run scamp
./create_scamp.sh /home/obsastro2/20190207 pointing021-r OFC
./create_stats_table.sh /home/obsastro2/20190207 pointing021-r OFC headers
./create_absphotom_coadd.sh /home/obsastro2/20190207 pointing021-r

# subtract median of sky
./create_skysubconst_clean.sh /home/obsastro2/20190207 pointing021-r
./parallel_manager.sh create_skysubconst_para.sh /home/obsastro2/20190207 pointing021-r OFC ALL

# coaddition
./prepare_coadd_swarp.sh /home/obsastro2/20190207 pointing021-r OFC.sub
./parallel_manager.sh resample_coadd_swarp_para.sh /home/obsastro2/20190207 pointing021-r OFC.sub
./perform_coadd_swarp.sh /home/obsastro2/20190207 pointing021-r
./update_coadd_header.sh /home/obsastro2/20190207 pointing021-r OFC

'''
