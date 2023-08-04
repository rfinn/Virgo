#!/usr/bin/env python


"""
GOAL:
* subtract the continuum using color-correction from legacy images

USEAGE:

python subtract_continuum.py dirname

this assumes the file naming convention used in the Virgo Filament Survey, with 

dirname-R.fits
dirname-Ha.fits

subdirectory called
legacy/

that has g and r images that have been reprojected onto the halpha image pixel scale


PROCEDURE:
* create a g-r image
   2.5 log10(flux_r/flux_g)

* reproject onto halpha image

* mask r-band image

* calculate rms in r-band outside the masked regions

* apply color correction to continuum-subtracted images
  

FILTER TRANSFORMATIONS FROM MATTEO

---------------------------------------------------------------------
Number of stars meeting the color and brightness cuts: 54594
---------------------------------------------------------------------
Best fit linear    Ha4 - KPHr = -0.0465 * (PS1_g-PS1_r) + 0.0012
Best fit quadratic Ha4 - KPHr = -0.0156*(PS1_g-PS1_r)^2 + -0.0262*(PS1_g-PS1_r) + -0.0048
---------------------------------------------------------------------
Initial bias/scatter:             -0.0265/0.0158
Corr linear fit              :    0.0000/0.0130
Corr quadratic fit           :    0.0000/0.0130
---------------------------------------------------------------------

---------------------------------------------------------------------
Number of stars meeting the color and brightness cuts: 54692
---------------------------------------------------------------------
Best fit linear    Ha4 - KPSr = -0.1804 * (PS1_g-PS1_r) + 0.0158
Best fit quadratic Ha4 - KPSr = -0.0066*(PS1_g-PS1_r)^2 + -0.1718*(PS1_g-PS1_r) + 0.0133
---------------------------------------------------------------------
Initial bias/scatter:             -0.0916/0.0375
Corr linear fit              :    0.0000/0.0138
Corr quadratic fit           :    0.0000/0.0138
---------------------------------------------------------------------

---------------------------------------------------------------------
Number of stars meeting the color and brightness cuts: 54712
---------------------------------------------------------------------
Best fit linear    Intha - INTSr = -0.2334 * (PS1_g-PS1_r) + 0.0711
Best fit quadratic Intha - INTSr = 0.0160*(PS1_g-PS1_r)^2 + -0.2541*(PS1_g-PS1_r) + 0.0772
---------------------------------------------------------------------
Initial bias/scatter:             -0.0682/0.0491
Corr linear fit              :    0.0000/0.0193
Corr quadratic fit           :    0.0000/0.0192
——————————————————————————————————
"""
import sys
import os
from astropy.io import fits
from astropy import stats, convolution
from astropy import wcs

import glob
from reproject import reproject_interp

# Halpha filter width in angstrom
filter_width_AA = {'BOK':80.48,'HDI':80.48,'INT':95,'MOS':80.48,'INT6657':80}

filter_lambda_c_AA = {'BOK':6620.52,'HDI':6620.52,'INT':6568,'MOS':6620.52,'INT6657':6657}
52

def filter_transformation(telescope,rfilter, gr_col):
    if telescope == 'BOK' | ((telescope == 'HDI' and (rfilter == 'r'):
        #Ha4_KPSr = -0.1804 * (gr_col) + 0.0158
        ha_r = -0.1804 * (gr_col) + 0.0158
    elif telescope == 'INT':
        #Intha_INTSr = -0.2334 * (gr_col) + 0.0711
        ha_r = -0.2334 * (gr_col) + 0.0711
    elif telescope == 'MOS' | ((telescope == 'HDI') and (rfilter == 'R')):
        #Ha4_KPHr = -0.0465 * (gr_col) + 0.0012
        ha_r = -0.0465 * (gr_col) + 0.0012
            
def get_gr(gfile,rfile,mask=None):
    
    """ take g and r filenames, return g-r data and save g-r color image """
    g = fits.open(gfile)
    r = fits.open(rfile)
    data_g = g[0].data
    data_r = r[0].data
    g.close()
    r.close()
    
    ###
    # the following is from Matteo Fossati
    ###
    
    # get noise in the image    
    #stat is a tuple of mean, median, sigma
    print('Computing median values for g and r images')
    stat_r = stats.sigma_clipped_stats(data_r,mask=mask)
    print('Subtracting {0:3.2f} from r-band image'.format(stat_r[1]))
    data_r -= stat_r[1]
    stat_g = stats.sigma_clipped_stats(data_g,mask=mask)
    print('Subtracting {0:3.2f} from g-band image'.format(stat_g[1]))
    data_g -= stat_g[1]

    # create a mask, where SNR > 10    
    usemask = (data_r>10*stat_r[2])    
    gr_col = -2.5*np.log10(data_g/data_r)
    print('Smoothing images for color calculation')
    gr_col = convolution.convolve_fft(gr_col, convolution.Box2DKernel(20), allow_huge=True, nan_treatment='interpolate')

    gr_col[np.logical_not(usemask)] = np.nan
    # save gr color image
    hdu = fits.PrimaryHDU(gr_col, header=r[0].header)
    outimage = rfile.replace('r.fits','gr.fits')
    hdu.writeto(outimage, overwrite=True)
    hdu.close()
    return gr_col

def subtract_continuum(Rfile, Hfile, gfile, rfile, mask=None,overwrite=False):
    """
    reproject infile to reffile image

    PARAMS:
    Rfile : r-band image taken with halpha, to be used for continuum
    Hfile : halpha image filename
    gfile : g-band filename to be used for calculating g-r color (legacy image)
    rfile : r-band filename to be used for calculating g-r color (legacy image)
    outname: output name 


    RETURN:

    
    """
    outname = Hfile.replace('Ha.fits','CS-gr.fits')
    if os.path.exists(outname) & (not overwrite):
        print("continuum-subtracted image exists - not redoing it")
        return

    outimage = rfile.replace('r.fits','gr.fits')
    if os.path.exists(outimage):
        print("found g-r image.  not remaking this")
        hdu = fits.open(outimage)
        gr_col = hdu[0].data
        hdu.close()
    else:
        gr_col = get_gr(gfile,rfile,mask=mask)
    usemask = gr_col == np.nan
    fileroot = Rfile.replace('-R.fits','')        
    hhdu = fits.open(Hfile)
    rhdu = fits.open(Rfile)
    rZP = rhdu[0].header['PHOTZP']
    hZP = hhdu[0].header['PHOTZP']
    # get filter names
    rfilter = rhdu[0].header['FILTER']
    hfilter = hhdu[0].header['FILTER']    

    tels = ['BOK','INT','HDI','MOS']
    for t in tels:
        if t in dirname:
            telescope = t
            print('telescope = ',t,dirname)
            break
    ##
    # The following is from Matteo Fossati
    ##

    print('Generate NET image')

    # DONE: TODO - subtract sky from r-band image
    # DONE: TODO - subtract sky from Halpha image
    print('Computing median values for r and halpha images')
    stat_r = stats.sigma_clipped_stats(rhdu[0].data,mask=mask)
    print('Subtracting {0:3.2f} from r-band image'.format(stat_r[1]))
    data_r = rhd[0].data - stat_r[1]
    stat_h = stats.sigma_clipped_stats(hhdu[0].data,mask=mask)
    print('Subtracting {0:3.2f} from halpha image'.format(stat_h[1]))
    data_h = hhdu[0].data - stat_g[1]

    # Generate the r band mag image and the the r band calibrated to Halpha wave
    # This works only for positive flux pixels. Take this into account
    # Best fit with offset    CFHT_Ha = CFHT_r + -0.1713+/-0.000017 * (CFHT_gs-CFHT_r) + 0.0717+/-0.003070

    # DONE: TODO - change ZP - get this from image header
    mag_r = -2.5*np.log10(data_r) + rZP

    # DONE: TODO - change the color conversion - need to use different conversion for each Halpha/r combo
    mag_r_to_Ha = mag_r + filter_transformation(telescope,rfilter, gr_col)

    #Back to calibrated flux units
    data_r_to_Ha = np.copy(data_r)
    data_r_to_Ha = convolution.convolve_fft(data_r, convolution.Box2DKernel(5), allow_huge=True, nan_treatment='interpolate')
    # DONE: TODO - change ZP from 30 to value in image header
    data_r_to_Ha[usemask] = 10**(-0.4*(mag_r_to_Ha[usemask]-rZP))

    #Go to cgs units
    fnu_NB  = 3.631E3*data_NB*1E-12
    # DONE: TODO - change the filter EWs - need a dictionary for each Halpha filter
    flam_NB = 2.99792458E-5*fnu_NB/(filter_lambda_c_AA(telescope)**2) *1E18

    cnu_NB  = 3.631E3*data_r_to_Ha*1E-12
    clam_NB = 2.99792458E-5*cnu_NB/(filter_lambda_c_AA**2) *1E18 # change central wavelength

    # DONE: TODO - change width of the filter
    flam_net = filter_width_AA(telescope)*(flam_NB-clam_NB) #106 is the width of the filter

    # TODO - change output image name
    hdu = fits.PrimaryHDU(flam_NB, header=hhdu[0].header)
    hdu.writeto(fileroot+'_net_new.fits', overwrite=True) #NB image in F_lambda units, before
    # TODO - change output image name
    hdu = fits.PrimaryHDU(clam_NB, header=hhdu[0].header)
    hdu.writeto(fileroot+'_cont_new.fits', overwrite=True)

    #Calculate clipped statistic
    #stat is a tuple of mean, median, sigma
    stat = stats.sigma_clipped_stats(flam_net,mask=mask)
    flam_net -= stat[1]

    print('Unbinned SB limit 1sigma {0:3.1f} e-18'.format(stat[2]/(pscale_NB[0]**2)))

    # DONE: TODO - change output image name
    hdu = fits.PrimaryHDU(flam_net, header=hhdu[0].header)
    hdu.writeto(fileroot+'_net_flux.fits', overwrite=True)

    # DONE: TODO - change pixel scale
    sblam_net = flam_net/(pscale_NB[0]**2)
    
    # DONE: TODO - change output image name
    hdu = fits.PrimaryHDU(sblam_net, header=hhdu[0].header)
    hdu.writeto(fileroot+'_net_sb.fits', overwrite=True)


    print('Smoothing net image')

    flam_net_smooth = convolution.convolve_fft(flam_net, convolution.Box2DKernel(15), allow_huge=True, nan_treatment='interpolate')

    hdu = fits.PrimaryHDU(flam_net_smooth, header=hhdu[0].header)
    # DONE: TODO - change output image name
    hdu.writeto(fileroot+'_net_smooth.fits', overwrite=True)

    stat_sm = stats.sigma_clipped_stats(flam_net_smooth,mask=mask)
    # TODO - add this to image header

    print('Smoothed {1}x{1} SB limit 1sigma {0:3.1f} e-18'.format(stat_sm[2]/(pscale_NB[0]**2), 15))

    # close hdu files
    hhdu.close()
    rhdu.close()

 

if __name__ == '__main__':
    dirname = sys.argv[1]

    t = dirname.split('-')
    prefix = t[0]+'-'+t[1]
    # define the file names
    Rfile = os.path.join(dirname,dirname+'-R.fits') # r-band image taken with same telescope as halpha
    Hfile = os.path.join(dirname,dirname+'-Ha.fits')  # halpha image

    rfiles = glob.glob(os.path.join(dirname,'legacy',prefix+'*r-ha.fits'))
    if len(rfiles) < 1:
        print("problem getting r-ha.fits legacy image")
        sys.exit()
    else:
        rfile = rfiles[0] # legacy r-band image
        
    # legacy g-band image, shifted to match halpha footprint and pixel scale
    gfiles = glob.glob(os.path.join(dirname,'legacy',prefix+'*g-ha.fits'))
    if len(gfiles) < 1:
        print("problem getting g-ha.fits legacy image")
        sys.exit()
    else:
        gfile = gfiles[0] # legacy r-band image

    maskfile = Rfile.replace('-R.fits','-R-mask.fits')
    if not os.path.exists(maskfile):
        print("WARNING: no mask found")
        mask = None
    else:
        mask = fits.getdata(maskfile)
        mask = mask > 0
    subtract_continuum(Rfile, Hfile, gfile, rfile,mask=None,overwrite=False)
