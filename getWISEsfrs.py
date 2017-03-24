from astropy.io import fits
import sys
sys.path.append('./CharyElbaz/')
import chary_elbaz_24um as chary

catpath = '/Users/rfinn/github/Virgo/tables/'



wise = fits.getdata(catpath+'WISE_virgo.fits')
nsa = fits.getdata(catpath+'VirgoCatalog.fits')

# galaxies with reliable W4 flux
wflag =  (wise.W4MPRO > 0.1) & (wise.W4SNR > 2.)

# set up arrays to hold flux, sfr, Lir
w4_flux_Jy = np.zeros(len(wise),'f')
ce_sfr = np.zeros(len(wise),'f')
ce_lir = np.zeros(len(wise),'f')

# zeropoint of magnitude scale from
# http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#Summary
F0 = 8.363

# calculate flux in Jy for galaxies with W4 detections
w4_flux_Jy[w4_flux_Jy] = F0*10**(-1*wise.W4MPRO[w4_flux+Jy]/2.5)

# calculate chary & elbaz sfr and Lir

ce_lir[wflag],ce_sfr[wflag]=chary.chary_elbaz_24um(nsa.Z[wflag],w4_flux_Jy[wflag]*1.e6)

