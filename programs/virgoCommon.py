import os


#Get current path so program can tell if this is being run on Kelly's or Rose's computer
mypath=os.getcwd()
if mypath.find('rfinn') > -1:
    print("Running on Rose's computer")
    #agcfile='/Users/rfinn/idl/programs/idl_alfa/agctotal.sav'
    gitpath='/Users/rfinn/github/'
    nsapath = '/Users/rfinn/research/NSA/'
    gswlpath = '/Users/rfinn/Dropbox/Research/GSWLC/'
    agcpath = '/Users/rfinn/research/AGC/'
    outfile_prefix = '/Users/rfinn/Dropbox/Research/Virgo/finding-charts/'
    tablepath = '/Users/rfinn/github/Virgo/tables/'

elif mypath.find('kelly') > -1:
    print("Running on Kellys's computer")
    gitpath='/Users/kellywhalen/Github/'
    nsapath = '/Users/kellywhalen/RESEARCH/NSA_table/'
    gswlpath = '/Users/kellywhalen/RESEARCH/GSWLC/'
    agcpath = '/Users/kellywhalen/RESEARCH/AGC/'

elif mypath.find('grudnick') > -1:
    print("Running on Greg's computer")
    gitpath='/Users/grudnick/Work/Virgo_outskirts/Rfinn_github/'
    ###dummy variables as I don't have the full NSA, GSWL, or AGC catalogs
    nsapath = '/Users/grudnick/Temp/'
    gswlpath = '/Users/grudnick/Temp/'
    agcpath = '/Users/grudnick/Temp/'
    outfile_prefix = '/Users/grudnick/Dropbox/Virgo_filaments/finding-charts/'
    tablepath = '/Users/grudnick/Work/Virgo_outskirts/Rfinn_github/Virgo/tables/'

    
tablepath = gitpath+'Virgo/tables/'
galex_file = gitpath+'Virgo/tables/GALEX-WISE-allsky_virgo.fits'
nsa_file = gitpath+'Virgo/tables/VirgoCatalog.fits'
mass_file = gitpath+'Virgo/tables/StellarMasstoNSA_virgo.fits'
co_targets = gitpath+'Virgo/tables/CO-HI_virgo.fits'
wise_file = gitpath+'Virgo/tables/WISE_virgo.fits'
full_nsa = nsapath + 'nsa_v0_1_2.fits'

galex_file = gitpath+'Virgo/tables/GALEX-WISE-allsky_virgo.fits'
nsa_file = gitpath+'Virgo/tables/nsa.virgo.fits'
mass_file = gitpath+'Virgo/tables/nsa_mstar.virgo.fits'
co_targets = gitpath+'Virgo/tables/nsa_CO-Gianluca.virgo.fits'
wise_file = gitpath+'Virgo/tables/nsa_wise.virgo.fits'
halpha_file = gitpath+'Virgo/tables/nsa_Halpha.virgo.fits'
full_nsa = nsapath + 'nsa_v0_1_2.fits'

