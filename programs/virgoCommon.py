import os
from astropy.table import Table,join,Column
import numpy as np

#Get current path so program can tell if this is being run on Kelly's or Rose's computer
mypath=os.getcwd()
if mypath.find('rfinn') > -1:
    homedir = os.getenv("HOME")
    print("Running on Rose's computer")
    #agcfile='/Users/rfinn/idl/programs/idl_alfa/agctotal.sav'
    gitpath=homedir+'/github/'
    nsapath = homedir+'/research/NSA/'
    gswlpath = homedir+'/research/GSWLC/'
    agcpath = homedir+'/research/AGC/'
    outfile_directory = homedir+'/research/Virgo/finding-charts/'
    tablepath = homedir+'/github/Virgo/tables/'

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
    outfile_directory = '/Users/grudnick/Dropbox/Virgo_filaments/finding-charts/'
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


def myjoinleft(table1,table2,keys=None):
    table1 = Table(table1)
    table2 = Table(table2)
    table1.add_column(Column(np.arange(len(table1))),name='originalOrder')
    temp = join(table1,table2,keys=keys,join_type='left')
    table1.remove_column('originalOrder')
    temp = temp[np.argsort(temp['originalOrder'])]                           
    temp.remove_column('originalOrder')    
    # join returns a table that is sorted by the key columns
    # the following command gets table back into its original order
    return temp

