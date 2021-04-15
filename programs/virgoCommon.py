import os
from astropy.table import Table,join,Column
import numpy as np
import glob
from matplotlib import pyplot as plt
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
from matplotlib.patches import Rectangle

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


# super galactic coordinate ranges of filaments
'''
filaments = ['LeoII_A_Filament','LeoII_B_Filament','Leo_Minor_Filament',\
             'NGC5353_4_Filament','Filament_3','Filament_3b',\
             'VirgoIII_Filament','VirgoIII_Filament_Extension','W-M_Sheet',\
             'Virgo_Northern_Filament','Ursa_Major_Cloud']


SGX = [[0.43   , 9.18],   
       [2.42   ,13.35],    
       [0.59   , 5.64],    
       [-12.17 , 8.88],      
       [2.15   , 6.04],    
       [5.84   ,10.83],    
       [-10.94 ,-5.51],     
       [-5.48  ,-0.99],     
       [-9.28  ,-3.41],    
       [13.98  ,18.77],    
       [0.50   , 8.74]]


SGY = [[12.25 ,14.47],
       [13.65 ,13.97],
       [4.06  ,6.30 ],
       [25.40 ,27.97],
       [12.89 ,38.31],
       [18.28 ,21.53],
       [11.93 ,17.32],
       [11.02 ,17.01],
       [20.30 ,23.27],
       [16.71 ,21.50],
       [2.67  ,13.95]]

SGZ = [[-14.26 ,-6.73],
       [-8.37  ,-4.30],
       [-2.83  ,-1.77],
       [0.28   , 9.56],
       [-4.92  ,-1.89],
       [-7.01  ,-5.41],
       [3.40   ,11.09],
       [9.47   ,33.14],
       [-2.67  ,-1.91],
       [14.90  ,22.51],
       [0.12   , 1.37]]
'''
# from GL's CO paper
filaments= ['LeoII_A_Filament','LeoII_B_Filament','Leo_Minor_Filament',\
            'VirgoIII_Filament','NGC5353_4_Filament', 'Virgo_Draco_Filament',\
            'Virgo_Serpens_Filament', 'Virgo_Coma_Berenices_Filament',\
            'Leo_Minor_B_Filament','Ursa_Major_Cloud',\
            'W-M_Sheet','Canes_Venatici_Filament']

#names = [‘Leo Minor’, ‘Ursa Major Cloud’, ‘Canes Venatici’,  ‘LeoII B’, ‘LeoII A’, ‘Virgo III’, ‘Leo Minor B’, ‘W-M Sheet’, ‘Bootes’, ‘Coma Berenices’, ‘NGC5353/4 ’, ‘Serpens’, ‘Draco’, ] 

SGX =[ [-1,12],
       [0,16],
       [-0.5,7],
       [-11.5,-3],
       [-15,15],
       [13,20.5],
       [-9,1],
       [0,16],
       [4.5,12],
       [-1,10],
       [-14,-1],
       [0,5]]


SGY = [[8,16],
       [10,15.5],
       [3,10],
       [8,19],
       [16,27],
       [14,26],
       [8,20],
       [7.5,27.5],
       [17,23],
       [1,16],
       [15,26],
       [6,14]]

SGZ = [[-16,-3],
       [-11,-2],
       [-6,0],
       [2,13],
       [-2,15],
       [14,24],
       [10,35],
       [-6,0],
       [-8,-4.5],
       [-1,2.5],
       [-4,0],
       [1,5]]

filaments.append('Virgo_Bootes_Filament')
SGX.append([5,15])
SGZ.append([0,10])
SGY.append([14,18])


SGXrange = dict((a,b) for a,b in zip(filaments,SGX))
SGYrange = dict((a,b) for a,b in zip(filaments,SGY))
SGZrange = dict((a,b) for a,b in zip(filaments,SGZ))

# from paper
fil_lengths = {'LeoII_A_Filament':13.93,\
               'LeoII_B_Filament':12.67,\
               'Leo_Minor_Filament':7.63,\
               'VirgoIII_Filament':11.72,\
               'NGC5353_4_Filament':24.01,\
               'Virgo_Draco_Filament':12.31,\
               'Virgo_Serpens_Filament':25.63,\
               'Virgo_Coma_Berenices_Filament':26.27,\
               'Leo_Minor_B_Filament':7.95,\
               'Ursa_Major_Cloud':15.44,\
               'W-M_Sheet':8.45,\
               'Canes_Venatici_Filament':10.53,\
               'Virgo_Bootes_Filament':11.83}# just making up a number


def plot_spines():
    sfiles = glob.glob(homedir+'/research/Virgo/tables-north/spines/filament*.fits')
    ncolor = 0
    for i,f in enumerate(sfiles):
        spine  = Table.read(f)
        plt.plot(spine['ra'],spine['dec'],c=mycolors[ncolor],label=os.path.basename(f).replace('filament_spine_','').replace('.fits','').replace('_Filament',''),lw=3)
        ncolor += 1
        if ncolor > len(mycolors)-1:
            ncolor = 0
