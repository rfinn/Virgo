#!/usr/bin/env python

'''
GOAL:
* download spreadsheet from google, save as excel
* read in excel sheets
* combine into one csv/fits table
* write out table
  - called virgo_check_sample_by_eye.csv

USAGE:
* first download latest version of virgo_check_by_eye from google spreadsheet, save as excel file
* move to ~/research/Virgo/google-tables

mv ~/Downloads/virgo_check_sample_by_eye.xlsx virgo_check_sample_by_eye.finished.xlsx


* run this program 

python ~/github/Virgo/programs/collate_check_by_eye_results.py

or from within ipython -pylab

run ~/github/Virgo/programs/collate_check_by_eye_results.py

'''

import pandas as pd
import os
from matplotlib import pyplot as plt
from astropy.table import Table
import numpy as np

homedir = os.getenv("HOME")
# download spreadsheet from google, save as excel
tablepath =  homedir+'/research/Virgo/google-tables/'
results_spreadsheet = tablepath+'virgo_check_sample_by_eye_v1.finished.xlsx'
# read in excel sheets

# combine into one csv file



class mydataframe:
    def __init__(self,spreadsheet):
        self.df = pd.concat(pd.read_excel(spreadsheet, sheet_name=None),ignore_index=True)
        self.print_statistics()
        self.print_statistics_north()
        #self.classification_hist()
        self.sheet_name = spreadsheet

    def print_statistics(self):
        print('number of objects with class=1 = ',sum(self.df['class'] == 1))
        flag = (self.df['class'] == 0) | (self.df['class'] == 2) | (self.df['class'] == 4)
        nremoved = sum(flag)
        print('number of objects to be removed (class=0, 2, 4) = ',sum(flag))
        class_individual = [0,1,2,3,4,5,6,7,8,9,16]
        for i in class_individual:
            print('number of objects with class',i,' = ',sum(self.df['class'] == i))
        print('percent of sample removed = %.1f'%(nremoved/len(self.df)*100))            
    def print_statistics_north(self):
        print('DEC > -1 galaxies only')
        dflag = self.df['DEC'] > -1
        print('number of objects with class=1 = ',sum((self.df['class'] == 1) & dflag))
        flag = ((self.df['class'] == 0) | (self.df['class'] == 2) | (self.df['class'] == 4) )  &dflag
        nremoved = sum(flag)
        print('number of objects to be removed (class=0, 2, 4) = ',sum(flag))
        class_individual = [0,1,2,3,4,5,6,7,8,9,16]
        for i in class_individual:
              print('number of objects with class',i,' = ',sum((self.df['class'] == i)& dflag))
        print('percent of sample removed = %.1f'%(nremoved/sum(dflag)*100))
    def print_id9_for_gianluca(self):
        flag = (self.df['class'] == 9)
        ids = self.df['galnumber'][flag]
        #hl = self.df['HL'][flag]
        print('problem with GL catalog')
        print(ids)
    def classification_hist(self):
        plt.figure()
        plt.hist(self.df['class'])
        plt.gca().set_yscale('log')
        plt.xlabel('Classification')
        plt.show()
        
    def write_table(self):
        t2 = Table.from_pandas(self.df)
        outfile = self.sheet_name.split('.')[0]+'.csv'
        t2.write(outfile,format='csv',overwrite=True)
        pass
    
if __name__ == '__main__':
    df = mydataframe(results_spreadsheet)
    df.write_table()
