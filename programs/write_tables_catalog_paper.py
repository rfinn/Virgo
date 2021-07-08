#!/usr/bin/env python
'''
GOAL:
* write out tables for Virgo catalog paper
* write out fits tables for full tables

OUTPUT:
* table1.tex
  - table 1 from 

* table2.tex 
  - table 2 from 

* fits versions of table1 and table2
  - .table1.DATE.fits
  - 

USAGE:

python writetables.py

UPDATES:

'''


import numpy as np
import os
import shutil
from astropy.io import fits, ascii
from astropy.table import Table
from datetime import datetime
from readtables import vtables

homedir = os.getenv('HOME')
tablepath = homedir+'/research/APPSS/tables/'
tablepath = homedir+'/research/Virgo/tables-north/v1/'
latextablepath = homedir+'/research/Virgo/papers/catalog_paper/'


class latextable(vtables):
    def clean_arrays(self):
        '''
        remove bogus values from SFR estimates and other arrays with null values

        '''
        pass
    def print_cat_table(self,nlines=10,filename=None,papertableflag=True,startindex=3000):
        '''write out latex version of table 1 '''
        if filename is None:
            fname=latextablepath+'cat_table.tex'
        else:
            fname = filename 
        outfile = open(fname,'w')
        print(fname)
        ### ADD FLAGS FOR CO, HI, ALFALFA
        outfile.write('\\begin{sidewaystable*}%[ptbh!]\n')
        outfile.write('\\begin{center}\n')
        outfile.write('\\scriptsize\n')
        outfile.write('\\setlength\\tabcolsep{3.0pt} \n')
        outfile.write('\\tablenum{4} \n')
        outfile.write('\\caption{Main Catalog with Cross IDs\label{tab:main}  } \n')
        outfile.write('\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n')
        outfile.write('\\hline \n')
        outfile.write('\\toprule \n')
        outfile.write('VFID   & RA &	DEC &	$v_{r}$ & $v_{\\rm cosmic}$ &  $v_{\\rm model}$  & HL~name & NSAID V0 & NSAID V1 & AGC Name & NED Name & CO  & HL & NSA & NSAV0  & A100  \\\\ \n')
        outfile.write('& (deg, J2000) & (deg, J2000) & $\\rm km~s^{-1}$ & $\\rm km~s^{-1}$ & $\\rm km~s^{-1}$ & & &  & & && & & &\\\\ \n')
        outfile.write('(1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12)& (13) & (14) & (15) &(16) \\\\ \n')
        outfile.write('\\midrule \n')
        outfile.write('\\hline \n')
        for i in range(startindex,nlines+startindex): # print first N lines of data
            # vfid ra dec z
            # replace ids with values of zero as nodata
            if self.main['AGC'][i] == 0:
                agc = '\\nodata'
            else:
                agc = self.main['AGC'][i] 
            if self.main['NSAIDV0'][i] == 0:
                nsaidv0 = '\\nodata'
            else:
                nsaidv0 = self.main['NSAIDV0'][i] 
            if self.main['NSAID'][i] == 0:
                nsaid = '\\nodata'
            else:
                nsaid = self.main['NSAID'][i] 
                
            format_s = '{0:s} & {1:9.6f} &{2:9.5f} & {3:.0f} & {4:.0f} & {5:.0f} &{6}& {7}& {8} & {9}& {10}& {11}&{12}&{13}&{14}&{15}\\\\ \n'
            s = format_s.format(self.main['VFID'][i],self.main['RA'][i],self.main['DEC'][i],\
                                self.main['vr'][i],self.env['Vcosmic'][i],self.env['Vmodel'][i],\
                                self.main['objname'][i],nsaidv0,nsaid,\
                                agc,self.main['NEDname'][i], self.main['COflag'][i],\
                                self.main['HLflag'][i],self.main['NSAflag'][i],self.main['NSAV0flag'][i],\
                                self.main['A100flag'][i])
            if papertableflag:
                # replace nans with \\nodata
                s=s.replace('nan','\\nodata')
                #s=s.replace('0','\\nodata')                

            outfile.write(s)

        outfile.write('\\bottomrule \n')
        outfile.write('\\hline \n')
        outfile.write('\\end{tabular} \n')
        outfile.write('\\end{center} \n')
        outfile.write('\\tablecomments{This table is published in its entirety in machine-readable format.  A portion is shown here for guidance regarding its form and content.}')        
        outfile.write('\\end{sidewaystable*} \n')
        outfile.close()
    def print_env_table(self,nlines=10,filename=None,papertableflag=True,startindex=3000):
        '''
        write out latex version of table 2

        changes from old verion:
        - remove Cluver stellar mass (2)
        - add error SFR22
        - add err SFRUV
        - add err SFRUVcorr

        net gain of one column
        '''
        if filename is None:
            fname=latextablepath+'env_table.tex'
        else:
            fname = filename 
        outfile = open(fname,'w')
        print(fname)
        # ADD 2d MEMBER TO FILAMENTS
        # PUT CLOSEST FILAMENT FIRST
        outfile.write('\\begin{sidewaystable*}%[ptbh!]\n')
        outfile.write('\\begin{center}\n')
        outfile.write('\\scriptsize\n')
        outfile.write('\\setlength\\tabcolsep{3.0pt} \n')
        outfile.write('\\tablenum{5} \n')
        outfile.write('\\caption{Environmental Properties of Catalog Galaxies\label{tab:environment}  } \n')
        outfile.write('\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n')
        outfile.write('\\hline \n')
        outfile.write('\\toprule \n')
        columns = ['VFID','SGX', 'SGY', 'SGZ',\
                   'n$_{5,2D}$', 'err(n$_{5,2D}$)' ,'n$_{5,3D}$','err(n$_{5,3D}$)' ,\
                   'Nearest Filament',\
                   '$\\rm D_{Filament}~2D$',\
                   '$\\rm D_{Filament}~3D$',\
                   'Filament Memb',\
                   'Group','Cluster','Pure Field']
        

        units = ['','$h^{-1}$~Mpc', '$h^{-1}$~Mpc', '$h^{-1}$~Mpc', \
                 '$h^{2}$~Mpc$^{-2}$', '$h^{2}$~Mpc$^{-2}$',\
                 '$h^{3}$~Mpc$^{-3}$', '$h^{3}$~Mpc$^{-3}$', '',\
                 '$h^{-1}$~Mpc','$h^{-1}$~Mpc',\
                 '',\
                 '','','']
        # build the latex for header and units
        outstring = ""
        for i in range(len(columns)):
            if i == len(columns)-1:
                outstring += "{} \\\\ \n"
            else:
                outstring += "{} &"

        # combine latex formatting with the column names 
        header = outstring.format(*columns)
        outfile.write(header)
        
        # combine latex formatting with the column names 
        latexunits = outstring.format(*units)
        outfile.write(latexunits)

        # column numbers
        outstring = ""
        for i in range(len(columns)):
            if i == len(columns)-1:
                outstring += "({}) \\\\ \n".format(i+1)
            else:
                outstring += "({}) &".format(i+1)
        outfile.write(outstring)
                                                   
        outfile.write('\\midrule \n')
        outfile.write('\\hline \n')
                                                   
        # build the data string 
        outstring = ""
        for i in range(len(columns)):
            if ((i > 0) & (i < 8)) | (i == 9) | (i == 10):
                outstring += "{:.1f} & "
            # don't add & after the last column
            # also add line end in latex (\\\\) and newline in text file (\n)
            elif i == len(columns)-1:
                outstring += "{} \\\\ \n"
            else:
                outstring += "{} &"

        for i in range(startindex,startindex+nlines): # print first N lines of data
            # vfid ra dec z
            # replace ids with values of zero as nodata

            #print(self.main['VFID'][i],self.env['Vcosmic'][i],\
            #                    self.fil['SGX'][i],self.fil['SGY'][i],self.fil['SGZ'][i],\
            #                    self.env['n5th_2D'][i],self.env['n5th_2D_err'][i],\
            #                    self.env['n5th'][i],self.env['n5th_err'][i],\
            #                    self.fil['filament'][i])

            #columns = ['VFID','SGX', 'SGY', 'SGZ',\
            #       'n5th\_2D', 'err' ,'n5th\_3D','err' ,\
            #       'dist\_3D','err', \
            #       'Filament', 'Filament Memb',\
            #       'Group','Cluster','Pure Field'
            #       ]

            s = outstring.format(self.main['VFID'][i],\
                                 self.fil['SGX'][i],self.fil['SGY'][i],self.fil['SGZ'][i],\
                                 self.env['n5th_2D'][i],self.env['n5th_2D_err'][i],\
                                 self.env['n5th'][i],self.env['n5th_err'][i],\
                                 self.fil['filament'][i].replace('_Filament','').replace('Virgo_','').replace('_','\_'),\
                                 self.fil['filament_dist_2D'][i],\
                                 self.fil['filament_dist_3D'][i],\
                                 int(self.fil['filament_member'][i]),\
                                 int(self.env['flag_gr'][i]),int(self.env['flag_clus'][i]),\
                                 int(self.env['flag_pf'][i]))
                                 
                                
            if papertableflag:
                # replace nans with \\nodata
                s=s.replace('nan','\\nodata')
                #s=s.replace('0','\\nodata')                

            outfile.write(s)

        outfile.write('\\bottomrule \n')
        outfile.write('\\hline \n')
        outfile.write('\\end{tabular} \n')
        outfile.write('\\end{center} \n')
        outfile.write('\\tablecomments{This table is published in its entirety in machine-readable format.  A portion is shown here for guidance regarding its form and content.}')        
        outfile.write('\\end{sidewaystable*} \n')
        outfile.close()

        
    def write_full_latex_tables(self):
        ''' write full verion of catalog and environment tables in latex format '''
        dateTimeObj = datetime.now()
        myDate = dateTimeObj.strftime("%d-%b-%Y")
        cat_tab_filename = latextablepath+'castignani2021_cat_table_full.'+myDate+'.tex'
        env_tab_filename = latextablepath+'castignani2021_env_table_full.'+myDate+'.tex'        
        print(cat_tab_filename)
        self.print_cat_table(nlines=len(self.main),filename = cat_tab_filename,startindex=0)




        self.print_env_table(nlines=len(self.main),filename = env_tab_filename,startindex=0)        

     
        pass

    def write_full_fits_tables(self):
        ''' write full version of catalog and environment tables in fits format '''
        dateTimeObj = datetime.now()
        myDate = dateTimeObj.strftime("%d-%b-%Y")
        cat_tab_filename = tablepath+'castignani2021_cat_table_full.'+myDate+'.fits'
        env_tab_filename = tablepath+'castignani2021_env_table_full.'+myDate+'.fits'        

        ####################
        ### CATALOG TABLE
        ####################
        print(cat_tab_filename)        
        colnames = ['VFID','RA','DEC',\
                    'v_r','v_cosmic','v_model ',\
                    'HL_name','NSAID_V0','NSAID_V1','AGC_Name','NED_Name',\
                    'CO','HL','NSA','NSAV0','A100']
        colunits = ['','deg','deg',\
                    'km/s','km/s','km/s',\
                    '','','','','',\
                    '','','','','']
        cat_tab = Table([self.main['VFID'],self.main['RA'],self.main['DEC'],\
                         self.main['vr'],self.env['Vcosmic'],self.env['Vmodel'],\
                         self.main['objname'],self.main['NSAIDV0'],self.main['NSAID'],\
                         self.main['AGC'],self.main['NEDname'], \
                         self.main['COflag'],self.main['HLflag'],\
                         self.main['NSAflag'],self.main['NSAV0flag'],\
                         self.main['A100flag']],names=colnames)
        cat_tab.write(cat_tab_filename,format='fits',overwrite=True)

        ######################
        ### ENVIRONMENT TABLE
        ######################        

        print(env_tab_filename)        
        colnames = ['VFID','SGX', 'SGY', 'SGZ',\
                   'n5_2D', 'err_n5_2D' ,'n5_3D','err_n5_3D',\
                   'Nearest_Filament',\
                   'D_Filament_2D',\
                   'D_Filament_3D',\
                   'Filament_Memb',\
                   'Group','Cluster','Pure_Field']
        filament = []
        for f in self.fil['filament']:
            filament.append(f.replace('_Filament','').replace('Virgo_',''))
        env_tab = Table([self.main['VFID'],\
                         self.fil['SGX'],self.fil['SGY'],self.fil['SGZ'],\
                         self.env['n5th_2D'],self.env['n5th_2D_err'],\
                         self.env['n5th'],self.env['n5th_err'],\
                         filament,\
                         self.fil['filament_dist_2D'],\
                         self.fil['filament_dist_3D'],\
                         self.fil['filament_member'],\
                         np.array(self.env['flag_gr'],'i'),\
                         np.array(self.env['flag_clus'],'i'),\
                         np.array(self.env['flag_pf'],'i')],\
                        names=colnames)
        env_tab.write(env_tab_filename,format='fits',overwrite=True)
        

        pass
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v1/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_north_v1_', help = 'prefix for tables; default is vf_north_v1')                               
    args = parser.parse_args()
    
    v = latextable(args.tabledir,args.tableprefix)
    v.read_all()

    
