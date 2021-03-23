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
  - Durbala2020.table2.DATE.fits

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
    def print_table4(self,nlines=10,filename=None,papertableflag=True,startindex=3000):
        '''write out latex version of table 1 '''
        if filename is None:
            fname=latextablepath+'table4.tex'
        else:
            fname = filename 
        outfile = open(fname,'w')
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
        outfile.write('VFID   & RA &	DEC &	V$_{helio}$ & V$_{cosmic}$ &  V$_{model}$  & HL~name & NSAID V0 & NSAID V1 & AGC Name & NED Name & CO  & HL & NSA & NSAV0  & A100  \\\\ \n')
        outfile.write('& J2000 & J2000 & $\\rm km~s^{-1}$ & $\\rm km~s^{-1}$ & $\\rm km~s^{-1}$ & & &  & & && & & &\\\\ \n')
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
    def print_table5(self,nlines=10,filename=None,papertableflag=True,startindex=3000):
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
            fname=latextablepath+'table5.tex'
        else:
            fname = filename 
        outfile = open(fname,'w')

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
                   'n5th\_2D', 'err' ,'n5th\_3D','err' ,\
                   'Nearest Filament',\
                   '$\\rm D_{Filament}~2D$',\
                   '$\\rm D_{Filament}~3D$',\
                   'Filament Memb',\
                   'Group','Cluster','Pure Field']
        

        units = ['','Mpc', 'Mpc', 'Mpc', \
                 '', '', '', '', '',\
                 'Mpc','Mpc',\
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
                                 self.fil['filament'][i].replace('_','\_'),\
                                 self.fil['filament_dist_2D'][i],\
                                 self.fil['filament_dist_3D'][i],\
                                 self.fil['filament_member'][i],\
                                 int(self.env['flag_gro'][i]),int(self.env['flag_clus'][i]),\
                                 int(self.env['flag_fie'][i]))
                                 
                                
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

        
    def write_full_tables(self):
        # table to assist with UAT summer research
        # will add columns that should be useful for lots of projects
        #
        # basically all of table 2:
        # a100
        # GSWLC mass, sfr
        # logMstarTaylor, err
        # logMstarMcGaugh, err
        # logSFR22, err
        # logSFRNUV, err
        # log SFR_NUVIR, err
        # logMHI, err
        # g,i
        #
        # plus RA, DEC, vhelio, D(err)
        dateTimeObj = datetime.now()
        myDate = dateTimeObj.strftime("%d-%b-%Y")

        tab1 = Table([self.tab['AGC'],self.tab['sdssPhotFlag'],self.tab['objID_1'],\
                      self.tab['RAdeg_Use'],self.tab['DECdeg_Use'],self.tab['Vhelio'],\
                      self.tab['Dist'],self.tab['sigDist'],\
                      self.tab['extinction_g'],self.tab['extinction_i'],\
                      self.tab['expAB_r'],self.expAB_r_err,\
                      self.tab['cModelMag_i'],self.tab['cModelMagErr_i']],
                     names=['AGC','sdssPhotFlag','sdss_objid',\
                            'RA','DEC','Vhelio',\
                            'Dist','sigDist',\
                            'extinction_g','extinction_i','expAB_r','expAB_r_err',\
                            'cModelMag_i','cModelMagErr_i'])
        tab1.write(tablepath+'durbala2020-table1.'+myDate+'.fits',format='fits',overwrite=True)


        tab2 = Table([self.tab['AGC'],self.tab['gamma_g'],self.tab['gamma_i'],self.tab['absMag_i_corr'],\
                      self.absMag_i_corr_err,self.tab['gmi_corr'],self.gmi_corr_err,\
                      self.tab['logMstarTaylor'],self.logMstarTaylor_err,\
                      self.tab['logMstarMcGaugh'],self.logMstarMcGaugh_err,\
                      self.gsw_mstar,self.gsw_mstar_err,\
                      self.sfr22,self.sfr22_err, self.sfrnuvir,self.sfrnuvir_err,\
                      self.gsw_sfr,self.gsw_sfr_err,\
                      self.tab['logMH'],self.tab['siglogMH']], \
                     names=['AGC','gamma_g','gamma_i','absMag_i_corr',\
                            'absMag_i_corr_err','gmi_corr','gmi_corr_err',\
                            'logMstarTaylor','logMstarTaylor_err',\
                            'logMstarMcGaugh','logMstarMcGaugh_err',\
                            'logMstarGSWLC','logMstarGSWLC_err',\
                            'logSFR22','logSFR22_err', 'logSFRNUVIR','logSFRNUVIR_err',\
                            'logSFRGSWLC','logSFRGSWLC_err',\
                            'logMH','logMH_err'])

        tab2.write(tablepath+'durbala2020-table2.'+myDate+'.fits',format='fits',overwrite=True)

        # machine readable tables for AAS journal
        #tab1.write(latextablepath+'durbala2020-table1.'+myDate+'.csv',format='csv',overwrite=True)        
        #tab2.write(latextablepath+'durbala2020-table2.'+myDate+'.csv',format='csv',overwrite=True)
        #ascii.write(tab1,latextablepath+'durbala2020-table1.'+myDate+'.txt',format='cds',overwrite=True)        
        #ascii.write(tab2,latextablepath+'durbala2020-table2.'+myDate+'.txt',format='cds',overwrite=True)

        # write out full latex tables for AAS journal
        self.print_table1(nlines=len(self.tab),filename=latextablepath+'table1_long.'+myDate+'.tex',papertableflag=False)
        self.print_table2(nlines=len(self.tab),filename=latextablepath+'table2_long.'+myDate+'.tex',papertableflag=False)        
        shutil.copy(latextablepath+'table1_long.'+myDate+'.tex',latextablepath+'table1_long.tex')
        shutil.copy(latextablepath+'table2_long.'+myDate+'.tex',latextablepath+'table2_long.tex')        
        pass
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description ='Read in all virgo filament tables')
    parser.add_argument('--tabledir', dest = 'tabledir', default = '/home/rfinn/research/Virgo/tables-north/v1/', help = 'directory where tables are stored')
    parser.add_argument('--tableprefix', dest = 'tableprefix', default = 'vf_north_v1_', help = 'prefix for tables; default is vf_north_v1')                               
    args = parser.parse_args()
    
    v = latextable(args.tabledir,args.tableprefix)
    v.read_all()

    
