#!/usr/bin/env python
'''
GOAL:
* write out tables for A100-SDSS paper (20 rows)
* write out fits tables for full A100 sample

OUTPUT:
* tab_catalog.tex
  - table X from Castignani+2021

* tab_environment.tex 
  - table X+1 from Castignani+2021

* fits versions of table1 and table2

USAGE:

python write_tables_catalog_paper.py

UPDATES:
* 2020-08-06
  - making empty array elements nan instead of, e.g. -99

'''


import numpy as np
import os
import shutil
from astropy.io import fits, ascii
from astropy.table import Table
from datetime import datetime

homedir = os.getenv('HOME')
tablepath = homedir+'/research/APPSS/tables/'
tablepath = homedir+'/research/Virgo/tables-north/v1/'
latextablepath = homedir+'/research/Virgo/papers/catalog_paper/'


class latextable():
    def __init__(self):
        # read in virgo tables
        
        self.tab = fits.getdata(tablepath+'a100-sdss-wise.fits')
        # the next table has the NSA data appended to the end columns
        # it contains the full NSA, but the first rows should be line-matched to the a100-sdss-wise

        ### THIS IS THE WRONG TABLE
        ### SHOULD BE USING a100-sdss-wise-nsa-gswlcA2.fits
        #self.tab2 = fits.getdata(tablepath+'full-a100-sdss-wise-nsa-gswlcA2.fits')
        self.tab2 = fits.getdata(tablepath+'a100-sdss-wise-nsa-gswlcA2.fits')        

        
        self.tab2 = self.tab2[0:len(self.tab)] # trim table to keep only A100 rows
        print('length of tab2 = ',len(self.tab2),' should be 31502')
        self.agcdict2=dict((a,b) for a,b in zip(self.tab2['AGC'],np.arange(len(self.tab2['AGC']))))
        # calculate distance

        # dist
        self.dist_Virgo
        pass
    def clean_arrays(self):
        '''
        remove bogus values from SFR estimates and other arrays with null values

        '''
        pass
    def print_table1(self,nlines=10,filename=None,papertableflag=True):
        '''write out latex version of table 1 '''
        if filename is None:
            fname=latextablepath+'catalog1.tex'
        else:
            fname = filename 
        outfile = open(fname,'w')
        
        outfile.write('\\begin{table*}%[ptbh!]\n')
        outfile.write('\\begin{center}\n')
        outfile.write('\\scriptsize\n')
        outfile.write('\\setlength\\tabcolsep{3.0pt} \n')
        #outfile.write('\\tablenum{1} \n')
        outfile.write('\\caption{Basic Optical Properties of Cross-listed objects in the ALFALFA-SDSS Catalog.\label{tab:catalog1}  } \n')
        outfile.write('\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n')
        outfile.write('\\hline \n')
        outfile.write('\\toprule \n')
        outfile.write('AGC &	Flag &	SDSS objID  & RA &	DEC &	V$_{r}$ &	D &	$\sigma_D$  &	Ext$_g$	& Ext$_i$	& expAB$_r$  & $\sigma_{expAB_r}$ &	cmodel$_i$ & $\sigma_{cmodel_i}$ \\\\ \n')
        outfile.write('& & & J2000 & J2000 & $km~s^{-1}$ & Mpc & Mpc & mag & mag & & & mag & mag \\\\ \n')
        outfile.write('(1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14)  \\\\ \n')
        outfile.write('\\midrule \n')
        outfile.write('\\hline \n')
        for i in range(nlines): # print first N lines of data
            # vfid ra dec z 
            
            s = '{0:d} & {1:d} & {2:d} & {3:9.6f}&{4:9.5f} & {5:d} & {6:.1f} & {7:.1f} &  {8:.2f} & {9:.2f}& {10:.2f}& {11:.2f}& {12:.2f} &{13:.2f}\\\\ \n'.format(self.tab['AGC'][i],self.tab['sdssPhotFlag'][i],self.tab['objID_1'][i],self.tab['RAdeg_Use'][i],self.tab['DECdeg_Use'][i],self.tab['Vhelio'][i],self.tab['Dist'][i],self.tab['sigDist'][i],self.tab['extinction_g'][i],self.tab['extinction_i'][i],self.tab['expAB_r'][i],self.expAB_r_err[i],self.tab['cModelMag_i'][i],self.tab['cModelMagErr_i'][i])
            if papertableflag:
                # replace nans with \\nodata
                s=s.replace('nan','\\nodata')

            outfile.write(s)

        outfile.write('\\bottomrule \n')
        outfile.write('\\hline \n')
        outfile.write('\\end{tabular} \n')
        outfile.write('\\end{center} \n')
        outfile.write('\\tablecomments{Table 1 is published in its entirety in the machine-readable format.  A portion is shown here for guidance regarding its form and content.}')        
        outfile.write('\\end{table*} \n')
        outfile.close()
    def print_table2(self,nlines=10,filename=None,papertableflag=True):
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
            fname=latextablepath+'catalog2.tex'
        else:
            fname = filename 
        outfile = open(fname,'w')
        outfile.write('\\begin{sidewaystable*}%[ptbh!]%\\tiny \n')
        outfile.write('\\begin{center}\n')
        outfile.write('\\tablewidth{0.5\\textwidth} \n')
        outfile.write('\\scriptsize \n')
        outfile.write('%\\footnotesize \n')
        outfile.write('\\setlength\\tabcolsep{1.0pt} \n')
        #outfile.write('\\tablenum{2}\n')
        outfile.write('\\caption{Derived Properties of Cross-listed objects in the ALFALFA-SDSS Catalog.\\label{tab:catalog2}  }\n')
        outfile.write('\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n')
        outfile.write('\\hline\n')
        outfile.write('\\toprule\n')
        #outfile.write('AGC & $\\gamma_g$ & $\\sigma_{\\gamma_g}$ & $\\gamma_i$ & $\\sigma_{\\gamma_i}$ & M$_{icorr}$ &	$\\sigma_{M_{icorr}}$ &	(g-i)$_{corr}$	& $\\sigma_{(g-i)_{corr}}$ & log M$_{\\star, Taylor}$ &	$\\sigma_{log M_{\\star,Taylor}}$  & log M$_{\\star, McGaugh}$ &	$\\sigma_{log M_{\\star, McGaugh}}$ & SFR$_{22}$ & $\\sigma_{log SFR_{22}}$ & ${SFR_{NUV}}$ & $\\sigma_{log SFR_{NUV}}$ & SFR$_{NUVIR}$ & $\\sigma_{log SFR_{UVIR}}$ & M$_{HI}$ & $\\sigma_{M_{HI}}$  \\\\\n')
        outfile.write('AGC & $\\gamma_g$ &  $\\gamma_i$  & M$_{icorr}$ &	$\\rm \\sigma_{M_{icorr}}$ &	(g-i)$_{corr}$	& $\\sigma_{(g-i)_{corr}}$ & log M$_{\\star}$ &	$\\rm \\sigma_{log M_{\\star}}$  & log M$_{\\star}$ &	$\\rm \\sigma_{log M_{\\star}}$ & log M$_{\\star}$& $\\rm \\sigma_{log M_{\\star}}$ & logSFR$_{22}$ & $\\rm \\sigma_{log SFR_{22}}$  & logSFR$\\rm _{NUVcor}$ & $\\rm \\sigma_{log SFR_{NUVcor}}$  & logSFR& $\\rm \\sigma_{logSFR}$ & M$_{HI}$ & $\\rm \\sigma_{M_{HI}}$  \\\\\n')
        outfile.write('&   & &  & &	& & Taylor & Taylor  & McGaugh & McGaugh &  GSWLC &GSWLC& & &  & & GSWLC & GSWLC &  &   \\\\\n')
        #outfile.write(' & mag & mag & mag & mag & mag & mag & $log(M_\\odot)$ & $log(M_\\odot)$ & $log(M_\\odot)$ & $log(M_\\odot)$& $log(M_\\odot)$& $log(M_\\odot)$ & $\\rm log(M_\\odot~yr^{-1})$ & $\\rm log(M_\\odot yr^{-1})$ & $\\rm log(M_\\odot~yr^{-1})$  & $\\rm log(M_\\odot) yr^{-1}$ &  $log(M_\\odot~yr^{-1})$&  $log(M_\\odot~yr^{-1})$ & $log(M_\\odot)$ & $log(M_\\odot)$    \\\\\n')

        ## removing units from error columns to decrease table width
        outfile.write(' & mag & mag & mag & mag & mag & mag & $log(M_\\odot)$ &  & $log(M_\\odot)$ & & $log(M_\\odot)$&  & $\\rm log(M_\\odot~yr^{-1})$ &  & $\\rm log(M_\\odot~yr^{-1})$  &  &  $log(M_\\odot~yr^{-1})$&  & $log(M_\\odot)$ &     \\\\\n')        
        outfile.write('(1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) & (14) & (15) & (16) & (17) & (18) & (19) &(20) &(21)  \\\\\n')
        outfile.write('\\midrule\n')
        outfile.write('\\hline\n')
        for i in range(nlines):
            try:
                j = self.agcdict2[self.tab['AGC'][i]]
            except KeyError:
                print(self.tab['AGC'][i])
            #print(self.tab['AGC'][i],j,self.tab2['AGC'][j])
            s=' {0:d} & {1:.2f} & {2:.2f} & {3:.2f} & {4:.2f}& {5:.2f}  & {6:.2f} & {7:.2f} & {8:.2f} & {9:.2f}& {10:.2f}&{11:.2f} &{12:.2f} &{13:.2f} &{14:.2f} &{15:.2f} &{16:.2f}&{17:.2f}&{18:.2f}&{19:.2f}&{20:.2f}  \\\\ \n'.format(self.tab['AGC'][i],self.tab['gamma_g'][i],self.tab['gamma_i'][i],self.tab['absMag_i_corr'][i],self.absMag_i_corr_err[i],self.tab['gmi_corr'][i],self.gmi_corr_err[i],self.tab['logMstarTaylor'][i],self.logMstarTaylor_err[i],self.tab['logMstarMcGaugh'][i],self.logMstarMcGaugh_err[i],self.gsw_mstar[i],self.gsw_mstar_err[i],self.sfr22[i],self.sfr22_err[i], self.sfrnuvir[i],self.sfrnuvir_err[i],self.gsw_sfr[i],self.gsw_sfr_err[i],self.tab['logMH'][i],self.tab['siglogMH'][i])
            if papertableflag:
                # replace nans with \\nodata
                s=s.replace('nan','\\nodata')
            outfile.write(s)
        outfile.write('\\bottomrule\n')
        outfile.write('\\hline\n')
        outfile.write('\\end{tabular}\n')
        outfile.write('\\end{center} \n')
        outfile.write('\\tablecomments{Table 2 is published in its entirety in the machine-readable format.  A portion is shown here for guidance regarding its form and content.}')
        
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
    t = latextable()
    t.calculate_errors()
    t.clean_arrays()
    t.print_table1()
    t.print_table2()
    t.write_full_tables()
