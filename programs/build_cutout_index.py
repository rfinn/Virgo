#!/usr/bin/env python

'''
GOAL:
* create index web page that points to html pages for all cutouts 

USAGE:
* run from html-dev/cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)

'''

import os
import numpy as np
import glob

from astropy.io import fits
homedir = os.getenv("HOME")

###########################################################
####  FUNCTIONS
###########################################################

    
def write_coadd_prop_table(html,filter,zp,fwhm_arcsec):
    html.write('<h3>Image Characteristics</h3>\n')
    html.write('<table>\n')
    html.write('<tr>')
    html.write('<th ">Filter</th>\n')
    html.write('<th ">ZP<br />mag</th>\n')
    html.write('<th ">PSF FWHM <br />arcsec</th>\n')
    html.write('</tr>\n')
    html.write('<tr><td>{}</td><td>{:.2f}</td><td>{:.2f}</td>\n'.format(filter,zp,fwhm_arcsec))
    html.write('</tr>\n')
    html.write('</table>\n')

def write_table(html,images,labels):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for i in images:
        html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
    html.write('</tr>\n')            
    html.write('</table>\n')

def write_text_table(html,labels,data):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for d in data:
        html.write('<td>{}</td>'.format(d))
    html.write('</tr>\n')            
    html.write('</table>\n')


def legacy_link(ra,dec):
    return "https://www.legacysurvey.org/viewer?ra={:.4f}&dec={:.4f}&layer=ls-dr9&zoom=13".format(ra,dec)    
###########################################################
####  CLASSES
###########################################################
    
class build_html_cutout():

    def __init__(self,gallist,outdir,co=False):
        ''' pass in instance of cutout_dir class and output directory '''

        self.coflag = co
        if self.coflag:
            outfile = os.path.join(outdir,'indexco.html')
        else:
            outfile = os.path.join(outdir,'index.html')
        self.outdir = outdir

        #print('inside build html')
        #print('coutdir = ',coutdir)
        #print('outfile = ',outfile)        
        self.html = open(outfile,'w')
        self.galnames = gallist
        self.build_html()
        
    def build_html(self):
        self.write_header()
        self.write_gal_list()
        #self.write_navigation_links()
        self.close_html()
    def write_header(self):
        self.html.write('<html><body>\n')
        self.html.write('<style type="text/css">\n')
        self.html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        self.html.write('</style>\n')
        self.html.write('<h1>VF Galaxies</h1>\n')        
    def write_gal_list(self):
        self.html.write('<table width="90%">\n')
        self.html.write('<tr>')
        colnames = ['Index','VFID','Legacy Image<br> and link','Cutouts and <br>Analysis  ','RA','DEC','vr <br> (km/s)','logMstar','logSFR','logsSFR','CO Flag','ALFALFA','Filament member','Nearest Filament']
        for i,l in enumerate(colnames):
            if i == 1:
                colspan=2
            else:
                colspan=1
            self.html.write('<th colspan="{}">{}</th>'.format(colspan,l))
        self.html.write('</tr></p>\n')            

        vfindices = np.arange(len(vfmain))
        # write one row for each galaxy
        galindex = 1
        ids = []
        for i,g in enumerate(self.galnames):
            #print(g)
            vfid = g.split('-')[0]
            vfindex = vfindices[vfmain['VFID'] == vfid][0]
            ra = vfmain['RA'][vfindex]
            dec = vfmain['DEC'][vfindex]
            if self.coflag & ~vfmain['COflag'][vfindex]:
                continue
            ids.append(vfmain['VFID'][vfindex])
            #print(vfindex)
            # get legacy jpg name

            jpg_path = os.path.join(self.outdir,g)
            search_path = os.path.join(jpg_path,'*legacy*.jpg')
            #print(search_path)
            #legacy_jpg = glob.glob(search_path)[0]            
            try:
                #print()
                #print("looking for legacy image")
                #print(glob.glob(search_path))
                legacy_jpg = glob.glob(search_path)[0]
                legacy_flag = True                
            except:
                legacy_flag = False
                legacy_jpg = None
                print('WARNING: no legacy image for ',g)
                #print('\t Skipping galaxy for now')
                #continue
            self.html.write('<tr>')
            self.html.write('<td>{}</td>'.format(galindex))
            self.html.write('<td>{}</td>'.format(vfid))            
            htmlpage = "{}/{}.html".format(g,g)
            if legacy_flag:
                relative_path_2legacyjpg = '{}/{}'.format(g,os.path.basename(legacy_jpg))
                self.html.write('<td colspan="2"><a href="{0}" target="_blank"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(legacy_link(ra,dec),relative_path_2legacyjpg))

                
                
            else:
                self.html.write('<td colspan="2"><a href="{0}" target="_blank">Missing</td>'.format(legacy_link(ra,dec)))
            self.html.write('<td><a href="{}">{}</td>'.format(htmlpage,g))                
            self.html.write('<td>{:.3f}</td>'.format(ra))
            self.html.write('<td>{:.3f}</td>'.format(dec))
            self.html.write('<td>{:.0f}</td>'.format(vfmain['vr'][vfindex]))
            self.html.write(f"<td>{vfmagphys['logMstar'][vfindex]:.2f}</td>")
            self.html.write(f"<td>{vfmagphys['logSFR'][vfindex]:.2f}</td>")
            self.html.write(f"<td>{vfmagphys['logsSFR'][vfindex]:.2f}</td>")            
            if vfmain['COflag'][vfindex]:
                text='CO'
            else:
                text='-'
            self.html.write('<td>{}</td>'.format(text))
            if vfmain['A100flag'][vfindex]:
                text='HI'
            else:
                text='-'
            self.html.write('<td>{}</td>'.format(text))
            if vffil['filament_member'][vfindex]:
                text='memb'
            else:
                text='-'
            self.html.write('<td>{}</td>'.format(text))
            self.html.write('<td>{}</td>'.format(vffil['filament'][vfindex]))            

            self.html.write('</tr>\n')
            galindex += 1
        self.html.write('</table>\n')
        if self.coflag:
            print('')
            print('################  CO STATS  ###############')

        else:
            print('')
            print('################  HALPHA STATS  ###############')
        print('\t number of images = {}'.format(len(ids)))
        print('\t number of unique targets = {}'.format(len(set(ids))))
    
    def close_html(self):
        self.html.close()
# wrap

if __name__ == '__main__':
    # work through coadd directory
    #global vmain

    VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
    VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'    
    vfmain = fits.getdata(VFMAIN_PATH)

    # updating for v2 catalogs
    VFFIL_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main_filament_membership_allgalaxies.fits'
    VFFIL_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_environment.fits'    
    vffil = fits.getdata(VFFIL_PATH)

    VFMAGPHYS_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_magphys_10-Jul-2023.fits'    
    vfmagphys = fits.getdata(VFMAGPHYS_PATH)
    
    outdir = homedir+'/research/Virgo/html-dev/cutouts/'
    outdir = '/data-pool/Halpha/html_dev/cutouts/'    

    # this should contain a list of all the galaxy folders
    flist1 = os.listdir(outdir)
    flist1.sort()
    galnames=[]
    for i,subdir in enumerate(flist1): # loop through list
        

        #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
        if (os.path.isdir(subdir)) & (subdir.startswith('VF')):
            #print('adding ',subdir)
            galnames.append(subdir)
    print('number of subdirectories = ',len(galnames))
    h = build_html_cutout(galnames,outdir)
    h = build_html_cutout(galnames,outdir,co=True)
