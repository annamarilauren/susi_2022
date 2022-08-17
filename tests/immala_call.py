# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 12:10:20 2022

@author: alauren
"""

import numpy as np
import datetime
from susi_utils import  get_motti, read_FMI_weather, get_mese_input
from susi_para import get_susi_para
from susi_main import Susi
import susi_io
import glob
import os




def run_immala(s, folderName, susiPath, wpath, wdata, mottipath, age, strip_width, sfc, startyr, endyr ):

    

    mottifile = {'path':mottipath,
          'dominant':{1: s + '.xlsx'},
          'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
          'under':{0:'susi_motti_input_lyr_2.xlsx'}} 

    start_date = datetime.datetime(startyr, 1, 1)
    end_date=datetime.datetime(endyr,12,31)
    start_yr = start_date.year; end_yr = end_date.year
    yrs = (end_date - start_date).days/365.25
    length = (end_date - start_date).days +1
    
    sarkaSim = strip_width
    n = int(sarkaSim / 2)
    sfc =  np.ones(n, dtype = int)*int(sfc)                                                                         #soil fertility class
    
    ageSim = {'dominant':np.ones(n)*age ,
              'subdominant': 0*np.ones(n),
              'under': 0*np.ones(n)}                                                         # age of the stand in each node
        

    site = 'develop_scens'
    
    forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
    
    
    wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                              folderName=folderName, hdomSim=None,  
                                                                              ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                              susiPath=susiPath,
                                                                              n=n)
    spara['sarkaSim'] = sarkaSim
    spara['n'] = n
    
    outpara['netcdf'] = s + '.nc'
    
    spara['scenario name']=[s+'_30', s+'_50', s+'_70', s+'_90']
    
    susi = Susi()
    susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                            mottifile=mottifile, peat= 'other', photosite='All data', 
                            folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath)



# mottipath =  r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Immala/susi_inputs/'
# ifiles = glob.glob(mottipath +"*.xlsx")
# wdata  = c

# sites=[]
# for ifile in ifiles:    
#     _, filename = os.path.split(ifile)
#     fn, _ = filename.split(sep='.')
#     sites.append(fn)

# print (sites)

# #sites                           Susi-input (motti-simuloinneista) tiedostojen nimet
# ages = [50, 60, 70]             # iät vuosina
# strip_widths =[40, 40, 50]            # sarkaleveydet m
# sfcs = [3, 3, 4]                #site fertility classes

# folderName=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Immala/susi_22_out/'
# susiPath = r'C:/Users/alauren/Documents/Susi_10/'
# wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/'
# startyr = 2005
# endyr = 2019
                                    
# for s, age, strip_width, sfc in zip(sites, ages, strip_widths, sfcs):
#     print (s, mottipath, age, strip_width, sfc)
#     run_immala(s, folderName, susiPath, wpath, wdata, mottipath, age, strip_width, sfc, startyr, endyr )
#     import sys; sys.exit()

# #%%
# from figures import *
# folderName=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Immala/susi_22_out/'

# ff = folderName + sites[0] + '.nc'
# print (ff)
# scen = 3
# hydrology(ff, scen)
# stand(ff, scen)
# # mass(ff, scen)
# # carbon(ff, scen)        
# compare_scens(ff)

# #%%
# stand(ff, 3)

# print (sites[0])