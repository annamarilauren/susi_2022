# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import pandas as pd
import datetime
from susi_utils import  get_motti, read_FMI_weather, get_mese_input
from susi_para import get_susi_para
from susi_main import run_susi
import susi_io


#***************** local call for SUSI*****************************************************
folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/' #'sensitivity/'
susiPath = r'C:/Users/alauren/Documents/Susi_9/'
wpath = r'C:/Users/alauren/Documents/Susi_9/'
#mottipath =  r'C:/Users/alauren/Documents/Susi_9/'
#mottipath =  r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/'

#mf='motti viitasaari_mtkg.xls'
#mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/',
#              'dominant':{1: 'susi_motti_input_lyr_0.xlsx', 2:'susi_motti_input_lyr_1.xlsx'},
#              'subdominant':{1:'susi_motti_input_lyr_1.xlsx'},
#              'under':{1:'susi_motti_input_lyr_2.xlsx'}} 

mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/',
              'dominant':{1: 'motti viitasaari_mtkg.xls'},
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 

wdata='parkano_weather.csv'

start_date = datetime.datetime(2000,1,1)
end_date=datetime.datetime(2009,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 45.*np.ones(n),
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                                         # age of the stand in each node

sfc =  np.ones(n, dtype=int)*3                                                                        # site fertility class

#ageSim['dominant'][int(n/2):] = 2.

#ageSim[4:-4] = 2.
site = 'develop_scens'

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=susiPath,
                                                                          n=n)
#spara['canopylayers']['dominant'][int(n/2):] = 2                                                                        
#spara['canopylayers']['subdominant'][:int(n/2)] = 1                                                                        

 
run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath)
    
          
             
#%%
from figures import *
ff = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/susi.nc'
scen = 0
hydrology(ff, scen)
stand(ff, scen)
mass(ff, scen)
carbon(ff, scen)
