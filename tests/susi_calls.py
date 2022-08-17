# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import datetime
from susi.susi_utils import read_FMI_weather
from inputs.susi_para import get_susi_para
from susi.susi_main import Susi


#***************** local call for SUSI*****************************************************
folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/' #'sensitivity/'
wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/inputs/'


mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/inputs/',
              'dominant':{1: 'motti viitasaari_mtkg.xls'},
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 

wdata='parkano_weather.csv'

start_date = datetime.datetime(2000,1,1)
end_date=datetime.datetime(2002,12,31)
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
                                                                          n=n)
#spara['canopylayers']['dominant'][int(n/2):] = 2                                                                        
#spara['canopylayers']['subdominant'][:int(n/2)] = 1                                                                        

susi = Susi()
 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    
          
             
#%%
from susi.figures import *
ff = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/susi.nc'
# scen = 1
# hydrology(ff, scen)
# stand(ff, scen)
# mass(ff, scen)
# carbon(ff, scen)
# nutrient_balance(ff, 'N', scen)
# nutrient_balance(ff, 'P', scen)
# nutrient_balance(ff, 'K', scen)
# compare_1(ff, [0,3])
compare_scens(ff)