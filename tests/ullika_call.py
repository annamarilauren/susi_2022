# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 09:49:39 2022

@author: alauren
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:52:35 2022

@author: alauren
"""

import numpy as np
import pandas as pd
import datetime
from susi.susi_utils import  get_motti, read_FMI_weather, get_mese_input
from inputs.susi_para import get_susi_para
from susi.susi_main import Susi
import susi.susi_io

#***************** local call for SUSI*****************************************************
folderName=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/'

wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'
wdata='katila_weather.csv'


mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/',
              'dominant':{1: 'ullika_1.xlsx', 2: 'ullika_2.xlsx', 3: 'ullika_3.xlsx' },
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 


start_date = datetime.datetime(2005,1,1)
end_date=datetime.datetime(2014,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 80.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 120.*np.ones(n),
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                                         # age of the stand in each node

sfc =  np.ones(n, dtype=int)*5                                                                        # site fertility class

#ageSim['dominant'][int(n/2):] = 2.

#ageSim[4:-4] = 2.
#site = 'develop_scens'
site = 'ullika'


forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=None,
                                                                          n=n)
spara['canopylayers']['dominant'][int(n/3):int(n/3*2)] = 2
spara['canopylayers']['dominant'][int(n/3*2):] = 3
outpara['netcdf'] = 'ullika' + '.nc'

susi = Susi() 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=None)
    
          
#%%
from figures import *
ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/ullika.nc'
scen = 1
hydrology(ff, scen)
stand(ff, scen)
mass(ff, scen)
carbon(ff, scen)
             
scens = [0,1]
compare_1(ff, scens)

