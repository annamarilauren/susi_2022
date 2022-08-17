# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 08:05:30 2022

@author: alauren
"""
import numpy as np

class Methane():
    def __init__(self, ncols, yrs):
        self.ncols = ncols
        self.yrs = yrs
        
    def run_ch4_yr(self, yr, dfwt):
        """
        Leena Stenberg 2021
        Input:
            yr: year
            dfwt: dataframe of water table depths (m, negative down)
        Output: 
            CH4 flux, node-wise kg CH4 ha-1 year-1
            CH4 flux, mean kg CH4 ha-1 year-1 over computation nodes
            CH4 flux, mean kg CO2-eq. ha-1 year-1 over computation nodes
        """
        
        wt = dfwt[str(yr)+'-05-01':str(yr)+'-10-31'].mean().values*-100.     # convert to cm positive down
        # wt = wt[1:-1] # exclude ditch nodes
        
        CH4 = (-0.378 + 12.3*np.exp(-0.121*wt))*10. # Ojanen et al. 2010, Fig. 6, convert to kg CH4 /ha/year
        
        CH4mean = np.mean(CH4, axis=0)
    
        return CH4, CH4mean, CH4mean*54. # convert to kg CO2-eq./ha/year, SGWP100 coefficient = 54 (Ojanen et al. 2021)

