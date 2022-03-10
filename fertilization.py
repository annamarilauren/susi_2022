# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 19:31:27 2022

@author: alauren
"""
import numpy as np

class Fertilization():
    def __init__(self, spara):
        self.spara = spara        
        self.ncols = self.spara['n']
        
    def ph_effect(self, yr):
        tfert = yr - self.spara['fertilization']['application year']
        pH_increment = self.spara['fertilization']['pH_increment']*np.exp(-0.1*tfert)
        return pH_increment

    def nutrient_release(self, yr):        
        self.release={'N':0.0, 'P':0.0, 'K':0.0}                                  
        tfert = yr - self.spara['fertilization']['application year']
        if tfert >= 0:
            for nutr in ['N', 'P', 'K']:
                nut_efficiency = self.spara['fertilization'][nutr]['eff']       # nutrient use efficiency
                dose = self.spara['fertilization'][nutr]['dose']                # dose of compound in fertilizer 
                decay_k = self.spara['fertilization'][nutr]['decay_k']           # decay rate in 1/year in fertilizer     
                self.release[nutr] = (dose*np.exp(-decay_k*tfert)\
                                   - dose*np.exp(-decay_k*(tfert+1)))*nut_efficiency\
                                    * np.ones(self.ncols)  # nutrient release from the fertilizer
        
        return self.release