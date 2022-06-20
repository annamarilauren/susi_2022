# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 08:47:45 2022

@author: alauren
"""

from strip import StripHydrology
import numpy as np
import matplotlib.pylab as plt


spara ={ 

    'develop_scens':{
    'species': 'Pine', 'sfc':3, 'sfc_specification': 1,
    'hdom':20, 'vol':150, 'age':15, 'smc': 'Peatland',
    'nLyrs':100, 'dzLyr': 0.05, 'L': 40, 'n':20, 
    'ditch depth west': [-0.5],   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
    'ditch depth east': [-0.5],
    'ditch depth 20y west': [-0.5],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    'ditch depth 20y east': [-0.5],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
    'scenario name': ['Control'], #kasvunlisaykset
    'initial h': -0.2, 'slope': 0.0, 
    'peat type':['A','A'], 
    'peat type bottom':['A'],'anisotropy':10.,
    'vonP': True,
    'vonP top':  [4,4], 
    'vonP bottom': 4,
    'bd top':None, 'bd bottom': 0.16,
    'peatN':1, 'peatP':1, 'peatK':1,
    'enable_peattop': True, 'enable_peatmiddle': True,
    'enable_peatbottom': True,
    'depoN': 4.0, 'depoP':0.1, 'depoK':1.0,
    'fertilization': {
            'application year': 2201,
            'N':{'dose': 0.0, 'decay_k': 0.5, 'eff': 1.0},                              # fertilization dose in kg ha-1, decay_k in yr-1
            'P':{'dose': 45.0, 'decay_k': 0.2, 'eff': 1.0},
            'K':{'dose': 100.0, 'decay_k': 0.3, 'eff': 1.0},
            'pH_increment':1.0},  

    'canopylayers': {'dominant': np.ones(20, dtype=int),
                     'subdominant': np.zeros(20, dtype=int),
                     'under': np.zeros(20, dtype=int)}    

        }}

wt = np.arange(0,5,0.01) * -1

stp4 = StripHydrology(spara['develop_scens'])
storage = stp4.dwtToSto(wt)
S4 =  np.gradient(storage) / np.gradient(wt) 

spara['develop_scens']['vonP top'] = [6,6]
spara['develop_scens']['vonP bottom'] = 6
stp6 = StripHydrology(spara['develop_scens'])
storage = stp6.dwtToSto(wt)
S6 = np.gradient(storage) / np.gradient(wt) 


spara['develop_scens']['vonP top'] = [8,8]
spara['develop_scens']['vonP bottom'] = 8
stp8 = StripHydrology(spara['develop_scens'])
storage = stp8.dwtToSto(wt)
S8 = np.gradient(storage) / np.gradient(wt) 


fig = plt.figure(num='s', figsize=(15,8))
plt.plot(S4, wt, 'b-')
plt.plot(S6, wt, 'r-')
plt.plot(S8, wt, 'g-')
