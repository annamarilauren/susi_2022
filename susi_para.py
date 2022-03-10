# -*- coding: utf-8 -*-
"""
Created on Mon May 21 18:40:35 2018

@author: lauren
"""
import numpy as np


def get_susi_para(wlocation=None, peat=None, photosite='All data', folderName=None, hdomSim=None, volSim=None, 
                  ageSim=None, sarkaSim=None, sfc=None, susiPath = None, ddwest=None, ddeast=None, n=None, bd=None,
                  peatN=None, peatP=None, peatK=None):
    #********** Stand parameters and weather forcing*******************
    #--------------Weather variables 10 km x 10 km grid   
    if susiPath is None: susiPath=""
    wpara ={

        'undefined': {
        'infolder': susiPath + '\\wfiles\\',
        'infile_d':'Tammela_weather_1.csv',
        'start_yr': 1980, 'end_yr': 1984, 
        'description': 'Undefined, Finland',
        'lat': 65.00, 'lon': 25.00},

        }

    cpara = {'dt': 86400.0,
            'flow' : { # flow field
                     'zmeas': 2.0,
                     'zground': 0.5,
                     'zo_ground': 0.01
                     },
            'interc': { # interception
                        'wmax': 0.5, # storage capacity for rain (mm/LAI)
                        'wmaxsnow': 4.0, # storage capacity for snow (mm/LAI),
                        },
            'snow': {
                    # degree-day snow model
                    'kmelt': 2.8934e-05, # melt coefficient in open (mm/s)
                    'kfreeze': 5.79e-6, # freezing coefficient (mm/s)
                    'r': 0.05 # maximum fraction of liquid in snow (-)
                    },

            'physpara': {
                        # canopy conductance
                        'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                        'g1_conif': 2.1, # stomatal parameter, conifers
                        'g1_decid': 3.5, # stomatal parameter, deciduous
                        'q50': 50.0, # light response parameter (Wm-2)
                        'kp': 0.6, # light attenuation parameter (-)
                        'rw': 0.20, # critical value for REW (-),
                        'rwmin': 0.02, # minimum relative conductance (-)
                        # soil evaporation
                        'gsoil': 1e-2 # soil surface conductance if soil is fully wet (m/s)
                        },
            'phenopara': {
                        #seasonal cycle of physiology: smax [degC], tau[d], xo[degC],fmin[-](residual photocapasity)
                        'smax': 18.5, # degC
                        'tau': 13.0, # days
                        'xo': -4.0, # degC
                        'fmin': 0.05, # minimum photosynthetic capacity in winter (-)
                        },
            'state': {
                       'lai_conif': 3.0, # conifer 1-sided LAI (m2 m-2)
                       'lai_decid_max': 0.01, # maximum annual deciduous 1-sided LAI (m2 m-2): 
                       'hc': 16.0, # canopy height (m)
                       'cf': 0.7, # canopy closure fraction (-)
                       #initial state of canopy storage [mm] and snow water equivalent [mm]
                       'w': 0.0, # canopy storage mm
                       'swe': 0.0, # snow water equivalent mm
                       }
            }
    org_para = {
           'org_depth': 0.04, # depth of organic top layer (m)
           'org_poros': 0.9, # porosity (-)
           'org_fc': 0.3, # field capacity (-)
           'org_rw': 0.24, # critical vol. moisture content (-) for decreasing phase in Ef
           'pond_storage_max': 0.01, # max ponding allowed (m)
           #initial states
           'org_sat': 1.0, # organic top layer saturation ratio (-)
           'pond_storage': 0.0 # pond storage
            }
        
    # Hannun parametrit
    #------------ Soil and stand parameters ----------------------------------
    spara ={ 

        'develop_scens':{
        'species': 'Pine', 'sfc':sfc, 'sfc_specification': 1,
        'hdom':hdomSim, 'vol':volSim, 'age':ageSim, 'smc': 'Peatland',
        'nLyrs':30, 'dzLyr': 0.05, 'L': sarkaSim, 'n':n, 
        'ditch depth west': [-0.3],   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
        'ditch depth east': [-0.3],
        'ditch depth 20y west': [-0.3],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
        'ditch depth 20y east': [-0.3],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
        'scenario name': ['Control'], #kasvunlisaykset
        'initial h': -0.2, 'slope': 2.0, 
        'peat type':['A','A','A','A','A','A','A','A'], 
        'peat type bottom':['A'],'anisotropy':10.,
        'vonP': True,
        'vonP top':  [2,2,2,3,4,5,6,6], 
        'vonP bottom': 8,
        'bd top':None, 'bd bottom': 0.16,
        'peatN':peatN, 'peatP':peatP, 'peatK':peatK,
        'enable_peattop': True, 'enable_peatmiddle': True,
        'enable_peatbottom': True,
        'depoN': 4.0, 'depoP':0.1, 'depoK':1.0,
        'fertilization': {
                'application year': 2201,
                'N':{'dose': 0.0, 'decay_k': 0.5, 'eff': 1.0},                              # fertilization dose in kg ha-1, decay_k in yr-1
                'P':{'dose': 45.0, 'decay_k': 0.2, 'eff': 1.0},
                'K':{'dose': 100.0, 'decay_k': 0.3, 'eff': 1.0},
                'pH_increment':1.0},  

        'canopylayers': {'dominant': np.ones((int(n)), dtype=int),
                         'subdominant': np.zeros((int(n)), dtype=int),
                         'under': np.zeros((int(n)), dtype=int)}    

            },

        'krycklan':{
        'species': 'Spruce', 'sfc':sfc, 'sfc_specification': 1,
        'hdom':hdomSim, 'vol':volSim, 'age':ageSim, 'smc': 'Peatland',
        'nLyrs':30, 'dzLyr': 0.05, 'L': sarkaSim, 'n':n, 
        'ditch depth west': [-0.2, -0.9],   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
        'ditch depth east': [-0.2, -0.9],
        'ditch depth 20y west': [-0.2, -0.9],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
        'ditch depth 20y east': [-0.2, -0.9],                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
        'scenario name': ['Control', 'DNM90'], #kasvunlisaykset
        'initial h': -0.2, 'slope': 3.0, 
        'peat type':['A','A','A','A','A','A','A','A'], 
        'peat type bottom':['A'],'anisotropy':10.,
        'vonP': True,
        'vonP top':  [2,2,2,2,4,5,5,5], 
        'vonP bottom': 5,
        'bd top':None, 'bd bottom': 0.16,
        'peatN':1.2, 'peatP':0.12, 'peatK':0.07,                               #peat nutrient contents in gravimetric %
        'enable_peattop': True, 'enable_peatmiddle': False,
        'enable_peatbottom': False,
        'depoN': 4.0, 'depoP':0.1, 'depoK':1.0,
        'fertilization': {
                'application year': 2201,
                'N':{'dose': 0.0, 'decay_k': 0.5, 'eff': 1.0},                              # fertilization dose in kg ha-1, decay_k in yr-1
                'P':{'dose': 45.0, 'decay_k': 0.2, 'eff': 1.0},
                'K':{'dose': 100.0, 'decay_k': 0.3, 'eff': 1.0},
                'pH_increment':1.0},  

        'canopylayers': {'dominant': np.ones((int(n)), dtype=int),
                         'subdominant': np.zeros((int(n)), dtype=int),
                         'under': np.zeros((int(n)), dtype=int)}    

            },

        }
    #------------  Output parameters -------------------------------------------------    
    outpara ={
        'outfolder':folderName, 
        'netcdf': 'susi.nc',
        'startday': 1, 'startmonth':7, # Päivä, josta keskiarvojen laskenta alkaa
        'endday':31, 'endmonth':8, # Päivä, johon keskiarvojen laskenta loppuu
        #'figs': True, 'to_file':True, 'static stand':False, 'hydfig':True, 'DOCfig':False, 
        }    
    photopara = {
              'All data':
                  {'beta':0.513,
                   'gamma':0.0196,
                   'kappa': -0.389,
                   'tau':7.2,
                   'X0': -4.0,
                   'Smax':17.3,
                   'alfa': 1.,
                   'nu': 5.
                   },
              'Sodankyla':
                  {'beta':0.831,
                   'gamma':0.065,
                   'kappa': -0.150,
                   'tau':10.2,
                   'X0': -0.9,
                   'Smax':16.4,
                   'alfa': 1.,
                   'nu': 5.
                   },
              'Hyytiala':
                  {'beta':0.504,
                   'gamma':0.0303,
                   'kappa': -0.235,
                   'tau':11.1,
                   'X0': -3.1,
                   'Smax':17.3,
                   'alfa': 1.,
                   'nu': 5.
                   },
              'Norunda':
                  {'beta':0.500,
                   'gamma':0.0220,
                   'kappa': -0.391,
                   'tau':5.7,
                   'X0': -4.0,
                   'Smax':17.6,
                   'alfa':1.062,
                   'nu': 11.27,
                   },
              'Tharandt':
                  {'beta':0.742,
                   'gamma':0.0267,
                   'kappa': -0.512,
                   'tau':1.8,
                   'X0': -5.2,
                   'Smax':18.5,
                   'alfa':1.002,
                   'nu': 442.
                   },
              'Bray':
                  {'beta':0.459,
                   'gamma':-0.000669,
                   'kappa': -0.560,
                   'tau':2.6,
                   'X0': -17.6,
                   'Smax':45.0,
                   'alfa':0.843,
                   'nu': 2.756,
                   },
                           }
    #----------- Arrange and make coherent------
    #cpara['lat']=wpara[wlocation]['lat']; cpara['lon']=wpara[wlocation]['lon']

    o_w = wpara[wlocation] if wlocation is not None  else wpara 
    o_s = spara[peat] if peat is not None else spara
    o_p = photopara[photosite] if photosite is not None else photopara  
    return o_w, cpara, org_para, o_s, outpara,o_p
