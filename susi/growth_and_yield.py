# -*- coding: utf-8 -*-
"""
This module applies growth and yield model by Pukkala (2021) to produce 
simmilar input array to Motti-simulator.
This is aimed to be an intermediate phase to fully integrated growth and yield and allometry module
"""

import numpy as np
import pandas as pd
import xlsxwriter

from scipy.optimize import minimize
from susi.stem_curve import StemCurve

class Growth_and_Yield_Table():
   
    def __init__(self,G, N, Dg, D, Hg, Hdom, DDY, species, fertility_class,W = None, Nd = None):
        self.G = G                           # basal area, m2/ha
        self.N = N                           # stem number, /ha
        self.Dg = Dg                         # mean diameter, cm
        self.D = D                           # diameter classes, array: 1...50 cm
        self.W = W                           # Weibull distribution, sum = 1
        self.Hg = Hg                         # mean height, m
        self.Hdom = Hdom                     # dominant height, m
        self.DDY = DDY                       # temperature sum, degree days
        self.species = species               #'pine' tree species as string
        self.fertility_class = fertility_class   # as integer 1...6
        
        #herb = 0                            # fertility class dummy variables
        #mesic = 0
        #subxeric = 1
        #xeric = 0
        self.peat = 1
        self.dxpendula = 0
        self.aspen = 0
        self.birch = 0 
        sp ={'pine':1, 'spruce':2, 'birch':3}
        self.speciesid = sp[species]                   # species id as integer              

        
        if Nd is not None:
            self.Nd = Nd                     # Number of trees in diameter classes
        else:
            b,c =  self.recover_weibull(self.Dg,self.N,self.G)           # Recover Weibull parameters using known basal area
            self.W = c / b * (D/b)**(c-1) * np.exp(-(D/b)**c)
            self.Nd = self.W * self.N
    
    def naslund_to_hdom(self, D, Nd, species):
        """
        Calculates individual tree height and the subsequent dominant height using diameter only. 
        Hdom is then used in the actual tree height calculation.
        Siipilehto & Kangas 2015, mean model.
        D - diameter classes 1...50
        species - tree species as string
        """
        para = {'pine': {'m':2, 'b0': 1.195,'b1': 0.242, 'sy': 0.238},
               'spruce':{'m':3, 'b0': 1.635, 'b1': 0.330,  'sy': 0.227 },
               'birch': {'m':2, 'b0': 0.898, 'b1': 0.242, 'sy': 0.169 }}
        p=para[species]  
        b0, b1, sy = p['b0'],   p['b1'],  p['sy']                                                                             # light species m=2, shade tolerant m=3
        if species == 'spruce':
           h =  np.power(D, p['m']) / np.power((b0 + b1*D), p['m']) + (6*np.power(D, p['m']) / np.power((b0 + b1*D), 5))*sy + 1.3
        else:
           h =  np.power(D, p['m']) / np.power((b0 + b1*D), p['m']) + (3*np.power(D, p['m']) / np.power((b0 + b1*D), 4))*sy + 1.3
        
        h = np.square(D) / np.square(b0+b1*D) + (3*np.power(D,3)/ np.power((b0+b1*D),4)) *sy +1.3
        Nd = Nd.values
        tmp = np.cumsum(Nd)[-1]-np.cumsum(Nd)                   # cumulative number of trees from thickest to thinnest dimater
        ix = np.where(tmp<100)[0]  
        #print (np.shape(ix))
        #print (h[ix], ix[0], Nd[ix])
        #print (np.average(h[ix],weights=Nd[ix]))
        #import matplotlib.pylab as plt
        #plt.plot(ix[0], tmp[ix])
        #plt.plot(ix[0], h[ix])
        #import sys; sys.exit()
        return np.around(np.average(h[ix],weights=Nd[ix]), decimals=2)   # weighed average of thick diameters (weighed with number of trees)

     

    def naslund_height(self, DDY, HDOM, D, G, species):                                      # individual tree height, arguments: temperature sum, dominant height, basal area, stem diameter
       """
       Siipilehto, J. & Kangas, A. 2015. Näslundin pituuskäyrä ja siihen perustuvia malleja läpimitan ja
       pituuden välisestä riippuvuudesta suomalaisissa talousmetsissä. Metsätieteen aikakauskirja 4/2015:
       215–236.
       Siipilehto edition 2022, birch:
       input: 
           DDY  - temperature sum, degree days, float
           HDOM  -  dominant height of the stand, m
           D  -  diameter classes, array 1...50, integer
           G  -  basal area m2 ha-1
           species - species as string
        return:
            tree height, m, array of floats, dimension: to all diameter classes
       """
       para = {'pine': {'m':2, 'k0': 7.136,'k1': -0.686, 'k2': -0.273, 'k3': 0.00139 , 'k4':- 0.329, 'k5':0.0,
                                      'm0': 0.530 ,'m1': 0.0136 , 'm2': -0.128, 'm3':0.0, 'varf': 0.944},
               'spruce':{'m':3,'k0': 9.833 ,'k1': -1.185, 'k2': 0.0, 'k3': 0.00176 , 'k4': -0.188, 'k5': 0.0,
                                      'm0': 0.710 ,'m1': 0.0  , 'm2': -0.133, 'm3':0.0, 'varf':0.813},
               'birch':{'m':2,'k0': 4.206 ,'k1':-0.708 , 'k2': 0.0 , 'k3': 9.15*10**(-4) , 'k4': 0.0, 'k5':0.098 ,
                                      'm0': 0.474 ,'m1':0.0 , 'm2': 0.0, 'm3':-0.104, 'varf': 1.056 }}
       p=para[species]                                                                                 # light species m=2, shade tolerant m=3
    
       b0 =  p['k0'] + p['k1']* np.log(DDY) +p['k2'] * np.log(G + 1) + \
               p['k3']* np.square(HDOM) + p['k4'] * np.log(HDOM - 1) + p['k5'] * np.log(HDOM - 4)
    
       b1 = p['m0'] + p['m1']* np.log(G + 1) + p['m2'] * np.log(HDOM-1) + p['m3'] * np.log(HDOM-4)
    
       if species == 'spruce':
          h =  np.power(D, p['m']) / np.power(b0 + b1*D, p['m']) + (6*np.power(D, p['m']) / np.power(b0 + b1*D, 5))*p['varf']*1000/DDY + 1.3
       else:
          h =  np.power(D, p['m']) / np.power(b0 + b1*D, p['m']) + (3*np.power(D, p['m']) / np.power(b0 + b1*D, 4))*p['varf']*1000/DDY + 1.3
       return h       # tree height in m


    def get_biomass(self, D, h, Nd, species):                                                       # generic biomass equation for all biomass components
        """
        Repola 2009, Biomass equations for Scots pine and Norway spruce in Finland. Silva Fennica 43(4): 625–647.
        Repola 2008, Biomass equations for birch in Finland. Silva Fennica 42(4): 605–624.
        input:
            D  -  diameter classes, unit cm, array of integers, 1...50
            h  -  tree height, m, in all diameter classes
            Nd  - number of trees in each diameter class, array of floats
        return:
            out - stand biomass, kg ha-1, in each biomass component
            outs  - tree biomass, kg/tree, in each biomass component
        """
    
        #----- Model parameters  by species and biomass components -----------------
        para = {'pine':{
                       'ag':{'Eq':0,'b0':-3.198, 'b1': 9.574, 'b2': 3.241, 'gamma':12.0, 'eta':20.0},  #  Repola 2009, Table 5 page 632
                     'bark':{'Eq':1,'b0':-4.548, 'b1': 7.997 , 'b2': 0.357, 'gamma':12.0, 'eta':0.0},
              'branch_dead':{'Eq':0,'b0':-5.201, 'b1': 10.574,'b2': 0.0  , 'gamma':16.0, 'eta':0.0},
            'branch_living':{'Eq':0,'b0':-6.162, 'b1': 15.075, 'b2':-2.618,'gamma':12.0, 'eta':12.0},
                     'leaf':{'Eq':0,'b0':-6.303, 'b1': 14.472, 'b2':-3.976,'gamma':6.0, 'eta':1.0},
                    'roots':{'Eq':0,'b0':-5.550, 'b1': 13.408, 'b2': 0.0,  'gamma':15.0, 'eta':1.0},
                 'stemwood':{'Eq':0,'b0':-3.721, 'b1': 8.103, 'b2': 5.066, 'gamma':14.0, 'eta':12.0},
                    'stump':{'Eq':0,'b0':-6.753, 'b1': 12.681, 'b2': 0.0 , 'gamma':12.0, 'eta':1.0}
                       },
                'spruce':{
                       'ag':{'Eq':1,'b0':-1.808, 'b1': 9.482, 'b2': 0.469, 'gamma':20.0, 'eta':0.0},  #  Repola 2009, Table 7 page 634
                     'bark':{'Eq':1,'b0':-4.548, 'b1': 9.448, 'b2': 0.436, 'gamma':18.0, 'eta':0.0},
              'branch_dead':{'Eq':1,'b0':-4.850, 'b1': 7.702, 'b2':0.513 , 'gamma':18.0, 'eta':0.0},
            'branch_living':{'Eq':0,'b0':-4.214, 'b1':14.508 , 'b2': -3.277, 'gamma':13.0, 'eta':5.0},
                     'leaf':{'Eq':0,'b0':-2.994, 'b1': 12.251, 'b2': -3.415, 'gamma':10.0, 'eta':1.0},
                    'roots':{'Eq':0,'b0':-2.294, 'b1': 10.646, 'b2': 0.0, 'gamma':24.0, 'eta':0.0},
                 'stemwood':{'Eq':2,'b0':-3.555, 'b1': 8.042, 'b2': 0.869, 'b3':0.015 , 'gamma':14.0, 'eta':0.0},
                    'stump':{'Eq':0,'b0':-3.964, 'b1': 11.730, 'b2': 0.0, 'gamma':26.0, 'eta':0.0}
                       },
                'birch':{
                       'ag':{'Eq':0,'b0':-3.654, 'b1': 10.582, 'b2': 3.018, 'gamma':12.0, 'eta':22.0},  #  Repola 2008, Table 4 page 613
                     'bark':{'Eq':0,'b0':-5.401, 'b1':10.061, 'b2': 2.657,  'gamma':12.0, 'eta':20.0},
              'branch_dead':{'Eq':0,'b0':-8.335, 'b1': 12.402,'b2':0.0   ,  'gamma':16.0, 'eta':0.0},
            'branch_living':{'Eq':0,'b0':-4.152, 'b1': 15.874, 'b2': -4.407,'gamma':16.0, 'eta':10.0},
                     'leaf':{'Eq':0,'b0':-29.566, 'b1': 33.372, 'b2': 0.0,  'gamma':2.0, 'eta':0.0},
                    'roots':{'Eq':1,'b0':-3.223, 'b1': 6.497, 'b2': 1.033,  'gamma':22.0, 'eta':0.0},
                 'stemwood':{'Eq':1,'b0':-4.879, 'b1': 9.651, 'b2':1.012 ,  'gamma':12.0, 'eta':0.0},
                    'stump':{'Eq':0,'b0':-3.574, 'b1':11.304 , 'b2': 0.0,   'gamma':26.0, 'eta':1.0}
                       },
                }
    
        icomponents = ['leaf', 'ag', 'branch_dead', 'branch_living', 'stump', 'roots']                     # Selected biomass components
        out = {}                                                                    # total biomass over the diameter distribution
        outs = {}                                                                   # biomass in each diameter class
        dski=2.+1.25*D                                                              # approximate stump diameter
        for component in icomponents:
            out[component] = {}
            if para[species][component]['Eq']==0:
                lnmass = para[species][component]['b0'] + \
                        para[species][component]['b1'] *(dski)/(dski+ para[species][component]['gamma']) + \
                        para[species][component]['b2']*(h/(h+para[species][component]['eta']))                 # biomass equation
            elif para[species][component]['Eq']==1:
                lnmass = para[species][component]['b0'] + \
                        para[species][component]['b1'] *(dski)/(dski+ para[species][component]['gamma']) + \
                        para[species][component]['b2']*(np.log(h))
            elif para[species][component]['Eq']==2:
                lnmass = para[species][component]['b0'] + \
                        para[species][component]['b1'] *(dski)/(dski+ para[species][component]['gamma']) + \
                        para[species][component]['b2']*(np.log(h)) + \
                        para[species][component]['b3']*h
    
            mi = np.exp(lnmass)*Nd                                                  # return from ln to linear
            out[component] = sum(mi)                                                # this produces kg/ha
            outs[component] =  np.exp(lnmass)                                       # biomass / tree kg
        return out, outs                                                            #

    def getAssortmentVolumes(self, D, h, speciesid):
        """
        input:
            D - diameter classes, cm, array 1...50
            h - stem height,m, in each diameter class
            speciesid in each diameter class 1 - pine, 2 - spruce, 3 - deciduous
        """
        assortmentVolumes = {'log':np.zeros(len(D)), 
                             'pulp':np.zeros(len(D)), 
                             'residual':np.zeros(len(D)), 
                             'total':np.zeros(len(D))}
        speciesid = np.clip(speciesid,1,3)                                     # all speciesid >3 are changed to 3 "other deciduous"
        for i, (d, hh, sp) in enumerate(zip(D.values,h.values,speciesid.values)):
            assortments = StemCurve().predictAssortmentVolumes(d, hh, sp)
            assortmentVolumes['log'][i] = np.around(assortments['log'], decimals=2)
            assortmentVolumes['pulp'][i] = np.around(assortments['pulp'], decimals=2)
            assortmentVolumes['residual'][i] = np.around(assortments['residual'], decimals=2)
            assortmentVolumes['total'][i] = np.around(assortments['total'], decimals=2)
                
        return assortmentVolumes


# Survival model (Pukkala et al. 2021)
    def get_survival(self, speciesid, D, BAL_Total, BAL_Pine, BAL_Spruce, BAL_S_B, Peat, Aspen, Birch):
        
        # Parameters (Pukkala et al. 2021):
        S_param = { 'Intercept':[1.41223, 5.01677, 1.60895], 'sqrt_d':[1.8852, 0.36902, 0.71578], 'd':[-0.21317, -0.07504, -0.08236], \
                    'BAL_Total':[-0.25637, 0, 0], 'BAL_Pine':[0, 0, -0.04814], 'BAL_Spruce':[0, -0.2319, 0], 'BAL_Spruce_Broadleaf':[0, 0, -0.13481], \
                    'Peat':[-0.39878, -0.47361, -0.31789], 'Aspen':[0, 0, 0.56311], 'Birch':[0, 0, 1.40145] }
    
        fx = np.zeros(len(D))
        for i in range(len(D)):
            # Reclassify species code to represent species-index 'spi': 0=pine, 1=spruce, 2=other
            if speciesid[i] <= 2:
                spi = speciesid[i] - 1
            else:
                spi = 2
            fx[i] = S_param['Intercept'][spi] + S_param['sqrt_d'][spi] * np.sqrt(D[i]) + S_param['d'][spi] * D[i] + \
                    S_param['BAL_Total'][spi] * (BAL_Total[i]/np.sqrt(D[i]+1)) + S_param['BAL_Pine'][spi] * (BAL_Pine[i]/np.sqrt(D[i]+1)) + \
                    S_param['BAL_Spruce'][spi] * (BAL_Spruce[i]/np.sqrt(D[i]+1)) + S_param['BAL_Spruce_Broadleaf'][spi] * (BAL_S_B[i]/np.sqrt(D[i]+1)) + \
                    S_param['Peat'][spi] * Peat[i] + S_param['Aspen'][spi] * Aspen[i] + S_param['Birch'][spi] * Birch[i]
        survival = 1 / (1 + np.exp(-fx))
        return survival
    
    # Predict survival rate:
    def predict_survival_5_years(self, Nd, speciesid, D, Peat=1):
        """
        Pukkala et al. 2021. https://doi.org/10.1093/forestry/cpab008
        input:
            Nd - number of trees in each diameter class, array
            speciesid - tree species as integer value
            D - diameter array 1...50 cm
        """
        Peat = np.repeat(Peat, len(D))
        Aspen = np.zeros(len(D))
        Birch = np.zeros(len(D))
        for i in range(len(D)):
            if speciesid[i]==5:
                Aspen[i] = 1
            if (speciesid[i]==3) | (speciesid[i]==4):
                Birch[i] = 1
    
        G_DBH_class = Nd * np.pi*(D/2)**2 / 10000
        G_plot = np.repeat(np.sum(G_DBH_class), len(D))
    
        BAL_Total = np.zeros(len(D))
        BAL_Pine = np.zeros(len(D))
        BAL_Spruce = np.zeros(len(D))
        BAL_S_B = np.zeros(len(D))
        for i in range(len(D)):
            BAL_Total[i] = np.sum(G_DBH_class[D > D[i]])
            BAL_Pine[i] = np.sum(G_DBH_class[(D > D[i]) & (speciesid ==1)])
            BAL_Spruce[i] = np.sum(G_DBH_class[(D > D[i]) & (speciesid == 2)])
            BAL_S_B[i] = np.sum(G_DBH_class[(D > D[i]) & (speciesid != 1)])
    
        return np.around(self.get_survival(speciesid, D, BAL_Total, BAL_Pine, BAL_Spruce, BAL_S_B, Peat, Aspen, Birch), decimals=2)
    
    

    # Diameter increment model (Pukkala et al. 2021):
    def get_diameter_increment(self, initialDiameter, speciesid, G_plot, BAL_Total, BAL_Spruce, BAL_S_B, TS, SiteType, Peat, Pendula_or_Aspen):
        
        # Parameters (Pukkala et al. 2021):
        D_param = { 'Intercept':[-7.1552, -12.7527, -8.6306], 'sqrt_d':[0.4415, 0.1693, 0.5097], 'd':[-0.0685, -0.0301, -0.0829], \
                    'ln_G_1':[-0.2027, -0.1875, -0.3864], 'BAL_Total':[-0.1236, -0.0563, 0], 'BAL_Spruce':[0, -0.0870, 0], \
                    'BAL_Spruce_Broadleaf':[0, 0, -0.0545], 'ln_TS':[1.1198, 1.9747, 1.3163], 'Peat':[-0.2425, 0, 0], 'd_Pendula_or_Aspen':[0, 0, 0.0253], \
                    'Fertility':{'Herb-rich':[0.1438,0.2688,0.2566], 'Mesic':[0,0,0], 'Sub-xeric':[-0.1754,-0.2145,-0.2256], 'Xeric':[-0.5163,-0.6179,-0.3237]} }
    
        ln_D_increment = np.zeros(len(initialDiameter))
        for i in range(len(initialDiameter)):
            # Reclassify species code to represent species-index 'spi': 0=pine, 1=spruce, 2=other
            if speciesid[i] <= 2:
                spi = speciesid[i] - 1
            else:
                spi = 2
            ln_D_increment[i] = D_param['Intercept'][spi] + D_param['sqrt_d'][spi] * np.sqrt(initialDiameter[i]) + D_param['d'][spi] * initialDiameter[i] + \
                                D_param['ln_G_1'][spi] * np.log(G_plot[i]+1) + D_param['BAL_Total'][spi] * (BAL_Total[i]/np.sqrt(initialDiameter[i]+1)) + \
                                D_param['BAL_Spruce'][spi] * (BAL_Spruce[i]/np.sqrt(initialDiameter[i]+1)) + \
                                D_param['BAL_Spruce_Broadleaf'][spi] * (BAL_S_B[i]/np.sqrt(initialDiameter[i]+1)) + D_param['ln_TS'][spi] * np.log(TS[i]) + \
                                D_param['Fertility'][SiteType[i]][spi] + D_param['Peat'][spi] * Peat[i] + \
                                D_param['d_Pendula_or_Aspen'][spi] * Pendula_or_Aspen[i] * initialDiameter[i]
        return np.exp(ln_D_increment)
    
    # Predict diameter increment for next 5 years:
    def predict_diameter_increment_5_years(self, Nd, speciesid, D, fertilityClass, temperatureSum=1200, peat=1):
        """
        Timo Pukkala and others, Self-learning growth simulator for modelling forest stand dynamics in
        changing conditions, Forestry: An International Journal of Forest Research, Volume 94, Issue 3, 
        July 2021, Pages 333–346, https://doi.org/10.1093/forestry/cpab008
        input:
            Nd - number of trees in each diameter class, array
            speciesid - tree species as integer
            D - diamtere class array 1...50 cm
            fertilityClass - fertility class expressed as string
            temperatureSum - temperature sum as degree-days
            peat - dummy variable assigning the calculations to peatlands
        """

        siteTypeDict = {1: 'Herb-rich', 2: 'Herb-rich', 3: 'Mesic', 4: 'Sub-xeric', 5: 'Xeric', 6: 'Xeric'}
        siteType = np.repeat(siteTypeDict[fertilityClass], len(D))
        TS = np.repeat(temperatureSum, len(D))
        peat = np.repeat(peat, len(D))
    
        Pendula_or_Aspen = np.zeros(len(D))
        for i in range(len(D)):
            if (speciesid[i]==3) | (speciesid[i]==5):
                Pendula_or_Aspen[i] = 1
    
        G_DBH_class = Nd * np.pi*(D/2)**2 / 10000
        G_plot = np.repeat(np.sum(G_DBH_class), len(D))
    
        BAL_Total = np.zeros(len(D))
        BAL_Spruce = np.zeros(len(D))
        BAL_S_B = np.zeros(len(D))
        for i in range(len(D)):
            BAL_Total[i] = np.sum(G_DBH_class[D > D[i]])
            BAL_Spruce[i] = np.sum(G_DBH_class[(D > D[i]) & (speciesid == 2)])
            BAL_S_B[i] = np.sum(G_DBH_class[(D > D[i]) & (speciesid != 1)])
    
        return np.around(self.get_diameter_increment(D, speciesid, G_plot, BAL_Total, BAL_Spruce, BAL_S_B, TS, siteType, peat, Pendula_or_Aspen), decimals=2)

    def dominant_height(self, Nd, h):
        """
        input:
            Nd  - number of trees in each diameter class
            h  - height of trees, m, in each diameter class
        output:
            dominant height, m, mean diameter of 100 thickest trees
        """
        tmp = np.cumsum(Nd.values)[-1]-np.cumsum(Nd.values)                   # cumulative number of trees from thickest to thinnest dimater
        ix = np.where(tmp<100)                                                # find 100 thickest trees
        return np.around(sum((h.values[ix] * Nd.values[ix]))/sum(Nd.values[ix]), decimals=2)   # weighed average of thick diameters (weighed with number of trees)

    
    def forest_data_to_susi(self, distrib, species, age_ini = 35):
        """
        input:
            distrib - dictionary containing the stand development variables
            species - tree species as string
        """
        keysN = ['N0', 'N5', 'N10', 'N15', 'N20', 'N25', 'N30', 'N35', 'N40', 'N45', 'N50']
        keysD = ['D0', 'D5', 'D10', 'D15', 'D20', 'D25', 'D30', 'D35', 'D40', 'D45', 'D50']
        keysH = ['H0', 'H5', 'H10', 'H15', 'H20', 'H25', 'H30', 'H35', 'H40', 'H45', 'H50']
        keysV = ['VOL0', 'VOL5', 'VOL10', 'VOL15', 'VOL20', 'VOL25', 'VOL30', 'VOL35', 'VOL40', 'VOL45', 'VOL50']
        keysL = ['LOG0', 'LOG5', 'LOG10', 'LOG15', 'LOG20', 'LOG25', 'LOG30', 'LOG35', 'LOG40', 'LOG45', 'LOG50']
        keysP = ['PULP0', 'PULP5', 'PULP10', 'PULP15', 'PULP20', 'PULP25', 'PULP30', 'PULP35', 'PULP40', 'PULP45', 'PULP50']
        biom_components = ['leaf', 'ag', 'branch_dead', 'branch_living', 'stump', 'roots']
        biom_out = ['lehdet', 'ag', 'kuolleet oksat', 'elävät oksat', 'Kannot', 'Juuret >2mm']

        susi_input = pd.DataFrame(data = {'Kasvatus' : np.repeat(1, 11), 'Vuosi' : np.arange(0,51,5)})
        
        susi_input['Ikä'] = np.arange(age_ini,age_ini+51,5)        
        susi_input['N']   = [int(np.sum(distrib[kn])) for kn in keysN]    
        susi_input['PPA'] = [np.sum(distrib[kn] * np.pi * (distrib[kd]/2/100)**2) for kn,kd in zip(keysN, keysD)] 
        susi_input['Hg'] = [np.sum(distrib[kn] * distrib[kh] * distrib[kd]**2) / np.sum(distrib[kn] * distrib[kd]**2) for kn,kh,kd in zip(keysN, keysH, keysD)]     
        susi_input['Dg'] = [np.sum(distrib[kn] * distrib[kd]**3) / np.sum(distrib[kn] * distrib[kd]**2) for kn,kd in zip(keysN, keysD)]     
        susi_input['Hdom'] = [self.dominant_height(distrib[kn], distrib[kh]) for kn,kh in zip(keysN, keysH)]
        susi_input['Tilavuus'] = [np.sum(distrib[kn] * distrib[kv] / 1000) for kn,kv in zip(keysN, keysV)]
        susi_input['Tukki'] = [np.sum(distrib[kn] * distrib[kl] / 1000) for kn,kl in zip(keysN,keysL)]
        susi_input['Kuitu'] = [np.sum(distrib[kn] * distrib[kp] / 1000) for kn,kp in zip(keysN, keysP)]
        susi_input['Hukka'] = susi_input['Tilavuus'] - susi_input['Tukki'] - susi_input['Kuitu']
 
        mort_share = -np.gradient(susi_input['N'])/susi_input['N']        
        vs= susi_input['Tilavuus'].values
        susi_input['Tuotos'] = np.cumsum(mort_share * (vs + np.gradient(vs)/2)) + susi_input['Tilavuus']
        susi_input['Kuolleisuus'] = np.cumsum(mort_share * (vs + np.gradient(vs)/2))
    
        for bo,bc in zip(biom_out, biom_components):
            susi_input[bo] = [self.get_biomass(distrib[kd], distrib[kh], distrib[kn], species)[0][bc] for kd,kh,kn in zip(keysD, keysH, keysN)]
            susi_input[bo] = susi_input[bo] / 1000.
    
        susi_input['runko(hukka)'] = susi_input['Hukka'] *420. /1000.  # vol times density of wood -> to tons/ha
        susi_input['runko(aines)'] = susi_input['ag'] - susi_input['elävät oksat'] - susi_input['kuolleet oksat'] - susi_input['Kannot']       
        susi_input['Hienojuuret'] = susi_input['Juuret >2mm'] * 0.02
        
        susi_input = susi_input[['Kasvatus', 'Vuosi', 'Ikä', 'N', 'PPA', 'Hg', 'Dg', 'Hdom', 'Tilavuus', 
                                'Tukki', 'Kuitu', 'Hukka','Tuotos', 'Kuolleisuus', 'runko(aines)', 
                                'runko(hukka)', 'elävät oksat', 'kuolleet oksat', 'lehdet', 
                                'Kannot',  'Juuret >2mm', 'Hienojuuret']]
    
        return susi_input

    def create_development(self, species, D, N, W, G,Hdom, DDY, fertility_class):
        
        keysN = ['N0', 'N5', 'N10', 'N15', 'N20', 'N25', 'N30', 'N35', 'N40', 'N45', 'N50']
        keysD = ['D0', 'D5', 'D10', 'D15', 'D20', 'D25', 'D30', 'D35', 'D40', 'D45', 'D50']
        keysH = ['H0', 'H5', 'H10', 'H15', 'H20', 'H25', 'H30', 'H35', 'H40', 'H45', 'H50']
        keysV = ['VOL0', 'VOL5', 'VOL10', 'VOL15', 'VOL20', 'VOL25', 'VOL30', 'VOL35', 'VOL40', 'VOL45', 'VOL50']
        keysL = ['LOG0', 'LOG5', 'LOG10', 'LOG15', 'LOG20', 'LOG25', 'LOG30', 'LOG35', 'LOG40', 'LOG45', 'LOG50']
        keysP = ['PULP0', 'PULP5', 'PULP10', 'PULP15', 'PULP20', 'PULP25', 'PULP30', 'PULP35', 'PULP40', 'PULP45', 'PULP50']

        if species == "pine":
            sp_code = 1
        if species == "spruce":
            sp_code = 2
        
        DBH_classes = pd.DataFrame(data = {'speciesid' : sp_code, 'D0' : D, 'N0' : N*W})
        
        DBH_classes['H0'] = self.naslund_height(DDY, Hdom, D, G, species)
        
        #Initial values        
        # ['log', 'pulp', 'residual', 'total']
        assortments = self.getAssortmentVolumes(DBH_classes['D0'], DBH_classes['H0'], DBH_classes['speciesid'])
        DBH_classes['VOL0'] = assortments['total'] 
        DBH_classes['LOG0'] = assortments['log'] 
        DBH_classes['PULP0'] = assortments['pulp']  

        # nn and nd refer to previous timestep N nd D, kn,kd... to the ongoing timestep        
        for nn,nd,kn,kd,kh,kv,kl,kp in zip(keysN[:-1], keysD[:-1],keysN[1:],keysD[1:],keysH[1:], keysV[1:],keysL[1:],keysP[1:]):            
            DBH_classes[kn] = np.around(DBH_classes[nn] * self.predict_survival_5_years(DBH_classes[nn], DBH_classes['speciesid'], DBH_classes[nd]), decimals=2)
            DBH_classes[kd] = DBH_classes[nd] + self.predict_diameter_increment_5_years(DBH_classes[nn], DBH_classes['speciesid'], DBH_classes[nd], fertility_class, DDY)
            Hdom = self.naslund_to_hdom(DBH_classes[kd], DBH_classes[kn], species)
            
            import matplotlib.pyplot as plt
            plt.plot(DBH_classes[kd], DBH_classes[kn], label=kn)
            plt.legend(loc='upper right')
            print ('mean para', kh, Hdom)
            
            DBH_classes[kh] = self.naslund_height(DDY, Hdom, DBH_classes[kd], np.sum(DBH_classes[kn] * np.pi * np.square(DBH_classes[kd]/2./100.)), species)
            print ('specific', kh, self.dominant_height(DBH_classes[kn], DBH_classes[kh]))
            assortments = self.getAssortmentVolumes(DBH_classes[kd], DBH_classes[kh], DBH_classes['speciesid'])
            DBH_classes[kv] = assortments['total']
            DBH_classes[kl] = assortments['log']        
            DBH_classes[kp] = assortments['pulp']
            
            
        return DBH_classes

    def get_table(self, age_ini):
        DBH_classes =self.create_development(self.species, self.D, self.N, self.W, self.G, self.Hdom, self.DDY, self.fertility_class)
        susi_input = self.forest_data_to_susi(DBH_classes, self.species, age_ini)
        return susi_input

    
    def recover_weibull(self,d,N,G):
        def func(x,N,G):
            return ((sum((x[1]/x[0] * (self.D/x[0])**(x[1]-1) * np.exp(-(self.D/x[0])**x[1]))*N*np.pi*np.square(self.D/200.0))-G)**2)
        
        root = minimize(func,[d, 3.6], args = (N,G), bounds=((0.1, None),(0.1,None)))
      
        return root.x

"""
G = 16.7                   # basal area, m2/ha
N = 920                    # stem number, /ha
Dg =  20.0                  # mean diameter, cm
Hg = 16.6                  # mean height, m
Hdom = 18.7                # dominant height, m
DDY = 1100                 # temperature sum, degree days
species = 'pine'


#herb = 0
#mesic = 0
#subxeric = 1
#xeric = 0
peat = 1
dxpendula = 0
aspen = 0
birch = 0

FertilityClass = 4
fertility_class = 4


c = 2.051
b = 15.281

D = np.linspace(1,50, num=50)
W = c / b * (D/b)**(c-1) * np.exp(-(D/b)**c)
Nd = W * N


gy = Growth_and_Yield_Table(G, N, Dg, D, Hg, Hdom, DDY, species, fertility_class, W = W, Nd = Nd)
age_ini = 35

b,c =  gy.recover_weibull(Dg,N,G)
print(b,c)

susi_input = gy.get_table(age_ini)
susi_input.to_clipboard()


writer = pd.ExcelWriter(r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_input.xlsx', engine='xlsxwriter')
susi_input.to_excel(writer, sheet_name='StandData', index=False)

page2 = {'pine': {'idp': [1], 'puu': 'mänty'},
         'spruce':{'idp':[2], 'puu': 'kuusi'},
         'birch':{'idp':[3], 'puu': 'koivu'}
         }

dicpage2 = {'Kasvatus': 1, 'Vuosi': 0, 'Harvennus': 'Päätehakkuu', 
            'id Puulaji': page2[species]['idp'], 'Puulaji':page2[species]['puu']}
dfpage2 = pd.DataFrame(data=dicpage2) 

dfpage2.to_excel(writer, sheet_name='Kertymät')
writer.save()
"""