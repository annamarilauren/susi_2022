# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 17:37:43 2022

@author: alauren
"""

import numpy as np
from scipy.interpolate import interp1d
from susi.allometry import Allometry
from susi.susi_utils import  assimilation_yr

class Canopylayer():
    
    """
    UNITS: all units in /tree basis, except number of trees in the canopy layer, which is in /ha
    Canopylayer keeps track on the development of biomass components within the 
    different layers of canopy. Canopy layer is of single tree species and homegeneous in age. It is 
    initialized with growth and yield simulator outputfile that describes development
    of biomass components in time. Canopy layers is array shaped with dimensions 
    number of columns along the strip.
    """
    
    def __init__(self, name, nscens, yrs, ncols, nlyrs, sfc, agearr, mottipath, mottifile, ixs, photopara, nut_stat):
        self.name = name                                                       # name of the canopy layer e.g. 'dominant', 'subdominant', etc.
        self.nlyrs = nlyrs                                                     # number of different canopy layers along the strip
        self.ixs = ixs                                                         # indices for the location of the different canopy layers along the strip 
        self.ncols = ncols                                                     # number of columns in the strip
        self.agearr = agearr.copy()                                            # age of the canopy layer, yrs  
        self.photopara = photopara                                             # photosynthesis parameters for assimilation model (Mäkelä et al. 2008)
        self.nscens = nscens                                                   # number of scenarion in the simulation  
        self.yrs = yrs                                                         # number od years in the simulation
        self.remaining_share = np.ones(self.ncols)                             # share of remaining stems after thinning 0...1
        self.sfc = sfc.copy()
        self.tree_species = np.zeros(self.ncols, dtype=np.int8)                               # tree species 1 Scots pine, 2 Norway spruce
        
        # -------- Biomass interpolation functions-------------------
        self.allodic ={}                                                       # dictionary to contain all allometric functions
        for ncanopy in nlyrs:                                                  # numner of different allometric files along the strip in this canopy layer
            if ncanopy > 0:                                                    # zero indicates no tree in the layer
                self.sfc = int(np.median(self.sfc[self.ixs[ncanopy]]))         # site fertility class
                self.allodic[ncanopy] = Allometry()                            # allometry instance to the dictionary
                mfile = mottipath + mottifile[ncanopy]                         # mottifile where the allometry tables exist
                self.allodic[ncanopy].motti_development(mfile, self.sfc)       # run the allometry; interpolation functions in the instance
                self.tree_species[self.ixs[ncanopy]] = int(self.allodic[ncanopy].sp)
        self.initialize_domain(agearr, nut_stat)                               # create variables and set initial values
                
        print (self.name, 'initialized' )
        
    def initialize_domain(self, agearr, nut_stat):
        self.agearr = agearr.copy()
        self.remaining_share = np.ones(self.ncols)                             # share of remaining stems after thinning 0...1, in initialization should be one
        nlyrs = self.nlyrs                                                     # number of different canopy layers along the strip
        ixs = self.ixs                                                         # indices for the canopy layers along the strip
        
        ncols = self.ncols                                                     # number of columns along the strip
        self.stems = np.zeros(ncols, dtype=float)                              # stocking number of trees per ha

        self.basalarea = np.zeros(ncols, dtype=float)                          # basal area in the canopy layer m2/tree
        self.biomass = np.zeros(ncols, dtype=float)                            # total canopy layer biomass excluding leaves kg/tree
        self.n_demand = np.zeros(ncols, dtype=float)                           # N demand exluding leaves kg/tree/yr
        self.p_demand = np.zeros(ncols, dtype=float)                           # P demand exluding leaves kg/tree/yr 
        self.k_demand = np.zeros(ncols, dtype=float)                           # K demand exluding leaves kg/tree/yr
        self.hdom = np.zeros(ncols, dtype=float)                               # domainant height m in the canopy layer   
        self.leafarea = np.zeros(ncols, dtype=float)                           # leaf area in the canopy layer m2 m-2 tree-1
        self.lai_above = np.zeros(ncols, dtype=float)                          # one sided leaf area above this canopy layer, m2 m-2 tree-1 
        self.leafmass = np.zeros(ncols, dtype=float)                           # basic leaf mass kg/tree 
        self.logvolume = np.zeros(ncols, dtype=float)                          # volume of the saw logs m3/tree
        self.finerootlitter = np.zeros(ncols, dtype=float)                     # fineroot litter in the canopy layer kg/tree/yr
        self.n_finerootlitter = -np.zeros(ncols, dtype=float)                  # N in fineroot litter in the canopy layer kg/treea/yr 
        self.p_finerootlitter = np.zeros(ncols, dtype=float)                   # P in fineroot litter in the canopy layer kg/tree/yr      
        self.k_finerootlitter = np.zeros(ncols, dtype=float)                   # K in fineroot litter in the canopy layer kg/tree/yr  
        self.pulpvolume= np.zeros(ncols, dtype=float)                          # volume of pulpwood in the canopy layer m3/tree

        self.NPP= np.zeros(ncols, dtype=float)                                 # net primary production kg/tree/yr dry matter
        self.NPP_pot= np.zeros(ncols, dtype=float)                             # potential net primary production kg/tree/yr dry matter

        self.nonwoodylitter = np.zeros(ncols, dtype=float)                     # nonwoody litter kg/tree/yr
        self.n_nonwoodylitter = np.zeros(ncols, dtype=float)                   # N in nonwoody litter kg/tree/yr
        self.p_nonwoodylitter = np.zeros(ncols, dtype=float)                   # P in nonwoody litter kg/tree/yr
        self.k_nonwoodylitter = np.zeros(ncols, dtype=float)                   # K in nonwoody litter kg/tree/yr       
        self.species = np.zeros(ncols, dtype=int)                              # tree species: 1 Scots pine, 2: Norway spruce, 3: Betula pendula
        self.volume = np.zeros(ncols, dtype=float)                             # volume of the growing stock m3/tree
        self.volumegrowth = np.zeros(ncols, dtype=float)                       # volume growth of the growing stock m3/tree/yr
        self.woodylitter = np.zeros(ncols, dtype=float)                        # woody litter kg/tree/yr
        self.n_woodylitter = np.zeros(ncols, dtype=float)                      # N in woody litter kg/tree/yr
        self.p_woodylitter = np.zeros(ncols, dtype=float)                      # P in woody litter kg/tree/yr
        self.k_woodylitter = np.zeros(ncols, dtype=float)                      # K in woody litter kg/tree/yr
        self.yi = np.zeros(ncols, dtype=float)                                 # yied. here same as voolume
        
        self.nonwoody_lresid = np.zeros(ncols, dtype=float)                    # nonwoody logging residues kg/tree
        self.n_nonwoody_lresid = np.zeros(ncols, dtype=float)                  # N in nonwoody logging residues kg/tree      
        self.p_nonwoody_lresid = np.zeros(ncols, dtype=float)                  # P in nonwoody logging residues kg/tree  
        self.k_nonwoody_lresid = np.zeros(ncols, dtype=float)                  # K in nonwoody logging residues kg/tree  
        
        self.woody_lresid = np.zeros(ncols, dtype=float)                       # woody logging residues kg/tree
        self.n_woody_lresid = np.zeros(ncols, dtype=float)                     # N in woody logging residues kg/tree
        self.p_woody_lresid = np.zeros(ncols, dtype=float)                     # P in woody logging residues kg/tree
        self.k_woody_lresid = np.zeros(ncols, dtype=float)                     # K in woody logging residues kg/tree
        
        #----------- ARRANGE Leaf dynamics variables-----------------
        self.new_lmass = np.zeros(ncols, dtype=float)                          # leaf mass after dynamics computation kg/ha
        self.leaf_litter = np.zeros(ncols, dtype=float)                        # leaf litter kg/tree/yr
        self.C_consumption = np.zeros(ncols, dtype=float)                      # check the unit: C or mass? C consumption in leaf processes
        self.leafmax = np.zeros(ncols, dtype=float)                            # upper limit for leaf mass kg/tree
        self.leafmin = np.zeros(ncols, dtype=float)                            # lower limit for leaf mass kg/tree
        self.Nleafdemand = np.zeros(ncols, dtype=float)                        # N demand of leaf production kg/ha/tree
        self.Nleaf_litter = np.zeros(ncols, dtype=float)                       # N in litterfall kg/ha/tree
        self.N_leaf = np.zeros(ncols, dtype=float)                             # N in leaves kg/tree
        self.Pleafdemand = np.zeros(ncols, dtype=float)                        # P demand of leaf production kg/ha/tree
        self.Pleaf_litter= np.zeros(ncols, dtype=float)                        # P in litterfall kg/ha/tree
        self.P_leaf = np.zeros(ncols, dtype=float)                             # P in leaves kg/tree
        self.Kleafdemand = np.zeros(ncols, dtype=float)                        # K demand in leaf production kg/ha/tree 
        self.Kleaf_litter  = np.zeros(ncols, dtype=float)                      # K in litterfall kg/ha/tree
        self.K_leaf = np.zeros(ncols, dtype=float)                             # K in leaves kg/tree
        self.basNdemand = np.zeros(ncols, dtype=float)                         # basic N demand kg/tree, in table growth conditions, used in nutrient status calculation
        self.basPdemand = np.zeros(ncols, dtype=float)                         # basic P demand kg/tree, in table growth conditions, used in nutrient status calculation
        self.basKdemand = np.zeros(ncols, dtype=float)                         # basic K demand kg/tree, in table growth conditions, used in nutrient status calculation
        
        # ---------- Initial values from age -------------------------------------------------------
        for m in nlyrs:
            if m > 0:
                self.biomass[ixs[m]] = self.allodic[m].allometry_f['ageToBm'](self.agearr[ixs[m]])        # stem biomass from age kg/tree         
                self.stems[ixs[m]] = self.allodic[m].allometry_f['bmToStems'](self.biomass[ixs[m]])*self.remaining_share[ixs[m]]  # number of stems /ha from the biomass of tree       

                self.basalarea[ixs[m]] = self.allodic[m].allometry_f['ageToBa'](self.agearr[ixs[m]])      # stem basal area from age m2/tree  
                self.hdom[ixs[m]] = self.allodic[m].allometry_f['ageToHdom'](self.agearr[ixs[m]])         # dominant height m
                self.leafarea[ixs[m]] = self.allodic[m].allometry_f['bmToLAI'](self.biomass[ixs[m]])*nut_stat[ixs[m]]  # one sided LAI from the stem biomass M2/m2/tree
                self.leafmass[ixs[m]] = self.allodic[m].allometry_f['ageToLeaves'](self.agearr[ixs[m]])*nut_stat[ixs[m]]
                self.species[ixs[m]] = self.allodic[m].sp
                #self.volume[ixs[m]] = self.allodic[m].allometry_f['ageToVol'](self.agearr[ixs[m]])         # stem volume m3/tree
                self.volume[ixs[m]] = self.allodic[m].allometry_f['bmToVol'](self.biomass[ixs[m]])         # stem volume m3/tree
                
                self.basNdemand[ixs[m]] = self.allodic[m].allometry_f['bmToNLeafDemand'](self.biomass[ixs[m]])
                self.basPdemand[ixs[m]] = self.allodic[m].allometry_f['bmToPLeafDemand'](self.biomass[ixs[m]])
                self.basKdemand[ixs[m]] = self.allodic[m].allometry_f['bmToKLeafDemand'](self.biomass[ixs[m]])

    def update(self, bm):
        
        # bm in kg in a mean stem
        """ CHANGE all units here into /tree, Do we need remaining share?"""
        #------------Update all variables with new biomass--------------------------------------
        ixs = self.ixs
        for m in self.nlyrs:
            if m > 0:
                #print ('**********************')
                #print (np.round(np.mean(self.allodic[m].allometry_f['bmToVol'](bm[ixs[m]])*self.stems),2))
                #print (np.round(np.mean(self.allodic[m].allometry_f['ageToVol'](self.agearr[ixs[m]])*self.stems), 2))
                #print ('vol')
                # print (self.volume)
                # print ('n stems')
                # print (self.stems)
                #print ('**********************')
                
                self.stems[ixs[m]] = self.allodic[m].allometry_f['bmToStems'](bm[ixs[m]])*self.remaining_share[ixs[m]]
                
                self.basalarea[ixs[m]] = self.allodic[m].allometry_f['bmToBa'](bm[ixs[m]]) 
                self.biomass[ixs[m]] = bm[ixs[m]] 
                self.hdom[ixs[m]] = self.allodic[m].allometry_f['bmToHdom'](bm[ixs[m]]) 
                self.leafarea[ixs[m]] = self.allodic[m].allometry_f['bmToLAI'](bm[ixs[m]]) 
                self.leafmass[ixs[m]] = self.allodic[m].allometry_f['bmToLeafMass'](bm[ixs[m]]) 
                #self.volume[ixs[m]] = self.allodic[m].allometry_f['bmToYi'](bm[ixs[m]]) 
                self.volume[ixs[m]] = self.allodic[m].allometry_f['bmToVol'](bm[ixs[m]]) 
                self.n_demand[ixs[m]] = self.allodic[m].allometry_f['bmToNdemand'](bm[ixs[m]]) 
                self.p_demand[ixs[m]] = self.allodic[m].allometry_f['bmToPdemand'](bm[ixs[m]]) 
                self.k_demand[ixs[m]] = self.allodic[m].allometry_f['bmToKdemand'](bm[ixs[m]]) 
                self.logvolume[ixs[m]] = self.allodic[m].allometry_f['volToLogs'](self.volume[ixs[m]])         
                self.finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToFinerootLitter'](bm[ixs[m]]) 
                self.n_finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToNFineRootLitter'](bm[ixs[m]]) 
                self.p_finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToPFineRootLitter'](bm[ixs[m]]) 
                self.k_finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToKFineRootLitter'](bm[ixs[m]]) 

                self.nonwoodylitter[ixs[m]] = self.finerootlitter[ixs[m]] + self.leaf_litter[ixs[m]]           # nonwoody litter kg/tree/yr
                self.n_nonwoodylitter[ixs[m]] = self.n_finerootlitter[ixs[m]] + self.Nleaf_litter[ixs[m]]     # N in nonwoody litter kg/tree/yr
                self.p_nonwoodylitter[ixs[m]] = self.p_finerootlitter[ixs[m]] + self.Pleaf_litter[ixs[m]]     # P in nonwoody litter kg/tree/yr
                self.k_nonwoodylitter[ixs[m]] = self.k_finerootlitter[ixs[m]] + self.Kleaf_litter[ixs[m]]     # K in nonwoody litter kg/tree/yr

                self.pulpvolume[ixs[m]] = self.allodic[m].allometry_f['volToPulp'](self.volume[ixs[m]])                         
                self.woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToWoodyLitter'](bm[ixs[m]]) 
                self.n_woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToNWoodyLitter'](bm[ixs[m]]) 
                self.p_woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToPWoodyLitter'](bm[ixs[m]]) 
                self.k_woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToKWoodyLitter'](bm[ixs[m]]) 
                self.yi[ixs[m]] = self.allodic[m].allometry_f['bmToYi'](bm[ixs[m]]) 

                self.basNdemand[ixs[m]] = self.allodic[m].allometry_f['bmToNLeafDemand'](bm[ixs[m]])
                self.basPdemand[ixs[m]] = self.allodic[m].allometry_f['bmToPLeafDemand'](bm[ixs[m]])
                self.basKdemand[ixs[m]] = self.allodic[m].allometry_f['bmToKLeafDemand'](bm[ixs[m]])
                self.agearr[ixs[m]] = self.agearr[ixs[m]] + 1  
                # print ('vol')
                # print (self.volume)
                # print ('n stems')
                # print (self.stems)

    def assimilate(self, forc, wt, afp, previous_nut_stat, nut_stat, lai_above):
        """
        Calls photosynthesis model (Mäkelä et al. 2008, standwise model) and leaf dynamics model that
        accounts for leaf mass, longevity and nutrient contents. This is a canopy model instance
        and the leaf dynamics is solved for similar columns along the strip. 

        Parameters
        ----------
        forc : TYPE pandas dataframe
            DESCRIPTION. daily weather variables in year-long df
        wt : TYPE pandas dataframe 
            DESCRIPTION. simulated water tables along the strip, shape: days, ncols
        afp : TYPE pandas dataframe 
            DESCRIPTION.air-filled porosity in rooting zone, shape: days, ncols
        previous_nut_stat : TYPE array
            DESCRIPTION. nutrient status alonmg the strip in the previous year
        nut_stat : TYPE array 
            DESCRIPTION. nutrient staus along the strip in the current year
        lai_above: TYPE: array
            DESCRIPTION, leaf area above the canopy layer incoming unit: m2 m-2
        Returns
        -------
        None.

        """
        
        """ Change unit of all variables to /tree """
        """ assimilation_yr function operates in /ha basis,"""
        lai_above = lai_above*2                                                # lai_above is updated in stand object        
        self.NPP, self.NPP_pot = assimilation_yr(self.photopara, forc, wt, afp,
                                       self.leafarea*2 * self.stems, lai_above)     # double sided LAI required                     

        self.NPP = self.NPP * nut_stat / self.stems  * 1.1                     # returned back to tree basis unit
        self.NPP_pot = self.NPP_pot * nut_stat / self.stems * 1.1               # returned back to tree basis units
        
        bm_increment = self.NPP
        bm = self.biomass

        current_leafmass = self.leafmass
        
        ixs = self.ixs
        for m in self.nlyrs:
            if m > 0:
                self.new_lmass[ixs[m]], self.leaf_litter[ixs[m]], self.C_consumption[ixs[m]],\
                    self.leafmax[ixs[m]], self.leafmin[ixs[m]], self.Nleafdemand[ixs[m]],\
                    self.Nleaf_litter[ixs[m]],self.N_leaf[ixs[m]],\
                    self.Pleafdemand[ixs[m]], self.Pleaf_litter[ixs[m]],self.P_leaf[ixs[m]],\
                    self.Kleafdemand[ixs[m]], self.Kleaf_litter[ixs[m]],self.K_leaf[ixs[m]],\
                        self.leafarea[ixs[m]] = self.leaf_dynamics(bm[ixs[m]],\
                                                bm_increment[ixs[m]], current_leafmass[ixs[m]], 
                                                previous_nut_stat[ixs[m]], nut_stat[ixs[m]],\
                                                self.agearr[ixs[m]], self.allodic[m].allometry_f,\
                                                self.species[ixs[m]], printOpt=False)
                
                self.finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToFinerootLitter'](bm[ixs[m]]) 
                self.woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToWoodyLitter'](bm[ixs[m]]) 
        """
        if self.name=='dominant':
            print ('ooooooooooooooooooooooooo')
            print (self.name, np.mean(self.agearr))
            print (np.round(np.mean(self.NPP),2), 'npp' )
            print (np.round(np.mean(self.C_consumption),2), 'c cons')
            print (np.round(np.mean(self.finerootlitter), 2),'fr litter')
            print (np.round(np.mean(self.woodylitter),2),'woody l')
       """    
        
        
        delta_bm_noleaves = self.NPP - self.C_consumption - self.finerootlitter  - self.woodylitter 
        self.leafmass = self.new_lmass
        vol_ini =  self.volume.copy() 
        """
        if self.name=='dominant': 
            print (np.round(np.mean(self.biomass),2), 'biomass ini' )
            print (np.round(np.mean(self.volume*self.stems),2), 'volume ini' )
            print (np.round(np.mean(self.stems),2), 'stems ini' )
            print (np.round(np.mean(self.allodic[1].allometry_f['bmToVol'](bm)*self.stems),2))
            #print (np.round(np.mean(self.allodic[1].allometry_f['ageToVol'](self.agearr)*self.stems), 2))
        """
        self.update(self.biomass + np.maximum(delta_bm_noleaves,0.0))
        
        #if self.name=='dominant': print (np.round(np.mean(delta_bm_noleaves),2), 'delta no leaves' )                
        
        self.volumegrowth = self.volume - vol_ini
        """
        if self.name=='dominant':
            print (np.round(np.mean(self.volumegrowth*self.stems),2), 'volumegrowth')
            print (np.round(np.mean(self.biomass),2), 'biomass after' )
            print (np.round(np.mean(self.volume*self.stems),2), 'volume after' )
            print (np.round(np.mean(self.stems),2), 'stems after' )
        """
        
    def leaf_dynamics(self, bm, bm_increment, current_leafmass, previous_nut_stat,\
                      nut_stat, agenow, allometry_f, species, printOpt=False):
        """
        input:
             bm, current biomass, array, kg/tree without leaves
             bm_increment, array, total NPP kg/tree
             incoming columns are of single tree species canopy layers 
        """
        #******** Parameters *****************
        nuts = {'Pine':{
                    'Foliage':{'N': [10.0,20.0], 'P': [1.0, 2.2], 'K': [5.0, 6.5]}},
                'Spruce':{
                    'Foliage':{'N': [10.0,20.0], 'P': [1.0, 2.2], 'K': [3.5, 6.0]}},
                'Birch':{
                    'Foliage':{'N': [10.0,20.0], 'P': [1.0, 2.2], 'K': [3.5, 6.0]}},
                        }
                
        retrans = {'N': 0.69, 'P': 0.73, 'K':0.8}                              # Nieminen Helmisaari 1996 Tree Phys
        sla= {'Pine': 6.8, 'Spruce': 7.25, 'Birch':14.0}                       # specific leaf area Härkönen et al. 2015 BER 20, 181-195      
        species_codes = {1:'Pine', 2: 'Spruce', 3: 'Birch'}
        spe = species_codes[species[0]]                                        # Take the first item, all are same here
        longevityLeaves = {'Pine':[3.5, 2.5], 'Spruce':[6.0, 4.0], 'Birch':[1.0, 1.0]}                     # yrs, life span of leaves and fine roots    
    
        N_con = interp1d(np.array([0.66,1.5]), np.array(nuts[spe]['Foliage']['N']), fill_value=tuple(nuts[spe]['Foliage']['N']), bounds_error=False)
        P_con = interp1d(np.array([0.66,1.5]), np.array(nuts[spe]['Foliage']['P']), fill_value=tuple(nuts[spe]['Foliage']['P']), bounds_error=False)
        K_con = interp1d(np.array([0.66,1.5]), np.array(nuts[spe]['Foliage']['K']), fill_value=tuple(nuts[spe]['Foliage']['K']), bounds_error=False)
        
        longevity = interp1d(np.array([0.66,1.5]), np.array(longevityLeaves[spe]), fill_value= tuple(longevityLeaves[spe]), bounds_error=False)
    
        n = len(nut_stat)
        #************ Biomass and litter**************
        """ Change units to /tree here"""
        leafbase0 =  allometry_f['bmToLeafMass'](bm)                           # table growth leaf mass in the beginning of timestep
        leafbase1 = allometry_f['bmToLeafMass'](bm + bm_increment)             # table growth leaf mass in the end of timestep
        leafmass = leafbase1 * nut_stat                                        # actual leaf mass
    
        
        #-------------------------------------------------
        leafmax = leafbase1 * 1.5
        leafmin = leafbase1 / 1.5
        
        net_ch = leafbase1-leafbase0
        gr_demand = leafmass - current_leafmass
    
        max_ch = (leafmax-leafmin)/(2*longevity(nut_stat))
        
        allowed_ch = np.zeros(n, dtype=float)
        lmass_ch = np.zeros(n, dtype=float)
        
        ix1 =  np.where(gr_demand >= 0)
        ix0 =  np.where(gr_demand < 0)
        allowed_ch[ix1] = net_ch[ix1] + max_ch[ix1]
        allowed_ch[ix0] = net_ch[ix0] - max_ch[ix0]
        lmass_ch[ix1] = np.minimum(allowed_ch[ix1], gr_demand[ix1])
        lmass_ch[ix0] = np.maximum(allowed_ch[ix0], gr_demand[ix0])
    
        new_lmass = current_leafmass + lmass_ch
        new_lmass = np.minimum(leafmax, new_lmass)
        new_lmass = np.maximum(leafmin, new_lmass)
        leaf_litter = new_lmass / longevity(nut_stat)
        C_consumption = lmass_ch + leaf_litter
    
        N_net = np.maximum(np.zeros(n), N_con(nut_stat)/1000.*new_lmass - N_con(previous_nut_stat)/1000.*current_leafmass)  
        Nleaf_litter = leaf_litter * (1-retrans['N'])*N_con(nut_stat)/1000.   
        Ndemand = N_net + Nleaf_litter 
        N_leaf = N_con(nut_stat)/1000.*new_lmass
    
        P_net = np.maximum(np.zeros(n), P_con(nut_stat)/1000.*new_lmass - P_con(previous_nut_stat)/1000.*current_leafmass)  
        Pleaf_litter = leaf_litter * (1-retrans['P'])*P_con(nut_stat)/1000.   
        Pdemand = P_net + Pleaf_litter 
        P_leaf = P_con(nut_stat)/1000.*new_lmass
    
        K_net = np.maximum(np.zeros(n), K_con(nut_stat)/1000.*new_lmass - K_con(previous_nut_stat)/1000.*current_leafmass)  
        Kleaf_litter = leaf_litter * (1-retrans['K'])*K_con(nut_stat)/1000.   
        Kdemand = K_net + Kleaf_litter 
        K_leaf = K_con(nut_stat)/1000.*new_lmass

        LAI = new_lmass/10000.*sla[spe]
        
        if printOpt:
            print ('********************************************')
            print (leafmin,  leafbase0, leafbase1, leafmax)
            print ('net_change', net_ch)
            print ('demanded growth', gr_demand)
            print ('max_change', max_ch)
             
            print ('allowed change', allowed_ch )
            
            print ('leaf mass change', lmass_ch)
            print ('new leaf mass', new_lmass)
    
            print ('leaf_litter', leaf_litter)
            print ('basic consumption', C_consumption)
    
            print ('Ndemand ', N_net + Nleaf_litter)
            print ('Nlitter', Nleaf_litter)
    
            print ('Pdemand ', P_net + Pleaf_litter)
            print ('Plitter', Pleaf_litter)
    
            print ('Kdemand ', K_net + Kleaf_litter)
            print ('Klitter', Kleaf_litter)
            print ('Knet', K_net, K_con(nut_stat), K_con(previous_nut_stat))
    
            print('*******************************************')
            
            print ('nitrogen content,', N_con(nut_stat))
            print ('leaf longevity', longevity(nut_stat))
        return new_lmass, leaf_litter, C_consumption, leafmax, leafmin,\
            Ndemand, Nleaf_litter, N_leaf, Pdemand, Pleaf_litter, P_leaf, Kdemand,\
            Kleaf_litter, K_leaf, LAI 

    def cutting(self, yr, nut_stat, to_ba = 0.5 ):
        """Unit here /ha"""
        if to_ba < 1.0 :
            agearr = self.agearr
            ixs = self.ixs
            for m in self.nlyrs:
                if m > 0:

                    self.nonwoody_lresid[ixs[m]] = (self.new_lmass[ixs[m]] + self.allodic[m].allometry_f['bmToFineRoots'](self.biomass))*self.stems 
                    self.n_nonwoody_lresid[ixs[m]] = (self.N_leaf + self.allodic[m].allometry_f['bmToNFineRoots'](self.biomass)) * self.stems 
                    self.p_nonwoody_lresid[ixs[m]] = (self.P_leaf + self.allodic[m].allometry_f['bmToPFineRoots'](self.biomass)) * self.stems 
                    self.k_nonwoody_lresid[ixs[m]] = (self.K_leaf + self.allodic[m].allometry_f['bmToKFineRoots'](self.biomass)) * self.stems 
            
                    self.woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToWoodyLoggingResidues'](self.biomass)  * self.stems
                    self.n_woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToNWoodyLoggingResidues'](self.biomass)  * self.stems
                    self.p_woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToPWoodyLoggingResidues'](self.biomass)  * self.stems
                    self.k_woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToKWoodyLoggingResidues'](self.biomass)  * self.stems
        
                    agearr[ixs[m]] = np.ones(self.ncols)[ixs[m]]
                    self.initialize_domain(agearr, nut_stat)
                    print ('+        cutting in ' + self.name + ' year ' + str(yr)  )
        else:
            for m in self.nlyrs:
                if m > 0:
                    print ('******** Now cutting to: ' , to_ba)
                    print (self.name)
                    print ('basal area')
                    print (self.basalarea*self.stems)
                    print ('n stems')
                    print (self.stems)
                    print ('cut stems')
                    cut_stems = (1.0 - to_ba/(self.basalarea*self.stems)) * self.stems
                    print (cut_stems)
                    print ('remaining stems')
                    remaining_stems = (to_ba/(self.basalarea*self.stems))*self.stems
                    print (remaining_stems)
                    
                    self.remaining_share = to_ba/(self.allodic[m].allometry_f['bmToBa'](self.biomass)*self.stems)
                    
                    print ('remaining share')
                    print (self.remaining_share)
                    agearr = self.agearr
                    ixs = self.ixs
                    
                    print ('nonwoodyl')
                    print (self.new_lmass[ixs[m]]*cut_stems[ixs[m]])
                    print (self.allodic[m].allometry_f['bmToFineRoots'](self.biomass[ixs[m]]) * cut_stems[ixs[m]])
                    
                    self.nonwoody_lresid[ixs[m]] =  (self.new_lmass[ixs[m]]\
                        + self.allodic[m].allometry_f['bmToFineRoots'](self.biomass[ixs[m]])) * cut_stems[ixs[m]]
                    
                    self.n_nonwoody_lresid[ixs[m]] = (self.N_leaf[ixs[m]] + self.allodic[m].allometry_f['bmToNFineRoots'](self.biomass[ixs[m]])) * cut_stems[ixs[m]]
                    self.p_nonwoody_lresid[ixs[m]] = (self.P_leaf[ixs[m]] + self.allodic[m].allometry_f['bmToPFineRoots'](self.biomass[ixs[m]])) * cut_stems[ixs[m]]
                    self.k_nonwoody_lresid[ixs[m]] = (self.K_leaf[ixs[m]] + self.allodic[m].allometry_f['bmToKFineRoots'](self.biomass[ixs[m]])) * cut_stems[ixs[m]]
            
                    self.woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToWoodyLoggingResidues'](self.biomass[ixs[m]]) * cut_stems[ixs[m]]
                    self.n_woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToNWoodyLoggingResidues'](self.biomass[ixs[m]]) * cut_stems[ixs[m]]
                    self.p_woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToPWoodyLoggingResidues'](self.biomass[ixs[m]]) * cut_stems[ixs[m]]
                    self.k_woody_lresid[ixs[m]] = self.allodic[m].allometry_f['bmToKWoodyLoggingResidues'](self.biomass[ixs[m]]) * cut_stems[ixs[m]]

                    print ('nonwoodylogging resids after adding')
                    print (self.nonwoody_lresid[ixs[m]])
                    
                    print ('woody logging residues')
                    print (self.woody_lresid)
            
                    self.update(self.biomass)
                    """This update to stand or to main????? """
            
            
