# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 18:59:52 2022

@author: alauren
"""
import numpy as np
from scipy.interpolate import interp1d
from allometry import Allometry
from susi_utils import rew_drylimit, assimilation_yr

class Canopylayer():
    _registry=[]
    
    """
    Canopylayer keeps track on the development of biomass components within the 
    different layers of canopy. Canopy layer is of single tree species. It is 
    initialized with growth and yield simulator outputfile that describes development
    of biomass components in time. Canopy layers is array shaped with dimensions 
    number of columns along the strip.
    """
    
    def __init__(self, name, nscens, yrs, ncols, nlyrs, agearr, mottipath, mottifile, ixs, photopara):
        self._registry.append(self)                                            # locates canopy instances into a list
        self.name = name                                                       # name of the canopy layer e.g. 'dominant', 'subdominant', etc.
        self.nlyrs = nlyrs                                                     # number of different canopy layers along the strip
        self.ixs = ixs                                                         # indices for the location of the different canopy layers along the strip 
        self.ncols = ncols                                                     # number of columns in the strip
        self.agearr = agearr                                                   # age of the canopy layer, yrs  
        self.photopara = photopara
        self.nscens = nscens
        self.yrs = yrs
        
        # -------- Biomass interpolation functions-------------------
        self.allodic ={}
        for ncanopy in nlyrs:
            if ncanopy > 0:
                self.allodic[ncanopy] = Allometry()                            # allometry instance to the dictionary
                mfile = mottipath + mottifile[ncanopy]                         # mottifile where the allometry tables exist
                self.allodic[ncanopy].motti_development(mfile)                 # run the allometry; interpolation functions in the instance
        
        self.initialize_domain(agearr)                                         # create variables and set initial values

    def initialize_domain(self, agearr):
        nlyrs = self.nlyrs                                                     # number of different canopy layers along the strip
        ixs = self.ixs                                                         # indices for the canopy layers along the strip
        
        ncols = self.ncols                                                     # number of columns along the strip
        self.basalarea = np.zeros(ncols, dtype=float)                          # basal area in the canopy layer m2/ha
        self.biomass = np.zeros(ncols, dtype=float)                            # total canopy layer biomass excluding leaves
        self.n_demand = np.zeros(ncols, dtype=float)                           # N demand exluding leaves kg/ha/yr
        self.p_demand = np.zeros(ncols, dtype=float)                           # P demand exluding leaves kg/ha/yr 
        self.k_demand = np.zeros(ncols, dtype=float)                           # K demand exluding leaves kg/ha/yr
        self.hdom = np.zeros(ncols, dtype=float)                               # domainant height m in the canopy layer   
        self.leafarea = np.zeros(ncols, dtype=float)                           # leaf area in the canopy layer m2 m-2 
        self.lai_above = np.zeros(ncols, dtype=float)                          # one sided leaf area above this canopy layer, m2 m-2 
        self.leafmass = np.zeros(ncols, dtype=float)                           # basic leaf mass kg/ha
        self.logvolume = np.zeros(ncols, dtype=float)                          # volume of the saw logs m3/ha
        self.finerootlitter = np.zeros(ncols, dtype=float)                     # fineroot litter in the canopy layer kg/ha/yr
        self.n_finerootlitter = np.zeros(ncols, dtype=float)                   # N in fineroot litter in the canopy layer kg/ha/yr 
        self.p_finerootlitter = np.zeros(ncols, dtype=float)                   # P in fineroot litter in the canopy layer kg/ha/yr      
        self.k_finerootlitter = np.zeros(ncols, dtype=float)                   # K in fineroot litter in the canopy layer kg/ha/yr  
        self.pulpvolume= np.zeros(ncols, dtype=float)                          # volume of pulpwood in the canopy layer m3/ha

        self.NPP= np.zeros(ncols, dtype=float)                                 # net primary production kg/ha/yr dry matter
        self.NPP_pot= np.zeros(ncols, dtype=float)                             # potential net primary production kg/ha/yr dry matter

        self.nonwoodylitter = np.zeros(ncols, dtype=float)                     # nonwoody litter kg/ha/yr
        self.n_nonwoodylitter = np.zeros(ncols, dtype=float)                   # N in nonwoody litter kg/ha/yr
        self.p_nonwoodylitter = np.zeros(ncols, dtype=float)                   # P in nonwoody litter kg/ha/yr
        self.k_nonwoodylitter = np.zeros(ncols, dtype=float)                   # K in nonwoody litter kg/ha/yr
        
        self.species = np.zeros(ncols, dtype=int)                              # tree species: 1 Scots pine, 2: Norway spruce, 3: Betula pendula
        self.stems = np.zeros(ncols, dtype=float)                              # stocking number of trees per ha
        self.volume = np.zeros(ncols, dtype=float)                             # volume of the growing stock m3/ha
        self.woodylitter = np.zeros(ncols, dtype=float)                        # woody litter kg/ha/yr
        self.n_woodylitter = np.zeros(ncols, dtype=float)                      # N in woody litter kg/ha/yr
        self.p_woodylitter = np.zeros(ncols, dtype=float)                      # P in woody litter kg/ha/yr
        self.k_woodylitter = np.zeros(ncols, dtype=float)                      # K in woody litter kg/ha/yr
        self.yi = np.zeros(ncols, dtype=float)                                 # yied. here same as vcolume
        
        #----------- ARRANGE Leaf dynamics variables-----------------
        self.new_lmass = np.zeros(ncols, dtype=float)                          # leaf mass after dynamics computation kg/ha
        self.leaf_litter = np.zeros(ncols, dtype=float)                        # leaf litter kg/ha/yr
        self.C_consumption = np.zeros(ncols, dtype=float)                      # check the unit: C or mass? C consumption in leaf processes
        self.leafmax = np.zeros(ncols, dtype=float)                            # upper limit for leaf mass kg/ha
        self.leafmin = np.zeros(ncols, dtype=float)                            # lower limit for leaf mass kg/ha
        self.Nleafdemand = np.zeros(ncols, dtype=float)                        # N demand of leaf production kg/ha/yr
        self.Nleaf_litter = np.zeros(ncols, dtype=float)                       # N in litterfall kg/ha/yr
        self.N_leaf = np.zeros(ncols, dtype=float)                             # N in leaves kg/ha
        self.Pleafdemand = np.zeros(ncols, dtype=float)                        # P demand of leaf production kg/ha/yr
        self.Pleaf_litter= np.zeros(ncols, dtype=float)                        # P in litterfall kg/ha/yr
        self.P_leaf = np.zeros(ncols, dtype=float)                             # P in leaves kg/ha
        self.Kleafdemand = np.zeros(ncols, dtype=float)                        # K demand in leaf production kg/ha/yr 
        self.Kleaf_litter  = np.zeros(ncols, dtype=float)                      # K in litterfall kg/ha/yr
        self.K_leaf = np.zeros(ncols, dtype=float)                             # K in leaves kg/ha
        
        # ---------- Initial values from age -------------------------------------------------------
        for m in nlyrs:
            if m > 0:
                self.basalarea[ixs[m]] = self.allodic[m].allometry_f['ageToBa'](agearr[ixs[m]])
                self.biomass[ixs[m]] = self.allodic[m].allometry_f['ageToBm'](agearr[ixs[m]])
                self.hdom[ixs[m]] = self.allodic[m].allometry_f['ageToHdom'](agearr[ixs[m]])
                self.leafarea[ixs[m]] = self.allodic[m].allometry_f['bmToLAI'](self.biomass[ixs[m]])
                self.leafmass[ixs[m]] = self.allodic[m].allometry_f['ageToLeaves'](agearr[ixs[m]])
                self.species[ixs[m]] = self.allodic[m].sp
                self.stems[ixs[m]] = self.allodic[m].allometry_f['bmToStems'](self.biomass[ixs[m]])
                self.volume[ixs[m]] = self.allodic[m].allometry_f['ageToVol'](agearr[ixs[m]])
                

    def update(self, bm):
        #------------Update all variables with new biomass--------------------------------------
        ixs = self.ixs
        for m in self.nlyrs:
            if m > 0:
                self.basalarea[ixs[m]] = self.allodic[m].allometry_f['bmToBa'](bm[ixs[m]])
                self.biomass[ixs[m]] = bm[ixs[m]]
                self.hdom[ixs[m]] = self.allodic[m].allometry_f['bmToHdom'](bm[ixs[m]])
                self.leafarea[ixs[m]] = self.allodic[m].allometry_f['bmToLAI'](bm[ixs[m]])
                self.leafmass[ixs[m]] = self.allodic[m].allometry_f['bmToLeafMass'](bm[ixs[m]])
                self.stems[ixs[m]] = self.allodic[m].allometry_f['bmToStems'](bm[ixs[m]])
                self.volume[ixs[m]] = self.allodic[m].allometry_f['bmToYi'](bm[ixs[m]])
                self.n_demand[ixs[m]] = self.allodic[m].allometry_f['bmToNdemand'](bm[ixs[m]])
                self.p_demand[ixs[m]] = self.allodic[m].allometry_f['bmToPdemand'](bm[ixs[m]])
                self.k_demand[ixs[m]] = self.allodic[m].allometry_f['bmToKdemand'](bm[ixs[m]])
                self.logvolume[ixs[m]] = self.allodic[m].allometry_f['volToLogs'](self.volume[ixs[m]])        
                self.finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToFinerootLitter'](bm[ixs[m]])
                self.n_finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToNFineRootLitter'](bm[ixs[m]])
                self.p_finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToPFineRootLitter'](bm[ixs[m]])
                self.k_finerootlitter[ixs[m]] = self.allodic[m].allometry_f['bmToKFineRootLitter'](bm[ixs[m]])

                self.nonwoodylitter[ixs[m]] = self.finerootlitter[ixs[m]] + self.leaf_litter[ixs[m]]           # nonwoody litter kg/ha/yr
                self.n_nonwoodylitter[ixs[m]] = self.n_finerootlitter[ixs[m]] + self.Nleaf_litter[ixs[m]]     # N in nonwoody litter kg/ha/yr
                self.p_nonwoodylitter[ixs[m]] = self.p_finerootlitter[ixs[m]] + self.Pleaf_litter[ixs[m]]     # P in nonwoody litter kg/ha/yr
                self.k_nonwoodylitter[ixs[m]] = self.k_finerootlitter[ixs[m]] + self.Kleaf_litter[ixs[m]]     # K in nonwoody litter kg/ha/yr

                self.pulpvolume[ixs[m]] = self.allodic[m].allometry_f['volToPulp'](self.volume[ixs[m]])                        
                self.woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToWoodyLitter'](bm[ixs[m]])
                self.n_woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToNWoodyLitter'](bm[ixs[m]])
                self.p_woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToPWoodyLitter'](bm[ixs[m]])
                self.k_woodylitter[ixs[m]] = self.allodic[m].allometry_f['bmToKWoodyLitter'](bm[ixs[m]])
                self.yi[ixs[m]] = self.allodic[m].allometry_f['bmToYi'](bm[ixs[m]])

    def assimilate(self, forc, wt, afp, previous_nut_stat, nut_stat):
        
        self.NPP, self.NPP_pot = assimilation_yr(self.photopara, forc, wt, afp,
                                       self.leafarea*2, self.lai_above*2)  #double sided LAI required

        self.NPP = self.NPP * nut_stat
        self.NPP_pot = self.NPP_pot * nut_stat
        
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
                        self.leafarea[ixs[m]]  = self.leaf_dynamics(bm[ixs[m]],\
                                                bm_increment[ixs[m]], current_leafmass[ixs[m]], 
                                                previous_nut_stat[ixs[m]], nut_stat[ixs[m]],\
                                                self.agearr[ixs[m]], self.allodic[m].allometry_f,\
                                                self.species[ixs[m]], printOpt=False)
        
        delta_bm_noleaves = self.NPP - self.C_consumption - self.finerootlitter  - self.woodylitter 
        self.leafmass = self.new_lmass
        self.update(self.biomass + delta_bm_noleaves)
        

    def leaf_dynamics(self, bm, bm_increment, current_leafmass, previous_nut_stat,\
                      nut_stat, agenow, allometry_f, species, printOpt=False):
        """
        input:
             bm, current biomass, array, kg/ha without leaves
             bm_increment, array, total NPP 
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
        longevityLeaves = {'Pine':[4.0, 3.0], 'Spruce':[8.0, 3.5], 'Birch':[1.0, 1.0]}                     # yrs, life span of leaves and fine roots    
    
        N_con = interp1d(np.array([0.66,1.5]), np.array(nuts[spe]['Foliage']['N']), fill_value=tuple(nuts[spe]['Foliage']['N']), bounds_error=False)
        P_con = interp1d(np.array([0.66,1.5]), np.array(nuts[spe]['Foliage']['P']), fill_value=tuple(nuts[spe]['Foliage']['P']), bounds_error=False)
        K_con = interp1d(np.array([0.66,1.5]), np.array(nuts[spe]['Foliage']['K']), fill_value=tuple(nuts[spe]['Foliage']['K']), bounds_error=False)
        
        longevity = interp1d(np.array([0.66,1.5]), np.array(longevityLeaves[spe]), fill_value= tuple(longevityLeaves[spe]), bounds_error=False)
    
        n = len(nut_stat)
        #************ Biomass and litter**************
        leafbase0 =  allometry_f['bmToLeafMass'](bm)
        leafbase1 = allometry_f['bmToLeafMass'](bm + bm_increment)
        leafmass = leafbase1 * nut_stat
    
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

    def cutting(self, yr, method='all'):
        if method == 'all':
            agearr = self.agearr
            ixs = self.ixs
            for m in self.nlyrs:
                if m > 0:

                    self.nonwoodylitter[ixs[m]] += self.new_lmass[ixs[m]] + self.allodic[m].allometry_f['bmToFineRoots'](self.biomass)
                    self.n_nonwoodylitter[ixs[m]] += self.N_leaf + self.allodic[m].allometry_f['bmToNFineRoots'](self.biomass)
                    self.p_nonwoodylitter[ixs[m]] += self.P_leaf + self.allodic[m].allometry_f['bmToPFineRoots'](self.biomass)
                    self.k_nonwoodylitter[ixs[m]] += self.K_leaf + self.allodic[m].allometry_f['bmToKFineRoots'](self.biomass)
            
                    self.woodylitter[ixs[m]] += self.allodic[m].allometry_f['bmToWoodyLoggingResidues'](self.biomass)
                    self.n_woodylitter[ixs[m]] += self.allodic[m].allometry_f['bmToNWoodyLoggingResidues'](self.biomass)
                    self.p_woodylitter[ixs[m]] += self.allodic[m].allometry_f['bmToPWoodyLoggingResidues'](self.biomass)
                    self.k_woodylitter[ixs[m]] += self.allodic[m].allometry_f['bmToKWoodyLoggingResidues'](self.biomass)
        
                    agearr[ixs[m]] = np.ones(self.ncols)[ixs[m]]
                    self.initialize_domain(agearr)
                    print ('+        cutting in ' + self.name + ' year ' + str(yr)  )
        

class Stand():
    def __init__(self, nscens, yrs, canopylayers, ncols, agearr, mottifile, photopara):
        """
        Stand object composes of canopy layer objects and keeps track on the 
        stand-wise sums of the variables
        Stand is array-form and has dimensions of number of columns in the strip
        """
        self.ncols = ncols
        self.nscens = nscens
        self.yrs = yrs
        self.photopara = photopara
        ndominants = np.unique(canopylayers['dominant'])                       # number codes of different dominant layers, number refers to key in mottifiles dictionary, 0 implies no dominant layer
        nsubdominants = np.unique(canopylayers['subdominant'])                 # number codes of different subdominant layers, number refers to key in mottifiles dictionary, 0 implies no dominant layer
        nunder  = np.unique(canopylayers['under'])                             # number codes of different undermost layers, number refers to key in mottifiles dictionary, 0 implies no dominant layer  
        
        ixdominants ={}                                                        # location indices for dominant canopy layers, along the transect
        for m in ndominants:
            if m > 0: ixdominants[m] =  np.where(canopylayers['dominant'] == m)
        
        ixsubdominants = {}                                                    # location indices for subdominant canopy layers
        for m in nsubdominants:
            if m > 0: ixsubdominants[m] =  np.where(canopylayers['subdominant'] == m)
        
        ixunder = {}                                                           # location indices for undersmost canopy layer
        for m in nunder:
            if m > 0: ixunder[m] =  np.where(canopylayers['under'] == m)
        
        self.dominant = Canopylayer('dominant', nscens, yrs, ncols, ndominants, agearr['dominant'],\
                                    mottifile['path'], mottifile['dominant'], ixdominants, photopara)
        self.subdominant = Canopylayer('subdominant',  nscens, yrs, ncols, nsubdominants, agearr['subdominant'],\
                                       mottifile['path'], mottifile['subdominant'], ixsubdominants, photopara)
        self.under = Canopylayer('under',  nscens, yrs, ncols, nunder, agearr['under'], mottifile['path'], \
                                 mottifile['under'], ixunder, photopara)
        
        
        # ---------- create stand variables------------------------------------
        self.basalarea = np.zeros(ncols, dtype=float)                          # stand basal area m2/ha                          
        self.biomass = np.zeros(ncols, dtype=float)                            # stand dry biomass kg/ha
        self.n_demand = np.zeros(ncols, dtype=float)                           # stand N demand excluding leaves kg/ha/yr
        self.p_demand = np.zeros(ncols, dtype=float)                           # stand N demand excluding leaves kg/ha/yr
        self.k_demand = np.zeros(ncols, dtype=float)                           # stand N demand excluding leaves kg/ha/yr 
        self.hdom = np.zeros(ncols, dtype=float)                               # dominant height m
        self.leafarea = np.zeros(ncols, dtype=float)                           # one sided leaf area m2 m-2
        self.leafmass = np.zeros(ncols, dtype=float)                           # leaf dry biomass kg/ha
        self.logvolume = np.zeros(ncols, dtype=float)                          # saw log volume m3/ha
        self.finerootlitter = np.zeros(ncols, dtype=float)                     # fine root litter kg/ha/yr
        self.n_finerootlitter = np.zeros(ncols, dtype=float)                   # N in fine root litter kg/ha/yr
        self.p_finerootlitter = np.zeros(ncols, dtype=float)                   # P in fine root litter kg/ha/yr
        self.k_finerootlitter = np.zeros(ncols, dtype=float)                   # K in fine root litter kg/ha/yr
        self.nonwoodylitter = np.zeros(ncols, dtype=float)                        # woody litter kg/ha/yr
        self.n_nonwoodylitter = np.zeros(ncols, dtype=float)                      # N in woody litter kg/ha/yr
        self.p_nonwoodylitter = np.zeros(ncols, dtype=float)                      # P in woody litter kg/ha/yr
        self.k_nonwoodylitter = np.zeros(ncols, dtype=float)                      # K in woody litter kg/ha/yr   
        self.pulpvolume = np.zeros(ncols, dtype=float)                         # pulpwood volume m3/ha
        self.stems = np.zeros(ncols, dtype=float)                              # stocking, number of stems pcs/ha
        self.volume = np.zeros(ncols, dtype=float)                             # total volume of the growing stock m3/ha
        self.woodylitter = np.zeros(ncols, dtype=float)                        # woody litter kg/ha/yr
        self.n_woodylitter = np.zeros(ncols, dtype=float)                      # N in woody litter kg/ha/yr
        self.p_woodylitter = np.zeros(ncols, dtype=float)                      # P in woody litter kg/ha/yr
        self.k_woodylitter = np.zeros(ncols, dtype=float)                      # K in woody litter kg/ha/yr   
        self.yi = np.zeros(ncols, dtype=float)                                 # yield. here same as volume

        self.previous_nut_stat = np.ones(ncols)*0.75
        self.nut_stat = np.ones(ncols)*0.75
       
        for cl in Canopylayer._registry:                                       # sum standwise initial values from the canopy layers
            self.basalarea = self.basalarea + cl.basalarea
            self.biomass = self.biomass + cl.biomass
            self.hdom = np.maximum(self.hdom, cl.hdom)
            self.leafarea = self.leafarea + cl.leafarea
            self.leafmass = self.leafmass + cl.leafmass
            self.stems = self.stems + cl.stems
            self.volume = self.volume + cl.volume
    
    
    def reset_domain(self, agearr):
        # reset the stand and reinitialize with initial age for a new scenario                                            
        self.dominant.initialize_domain(agearr['dominant'])                      
        self.subdominant.initialize_domain(agearr['subdominant'])
        self.under.initialize_domain(agearr['under'])
        self.reset_vars()
        for cl in Canopylayer._registry:
            self.basalarea = self.basalarea + cl.basalarea
            self.biomass = self.biomass + cl.biomass
            self.hdom = np.maximum(self.hdom, cl.hdom)
            self.leafarea = self.leafarea + cl.leafarea
            self.leafmass = self.leafmass + cl.leafmass
            self.stems = self.stems + cl.stems
            self.volume = self.volume + cl.volume
        
        self.previous_nut_stat = np.ones(self.ncols)*0.75
        self.nut_stat = np.ones(self.ncols)*0.75

        self.out = {}
        
    def reset_vars(self):
            # returns variables summing the canopy layers to zero 
            self.basalarea = self.basalarea * 0.0 
            self.biomass = self.biomass * 0.0
            self.hdom = self.hdom * 0.0
            self.leafarea = self.leafarea * 0.0 
            self.leafmass = self.leafmass * 0.0
            self.stems = self.stems * 0.0
            self.volume = self.volume * 0.0      
            self.yi = self.yi * 0.0
            self.logvolume = self.logvolume * 0.0
            self.pulpvolume = self.pulpvolume  * 0.0
            self.finerootlitter = self.finerootlitter * 0.0 
            self.n_finerootlitter = self.n_finerootlitter  * 0.0 
            self.p_finerootlitter = self.p_finerootlitter  * 0.0
            self.k_finerootlitter = self.k_finerootlitter  * 0.0
            self.nonwoodylitter = self.nonwoodylitter * 0.0                        # woody litter kg/ha/yr
            self.n_nonwoodylitter = self.n_nonwoodylitter * 0.0                     # N in woody litter kg/ha/yr
            self.p_nonwoodylitter = self.p_nonwoodylitter * 0.0                     # P in woody litter kg/ha/yr
            self.k_nonwoodylitter = self.k_nonwoodylitter * 0.0                     # K in woody litter kg/ha/yr   
            self.woodylitter = self.woodylitter  * 0.0
            self.n_woodylitter = self.n_woodylitter  * 0.0
            self.p_woodylitter = self.p_woodylitter  * 0.0
            self.k_woodylitter = self.k_woodylitter  * 0.0
            self.n_demand = self.n_demand  * 0.0
            self.p_demand = self.p_demand  * 0.0
            self.k_demand = self.k_demand  * 0.0
            
    
    def update(self):
        # to be run after the self.layer.assimilate
        
        self.reset_vars()
        
        for cl in Canopylayer._registry:
            self.basalarea = self.basalarea + cl.basalarea
            self.biomass = self.biomass + cl.biomass
            self.hdom = np.maximum(self.hdom, cl.hdom)
            self.leafarea = self.leafarea + cl.leafarea
            self.leafmass = self.leafmass + cl.leafmass
            self.stems = self.stems + cl.stems
            self.volume = self.volume + cl.volume      
            self.yi = self.yi + cl.yi
            self.logvolume = self.logvolume + cl.logvolume
            self.pulpvolume = self.pulpvolume + cl.pulpvolume
            self.finerootlitter = self.finerootlitter + cl.finerootlitter 
            self.n_finerootlitter = self.n_finerootlitter + cl.n_finerootlitter 
            self.p_finerootlitter = self.p_finerootlitter + cl.p_finerootlitter
            self.k_finerootlitter = self.k_finerootlitter + cl.k_finerootlitter
            
            self.nonwoodylitter = self.nonwoodylitter +  cl.nonwoodylitter       # woody litter kg/ha/yr
            self.n_nonwoodylitter = self.n_nonwoodylitter + cl.n_nonwoodylitter  # N in woody litter kg/ha/yr
            self.p_nonwoodylitter = self.p_nonwoodylitter + cl.p_nonwoodylitter  # P in woody litter kg/ha/yr
            self.k_nonwoodylitter = self.k_nonwoodylitter + cl.k_nonwoodylitter  # K in woody litter kg/ha/yr   
            
            self.woodylitter = self.woodylitter + cl.woodylitter
            self.n_woodylitter = self.n_woodylitter + cl.n_woodylitter
            self.p_woodylitter = self.p_woodylitter + cl.p_woodylitter
            self.k_woodylitter = self.k_woodylitter + cl.k_woodylitter
            self.n_demand = self.n_demand + cl.n_demand + cl.Nleafdemand
            self.p_demand = self.p_demand + cl.p_demand+ cl.Pleafdemand
            self.k_demand = self.k_demand + cl.k_demand+ cl.Kleafdemand
            
    def assimilate(self, forc, wt, afp):
        #specieswise specific leaf area: 1 Scots pine, 2:Norway spruce, 3: Birc h
        # assimilates each canopy layer and updates allometric variables in canopy layers
        heightarray = np.vstack([self.dominant.hdom, self.subdominant.hdom, self.under.hdom])          # array of canopy layer heights
        h_order = np.argsort(heightarray*-1, axis=0)                                                   # sort along column in descending order, indices
        laiarray = np.vstack([self.dominant.leafarea, self.subdominant.leafarea, self.under.leafarea]) # array of leaf areas
        
        laiout = np.zeros((3,self.ncols))                                      # initialize temporary lai array
        lai_above = np.zeros((4,self.ncols))                                   # initialize the ablove-lai array (used later in assimilation)

        for layer in range(3):                                                 # loop through canopy layers
          order = h_order[layer]                                               # inddices in heght array (descending order, 0 for the tallest)
          col = np.arange(0,self.ncols,1)                                      # indices along the strip
          laiout[layer,:] = laiarray[order, col]                               # locate lai on the height order
        laiabove = np.cumsum(laiout, axis = 0)                                 # cumulative lai sum above
        for layer in range(3):
          order = h_order[layer]
          col = np.arange(0,self.ncols,1)
          lai_above[layer+1,:] = laiabove[order, col]
        
        self.dominant.lai_above = lai_above[0,:]
        self.subdominant.lai_above = lai_above[1,:]
        self.under.lai_above = lai_above[2,:]
        
        self.dominant.assimilate(forc, wt, afp, self.previous_nut_stat, self.nut_stat)                                  # npp, leaf dynamics and updating the canopylayers
        self.subdominant.assimilate(forc, wt, afp, self.previous_nut_stat, self.nut_stat)                               # npp, leaf dynamics and updating the canopylayers
        self.under.assimilate(forc, wt, afp, self.previous_nut_stat, self.nut_stat)                                     # npp, leaf dynamics and updating the canopylayers
      
        self.update()  
                                                                # updates the stand variables
        
    def update_nutrient_status(self, groundvegetation, N_supply, P_supply, K_supply):
        
        # gvshare = groundvegetation.gv_leafmass / (self.leafmass + groundvegetation.gv_leafmass )     # maximum share of ground vegetation from the nutrient uptake -> shared in proportion of green biomass
        # n_available = (1- gvshare) * N_supply
        # p_available = (1- gvshare) * P_supply
        # k_available = (1- gvshare) * K_supply
        
        # nstat = np.ones((3, self.ncols))
        # nstat[0,:] = n_available / self.n_demand
        # nstat[1,:] = p_available / self.p_demand
        # nstat[2,:] = k_available / self.k_demand

        self.previous_nut_stat = self.nut_stat.copy()
        
        nstat = np.ones((3, self.ncols))
        nstat[0,:] = N_supply / (self.n_demand + groundvegetation.nup)
        nstat[1,:] = P_supply / (self.p_demand + groundvegetation.pup)
        nstat[2,:] = K_supply / (self.k_demand + groundvegetation.kup)
        minnstat = np.min(nstat, axis=0)
        
        tau = 3.0
        for c in range(self.ncols):
            self.nut_stat[c] = self.nut_stat[c] + (minnstat[c] - self.nut_stat[c])/tau     
        
        
    def update_spara(self, spara):
        spara['vol'] = self.volume                                                # initial stand volume along the cross section m3/ha in each node
        spara['hdom'] = self.hdom                                                 # these are for printing purposes only
        return (spara)
    
