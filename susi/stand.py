# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 18:59:52 2022

@author: alauren
"""
import numpy as np
from susi.canopylayer import Canopylayer
            
class Stand():
    def __init__(self, nscens, yrs, canopylayers, ncols, sfc, agearr, mottifile, photopara):
        """
        ALL VARIABLES IN STAND OBJECT ARE IN ha AND kg -BASIS
        Creates canopy layer instances
        Stand object composes of canopy layer objects and keeps track on the 
        stand-wise sums of the variables
        Stand is array-form and has dimensions of number of columns in the strip
        Input:
            nscens , int, number of scenarios in the simulation
            yrs, int, number of years in the simulation
            canopylayers, dict in spara, contains integer arrays (len(ncols)) for each canopy layer pointing to specific Motti file
            ncols, int, number of columns along the strip 
            sfc, site fertility class
            agearr, dict of float arrays (len(ncols)) for stand age in the particular column and canopylayer
            mottifile, dict of dicts, telling the growth and yield (Motti files) in each canopy layer with key pointing to integer in the canopylayer dict 
            photopara - photosynthesis parameters used in the assimilation model
        """
        self.ncols = ncols                                                     # number of columns along the strip
        self.nscens = nscens                                                   # number of ditch depth scenarios in the simulation 
        self.yrs = yrs                                                         # number of years in the simulation
        self.photopara = photopara                                             # photosynthesis parameters for the assimilation function (Mäkelä et al. 2008)
        self.nut_stat = np.ones(ncols) #*0.5                                         # nutrient status, make this an argument

        
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
        
        self.dominant = Canopylayer('dominant', nscens, yrs, ncols, ndominants, sfc, agearr['dominant'],\
                                    mottifile['path'], mottifile['dominant'], ixdominants, photopara, self.nut_stat)
        self.subdominant = Canopylayer('subdominant',  nscens, yrs, ncols, nsubdominants, sfc, agearr['subdominant'],\
                                       mottifile['path'], mottifile['subdominant'], ixsubdominants, photopara, self.nut_stat)
        self.under = Canopylayer('under',  nscens, yrs, ncols, nunder, sfc, agearr['under'], mottifile['path'], \
                                 mottifile['under'], ixunder, photopara , self.nut_stat)
        self.clyrs = [self.dominant, self.subdominant, self.under]             # list of canopy layers, used later in loops
        
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
        self.nonwoodylitter = np.zeros(ncols, dtype=float)                     # non woody litter kg/ha/yr
        self.n_nonwoodylitter = np.zeros(ncols, dtype=float)                   # N in nonwoody litter kg/ha/yr
        self.p_nonwoodylitter = np.zeros(ncols, dtype=float)                   # P in nonwoody litter kg/ha/yr
        self.k_nonwoodylitter = np.zeros(ncols, dtype=float)                   # K in nonwoody litter kg/ha/yr   
        self.pulpvolume = np.zeros(ncols, dtype=float)                         # pulpwood volume m3/ha
        self.stems = np.zeros(ncols, dtype=float)                              # stocking, number of stems pcs/ha
        self.volume = np.zeros(ncols, dtype=float)                             # total volume of the growing stock m3/ha
        self.volumegrowth = np.zeros(ncols, dtype=float)                       # total volume growth of the growing stock m3/ha/yr      
        self.biomassgrowth = np.zeros(ncols, dtype=float)                      # total biomass growth of stand kg/ha/yr
        self.woodylitter = np.zeros(ncols, dtype=float)                        # woody litter kg/ha/yr
        self.n_woodylitter = np.zeros(ncols, dtype=float)                      # N in woody litter kg/ha/yr
        self.p_woodylitter = np.zeros(ncols, dtype=float)                      # P in woody litter kg/ha/yr
        self.k_woodylitter = np.zeros(ncols, dtype=float)                      # K in woody litter kg/ha/yr   
        self.yi = np.zeros(ncols, dtype=float)                                 # yield. here same as volume

        self.nonwoody_lresid = np.zeros(ncols, dtype=float)                    # nonwoody logging residues kg/ha
        self.n_nonwoody_lresid = np.zeros(ncols, dtype=float)                  # N in nonwoody logging residues kg/ha
        self.p_nonwoody_lresid = np.zeros(ncols, dtype=float)                  # P in nonwoody logging residues kg/ha
        self.k_nonwoody_lresid = np.zeros(ncols, dtype=float)                  # K in nonwoody logging residues kg/ha
        
        self.woody_lresid = np.zeros(ncols, dtype=float)                       # woody logging residues kg/ha
        self.n_woody_lresid = np.zeros(ncols, dtype=float)                     # N in woody logging residues kg/ha
        self.p_woody_lresid = np.zeros(ncols, dtype=float)                     # P in woody logging residues kg/ha
        self.k_woody_lresid = np.zeros(ncols, dtype=float)                     # K in woody logging residues kg/ha

        self.basNdemand = np.zeros(ncols, dtype=float)                         # basic N demand kg/tree, in table growth conditions, used in nutrient status calculation
        self.basPdemand = np.zeros(ncols, dtype=float)                         # basic P demand kg/tree, in table growth conditions, used in nutrient status calculation
        self.basKdemand = np.zeros(ncols, dtype=float)                         # basic K demand kg/tree, in table growth conditions, used in nutrient status calculation
        
        self.n_leaf_demand = np.zeros(ncols, dtype=float)                      # current leaf demand for N kg/ha
        self.p_leaf_demand = np.zeros(ncols, dtype=float)                      # current leaf demand for P kg/ha
        self.k_leaf_demand = np.zeros(ncols, dtype=float)                      # current leaf demand for K kg/ha
        

        self.previous_nut_stat = np.ones(ncols)                                # nutrient status in previous year
       
        """ ATTN cl in tree basis, convert to ha basis"""
        for cl in self.clyrs:                                                  # sum standwise initial values from the canopy layers
            self.basalarea = self.basalarea + cl.basalarea * cl.stems                     # Here: cl.basalarea in tree basis, conversion to ha basis    
            self.biomass = self.biomass + cl.biomass * cl.stems                           
            self.hdom = np.maximum(self.hdom, cl.hdom)
            self.leafarea = self.leafarea + cl.leafarea * cl.stems                           
            self.leafmass = self.leafmass + cl.leafmass * cl.stems
            self.stems = self.stems + cl.stems 
            self.volume = self.volume + cl.volume * cl.stems
            self.volumegrowth  = self.volumegrowth + cl.volumegrowth * cl.stems
            
        
    def reset_domain(self, agearr):
        """
        Resets stand domain, re-initializes canopy layer instances 
        by setting the initial values for the state variables
        Sums all canopy layers to gain initial values for the stand
        Parameters
        ----------
        agearr : TYPE dictionary of array of floats, len(ncols) 
            DESCRIPTION. dictionary of stand ages in the beginning of simulation in each canopy layer
        Returns
        -------
        None.

        """
        # reset the stand and reinitialize with initial age for a new scenario                                            
        self.previous_nut_stat = np.ones(self.ncols) 
        self.nut_stat = np.ones(self.ncols)#*0.5

        self.dominant.initialize_domain(agearr['dominant'], self.nut_stat)                      
        self.subdominant.initialize_domain(agearr['subdominant'], self.nut_stat)
        self.under.initialize_domain(agearr['under'], self.nut_stat)
        self.reset_vars()                                                      # sets all variables to zero
        
        """ ATTN cl in tree basis, convert to ha basis"""
        for cl in self.clyrs:                                                  # sum over the canopy layers to get standwise values
            self.basalarea = self.basalarea + cl.basalarea * cl.stems
            self.biomass = self.biomass + cl.biomass * cl.stems
            self.hdom = np.maximum(self.hdom, cl.hdom)
            self.leafarea = self.leafarea + cl.leafarea * cl.stems
            self.leafmass = self.leafmass + cl.leafmass * cl.stems
            self.stems = self.stems + cl.stems 
            self.volume = self.volume + cl.volume * cl.stems
        
        
    def reset_vars(self):
        """
        Stand variables are sums over canopy layer variables - therefore the sum variables
        are needed to set to zero before a new summation
 
        """
        self.basalarea = self.basalarea * 0.0 
        self.biomass = self.biomass * 0.0
        self.hdom = self.hdom * 0.0
        self.leafarea = self.leafarea * 0.0 
        self.leafmass = self.leafmass * 0.0
        self.stems = self.stems * 0.0
        self.volume = self.volume * 0.0   
        self.volumegrowth = self.volumegrowth * 0.0
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
        
        self.nonwoody_lresid = self.nonwoody_lresid * 0.0 
        self.n_nonwoody_lresid = self.n_nonwoody_lresid * 0.0
        self.p_nonwoody_lresid = self.p_nonwoody_lresid * 0.0
        self.k_nonwoody_lresid = self.k_nonwoody_lresid  * 0.0
        
        self.woody_lresid = self.woody_lresid  * 0.0 
        self.n_woody_lresid = self.n_woody_lresid  * 0.0 
        self.p_woody_lresid = self.p_woody_lresid  * 0.0 
        self.k_woody_lresid = self.k_woody_lresid  * 0.0

        self.basNdemand = self.basNdemand * 0.0                         # basic N demand kg/tree, in table growth conditions, used in nutrient status calculation
        self.basPdemand = self.basPdemand * 0.0                          # basic P demand kg/tree, in table growth conditions, used in nutrient status calculation
        self.basKdemand = self.basKdemand * 0.0                          # basic K demand kg/tree, in table growth conditions, used in nutrient status calculation

        self.n_leaf_demand = self.n_leaf_demand * 0.0
        self.p_leaf_demand = self.n_leaf_demand * 0.0
        self.k_leaf_demand = self.n_leaf_demand * 0.0
        
    def update(self):
        """
        ALL UNITS must be converted to ha BASIS
        Updates the stand variables by summing all the canopy layers
        Calls canopy layer instances
        Note: to be run after the self.layer.assimilate
        Returns
        -------
        None.

        """
        biomass_ini = self.biomass
        
        self.reset_vars()
        
        """ ATTN cl in tree basis, convert to ha basis"""
        for cl in self.clyrs:
            self.basalarea = self.basalarea + cl.basalarea * cl.stems
            self.biomass = self.biomass + cl.biomass * cl.stems
            self.hdom = np.maximum(self.hdom, cl.hdom)
            self.leafarea = self.leafarea + cl.leafarea * cl.stems
            self.leafmass = self.leafmass + cl.leafmass * cl.stems
            self.stems = self.stems + cl.stems 
            self.volume = self.volume + cl.volume * cl.stems      
            self.volumegrowth = self.volumegrowth + cl.volumegrowth * cl.stems
            self.yi = self.yi + cl.yi * cl.stems
            self.logvolume = self.logvolume + cl.logvolume * cl.stems
            self.pulpvolume = self.pulpvolume + cl.pulpvolume * cl.stems
            self.finerootlitter = self.finerootlitter + cl.finerootlitter * cl.stems 
            self.n_finerootlitter = self.n_finerootlitter + cl.n_finerootlitter * cl.stems 
            self.p_finerootlitter = self.p_finerootlitter + cl.p_finerootlitter * cl.stems
            self.k_finerootlitter = self.k_finerootlitter + cl.k_finerootlitter * cl.stems
            
            self.nonwoodylitter = self.nonwoodylitter +  cl.nonwoodylitter * cl.stems       # woody litter kg/ha/yr
            self.n_nonwoodylitter = self.n_nonwoodylitter + cl.n_nonwoodylitter * cl.stems # N in woody litter kg/ha/yr
            self.p_nonwoodylitter = self.p_nonwoodylitter + cl.p_nonwoodylitter * cl.stems # P in woody litter kg/ha/yr
            self.k_nonwoodylitter = self.k_nonwoodylitter + cl.k_nonwoodylitter * cl.stems # K in woody litter kg/ha/yr   
            
            
            self.woodylitter = self.woodylitter + cl.woodylitter * cl.stems
            self.n_woodylitter = self.n_woodylitter + cl.n_woodylitter * cl.stems
            self.p_woodylitter = self.p_woodylitter + cl.p_woodylitter * cl.stems
            self.k_woodylitter = self.k_woodylitter + cl.k_woodylitter * cl.stems

            #self.n_demand = self.n_demand + (cl.n_demand + cl.Nleafdemand) * cl.stems
            #self.p_demand = self.p_demand + (cl.p_demand+ cl.Pleafdemand) * cl.stems
            #self.k_demand = self.k_demand + (cl.k_demand+ cl.Kleafdemand) * cl.stems
            self.n_demand = self.n_demand + cl.n_demand * cl.stems
            self.p_demand = self.p_demand + cl.p_demand * cl.stems
            self.k_demand = self.k_demand + cl.k_demand * cl.stems

            self.basNdemand = self.basNdemand + cl.basNdemand * cl.stems  
            self.basPdemand = self.basPdemand + cl.basPdemand * cl.stems  
            self.basKdemand = self.basKdemand + cl.basKdemand * cl.stems  
            
            self.n_leaf_demand = self.n_leaf_demand + cl.Nleafdemand * cl.stems  
            self.p_leaf_demand = self.p_leaf_demand + cl.Pleafdemand * cl.stems  
            self.k_leaf_demand = self.k_leaf_demand + cl.Kleafdemand * cl.stems  
        
        self.biomassgrowth = self.biomass - biomass_ini
        
    def assimilate(self, forc, wt, afp):
        """
        Runs the photosyntheis function for all canopy layers 
        Calls canopy layer instances
        First it arranges the canopy layers into height order, and calculates the above leaf area
        
        Parameters
        ----------
        forc : TYPE   pandas dataframe
            DESCRIPTION. year-long measured daily weather variables
        wt : TYPE   pandas dataframe  
            DESCRIPTION. simulated water tables, shape: days, ncols
        afp : TYPE pandas dataframe 
            DESCRIPTION. air-filled porosity of rooting zone, shape days, ncols

        Returns
        -------
        None.

        """
        # specieswise specific leaf area: 1 Scots pine, 2:Norway spruce, 3: Birc h
        # assimilates each canopy layer and updates allometric variables in canopy layers
        heightarray = np.vstack([self.dominant.hdom, self.subdominant.hdom, self.under.hdom])          # array of canopy layer heights
        h_order = np.argsort(heightarray*-1, axis=0)                                                   # sort along column in descending order, indices
        laiarray = np.vstack([self.dominant.leafarea*self.dominant.stems,\
                              self.subdominant.leafarea*self.subdominant.stems,\
                              self.under.leafarea*self.under.stems])           # array of leaf areas converterd to m2/m2
        
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
        
        self.dominant.assimilate(forc, wt, afp, self.previous_nut_stat, self.nut_stat, lai_above[0,:])         # npp, leaf dynamics and updating the canopylayers
        self.subdominant.assimilate(forc, wt, afp, self.previous_nut_stat, self.nut_stat, lai_above[1,:])      # npp, leaf dynamics and updating the canopylayers
        self.under.assimilate(forc, wt, afp, self.previous_nut_stat, self.nut_stat, lai_above[2,:])            # npp, leaf dynamics and updating the canopylayers
        
        # updating call from the main program
        
    def update_nutrient_status(self, groundvegetation, N_supply, P_supply, K_supply):
        """
        Calculates nutrient status of the stand: supply/(stand demand + ground vegetation demand) 
        Change in nutrient status is delayed using time delay difference function

        Parameters
        ----------
        groundvegetation : TYPE instance of groundvegetation class
            DESCRIPTION. includes the nutrient demand for the ground vegetation
        N_supply : TYPE array len(ncols)
            DESCRIPTION. N supply from decomposition, atmospheric deposition and fertilization kg/ha/yr 
        P_supply : TYPE array len(ncols)
            DESCRIPTION. P supply from decomposition, atmospheric deposition and fertilization kg/ha/yr
        K_supply : TYPE array len(ncols)
            DESCRIPTION. K supply from decomposition, atmospheric deposition and fertilization kg/ha/yr

        Returns
        -------
        None.

        """
        self.previous_nut_stat = self.nut_stat.copy()
        
        nstat = np.ones((3, self.ncols))
        # nstat[0,:] = N_supply / (self.n_demand + groundvegetation.nup)
        # nstat[1,:] = P_supply / (self.p_demand + groundvegetation.pup)
        # nstat[2,:] = K_supply / (self.k_demand + groundvegetation.kup)
        nstat[0,:] = N_supply / (self.n_demand + self.basNdemand + groundvegetation.nup)
        nstat[1,:] = P_supply / (self.p_demand + self.basPdemand + groundvegetation.pup)
        nstat[2,:] = K_supply / (self.k_demand + self.basKdemand + groundvegetation.kup)

        minnstat = np.min(nstat, axis=0)
        
        tau = 3.0
        for c in range(self.ncols):
            self.nut_stat[c] = self.nut_stat[c] + (minnstat[c] - self.nut_stat[c])/tau     


        
    def update_spara(self, spara):
        spara['vol'] = self.volume                                                # initial stand volume along the cross section m3/ha in each node
        spara['hdom'] = self.hdom                                                 # these are for printing purposes only
        return (spara)
    
    def update_lresid(self):
        
        for cl in self.clyrs:
            self.nonwoody_lresid = self.nonwoody_lresid + cl.nonwoody_lresid #* cl.stems 
            self.n_nonwoody_lresid = self.n_nonwoody_lresid + cl.n_nonwoody_lresid #* cl.stems
            self.p_nonwoody_lresid = self.p_nonwoody_lresid + cl.p_nonwoody_lresid #* cl.stems
            self.k_nonwoody_lresid = self.k_nonwoody_lresid + cl.k_nonwoody_lresid #* cl.stems
            
            self.woody_lresid = self.woody_lresid  + cl.woody_lresid #* cl.stems
            self.n_woody_lresid = self.n_woody_lresid  + cl.n_woody_lresid #* cl.stems
            self.p_woody_lresid = self.p_woody_lresid  + cl.p_woody_lresid #* cl.stems
            self.k_woody_lresid = self.k_woody_lresid  + cl.k_woody_lresid #* cl.stems


    def reset_lresid(self):
        """
        Resets logging residue arrays after locating them to decomposition model

        Returns
        -------
        None.

        """
                    
        self.nonwoody_lresid = self.nonwoody_lresid * 0.0 
        self.n_nonwoody_lresid = self.n_nonwoody_lresid * 0.0
        self.p_nonwoody_lresid = self.p_nonwoody_lresid * 0.0
        self.k_nonwoody_lresid = self.k_nonwoody_lresid  * 0.0
        
        self.woody_lresid = self.woody_lresid  * 0.0 
        self.n_woody_lresid = self.n_woody_lresid  * 0.0 
        self.p_woody_lresid = self.p_woody_lresid  * 0.0 
        self.k_woody_lresid = self.k_woody_lresid  * 0.0

        for cl in self.clyrs:
            cl.nonwoody_lresid = cl.nonwoody_lresid * 0.0
            cl.n_nonwoody_lresid = cl.n_nonwoody_lresid * 0.0
            cl.p_nonwoody_lresid = cl.p_nonwoody_lresid * 0.0
            cl.k_nonwoody_lresid = cl.k_nonwoody_lresid * 0.0
            
            cl.woody_lresid = cl.woody_lresid * 0.0
            cl.n_woody_lresid = cl.n_woody_lresid * 0.0
            cl.p_woody_lresid = cl.p_woody_lresid * 0.0
            cl.k_woody_lresid = cl.k_woody_lresid * 0.0
            