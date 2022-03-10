# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 10:13:04 2022

@author: alauren
"""
from netCDF4 import Dataset 
from datetime import datetime
import numpy as np

class Outputs():
    def __init__(self, nscens, ncols, ndays, nyrs, nLyrs, fname):
        print('**** creating Susi netCDF4 file: ' + fname + ' ****')
        
        # create dataset & dimensions
        ff = fname 
        self.ncf = Dataset(ff, 'w')
        self.ncf.description = "Peatland simulator Susi results"
        self.ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.ncf.source = 'Susi v.2022,'

        self.ncf.createDimension('nscens', nscens)                                      # number of management scenarios
        self.ncf.createDimension('ncols', ncols)                                        # number of columns along the strip  
        self.nyrs = nyrs + 1                                                            # allow saving of the initial status      
        self.ncf.createDimension('nyrs', self.nyrs)                                          # number of years in the simulation
        self.ncf.createDimension('ndays', ndays)                                        # number of days in the simulation
        self.ncf.createDimension('nLyrs', nLyrs)                                        # numbre of peat layers in the vertical column
        

        # call as createVariable(varname,type,(dimensions))
        time = self.ncf.createVariable('time', 'f8', ('ndays',))
        time.units = "days since 0001-01-01 00:00:00.0"
        time.calendar = 'standard'
        
        print ('    + netCDF file created')
    def close(self):
        self.ncf.close()
        
    def initialize_stand(self):    
                
        basalarea = self.ncf.createVariable('/stand/basalarea','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        basalarea.units = 'stand basal area [m2/ha]'
        biomass = self.ncf.createVariable('/stand/biomass','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        biomass.units = 'stand dry biomass [kg/ha]'
        n_demand =  self.ncf.createVariable('/stand/n_demand','f4',('nscens','nyrs', 'ncols',))                          # stand N demand excluding leaves kg/ha/yr
        n_demand.units = 'stand N demand [kg/ha/yr]'
        p_demand =  self.ncf.createVariable('/stand/p_demand','f4',('nscens','nyrs', 'ncols',))                           # stand N demand excluding leaves kg/ha/yr
        p_demand.units = 'stand P demand [kg/ha/yr]'
        k_demand = self.ncf.createVariable('/stand/k_demand','f4',('nscens','nyrs', 'ncols',))                            # stand N demand excluding leaves kg/ha/yr 
        k_demand.units = 'stand K demand [kg/ha/yr]' 
        hdom =  self.ncf.createVariable('/stand/hdom','f4',('nscens','nyrs', 'ncols',))                               # dominant height m
        hdom.units = 'stand dominant height [m]'
        leafarea =  self.ncf.createVariable('/stand/leafarea','f4',('nscens','nyrs', 'ncols',))                           # one sided leaf area m2 m-2
        leafarea.units = 'stand one-sided leaf area [m2/m2]'
        leafmass =  self.ncf.createVariable('/stand/leafmass','f4',('nscens','nyrs', 'ncols',))                           # leaf dry biomass kg/ha
        leafmass.units = 'stand total leaf mass [kg/ha]'
        logvolume = self.ncf.createVariable('/stand/logvolume','f4',('nscens','nyrs', 'ncols',))                          # saw log volume m3/ha
        logvolume.units = 'stand volume of sawlogs [m3/ha]'
        finerootlitter = self.ncf.createVariable('/stand/finerootlitter','f4',('nscens','nyrs', 'ncols',))                      # fine root litter kg/ha/yr
        finerootlitter.units = 'stand fine root litterfall in dry biomass [kg/ha/yr] '
        n_finerootlitter = self.ncf.createVariable('/stand/n_finerootlitter','f4',('nscens','nyrs', 'ncols',))                    # N in fine root litter kg/ha/yr
        n_finerootlitter.units = 'N in stand fine root litterfall [kg/ha/yr] '
        p_finerootlitter = self.ncf.createVariable('/stand/p_finerootlitter','f4',('nscens','nyrs', 'ncols',))                    # P in fine root litter kg/ha/yr
        p_finerootlitter.units = 'P in stand fine root litterfall [kg/ha/yr] '
        k_finerootlitter =  self.ncf.createVariable('/stand/k_finerootlitter','f4',('nscens','nyrs', 'ncols',))                   # K in fine root litter kg/ha/yr
        k_finerootlitter.units = 'K in stand fine root litterfall [kg/ha/yr] '
        nonwoodylitter =   self.ncf.createVariable('/stand/nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                       # woody litter kg/ha/yr
        nonwoodylitter.units = 'stand non-woody litterfall: leaves + fineroots [kg/ha/yr]'
        n_nonwoodylitter =  self.ncf.createVariable('/stand/n_nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                      # N in woody litter kg/ha/yr
        n_nonwoodylitter.units = 'N in stand non-woody litterfall: leaves + fineroots [kg/ha/yr]'
        p_nonwoodylitter =  self.ncf.createVariable('/stand/p_nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                      # P in woody litter kg/ha/yr
        p_nonwoodylitter.units = 'P in stand non-woody litterfall: leaves + fineroots [kg/ha/yr]'
        k_nonwoodylitter =  self.ncf.createVariable('/stand/k_nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                      # K in woody litter kg/ha/yr   
        k_nonwoodylitter.units = 'K in stand non-woody litterfall: leaves + fineroots [kg/ha/yr]'    
        pulpvolume =  self.ncf.createVariable('/stand/pulpvolume','f4',('nscens','nyrs', 'ncols',))                         # pulpwood volume m3/ha
        pulpvolume.units = 'stand pulpwood volume [m3/ha]'
        stems =  self.ncf.createVariable('/stand/stems','f4',('nscens','nyrs', 'ncols',))                              # stocking, number of stems pcs/ha
        stems.units = 'stocking, number of stems [pcs/ha]'
        volume =  self.ncf.createVariable('/stand/volume','f4',('nscens','nyrs', 'ncols',))                             # total volume of the growing stock m3/ha
        volume.units ='stand volume [m3/ha]'
        woodylitter = self.ncf.createVariable('/stand/woodylitter','f4',('nscens','nyrs', 'ncols',))                   # woody litter kg/ha/yr
        woodylitter.units = 'stand woody litterfall [kg/ha/yr]'
        n_woodylitter =  self.ncf.createVariable('/stand/n_woodylitter','f4',('nscens','nyrs', 'ncols',))                      # N in woody litter kg/ha/yr
        n_woodylitter.units = 'N in stand woody litterfall [kg/ha/yr]'
        p_woodylitter = self.ncf.createVariable('/stand/p_woodylitter','f4',('nscens','nyrs', 'ncols',))                       # P in woody litter kg/ha/yr
        p_woodylitter.units = 'P in stand woody litterfall [kg/ha/yr]'
        k_woodylitter = self.ncf.createVariable('/stand/k_woodylitter','f4',('nscens','nyrs', 'ncols',))                       # K in woody litter kg/ha/yr   
        k_woodylitter.units = 'K stand woody litterfall [kg/ha/yr]'
        yi =  self.ncf.createVariable('/stand/yi','f4',('nscens','nyrs', 'ncols',))                            # yield. here same as volume
        yi.units = 'cumulative volume of growing stock, here same as volume [m3/ha] '
        previous_nut_stat = self.ncf.createVariable('/stand/previous_nut_stat','f4',('nscens','nyrs', 'ncols',)) 
        previous_nut_stat.units = 'nutrient supply/demand-ratio in previous year [-]'
        nut_stat = self.ncf.createVariable('/stand/nut_stat','f4',('nscens','nyrs', 'ncols',)) 
        nut_stat.units = 'nutrient supply/demand-ratio in ongoing year [-]'
    
    def initialize_canopy_layer(self, name):    
        basalarea = self.ncf.createVariable('/stand/' + name +'/basalarea','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        basalarea.units = 'canopy layer basal area [m2/ha]'
        biomass = self.ncf.createVariable('/stand/' + name +'/biomass','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        biomass.units = 'canopy layer dry biomass [kg/ha]'
        
        n_demand =  self.ncf.createVariable('/stand/' + name +'/n_demand','f4',('nscens','nyrs', 'ncols',))                          # stand N demand excluding leaves kg/ha/yr
        n_demand.units = 'canopy layer N demand [kg/ha/yr]'
        p_demand =  self.ncf.createVariable('/stand/' + name +'/p_demand','f4',('nscens','nyrs', 'ncols',))                           # stand N demand excluding leaves kg/ha/yr
        p_demand.units = 'canopy layer P demand [kg/ha/yr]'
        k_demand = self.ncf.createVariable('/stand/' + name +'/k_demand','f4',('nscens','nyrs', 'ncols',))                            # stand N demand excluding leaves kg/ha/yr 
        k_demand.units = 'canopy layer K demand [kg/ha/yr]' 
        hdom =  self.ncf.createVariable('/stand/' + name +'/hdom','f4',('nscens','nyrs', 'ncols',))                               # dominant height m
        hdom.units = 'canopy layer dominant height [m]'
        leafarea =  self.ncf.createVariable('/stand/' + name +'/leafarea','f4',('nscens','nyrs', 'ncols',))                           # one sided leaf area m2 m-2
        leafarea.units = 'canopy layer one-sided leaf area [m2/m2]'
        leafmass =  self.ncf.createVariable('/stand/' + name +'/leafmass','f4',('nscens','nyrs', 'ncols',))                           # leaf dry biomass kg/ha
        leafmass.units = 'canopy layer total leaf mass [kg/ha]'
        logvolume = self.ncf.createVariable('/stand/' + name +'/logvolume','f4',('nscens','nyrs', 'ncols',))                          # saw log volume m3/ha
        logvolume.units = 'canopy layer volume of sawlogs [m3/ha]'
        finerootlitter = self.ncf.createVariable('/stand/' + name +'/finerootlitter','f4',('nscens','nyrs', 'ncols',))                      # fine root litter kg/ha/yr
        finerootlitter.units = 'canopy layer fine root litterfall in dry biomass [kg/ha/yr] '
        n_finerootlitter = self.ncf.createVariable('/stand/' + name +'/n_finerootlitter','f4',('nscens','nyrs', 'ncols',))                    # N in fine root litter kg/ha/yr
        n_finerootlitter.units = 'N in canopy layer fine root litterfall [kg/ha/yr] '
        p_finerootlitter = self.ncf.createVariable('/stand/' + name +'/p_finerootlitter','f4',('nscens','nyrs', 'ncols',))                    # P in fine root litter kg/ha/yr
        p_finerootlitter.units = 'P in canopy layer fine root litterfall [kg/ha/yr] '
        k_finerootlitter =  self.ncf.createVariable('/stand/' + name +'/k_finerootlitter','f4',('nscens','nyrs', 'ncols',))                   # K in fine root litter kg/ha/yr
        k_finerootlitter.units = 'K in canopy layer fine root litterfall [kg/ha/yr] '
        nonwoodylitter =   self.ncf.createVariable('/stand/' + name +'/nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                       # woody litter kg/ha/yr
        nonwoodylitter.units = 'canopy layer non-woody litterfall: leaves + fineroots [kg/ha/yr]'
        n_nonwoodylitter =  self.ncf.createVariable('/stand/' + name +'/n_nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                      # N in woody litter kg/ha/yr
        n_nonwoodylitter.units = 'N in canopy layer non-woody litterfall: leaves + fineroots [kg/ha/yr]'
        p_nonwoodylitter =  self.ncf.createVariable('/stand/' + name +'/p_nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                      # P in woody litter kg/ha/yr
        p_nonwoodylitter.units = 'P in canopy layer non-woody litterfall: leaves + fineroots [kg/ha/yr]'
        k_nonwoodylitter =  self.ncf.createVariable('/stand/' + name +'/k_nonwoodylitter','f4',('nscens','nyrs', 'ncols',))                      # K in woody litter kg/ha/yr   
        k_nonwoodylitter.units = 'K in canopy layer non-woody litterfall: leaves + fineroots [kg/ha/yr]'    
        pulpvolume =  self.ncf.createVariable('/stand/' + name +'/pulpvolume','f4',('nscens','nyrs', 'ncols',))                         # pulpwood volume m3/ha
        pulpvolume.units = 'canopy layer pulpwood volume [m3/ha]'
        stems =  self.ncf.createVariable('/stand/' + name +'/stems','f4',('nscens','nyrs', 'ncols',))                              # stocking, number of stems pcs/ha
        stems.units = 'canopy layer stocking, number of stems [pcs/ha]'
        volume =  self.ncf.createVariable('/stand/' + name +'/volume','f4',('nscens','nyrs', 'ncols',))                             # total volume of the growing stock m3/ha
        volume.units ='canopy layer volume [m3/ha]'
        woodylitter = self.ncf.createVariable('/stand/' + name +'/woodylitter','f4',('nscens','nyrs', 'ncols',))                   # woody litter kg/ha/yr
        woodylitter.units = 'canopy layer woody litterfall [kg/ha/yr]'
        n_woodylitter =  self.ncf.createVariable('/stand/' + name +'/n_woodylitter','f4',('nscens','nyrs', 'ncols',))                      # N in woody litter kg/ha/yr
        n_woodylitter.units = 'N in canopy layer woody litterfall [kg/ha/yr]'
        p_woodylitter = self.ncf.createVariable('/stand/' + name +'/p_woodylitter','f4',('nscens','nyrs', 'ncols',))                       # P in woody litter kg/ha/yr
        p_woodylitter.units = 'P in canopy layer woody litterfall [kg/ha/yr]'
        k_woodylitter = self.ncf.createVariable('/stand/' + name +'/k_woodylitter','f4',('nscens','nyrs', 'ncols',))                       # K in woody litter kg/ha/yr   
        k_woodylitter.units = 'K canopy layer woody litterfall [kg/ha/yr]'
        yi =  self.ncf.createVariable('/stand/' + name +'/yi','f4',('nscens','nyrs', 'ncols',))                            # yield. here same as volume
        yi.units = 'canopy layer cumulative volume of growing stock, here same as volume [m3/ha] '       
        
        
        #--------Leaf dynamics variables---------------------------
        new_lmass = self.ncf.createVariable('/stand/' + name +'/new_lmass','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        new_lmass.units = 'new leaf mass in the canopy layer [kg/ha]'
        leaf_litter = self.ncf.createVariable('/stand/' + name +'/leaf_litter','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        leaf_litter.units = 'leaf litterfall in the canopy layer, dry biomass  [kg/ha/yr]'
        C_consumption = self.ncf.createVariable('/stand/' + name +'/C_consumption','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        C_consumption.units = 'C consumption in leaf growth and litterfall in the canopy layer  [kg/ha/yr]'
        leafmax = self.ncf.createVariable('/stand/' + name +'/leafmax','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        leafmax.units = 'upper limit for leaf mass in the canopy layer [kg/ha]'
        leafmin = self.ncf.createVariable('/stand/' + name +'/leafmin','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        leafmin.units = 'lower limit for leaf mass in the canopy layer  [kg/ha]'
        NPP = self.ncf.createVariable('/stand/' + name +'/NPP','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        NPP.units = 'net primary production in dry mass, subtracted by physical growth constraints [kg/ha/yr]'
        NPP_pot = self.ncf.createVariable('/stand/' + name +'/NPP_pot','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        NPP_pot.units = 'net primary production in dry mass, without  physical growth constraints [kg/ha/yr]'
        Nleafdemand = self.ncf.createVariable('/stand/' + name +'/Nleafdemand','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        Nleafdemand.units = 'N demand of leaf production [kg/ha/yr]'
        Nleaf_litter = self.ncf.createVariable('/stand/' + name +'/Nleaf_litter','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        Nleaf_litter.units = 'N in leaf litterfall [kg/ha/yr]'
        Pleafdemand = self.ncf.createVariable('/stand/' + name +'/Pleafdemand','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        Pleafdemand.units = 'P demand of leaf production [kg/ha/yr]'
        Pleaf_litter = self.ncf.createVariable('/stand/' + name +'/Pleaf_litter','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        Pleaf_litter.units = 'P in leaf litterfall [kg/ha/yr]'
        Kleafdemand = self.ncf.createVariable('/stand/' + name +'/Kleafdemand','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        Kleafdemand.units = 'K demand of leaf production [kg/ha/yr]'
        Kleaf_litter = self.ncf.createVariable('/stand/' + name +'/Kleaf_litter','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        Kleaf_litter.units = 'K in leaf litterfall [kg/ha/yr]'

    def initialize_gv(self):   
        
        gv_tot = self.ncf.createVariable('/groundvegetation/gv_tot','f4',('nscens','nyrs', 'ncols',))                           # stand basal area m2/ha                          
        gv_tot.units = 'total groundvegetation biomass in dry matter [kg/ha]'
        gv_field = self.ncf.createVariable('/groundvegetation/gv_field','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        gv_field.units = 'groud vegetation field layer dry biomass [kg/ha]'
        gv_bot = self.ncf.createVariable('/groundvegetation/gv_bot','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        gv_bot.units = 'groud vegetation bottom layer dry biomass [kg/ha]'
        gv_leafmass = self.ncf.createVariable('/groundvegetation/gv_leafmass','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        gv_leafmass.units = 'groud vegetation leaf mass dry biomass [kg/ha]'        
        ds_litterfall = self.ncf.createVariable('/groundvegetation/ds_litterfall','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        ds_litterfall.units = 'dwarf shrubs litterfall dry matter  [kg/ha/yr]'
        h_litterfall = self.ncf.createVariable('/groundvegetation/h_litterfall','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        h_litterfall.units = 'herbs and sedges litterfall dry matter  [kg/ha/yr]'    
        s_litterfall = self.ncf.createVariable('/groundvegetation/s_litterfall','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        s_litterfall.units = 'sphagnum mosses litterfall dry matter  [kg/ha/yr]'
        n_litter = self.ncf.createVariable('/groundvegetation/n_litter','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        n_litter.units = 'N in ground vegetation litterfall dry matter  [kg/ha/yr]'
        p_litter = self.ncf.createVariable('/groundvegetation/p_litter','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        p_litter.units = 'P in ground vegetation litterfall dry matter  [kg/ha/yr]'
        k_litter = self.ncf.createVariable('/groundvegetation/k_litter','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        k_litter.units = 'K in ground vegetation litterfall dry matter  [kg/ha/yr]'
        n_gv = self.ncf.createVariable('/groundvegetation/n_gv','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        n_gv.units = 'N in ground vegetation  [kg/ha]'
        p_gv = self.ncf.createVariable('/groundvegetation/p_gv','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        p_gv.units = 'P in ground vegetation  [kg/ha]'
        k_gv = self.ncf.createVariable('/groundvegetation/k_gv','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        k_gv.units = 'K in ground vegetation  [kg/ha]'
        nup = self.ncf.createVariable('/groundvegetation/nup','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        nup.units = 'N uptake ground vegetation  [kg/ha/yr]'
        pup = self.ncf.createVariable('/groundvegetation/pup','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        pup.units = 'P uptake ground vegetation  [kg/ha/yr]'
        kup = self.ncf.createVariable('/groundvegetation/kup','f4',('nscens','nyrs', 'ncols',))                             # stand dry biomass kg/ha
        kup.units = 'K uptake ground vegetation  [kg/ha/yr]'
        
    def initialize_esom(self, substance):

        L0L = self.ncf.createVariable('/esom/' + substance +'/L0L','f4',('nscens','nyrs', 'ncols',))                                                     
        L0L.units = 'L0L input of leaf and fine root litter [kg/m2/yr]'
        L0W = self.ncf.createVariable('/esom/' + substance +'/L0W','f4',('nscens','nyrs', 'ncols',))                                                     
        L0W.units = 'L0W input of branch and coarse root, i.e. woody litter [kg/m2/yr]'
        LL = self.ncf.createVariable('/esom/' + substance +'/LL','f4',('nscens','nyrs', 'ncols',))                                                    
        LL.units = 'storage of leaf litter [kg/m2]'
        LW = self.ncf.createVariable('/esom/' + substance +'/LW','f4',('nscens','nyrs', 'ncols',))                                                
        LW.units = 'storage of woody litter [kg/m2]'   
        FL = self.ncf.createVariable('/esom/' + substance +'/FL','f4',('nscens','nyrs', 'ncols',))                                                    
        FL.units = 'FL storage of leaf F material [kg/m2]'
        FW = self.ncf.createVariable('/esom/' + substance +'/FW','f4',('nscens','nyrs', 'ncols',))                                                   
        FW.units = 'FW storage of woody F material [kg/m2]' 
        H = self.ncf.createVariable('/esom/' + substance +'/H','f4',('nscens','nyrs', 'ncols',))                                             
        H.units = 'H storage of H material [kg/m2]'
        P1 = self.ncf.createVariable('/esom/' + substance +'/P1','f4',('nscens','nyrs', 'ncols',))                                                     
        P1.units = 'storage of top peat in dry matter [kg/m2]'
        P2 = self.ncf.createVariable('/esom/' + substance +'/P2','f4',('nscens','nyrs', 'ncols',))                                                  
        P2.units = 'storage of middle layer peat  [kg/m2]'
        P1 = self.ncf.createVariable('/esom/' + substance +'/P3','f4',('nscens','nyrs', 'ncols',))                                                    
        P1.units = 'storage of bottom peat [kg/m2]'
        out = self.ncf.createVariable('/esom/' + substance +'/out','f4',('nscens','nyrs', 'ncols',))
        out.units = 'annual output of substance [kg/ha/yr]'
        if substance == 'Mass':
            c_to_atm =  self.ncf.createVariable('/esom/' + substance +'/c_to_atm','f4',('nscens','nyrs', 'ncols',))
            c_to_atm.units = 'gaseous outflux of C in decomposition [kg/ha/yr]'
            co2 = self.ncf.createVariable('/esom/' + substance +'/co2','f4',('nscens','nyrs', 'ncols',))
            co2.units = 'co2 release [kg/ha/yr]'
            doc = self.ncf.createVariable('/esom/' + substance +'/doc','f4',('nscens','nyrs', 'ncols',))
            doc.units = 'doc release [kg/ha/yr]'
            lmwdoc = self.ncf.createVariable('/esom/' + substance +'/lmwdoc','f4',('nscens','nyrs', 'ncols',))
            lmwdoc.units = 'lmwdoc release [kg/ha/yr]'
            
    def initialize_cpy(self):
        
        interc = self.ncf.createVariable('/cpy/interc','f4',('nscens','ndays', 'ncols',))                          
        interc.units = 'interception storage [mm]'
        interc_yr = self.ncf.createVariable('/cpy/interc_yr','f4',('nscens','ndays', 'ncols',))                          
        interc_yr.units = 'mean interception storage [mm]'     
        evap = self.ncf.createVariable('/cpy/evap','f4',('nscens','ndays', 'ncols',))                          
        evap.units = 'evaporation [mm/day]'
        ET = self.ncf.createVariable('/cpy/ET','f4',('nscens','ndays', 'ncols',))                          
        ET.units = 'evapotranspiration [mm/day]'
        ET_yr = self.ncf.createVariable('/cpy/ET_yr','f4',('nscens','nyrs', 'ncols',))                          
        ET_yr.units = 'annual evapotranspiration [mm/yr]'
        transpi = self.ncf.createVariable('/cpy/transpi','f4',('nscens','ndays', 'ncols',))                          
        transpi.units = 'transpiration [mm/day]'
        transpi_yr = self.ncf.createVariable('/cpy/transpi_yr','f4',('nscens','nyrs', 'ncols',))                          
        transpi_yr.units = 'annual transpiration [mm/yr]'
        efloor = self.ncf.createVariable('/cpy/efloor','f4',('nscens','ndays', 'ncols',))                          
        efloor.units = 'evaporation from forest floor [mm/day]'
        efloor_yr = self.ncf.createVariable('/cpy/efloor_yr','f4',('nscens','nyrs', 'ncols',))                          
        efloor_yr.units = 'annual evaporation from forest floor [mm/yr]'
        SWE = self.ncf.createVariable('/cpy/SWE','f4',('nscens','ndays', 'ncols',))                          
        SWE.units = 'snow water equivalent [mm]'
        SWEmax = self.ncf.createVariable('/cpy/SWEmax','f4',('nscens','nyrs', 'ncols',))                          
        SWEmax.units = 'Annual maximum snow water equivalent [mm]'
        
    def initialize_strip(self, strip):
        elevation = self.ncf.createVariable('/strip/elevation','f4',('ncols',))
        elevation.units = 'soil surface elevation with respect to fixed datum, m'
        kmap = self.ncf.createVariable('/strip/kmap','f4',('ncols','nLyrs',))
        kmap.units = 'hydraulic conductivity map of the profile [m/s]'
        dwt = self.ncf.createVariable('/strip/dwt','f4',('nscens','ndays', 'ncols',))                          
        dwt.units = 'water table m below soil surface [m, negative down]'
        dwtyr = self.ncf.createVariable('/strip/dwtyr','f4',('nscens','nyrs', 'ncols',))                          
        dwtyr.units = 'annual water table m below soil surface [m, negative down]'
        dwtyr_growingseason = self.ncf.createVariable('/strip/dwtyr_growingseason','f4',('nscens','nyrs', 'ncols',))                          
        dwtyr_growingseason.units = 'annual growing season (May 1 - October 31) water table m below soil surface [m, negative down]'
        dwtyr_latesummer = self.ncf.createVariable('/strip/dwtyr_latesummer','f4',('nscens','nyrs', 'ncols',))                          
        dwtyr_latesummer.units = 'annual late summer (July 1 - August 31) water table m below soil surface [m, negative down]'     
        deltas = self.ncf.createVariable('/strip/deltas','f4',('nscens','nyrs', 'ncols',))                          
        deltas.units = 'annual flux through soil surface [m/yr]'
        H = self.ncf.createVariable('/strip/H','f4',('nscens','ndays', 'ncols',))                          
        H.units = 'height of water teble to common datum [m]'
        roff = self.ncf.createVariable('/strip/roff','f4',('nscens','ndays',))                          
        roff.units = 'total runoff including west and east ditches and surface runoff [m/day]'
        roffwest = self.ncf.createVariable('/strip/roffwest','f4',('nscens','ndays',))                          
        roffwest.units = 'runoff along the west ditch [m/day]'
        roffeast = self.ncf.createVariable('/strip/roffeast','f4',('nscens','ndays',))                          
        roffeast.units = 'runoff along the east ditch [m/day]'
        surfacerunoff = self.ncf.createVariable('/strip/surfacerunoff','f4',('nscens','ndays', 'ncols',))                          
        surfacerunoff.units = 'surface runoff generated from each column, has to be averaged through columns to be at the same unit than west east runoff [m/area/day]'

        residencetime = self.ncf.createVariable('/strip/residencetime','f4',('nscens','nyrs', 'ncols',))                          
        residencetime.units = 'residence time of water from column to ditch [days]'
        ixwest = self.ncf.createVariable('/strip/ixwest','f4',('nscens','nyrs', 'ncols',))                          
        ixwest.units = 'columns discharging to the west ditch, [index number of a column]'
        ixeast = self.ncf.createVariable('/strip/ixeast','f4',('nscens','nyrs', 'ncols',))                          
        ixeast.units = 'columns discharging to the east ditch, [index number of a column]'
        
        self.ncf['/strip/elevation'][:] = strip.ele
        self.ncf['/strip/kmap'][:, :] = strip.Kmap
        
    def initialize_temperature(self):
        T = self.ncf.createVariable('/temperature/T','f4',('nscens','ndays', 'nLyrs',))                          
        T.units = 'peat temperature at different depth of the profile [deg C]'

    def initialize_doc(self):
        
        DOC = self.ncf.createVariable('/doc/DOC','f4',('nscens','nyrs', 'ncols',))                          
        DOC.units = 'total DOC release [kg/ha/yr]'
        LMW = self.ncf.createVariable('/doc/LMW','f4',('nscens','nyrs', 'ncols',))                          
        LMW.units = 'light molecular weight DOC release [kg/ha/yr]'
        HMW = self.ncf.createVariable('/doc/HMW','f4',('nscens','nyrs', 'ncols',))                          
        HMW.units = 'high molecular weight DOC release [kg/ha/yr]'
 
    def initialize_methane(self):
        
        ch4 = self.ncf.createVariable('/methane/ch4','f4',('nscens','nyrs', 'ncols',))                          
        ch4.units = 'CH4 emission [kg/ha/yr]'
        ch4_in_co2 = self.ncf.createVariable('/methane/ch4_in_co2','f4',('nscens','nyrs', 'ncols',))                          
        ch4_in_co2.units = 'methane emissions in CO2 equivalents [kg/ha/yr]'

    def initialize_fertilization(self):
        
        fn = self.ncf.createVariable('/fertilization/n_release','f4',('nscens','nyrs', 'ncols',))                          
        fn.units = 'N release in fertilization [kg/ha/yr]'
        fp = self.ncf.createVariable('/fertilization/p_release','f4',('nscens','nyrs', 'ncols',))                          
        fp.units = 'P release in fertilization [kg/ha/yr]'
        fk = self.ncf.createVariable('/fertilization/k_release','f4',('nscens','nyrs', 'ncols',))                          
        fk.units = 'K release in fertilization [kg/ha/yr]'

    def initialize_export(self):
        hmwtoditch = self.ncf.createVariable('/export/hmwtoditch','f4',('nscens','nyrs', 'ncols',))
        hmwtoditch.units = 'hmwdoc reaching the ditch, reduced by the biodegradation [kg/ha/yr]'
        lmwtoditch = self.ncf.createVariable('/export/lmwtoditch','f4',('nscens','nyrs', 'ncols',))
        lmwtoditch.units = 'lmwdoc reaching the ditch, reduced by the biodegradation [kg/ha/yr]'

        hmwdoc_to_west = self.ncf.createVariable('/export/hmwdoc_to_west','f4',('nscens','nyrs',))                          
        hmwdoc_to_west.units = 'HMW doc export from a column reaching the west ditch [kg/ha/yr]'
        hmwdoc_to_east = self.ncf.createVariable('/export/hmwdoc_to_east','f4',('nscens','nyrs',))                          
        hmwdoc_to_east.units = 'HMW doc export from a column reaching the east ditch [kg/ha/yr]'

        lmwdoc_to_west = self.ncf.createVariable('/export/lmwdoc_to_west','f4',('nscens','nyrs',))                          
        lmwdoc_to_west.units = 'LMW doc export from a column reaching the west ditch [kg/ha/yr]'
        lmwdoc_to_east = self.ncf.createVariable('/export/lmwdoc_to_east','f4',('nscens','nyrs',))                          
        lmwdoc_to_east.units = 'LMW doc export from a column reaching the east ditch [kg/ha/yr]'
        

    def write_stand(self, scen, year, stand):
        self.ncf['stand']['basalarea'][scen, year, :] = stand.basalarea
        self.ncf['stand']['biomass'][scen, year, :] = stand.biomass
        self.ncf['stand']['n_demand'][scen, year, :] = stand.n_demand
        self.ncf['stand']['p_demand'][scen, year, :] = stand.p_demand
        self.ncf['stand']['k_demand'][scen, year, :] = stand.k_demand
        self.ncf['stand']['hdom'][scen, year, :] = stand.hdom
        self.ncf['stand']['leafarea'][scen, year, :] = stand.leafarea
        self.ncf['stand']['leafmass'][scen, year, :] = stand.leafmass
        self.ncf['stand']['logvolume'][scen, year, :] = stand.logvolume
        self.ncf['stand']['finerootlitter'][scen, year, :] = stand.finerootlitter
        self.ncf['stand']['n_finerootlitter'][scen, year, :] = stand.n_finerootlitter
        self.ncf['stand']['p_finerootlitter'][scen, year, :] = stand.p_finerootlitter
        self.ncf['stand']['k_finerootlitter'][scen, year, :] = stand.k_finerootlitter
        self.ncf['stand']['nonwoodylitter'][scen, year, :] = stand.nonwoodylitter   
        self.ncf['stand']['pulpvolume'][scen, year, :] = stand.pulpvolume
        self.ncf['stand']['stems'][scen, year, :] = stand.stems
        self.ncf['stand']['volume'][scen, year, :] = stand.volume
        self.ncf['stand']['woodylitter'][scen, year, :] = stand.woodylitter
        self.ncf['stand']['n_woodylitter'][scen, year, :] = stand.n_woodylitter
        self.ncf['stand']['p_woodylitter'][scen, year, :] = stand.p_woodylitter
        self.ncf['stand']['k_woodylitter'][scen, year, :] = stand.k_woodylitter
        self.ncf['stand']['yi'][scen, year, :] = stand.yi
        self.ncf['stand']['previous_nut_stat'][scen, year, :] = stand.previous_nut_stat
        self.ncf['stand']['nut_stat'][scen, year, :] = stand.nut_stat
        


    def write_canopy_layer(self, scen, year, name, layer):
        self.ncf['stand'][name]['basalarea'][scen, year, :] = layer.basalarea
        self.ncf['stand'][name]['biomass'][scen, year, :] = layer.biomass
        self.ncf['stand'][name]['n_demand'][scen, year, :] = layer.n_demand
        self.ncf['stand'][name]['p_demand'][scen, year, :] = layer.p_demand
        self.ncf['stand'][name]['k_demand'][scen, year, :] = layer.k_demand
        self.ncf['stand'][name]['hdom'][scen, year, :] = layer.hdom
        self.ncf['stand'][name]['leafarea'][scen, year, :] = layer.leafarea
        self.ncf['stand'][name]['leafmass'][scen, year, :] = layer.leafmass
        self.ncf['stand'][name]['logvolume'][scen, year, :] = layer.logvolume
        self.ncf['stand'][name]['finerootlitter'][scen, year, :] = layer.finerootlitter
        self.ncf['stand'][name]['n_finerootlitter'][scen, year, :] = layer.n_finerootlitter
        self.ncf['stand'][name]['p_finerootlitter'][scen, year, :] = layer.p_finerootlitter
        self.ncf['stand'][name]['k_finerootlitter'][scen, year, :] = layer.k_finerootlitter
        self.ncf['stand'][name]['nonwoodylitter'][scen, year, :] = layer.nonwoodylitter   
        self.ncf['stand'][name]['pulpvolume'][scen, year, :] = layer.pulpvolume
        self.ncf['stand'][name]['stems'][scen, year, :] = layer.stems
        self.ncf['stand'][name]['volume'][scen, year, :] = layer.volume
        self.ncf['stand'][name]['woodylitter'][scen, year, :] = layer.woodylitter
        self.ncf['stand'][name]['n_woodylitter'][scen, year, :] = layer.n_woodylitter
        self.ncf['stand'][name]['p_woodylitter'][scen, year, :] = layer.p_woodylitter
        self.ncf['stand'][name]['k_woodylitter'][scen, year, :] = layer.k_woodylitter
        self.ncf['stand'][name]['yi'][scen, year, :] = layer.yi
        
        self.ncf['stand'][name]['new_lmass'][scen, year, :] = layer.new_lmass
        self.ncf['stand'][name]['leaf_litter'][scen, year, :] = layer.leaf_litter
        self.ncf['stand'][name]['C_consumption'][scen, year, :] = layer.C_consumption
        self.ncf['stand'][name]['leafmax'][scen, year, :] = layer.leafmax
        self.ncf['stand'][name]['leafmin'][scen, year, :] = layer.leafmin
        
        self.ncf['stand'][name]['NPP'][scen, year, :] = layer.NPP
        self.ncf['stand'][name]['NPP_pot'][scen, year, :] = layer.NPP_pot
        self.ncf['stand'][name]['Nleafdemand'][scen, year, :] = layer.Nleafdemand
        self.ncf['stand'][name]['Nleaf_litter'][scen, year, :] = layer.Nleaf_litter
        self.ncf['stand'][name]['Pleafdemand'][scen, year, :] = layer.Pleafdemand
        self.ncf['stand'][name]['Pleaf_litter'][scen, year, :] = layer.Pleaf_litter
        self.ncf['stand'][name]['Kleafdemand'][scen, year, :] = layer.Kleafdemand
        self.ncf['stand'][name]['Kleaf_litter'][scen, year, :] = layer.Kleaf_litter
         
        
    def write_groundvegetation(self, scen, year, gv):
        self.ncf['groundvegetation']['gv_tot'][scen, year, :] = gv.gv_tot 
        self.ncf['groundvegetation']['gv_field'][scen, year, :] = gv.gv_field
        self.ncf['groundvegetation']['gv_bot'][scen, year, :] = gv.gv_bot
        self.ncf['groundvegetation']['gv_leafmass'][scen, year, :] = gv.gv_leafmass        
        self.ncf['groundvegetation']['ds_litterfall'][scen, year, :] = gv.ds_litterfall
        self.ncf['groundvegetation']['h_litterfall'][scen, year, :] = gv.h_litterfall
        self.ncf['groundvegetation']['s_litterfall'][scen, year, :] = gv.s_litterfall
        self.ncf['groundvegetation']['n_litter'][scen, year, :] = gv.n_litter
        self.ncf['groundvegetation']['p_litter'][scen, year, :] = gv.p_litter
        self.ncf['groundvegetation']['k_litter'][scen, year, :] = gv.k_litter       
        self.ncf['groundvegetation']['n_gv'][scen, year, :] = gv.n_gv
        self.ncf['groundvegetation']['p_gv'][scen, year, :] = gv.p_gv
        self.ncf['groundvegetation']['k_gv'][scen, year, :] = gv.k_gv
        self.ncf['groundvegetation']['nup'][scen, year, :] = gv.nup
        self.ncf['groundvegetation']['pup'][scen, year, :] = gv.pup
        self.ncf['groundvegetation']['kup'][scen, year, :] = gv.kup
        
    def write_esom(self, scen, year, substance, esom, inivals = False):
        if inivals:
            self.ncf['esom'][substance]['L0L'][scen, year, :] = np.zeros(esom.y)
            self.ncf['esom'][substance]['L0W'][scen, year, :] = np.zeros(esom.y)
        else:
            self.ncf['esom'][substance]['L0L'][scen, year, :] = esom.nonwoodylitter*10000.
            self.ncf['esom'][substance]['L0W'][scen, year, :] = esom.woodylitter*10000.
            
        self.ncf['esom'][substance]['LL'][scen, year, :] = esom.M[:,:, 2]*10000.
        self.ncf['esom'][substance]['LW'][scen, year, :] = esom.M[:,:, 3]*10000.
        self.ncf['esom'][substance]['FL'][scen, year, :] = esom.M[:,:, 4]*10000.   
        self.ncf['esom'][substance]['FW'][scen, year, :] = esom.M[:,:, 5]*10000.
        self.ncf['esom'][substance]['H'][scen, year, :] = esom.M[:,:, 6]*10000.
        self.ncf['esom'][substance]['P1'][scen, year, :] = esom.M[:,:, 7]*10000.
        self.ncf['esom'][substance]['P2'][scen, year, :] = esom.M[:,:, 8]*10000.
        self.ncf['esom'][substance]['P3'][scen, year, :] = esom.M[:,:, 9]*10000.
        if inivals:
            self.ncf['esom'][substance]['out'][scen, year, :] = np.zeros(esom.y)
        else: 
            self.ncf['esom'][substance]['out'][scen, year, :] = esom.out        
            if substance == 'Mass':
                mass_to_c = 0.5
                self.ncf['esom'][substance]['c_to_atm'][scen, year, :] = esom.out* 1/1.05 * mass_to_c  
                self.ncf['esom'][substance]['co2'][scen, year, :] = esom.out* 1/1.05 * mass_to_c * 44./12 
                self.ncf['esom'][substance]['doc'][scen, year, :] = esom.out * 1/1.05 * 0.05 * mass_to_c 
                self.ncf['esom'][substance]['lmwdoc'][scen, year, :] = esom.out * 1/1.05 * 0.05 * 0.04 * mass_to_c 
      
        

    def write_cpy(self, scen, start, days, yr, interc, evap, ET, transpi, efloor, SWE):
        
        self.ncf['cpy']['interc'][scen, start :start+days, :] = interc
        self.ncf['cpy']['evap'][scen, start :start+days, :] = evap
        self.ncf['cpy']['ET'][scen, start :start+days, :] = ET
        self.ncf['cpy']['transpi'][scen, start :start+days, :] = transpi
        self.ncf['cpy']['efloor'][scen, start :start+days, :] = efloor
        self.ncf['cpy']['SWE'][scen, start :start+days, :] = SWE
        self.ncf['cpy']['ET_yr'][scen, yr, :] = np.sum(ET[:,:], axis=0)
        self.ncf['cpy']['transpi_yr'][scen, yr, :] = np.sum(transpi[:,:], axis=0)
        self.ncf['cpy']['efloor_yr'][scen, yr, :] = np.sum(efloor[:,:], axis=0)
        self.ncf['cpy']['SWEmax'][scen, yr, :] = np.max(SWE[:,:], axis=0)
        self.ncf['cpy']['interc_yr'][scen, yr, :] = np.mean(interc[:,:], axis=0)


    def write_temperature(self, scen, start, days, T):
        self.ncf['temperature']['T'][scen, start :start+days, :] = T

    def write_strip(self, scen, start, days, yr, year, dfwt, stpout, outpara, stp):
                
        startdate = '-' + str(outpara['startmonth']) + '-' + str(outpara['startday'])  
        enddate = '-' +str(outpara['endmonth']) + '-' + str(outpara['endday'])  
        self.ncf['strip']['dwt'][scen, start :start+days, :] = stpout['dwts'][scen,start:start+days,:]
        self.ncf['strip']['dwtyr'][scen, year, :] = dfwt.mean(axis = 0)
        self.ncf['strip']['dwtyr_latesummer'][scen, year, :] = dfwt[str(yr)+startdate: str(yr)+enddate].mean()
        self.ncf['strip']['dwtyr_growingseason'][scen, year, :] = dfwt[str(yr)+'-05-01': str(yr)+'-10-31'].mean()      
        self.ncf['strip']['H'][scen, start :start+days, :] = stpout['hts'][scen,start:start+days,:]
        self.ncf['strip']['roff'][scen, start :start+days] = stpout['runoff'][scen,start:start+days] 
        self.ncf['strip']['roffwest'][scen, start :start+days] = stpout['runoffwest'][scen,start:start+days] 
        self.ncf['strip']['roffeast'][scen, start :start+days] = stpout['runoffeast'][scen,start:start+days] 
        self.ncf['strip']['surfacerunoff'][scen, start :start+days] = stpout['surfacerunoff'][scen,start:start+days,:] 
        self.ncf['strip']['deltas'][scen, year, :] = np.sum(stpout['deltas'][scen,start:start+days,:], axis=0)
        self.ncf['strip']['residencetime'][scen, year, :] = stp.residence_time
        
        self.ncf['strip']['ixwest'][scen, year, :] = np.zeros(stp.n) - 1.
        nwest = np.ravel(stp.ixwest)
        self.ncf['strip']['ixwest'][scen, year, : len(nwest)] = stp.ixwest 
        self.ncf['strip']['ixeast'][scen, year, :] = np.zeros(stp.n) - 1.
        neast = np.ravel(stp.ixeast)
        self.ncf['strip']['ixeast'][scen, year, :len(neast)] = stp.ixeast 
        
    
    def write_doc(self, scen, year, DOC, HMW):
        self.ncf['doc']['DOC'][scen, year, :] = DOC
        self.ncf['doc']['HMW'][scen, year, :] = HMW
        self.ncf['doc']['LMW'][scen, year, :] = DOC-HMW
    
    def write_methane(self, scen, year, ch4):
        self.ncf['methane']['ch4'][scen, year, :] = ch4
        self.ncf['methane']['ch4_in_co2'][scen, year, :] = ch4*54.0

    def write_fertilization(self, scen, year, ferti):
        self.ncf['fertilization']['n_release'][scen, year, :] = ferti.release['N']
        self.ncf['fertilization']['p_release'][scen, year, :] = ferti.release['P']
        self.ncf['fertilization']['k_release'][scen, year, :] = ferti.release['K']

    def write_export(self, scen, year, esmass):
        
        self.ncf['export']['hmwtoditch'][scen, year, :] = esmass.hmwtoditch
        self.ncf['export']['lmwtoditch'][scen, year, :] = esmass.lmwtoditch
        self.ncf['export']['hmwdoc_to_west'][scen, year] = esmass.hmw_to_west
        self.ncf['export']['hmwdoc_to_east'][scen, year] = esmass.hmw_to_east
        self.ncf['export']['lmwdoc_to_west'][scen, year] = esmass.lmw_to_west
        self.ncf['export']['lmwdoc_to_east'][scen, year] = esmass.lmw_to_east
        
