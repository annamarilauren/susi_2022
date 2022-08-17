# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:32:24 2022

@author: alauren
"""
import numpy as np
from pyproj import CRS, Transformer


class Gvegetation():
    def __init__(self, n, lat, lon, sfc, species):
        
        """
        INPUT: 
            n - number of computation nodes between ditches
            lat - latitude in [units]
            lon - longitude in [units]
            barea 0,1 - stand basal area [m2], start, end
            stems 0,1 - number of stems per ha, start, end
            yi 0,1 - stand volume as time series, start, end
            sp - dominant tree species - 1 pine, 2 spruce, 3 birch
            ts - temperature sum [dd]
            simtime  - simulation time in yrs
            sfc - site fertility class, 
            standAge in yrs
        SOURCE:
            Muukkonen & Makipaa, 2006. Bor.Env.Res. 11, 355-369.\n
        AUTHOR:
            Samuli Launiainen 18.06.2014, Modified for array operations by Ari Laurén 13.4.2020 \n
        NOTE:
             Multi-regression models not yet tested!
             In model equations independent variables named differently to M&M (2006): here x[0] = z1, x[1]=z2, ... x[7]=z8 and x[8]=z10\n
             \n
    
             Site nutrient fertility class (sfc) at mires:\n
                 1: herb-rich hw-spruce swamps, pine mires, fens, 
                 2: V.myrtillus / tall sedge spruce swamps, tall sedge pine fens, tall sedge fens,
                 3: Carex clobularis / V.vitis-idaea swamps, Carex globularis pine swamps, low sedge (oligotrophic) fens,
                 4: Low sedge, dwarf-shrub & cottongrass pine bogs, ombo-oligotrophic bogs,
                 5: S.fuscum pine bogs, ombotrophic and S.fuscum low sedge bogs.
             Drainage status x[8] at mires (Paavilainen & Paivanen, 1995):
                 1: undrained
                 2: Recently draines, slight effect on understory veg., no effect on stand
                 3: Transforming drained mires, clear effect on understory veg and stand
                 4: Transforming drained mires, veget. resembles upland forest site type, tree-stand forest-like.
      
        """
        
        #------------------------ GENERAL PARAMETERS -----------------------------------------------------------------------------------
        
        """
        Computed:
           - total biomass and bottom layer; field layer is gained as a difference of tot and bottom layer (more cohrent results)
           - N and P storage in the each pixel
           - annual use of N and P due to litterfall
        Muukkonen Mäkipää 2005 upland sites: field layer contains dwarf shrubs and (herbs + grasses), see Fig 1
            share     dwarf shrubs     herbs 
            - Pine       91%            9%
            - Spruce     71%            29%
            - broad l    38%            62%
        Peatland sites (assumption):
            share      dwarf shrubs    herbs
            - Pine bogs    90%          10%
            - Spruce mires 50%          50%
        Palviainen et al. 2005 Ecol Res (2005) 20: 652–660, Table 2
        Nutrient concentrations for
                                N              P           K
            - Dwarf shrubs      1.2%         1.0 mg/g     4.7 mg/g
            - herbs & grasses   1.8%         2.0 mg/g    15.1 mg/g
            - upland mosses     1.25%        1.4 mg/g     4.3 mg/g
        Nutrient concentrations for sphagna (FIND):
                                N              P     for N :(Bragazza et al Global Change Biology (2005) 11, 106–114, doi: 10.1111/j.1365-2486.2004.00886.x)
            - sphagnum          0.6%           1.4 mg/g     (Palviainen et al 2005)   
        Annual litterfall proportions from above-ground biomass (Mälkönen 1974, Tamm 1953):
            - Dwarf shrubs          0.2
            - herbs & grasses        1
            - mosses                0.3
            Tamm, C.O. 1953. Growth, yield and nutrition in carpets of a forest moss (Hylocomium splendens). Meddelanden från Statens Skogsforsknings Institute 43 (1): 1-140.
        
        We assume retranslocation of N and P away from senescing tissues before litterfall:
                                N           P
            - Dwarf shrubs     0.5         0.5
            - Herbs & grasses  0.5         0.5
            - mossess          0.0         0.0
        
        Turnover of total biomass including the belowground biomass is assumed to be 1.2 x above-ground biomass turnover
        
        """

        self.fl_share = {'description': 'share of dwarf shrubs (ds) and herbas & grasses (h) from field layer biomass, kg kg-1',
                    'pine_upland':{'ds': 0.91, 'h': 0.09}, 'spruce_upland':{'ds': 0.71, 'h': 0.29}, 
                    'broadleaved_upland':{'ds': 0.38, 'h': 0.62}, 'spruce_mire':{'ds': 0.50, 'h': 0.50}, 
                    'pine_bog':{'ds': 0.90, 'h': 0.10}}
        self.nut_con ={'description': 'nutrient concentration of dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit mg/g',
                  'ds':{'N':12.0, 'P':1.0, 'K': 4.7}, 'h':{'N':18.0, 'P':2.0, 'K': 15.1}, 'um':{'N':12.5, 'P':1.4, 'K':4.3}, 
                  's':{'N':6.0, 'P':1.4, 'K':4.3}}
        self.lit_share = {'description': 'share of living biomass that is lost as litter annually for dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit: kg kg-1',
                   'ds': 0.2, 'h': 0.5, 'um': 0.3, 's': 0.3}
        self.green_share = {'description': 'share of green mass from the total biomass of dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit kg/kg',
                       'ds': 0.2, 'h': 0.5, 'um': 0.3, 's': 0.3}
        self.retrans ={'description': 'share of nutrients retranslocated before litterfallfor dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit: kg kg-1',
                  'ds': {'N':0.69, 'P':0.73, 'K':0.87},'h': {'N':0.69, 'P':0.73, 'K':0.87}, 
                  'um': {'N':0.5, 'P':0.5, 'K':0.5},'s': {'N':0.5, 'P':0.5, 'K':0.5}}
        # retrans for pine {'N': 0.69, 'P': 0.73, 'K':0.87}
        #ATTN: changes 4.1.2021 fl_share pine_bog vs spruce_mire (vice versa)
        #check, and change back restranslocation for herbs 
        
        self.fl_to_total_turnover = 1.2    # converts the turnover of above-ground bionmass to total including root turnover
        self.fl_above_to_total = 1.7       # converts aboveground biomass to total biomass 
            
        # ------------------- SITE PARAMETERS-------------------------------------------------------
        self.n = n 
        self. sfc = sfc  
        self.dem = np.ones(n) * 80.  ## x2 surface elevation m asl

        #------------- classify and map pixels-------------------------------------------------------- 

        self.ix_spruce_mire = np.where(np.equal(species, 2))
        self.ix_pine_bog = np.where(np.equal(species, 1))
        self.ix_open_peat = np.where(np.equal(species, 4))
        self.drain_s = 4      # x8 drainage status, default value 4

        #---------------------------------------    
        inProj = CRS('epsg:3067')
        outProj = CRS('epsg:4326')
        transformer = Transformer.from_crs(inProj, outProj)
        self.latitude, self.longitude = transformer.transform(lon, lat)

        self.reset_domain()        

    def reset_domain(self):

        #--------- create output arrays -----------------------------
        self.gv_tot = np.zeros(self.n)                           # Ground vegetation mass kg ha-1
        self.gv_change = np.zeros(self.n)                        # biomass change during time step kg ha-1 yr-1
        self.gv_field = np.zeros(self.n)                         # Field layer vegetation mass
        self.gv_bot = np.zeros(self.n)                           # Bottom layer vegetation mass
        self.gv_leafmass = np.zeros(self.n)                      # Leaf mass in ground vegetation kg ha-1
        self.ds_litterfall = np.zeros(self.n)                    # dwarf shrub litterfall kg ha-1 yr-1
        self.h_litterfall = np.zeros(self.n)                     # herbs and grasses litterfall kg ha-1 yr-1
        self.s_litterfall = np.zeros(self.n)                     # sphagnum mosses litterfall kg ha-1 yr-1
        self.n_litter = np.zeros(self.n)                       # N uptake due to litterfall kg ha-1 yr-1
        self.p_litter = np.zeros(self.n)                       # P uptake due to litterfall kg ha-1 yr-1
        self.k_litter = np.zeros(self.n)                       # K uptake due to litterfall kg ha-1 yr-1
        self.n_gv = np.zeros(self.n)                             # N in ground vegetation kg ha-1
        self.p_gv = np.zeros(self.n)                             # P in ground vegetation kg ha-1
        self.k_gv = np.zeros(self.n)                             # K in ground vegetation kg ha-1

        
    def gv_biomass_and_nutrients(self,  ts, vol, Nstems, ba, age):   
        #--------------- nutrient contents in vegetation-----------------------
        gv_tot = np.zeros(self.n)                           # Ground vegetation mass kg ha-1
        gv_field = np.zeros(self.n)                         # Field layer vegetation mass
        gv_bot = np.zeros(self.n)                           # Bottom layer vegetation mass
        gv_leafmass = np.zeros(self.n)                      # Leaf mass in ground vegetation kg ha-1
        ds_litterfall = np.zeros(self.n)                    # dwarf shrub litterfall kg ha-1 yr-1
        h_litterfall = np.zeros(self.n)                     # herbs and grasses litterfall kg ha-1 yr-1
        s_litterfall = np.zeros(self.n)                     # sphagnum mosses litterfall kg ha-1 yr-1
        n_litter = np.zeros(self.n)                         # N uptake due to litterfall kg ha-1 yr-1
        p_litter = np.zeros(self.n)                         # P uptake due to litterfall kg ha-1 yr-1
        k_litter = np.zeros(self.n)                         # K uptake due to litterfall kg ha-1 yr-1
        n_gv = np.zeros(self.n)                             # N in ground vegetation kg ha-1
        p_gv = np.zeros(self.n)                             # P in ground vegetation kg ha-1
        k_gv = np.zeros(self.n)                             # K in ground vegetation kg ha-1


        """------ Ground vegetation models from Muukkonen & Mäkipää 2006 BER vol 11, Tables 6,7,8"""    
        
        #     ***************** Spruce mire ***************************************
        ix = self.ix_spruce_mire        
        gv_tot[ix] = np.square(35.52 +0.001*self.longitude*self.dem[ix] -1.1*self.drain_s**2 -2e-5*vol[ix]*Nstems[ix] \
                                +4e-5*Nstems[ix]*age[ix] +0.139*self.longitude*self.drain_s) -0.5 + 116.54 # Total, Eq.39, Table 9
        gv_bot[ix] =  np.square(-3.182 + 0.022*self.latitude*self.longitude +2e-4*self.dem[ix]*age[ix] \
                                -0.077*self.sfc[ix]*self.longitude -0.003*self.longitude*vol[ix] + 2e-4*np.square(vol[ix]))-0.5 + 98.10  #Bottom layer total, Eq. 35, Table 9
        gv_field[ix] =  np.square(23.24 -1.163*self.drain_s**2 +1.515*self.sfc[ix]*self.drain_s -2e-5*vol[ix]*Nstems[ix]\
                                +8e-5*ts*age[ix] +1e-5*Nstems[ix]*self.dem[ix])-0.5 +  162.58   #Field layer total, Eq. 37, Table 9
        
       # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix])        
        gv_field[ix] = np.maximum(gv_field[ix], gv_tot[ix] - gv_bot[ix])

        #annual litterfall rates
        #ATTN! vaihda nämä suoraan field layeriksi, poista tot ja bommomin kautta menevä yhteys
        ds_litterfall[ix] = self.fl_share['spruce_mire']['ds'] * (gv_tot[ix]-gv_bot[ix])*self.lit_share['ds'] * self.fl_to_total_turnover
        h_litterfall[ix]  = self.fl_share['spruce_mire']['h'] * (gv_tot[ix]-gv_bot[ix])*self.lit_share['h'] * self.fl_to_total_turnover
        s_litterfall[ix]  = gv_bot[ix]*self.lit_share['s']
        
        #ATTN! tarkista tämä, onko järkevä? Tee oma dictionary lehtimassalle
        gv_leafmass[ix]   = self.fl_share['spruce_mire']['ds'] * (gv_tot[ix]-gv_bot[ix]) * self.green_share['ds'] + \
                               self.fl_share['spruce_mire']['h']*(gv_tot[ix]-gv_bot[ix]) * self.green_share['h'] + \
                               gv_bot[ix] * self.green_share['s']
        
        
        n_gv[ix] = gv_field[ix] * self.fl_share['spruce_mire']['ds'] * self.nut_con['ds']['N']*1e-3 * self.fl_above_to_total \
                        + gv_field[ix] * self.fl_share['spruce_mire']['h'] * self.nut_con['h']['N']*1e-3 * self.fl_above_to_total \
                        + gv_bot[ix] * self.nut_con['s']['N']*1e-3
        p_gv[ix] = gv_field[ix] * self.fl_share['spruce_mire']['ds'] * self.nut_con['ds']['P']*1e-3 * self.fl_above_to_total \
                        + gv_field[ix] * self.fl_share['spruce_mire']['h'] * self.nut_con['h']['P']*1e-3 * self.fl_above_to_total \
                        + gv_bot[ix] *self.nut_con['s']['P']*1e-3
        k_gv[ix] = gv_field[ix] * self.fl_share['spruce_mire']['ds'] * self.nut_con['ds']['K']*1e-3 * self.fl_above_to_total \
                        + gv_field[ix] * self.fl_share['spruce_mire']['h'] *  self.nut_con['h']['K']*1e-3 * self.fl_above_to_total \
                        + gv_bot[ix] * self.nut_con['s']['K']*1e-3
        
        n_litter[ix] = ds_litterfall[ix] * self.nut_con['ds']['N']*1e-3 * (1.0 -self.retrans['ds']['N']) \
                        + h_litterfall[ix] * self.nut_con['h']['N']*1e-3 * (1.0 -self.retrans['h']['N']) \
                        + s_litterfall[ix] * self.nut_con['s']['N']*1e-3 * (1.0 -self.retrans['s']['N'])
        
        p_litter[ix] = ds_litterfall[ix] * self.nut_con['ds']['P']*1e-3 * (1.0 -self.retrans['ds']['P']) \
                        + h_litterfall[ix] * self.nut_con['h']['P']*1e-3 * (1.0 -self.retrans['h']['P']) \
                        + s_litterfall[ix] * self.nut_con['s']['P']*1e-3 * (1.0 -self.retrans['s']['P'])
        
        k_litter[ix] = ds_litterfall[ix] * self.nut_con['ds']['K']*1e-3 * (1.0 -self.retrans['ds']['K']) \
                        + h_litterfall[ix] * self.nut_con['h']['K']*1e-3 * (1.0 -self.retrans['h']['K']) \
                        + s_litterfall[ix] * self.nut_con['s']['K']*1e-3 * (1.0 -self.retrans['s']['K'])
        
       # ***************** Pine bogs ***************************************

        ix = self.ix_pine_bog            
        gv_tot[ix] =  np.square(50.098 + 0.005 * self.longitude*self.dem[ix] -1e-5 * vol[ix] * Nstems[ix] + 0.026 * self.sfc[ix] * age[ix] \
                      -1e-4 * self.dem[ix] * ts - 0.014 * vol[ix] * self.drain_s) - 0.5 + 167.40                #Total, Eq 45, Table 9           
        gv_bot[ix] =  np.square(31.809 + 0.008 * self.longitude * self.dem[ix] -3e-4 * Nstems[ix] * ba[ix] \
                                + 6e-5 * Nstems[ix] * age[ix] -0.188 * self.dem[ix]) -0.5 + 222.22                #Bottom layer total, Eq 41, Table 9
        gv_field[ix] =  np.square(48.12 - 1e-5 * ts**2 + 0.013 * self.sfc[ix] * age[ix] -0.04 * vol[ix] * self.drain_s \
                                + 0.026 * self.sfc[ix] * vol[ix]) - 0.5 +133.26                                        #Field layer total, Eq. 43, Table 9

        # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix])        
        gv_field[ix] = np.maximum(gv_field[ix], gv_tot[ix] - gv_bot[ix])

        # annual litterfall rates

        ds_litterfall[ix] = self.fl_share['pine_bog']['ds'] * (gv_tot[ix]-gv_bot[ix]) * self.lit_share['ds'] * self.fl_to_total_turnover
        h_litterfall[ix]  = self.fl_share['pine_bog']['h']  * (gv_tot[ix]-gv_bot[ix]) * self.lit_share['h'] * self.fl_to_total_turnover
        s_litterfall[ix]  = gv_bot[ix] * self.lit_share['s']
        gv_leafmass[ix]   = self.fl_share['pine_bog']['ds'] * (gv_tot[ix]-gv_bot[ix]) * self.green_share['ds'] + \
                            self.fl_share['pine_bog']['h'] * (gv_tot[ix]-gv_bot[ix]) * self.green_share['h'] + \
                            gv_bot[ix]*self.green_share['s']
      
        
        n_gv[ix] =  gv_field[ix] * self.fl_share['pine_bog']['ds'] * self.nut_con['ds']['N']*1e-3 * self.fl_above_to_total \
                        + gv_field[ix] * self.fl_share['pine_bog']['h'] * self.nut_con['h']['N']*1e-3 * self.fl_above_to_total \
                        + gv_bot[ix] * self.nut_con['s']['N']*1e-3
                        
        p_gv[ix] =  gv_field[ix] * self.fl_share['pine_bog']['ds'] * self.nut_con['ds']['P']*1e-3 * self.fl_above_to_total \
                        + gv_field[ix] * self.fl_share['pine_bog']['h'] * self.nut_con['h']['P']*1e-3 * self.fl_above_to_total \
                        + gv_bot[ix] * self.nut_con['s']['P']*1e-3
        
        k_gv[ix] =  gv_field[ix] * self.fl_share['pine_bog']['ds'] * self.nut_con['ds']['K']*1e-3 * self.fl_above_to_total \
                        + gv_field[ix] * self.fl_share['pine_bog']['h'] * self.nut_con['h']['K']*1e-3 * self.fl_above_to_total \
                        + gv_bot[ix] * self.nut_con['s']['K']*1e-3
        
        n_litter[ix] = ds_litterfall[ix] * self.nut_con['ds']['N']*1e-3 * (1.0 -self.retrans['ds']['N']) \
                        + h_litterfall[ix] * self.nut_con['h']['N']*1e-3 * (1.0 -self.retrans['h']['N']) \
                        + s_litterfall[ix] * self.nut_con['s']['N']*1e-3 * (1.0 -self.retrans['s']['N'])
        
        p_litter[ix] =  ds_litterfall[ix] * self.nut_con['ds']['P']*1e-3 * (1.0 -self.retrans['ds']['P']) \
                        + h_litterfall[ix] * self.nut_con['h']['P']*1e-3 * (1.0 -self.retrans['h']['P']) \
                        + s_litterfall[ix] * self.nut_con['s']['P']*1e-3 * (1.0 -self.retrans['s']['P'])
                        
        k_litter[ix] =  ds_litterfall[ix] * self.nut_con['ds']['K']*1e-3 * (1.0 -self.retrans['ds']['K']) \
                        + h_litterfall[ix] * self.nut_con['h']['K']*1e-3 * (1.0 -self.retrans['h']['K']) \
                        + s_litterfall[ix] * self.nut_con['s']['K']*1e-3 * (1.0 -self.retrans['s']['K'])

        
        #------------Change clear-cut areas: reduce to 1/3 of modelled ---------------------------------------------------
        to_cc = 0.33
        #ix_cc = np.where(np.logical_and(gisdata['age']<5.0, gisdata['smc']!=4))  #small stands excluding open peatlands
        ix_cc = np.where(age<5.0)
        n_gv[ix_cc] = n_gv[ix_cc] * to_cc 
        p_gv[ix_cc] = p_gv[ix_cc] * to_cc
        k_gv[ix_cc] = k_gv[ix_cc] * to_cc
        n_litter[ix_cc] = n_litter[ix_cc] * to_cc
        p_litter[ix_cc] = p_litter[ix_cc] * to_cc 
        k_litter[ix_cc] = k_litter[ix_cc] * to_cc 
        gv_tot[ix_cc] = gv_tot[ix_cc] * to_cc

        litterfall_gv = ds_litterfall + h_litterfall + s_litterfall
        
        return gv_tot, gv_field, gv_bot, n_gv, p_gv, k_gv, n_litter, p_litter, k_litter, \
                litterfall_gv, ds_litterfall, h_litterfall, s_litterfall, gv_leafmass

    def run(self, ba, stems, vol, sp, ts, age):
    
        #ATTN! convert barea, stems, yi, standAge from time series to list containing start and end state (adjustment to annual call)
        
        #---------------------------------------  
        gv_tot_ini = self.gv_tot 
        
        gv_tot, gv_field, gv_bot, n_gv, p_gv, k_gv, n_litter, p_litter, k_litter, \
                litterfall_gv, ds_litterfall, h_litterfall, s_litterfall,\
                    gv_leafmass = self.gv_biomass_and_nutrients(ts, vol, stems, ba,  age)
    
        
        # nutrient uptake due to net change in gv biomass, only positive values accepted, negative do not associate to nutrient uptake
        nup_net = np.where(n_gv - self.n_gv > 0.0, n_gv - self.n_gv, 0.0)
        pup_net = np.where(p_gv - self.p_gv > 0.0, p_gv - self.p_gv, 0.0)
        kup_net = np.where(k_gv - self.k_gv > 0.0, k_gv - self.k_gv, 0.0)
        
        self.nup = nup_net + n_litter        # total N uptake kg ha-1 simulation time (in yrs) -1
        self.pup = pup_net + p_litter        # total P uptake kg ha-1 simulation time (in yrs) -1
        self.kup = kup_net + k_litter        # total P uptake kg ha-1 simulation time (in yrs) -1
    
        self.n_gv = n_gv
        self.p_gv = p_gv
        self.k_gv = k_gv
        self.n_litter = n_litter 
        self.p_litter = p_litter
        self.k_litter = k_litter
        self.gv_tot = gv_tot 
        self.litterfall_gv = litterfall_gv
        self.gv_leafmass = gv_leafmass
        
        self.gv_tot = gv_tot                           # Ground vegetation mass kg ha-1
        self.gv_field = gv_field                         # Field layer vegetation mass
        self.gv_bot = gv_bot                           # Bottom layer vegetation mass
        self.ds_litterfall = ds_litterfall                    # dwarf shrub litterfall kg ha-1 yr-1
        self.h_litterfall = h_litterfall                     # herbs and grasses litterfall kg ha-1 yr-1
        self.s_litterfall = s_litterfall                     # sphagnum mosses litterfall kg ha-1 yr-1
        self.nonwoodylitter = litterfall_gv
        self.gv_change = gv_tot_ini - self.gv_tot
        
        
        return self.nup, self.pup, self.kup, self.n_litter, self.p_litter, self.k_litter, self.litterfall_gv, self.gv_leafmass

