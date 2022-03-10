# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 10:44:18 2022

@author: alauren
"""
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

class Allometry():
    def __init__(self):
        pass
    
    def get_motti(self, ifile, return_spe = False):
    #---read the Motti-simulation to be used as a basis for the Susi-simulation

        cnames=['yr', 'age', 'N', 'BA', 'Hg', 'Dg', 'hdom', 'vol', 'logs', 'pulp', 'loss', 'yield','mortality', 
                'stem', 'stemloss', 'branch_living', 'branch_dead', 'leaves', 'stump', 'roots_coarse', 'roots_fine']
        df = pd.read_excel(ifile, sheet_name=0, usecols=range(22), skiprows=1, header=None )
        df = df.drop([0], axis=1)
        df.columns=cnames
        cname=['idSpe']
        df2 = pd.read_excel(ifile, sheet_name=1, usecols=[4], skiprows=1, header=None )
        df2.columns = cname
        
        #---- find thinnings and add a small time to lines with the age to enable interpolation---------
        df = df.loc[df['age']!=0]    
        
        steps = np.array(np.diff(df['age']), dtype=float)
        idx = np.ravel(np.argwhere(steps<1.))+1
        df['age'][idx]=df['age'][idx]+5./365.
        
        if return_spe: 
            return df, df2['idSpe'][0]
        else:
            return df
    
    def motti_development(self, ifile):
        """
        Input:
            spara contains information about tree species in the stand
            Motti-input file name including the folder path
        Out:
            interpolation functions: 
                age in annual [yrs]  
                age [yrs] to variables: 
                    ageToHdom, [m] 
                    ageToBa, [m2/ha]
                    ageToVol, [m3/ha]
                    ageToYield, [m3/ha]
                    ageToBm [kg dry mass / ha]
                biomass [kg dry mass / ha] to variables:
                    bmToLeafMass, [kg/ha] 
                    bmToLAI, [m2/m2]
                    bmToHdom, [m]
                    bmToYi, [m3/ha]
                    bmToBa, [m2/ha]
                    bmToLitter, [kg/ha/day]
                    bmToStems [number/ha]
                volume or yield to variables:
                    yiToVol [m3]
                    yiToBm, [kg dry mass/ha]
                    volToLogs, [m3/ha]
                    volToPulp, [m3/ha]
                    sp    species
            Biomass models in Motti (Repola 2008, 2009) have been develped for mineral soils. In peatlands the 
            leaf mass is 35.5% lower. This bias is corrected in construction of the interpolation function
        Modifications needed:
            create new litter scheme: see Dec 21 esom model development
            biomass to nutrient interpolation functions: N, P, K
            nutrient to biomass interpolation functions
            litter: woody, nonwoody, locate to interpolation function 
            locate interpolation functions to dictionaty 
        """
        cnames=['yr', 'age', 'N', 'BA', 'Hg', 'Dg', 'hdom', 'vol', 'logs', 'pulp', 'loss', 'yield','mortality', 
                'stem', 'stemloss', 'branch_living', 'branch_dead', 'leaves', 'stump', 'roots_coarse', 'roots_fine']
        species_codes ={1:'Pine', 2:'Spruce', 3:'Birch'}
        df, sp = self.get_motti(ifile, return_spe=True)
        sp = sp if sp < 4 else 3
        spe = species_codes[sp]
    
    #----modify the data frame, include age = 0-------------------------
        row = np.zeros((np.shape(df)[1]), dtype=float)
        dfrow = pd.DataFrame([row])
        dfrow.columns = cnames
        df = pd.concat([dfrow, df], axis=0)
        df.loc[0, 'N'] = df.loc[1, 'N']                                             # modify stem number at 0 yrs 
        
        df[['stem', 'branch_living', 'branch_dead', 'leaves', 
            'stump', 'roots_coarse', 'roots_fine']] = df[['stem','branch_living', 
                                                          'branch_dead', 'leaves', 
                                                          'stump', 'roots_coarse', 
                                                          'roots_fine']] *1000.     # unit conversion for all biomass components tn/ha -> kg/ha
        a_arr = np.arange(0, max(df['age'].values), 1.)                             # stand age array from 0 to max in Motti simulation, time step year
    
        
        #---Nutrient concentrations in tree biomass components: Palviainen & Finer 2012 Eur J For Res 131: 945-964
        #*********** Parameters *********************************************
        #---Concentrations in mg/g
        nuts = {'Pine':{
                    'Foliage':{'N': 12.5, 'P': 1.3, 'K': 4.25},
                    'Stem':{'N':1.17, 'P':0.08, 'K':0.45}},
                'Spruce':{
                    'Foliage':{'N': 12.5, 'P': 1.3, 'K': 4.25},
                    'Stem':{'N':1.12, 'P':0.09, 'K':0.64}},
                'Birch':{
                    'Foliage':{'N': 12.5, 'P': 1.3, 'K': 4.25},
                    'Stem':{'N':1.51, 'P':0.15, 'K':0.58}},
                }
                
        retrans = {'N': 0.69, 'P': 0.73, 'K':0.8}                                   # Nieminen Helmisaari 1996 Tree Phys
        rho = {'Pine': 400., 'Spruce': 380., 'Birch':480.}                          # wood density kg/m3
        sla= {'Pine': 6.8, 'Spruce': 7.25, 'Birch':14.0}                            # one-sided specific leaf area Härkönen et al. 2015 BER 20, 181-195      
    
        longevityLeaves = {'Pine':2.5, 'Spruce':4., 'Birch':1.}                     # yrs, life span of leaves and fine roots    
        longevityFineRoots ={'Pine':1., 'Spruce':1., 'Birch':1.}                    # Yuan & Chen 2010, turnover 1.07 times per year    
        longevityBranch ={'Pine':22., 'Spruce':22., 'Birch':22.}                    # Pine Mäkinen 1999
        longevityCoarseRoots ={'Pine':22., 'Spruce':22., 'Birch':22.}               # assumption as branches
    
        #********** Interpolation data ****************************************
        df['leaves'] = df['leaves'] / 1.355                                         # adjusting to peatland sites (Data: hannu Hökkä 2022)
        df['leafarea'] = df['leaves'].values/10000. * sla[spe]                      # leaf area index m2 m-2
        df['stem_mass'] = df['yield'] * rho[spe]                                    # stem biomass
        df['stem_and_stump'] = df[['stem_mass', 'stump']].sum(axis=1)
        df['bm'] = df[['stem_mass', 'branch_living', 'branch_dead', 'leaves',
            'stump', 'roots_coarse', 'roots_fine']].sum(axis=1)
        df['bm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead',  
            'stump', 'roots_coarse', 'roots_fine']].sum(axis=1)
        df['N_leaves'] =  df[ 'leaves'] *nuts[spe]['Foliage']['N']/1000.   
        df['Nbm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead', 'stump', 'roots_coarse']].sum(axis=1) * nuts[spe]['Stem']['N'] /1000.\
                                     + df[['roots_fine']].sum(axis=1) *nuts[spe]['Foliage']['N']/1000.       
        df['P_leaves'] =  df[ 'leaves'] *nuts[spe]['Foliage']['P']/1000.   
        df['Pbm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead', 'stump', 'roots_coarse']].sum(axis=1) * nuts[spe]['Stem']['P'] /1000.\
                                     + df[['roots_fine']].sum(axis=1) *nuts[spe]['Foliage']['P']/1000.       
        df['K_leaves'] =  df[ 'leaves'] *nuts[spe]['Foliage']['K']/1000.   
        df['Kbm_noleaves'] =  df[['stem_mass','branch_living', 'branch_dead', 'stump', 'roots_coarse']].sum(axis=1) * nuts[spe]['Stem']['K'] /1000.\
                                     + df[['roots_fine']].sum(axis=1) *nuts[spe]['Foliage']['K']/1000.       

        df['woody_logging_residues'] = df[['branch_living', 'branch_dead', 'roots_coarse', 'stump']].sum(axis=1)
        df['N_woody_logging_residues'] = df['woody_logging_residues']*nuts[spe]['Stem']['N']/1000.
        df['P_woody_logging_residues'] = df['woody_logging_residues']*nuts[spe]['Stem']['P']/1000.
        df['K_woody_logging_residues'] = df['woody_logging_residues']*nuts[spe]['Stem']['K']/1000.

        df['N_fine_roots'] = df['roots_fine']*nuts[spe]['Foliage']['N'] / 1000.
        df['P_fine_roots'] = df['roots_fine']*nuts[spe]['Foliage']['P'] / 1000.
        df['K_fine_roots'] = df['roots_fine']*nuts[spe]['Foliage']['K'] / 1000.
        
        
        
        #********** Interpolation functions ******************************************
        ageToHdom = interp1d(df['age'].values, df['hdom'].values, fill_value=(df['hdom'].values[0], df['hdom'].values[-1]), bounds_error=False)
        ageToLAI= interp1d(df['age'].values, df['leafarea'].values,fill_value= (df['leafarea'].values[0], df['leafarea'].values[-1]), bounds_error=False)        
        ageToYield = interp1d(df['age'].values, df['yield'].values,fill_value=(df['yield'].values[0], df['yield'].values[-1]), bounds_error=False)
        ageToVol = interp1d(df['age'].values,df['vol'].values,fill_value=(df['vol'].values[0], df['vol'].values[-1]), bounds_error=False)
        ageToBa = interp1d(df['age'].values,df['BA'].values,fill_value=(df['BA'].values[0], df['BA'].values[-1]), bounds_error=False)
        ageToBm = interp1d(df['age'].values,df['bm'].values, fill_value=(df['bm'].values[0], df['bm'].values[-1]), bounds_error=False)
        ageToBmNoLeaves = interp1d(df['age'].values,df['bm_noleaves'].values, fill_value=(df['bm_noleaves'].values[0], df['bm_noleaves'].values[-1]), bounds_error=False)
        ageToStems = interp1d(df['age'].values,df['N'].values, fill_value=(df['N'].values[0], df['N'].values[-1]), bounds_error=False)
        ageToLeaves = interp1d(df['age'].values, df['leaves'].values, fill_value =(df['leaves'].values[0], df['leaves'].values[-1]), bounds_error=False)
        ageToFineRoots = interp1d(df['age'].values, df['roots_fine'].values, fill_value =(df['roots_fine'].values[0], df['roots_fine'].values[-1]), bounds_error=False)
        ageToBranchLiving = interp1d(df['age'].values, df['branch_living'].values, fill_value =(df['branch_living'].values[0], df['branch_living'].values[-1]), bounds_error=False)
        ageToBranchDead = interp1d(df['age'].values, df['branch_dead'].values, fill_value =(df['branch_dead'].values[0], df['branch_dead'].values[-1]), bounds_error=False)
        ageToCoarseRoots = interp1d(df['age'].values, df['roots_coarse'].values, fill_value =(df['roots_coarse'].values[0], df['roots_coarse'].values[-1]), bounds_error=False)
        ageToStemStump = interp1d(df['age'].values, df['stem_and_stump'].values, fill_value =(df['stem_and_stump'].values[0], df['stem_and_stump'].values[-1]), bounds_error=False)
        ageToNNoLeaves = interp1d(df['age'].values, df['Nbm_noleaves'].values, fill_value =(df['Nbm_noleaves'].values[0], df['Nbm_noleaves'].values[-1]), bounds_error=False)
        ageToPNoLeaves = interp1d(df['age'].values, df['Pbm_noleaves'].values, fill_value =(df['Pbm_noleaves'].values[0], df['Pbm_noleaves'].values[-1]), bounds_error=False)
        ageToKNoLeaves = interp1d(df['age'].values, df['Kbm_noleaves'].values, fill_value =(df['Kbm_noleaves'].values[0], df['Kbm_noleaves'].values[-1]), bounds_error=False)

            
        volToLogs = interp1d(df['vol'].values,df['logs'].values,fill_value=(df['logs'].values[0], df['logs'].values[-1]), bounds_error=False)  
        volToPulp = interp1d(df['vol'].values,df['pulp'].values,fill_value=(df['pulp'].values[0], df['pulp'].values[-1]), bounds_error=False)
        
        yiToVol = interp1d(df['yield'].values,df['vol'].values, fill_value=(df['vol'].values[0], df['vol'].values[-1]), bounds_error=False)
        yiToBm = interp1d(df['yield'].values, df['bm'].values, fill_value=(df['bm'].values[0], df['bm'].values[-1]), bounds_error=False)
    
        bmToYi= interp1d(df['bm_noleaves'].values,df['yield'].values, fill_value=(df['yield'].values[0], df['yield'].values[-1]), bounds_error=False)
        bmToBa= interp1d(df['bm_noleaves'].values,df['BA'].values, fill_value=(df['BA'].values[0], df['BA'].values[-1]), bounds_error=False)    
        bmToLeafMass = interp1d(df['bm_noleaves'].values, df['leaves'].values, fill_value=(df['leaves'].values[0], df['leaves'].values[-1]), bounds_error=False)
        bmToLAI = interp1d(df['bm_noleaves'].values, df['leaves'].values* sla[spe]/10000., fill_value=(df['leaves'].values[0]* sla[spe]/10000., df['leaves'].values[-1]* sla[spe]/10000.), bounds_error=False)
        bmToHdom = interp1d(df['bm_noleaves'].values,df['hdom'].values, fill_value=(df['hdom'].values[0], df['hdom'].values[-1]), bounds_error=False)
        bmToStems = interp1d(df['bm_noleaves'].values, df['N'].values, fill_value=(df['N'].values[0], df['N'].values[-1]), bounds_error=False)

 
        bmToFineRoots = interp1d(df['bm_noleaves'].values, df['roots_fine'].values, fill_value=(df['roots_fine'].values[0], df['roots_fine'].values[-1]), bounds_error=False)
        bmToNFineRoots = interp1d(df['bm_noleaves'].values, df['N_fine_roots'].values, fill_value=(df['N_fine_roots'].values[0], df['N_fine_roots'].values[-1]), bounds_error=False)
        bmToPFineRoots = interp1d(df['bm_noleaves'].values, df['P_fine_roots'].values, fill_value=(df['P_fine_roots'].values[0], df['P_fine_roots'].values[-1]), bounds_error=False)
        bmToKFineRoots = interp1d(df['bm_noleaves'].values, df['K_fine_roots'].values, fill_value=(df['K_fine_roots'].values[0], df['K_fine_roots'].values[-1]), bounds_error=False)


        bmToWoodyLoggingResidues = interp1d(df['bm_noleaves'].values, df['woody_logging_residues'].values, fill_value =(df['woody_logging_residues'].values[0], df['woody_logging_residues'].values[-1]), bounds_error=False)        
        bmToNWoodyLoggingResidues =  interp1d(df['bm_noleaves'].values, df['N_woody_logging_residues'].values, fill_value =(df['N_woody_logging_residues'].values[0], df['N_woody_logging_residues'].values[-1]), bounds_error=False)        
        bmToPWoodyLoggingResidues =  interp1d(df['bm_noleaves'].values, df['P_woody_logging_residues'].values, fill_value =(df['P_woody_logging_residues'].values[0], df['P_woody_logging_residues'].values[-1]), bounds_error=False)        
        bmToKWoodyLoggingResidues =  interp1d(df['bm_noleaves'].values, df['K_woody_logging_residues'].values, fill_value =(df['K_woody_logging_residues'].values[0], df['K_woody_logging_residues'].values[-1]), bounds_error=False)        
        
    
        #**********************************************************************
        # We need demand functions for mass, N,P,K: 
        #              net change + fineroot_litter + woody_litter + add foliage net demand litter from other function
        # Nutrient contents in Litterfall: 
        #              fineroot litter + woody_litter + mortality_fine_root + mortality_woody + foliage litter and net demand from another function
        # Arrange: 
        #              demand functions; mass, N, P, K 
        #              litter functions: mass, N, P, K
        #***********************************************************************
        
        """ Litter functions """
        #---- Litter arrays------------------
        fineroot_litter = ageToFineRoots(a_arr)/longevityFineRoots[spe]*np.gradient(a_arr)                  # unit kg / ha / yr 
        mortality_fineroot = -np.gradient(ageToStems(a_arr))/ageToStems(a_arr)*(ageToFineRoots(a_arr))      # unitkg / ha / yr 
        mortality_leaf = -np.gradient(ageToStems(a_arr))/ageToStems(a_arr)*(ageToLeaves(a_arr))
    
        woody_litter = ageToBranchLiving(a_arr)/longevityBranch[spe]*np.gradient(a_arr) \
                    + ageToBranchDead(a_arr)/longevityBranch[spe]*np.gradient(a_arr)  \
                    + ageToCoarseRoots(a_arr)/longevityCoarseRoots[spe]*np.gradient(a_arr)                  # litterfall kg/ha in timestep
    
        mortality_woody = -np.gradient(ageToStems(a_arr))/ageToStems(a_arr) *(ageToBranchDead(a_arr) + ageToBranchLiving(a_arr) 
                                                                         + ageToStemStump(a_arr) + ageToCoarseRoots(a_arr))    
        
        dbm = np.gradient(ageToBmNoLeaves(a_arr)) + fineroot_litter + woody_litter                          # biomass change without leaves kg/ha/yr   
        
        #---- Interpolation functions -----------------
        bmToDbm = interp1d(ageToBmNoLeaves(a_arr), dbm, fill_value=(dbm[0], dbm[-1]), bounds_error=False)   # from biomass to biomass change 
        bmToFinerootLitter = interp1d(ageToBmNoLeaves(a_arr), fineroot_litter, fill_value=(fineroot_litter[0], fineroot_litter[-1]), bounds_error=False )    
        bmToWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), woody_litter, fill_value=(woody_litter[0], woody_litter[-1]), bounds_error=False )
        bmToMortalityFineRoot = interp1d(ageToBmNoLeaves(a_arr), mortality_fineroot, fill_value=(mortality_fineroot[0], mortality_fineroot[-1]), bounds_error=False )
        bmToMortalityLeaves = interp1d(ageToBmNoLeaves(a_arr), mortality_leaf, fill_value=(mortality_leaf[0], mortality_leaf[-1]), bounds_error=False )
        bmToMortalityWoody =  interp1d(ageToBmNoLeaves(a_arr), mortality_woody, fill_value=(mortality_woody[0], mortality_woody[-1]), bounds_error=False )
    
    
        """ Demand functions """
        #-------- arrays--------------------------
        N_fineroot_litter = (1.0-retrans['N'])*nuts[spe]['Foliage']['N'] / 1000. * fineroot_litter + mortality_fineroot*nuts[spe]['Foliage']['N'] / 1000.
        P_fineroot_litter = (1.0-retrans['P'])*nuts[spe]['Foliage']['P'] / 1000. * fineroot_litter + mortality_fineroot*nuts[spe]['Foliage']['P'] / 1000.
        K_fineroot_litter = (1.0-retrans['K'])*nuts[spe]['Foliage']['K'] / 1000. * fineroot_litter + mortality_fineroot*nuts[spe]['Foliage']['K'] / 1000.
    
        N_woody_litter = (1.0-retrans['N'])*nuts[spe]['Stem']['N'] / 1000. * woody_litter + mortality_woody*nuts[spe]['Stem']['N'] / 1000.
        P_woody_litter = (1.0-retrans['P']) *nuts[spe]['Stem']['N'] / 1000.* woody_litter + mortality_woody*nuts[spe]['Stem']['P'] / 1000.
        K_woody_litter = (1.0-retrans['K']) *nuts[spe]['Stem']['N'] / 1000.* woody_litter + mortality_woody*nuts[spe]['Stem']['K'] / 1000.
    
        N_mortality_leaves = mortality_leaf *  nuts[spe]['Foliage']['N'] / 1000.
        P_mortality_leaves = mortality_leaf *  nuts[spe]['Foliage']['P'] / 1000.
        K_mortality_leaves = mortality_leaf *  nuts[spe]['Foliage']['K'] / 1000.
    
        N_mortality_fineroot = mortality_fineroot *  nuts[spe]['Foliage']['N'] / 1000.
        P_mortality_fineroot = mortality_fineroot *  nuts[spe]['Foliage']['P'] / 1000.
        K_mortality_fineroot = mortality_fineroot *  nuts[spe]['Foliage']['K'] / 1000.
    
        N_demand = np.gradient(ageToNNoLeaves(a_arr)) + (1.0-retrans['N'])*nuts[spe]['Foliage']['N'] / 1000. * fineroot_litter\
                                                     +  (1.0-retrans['N'])*nuts[spe]['Stem']['N'] / 1000. * woody_litter
    
        P_demand = np.gradient(ageToPNoLeaves(a_arr)) + (1.0-retrans['P'])*nuts[spe]['Foliage']['P'] / 1000. * fineroot_litter\
                                                     +  (1.0-retrans['P'])*nuts[spe]['Stem']['P'] / 1000. * woody_litter
    
        K_demand = np.gradient(ageToKNoLeaves(a_arr)) + (1.0-retrans['K'])*nuts[spe]['Foliage']['K'] / 1000. * fineroot_litter\
                                                     +  (1.0-retrans['K'])*nuts[spe]['Stem']['K'] / 1000. * woody_litter
    
        #---- Interpolation functions -----------------
        bmToNdemand = interp1d(ageToBmNoLeaves(a_arr), N_demand, fill_value=(N_demand[0], N_demand[-1]), bounds_error=False)   # from biomass to N demand No leaves here  
        bmToPdemand = interp1d(ageToBmNoLeaves(a_arr), P_demand, fill_value=(P_demand[0], P_demand[-1]), bounds_error=False)   # from biomass to P demand 
        bmToKdemand = interp1d(ageToBmNoLeaves(a_arr), K_demand, fill_value=(K_demand[0], K_demand[-1]), bounds_error=False)   # from biomass to K demand 
        
        bmToNFineRootLitter = interp1d(ageToBmNoLeaves(a_arr), N_fineroot_litter, fill_value=(N_fineroot_litter[0], N_fineroot_litter[-1]), bounds_error=False)
        bmToPFineRootLitter = interp1d(ageToBmNoLeaves(a_arr), P_fineroot_litter, fill_value=(P_fineroot_litter[0], P_fineroot_litter[-1]), bounds_error=False)
        bmToKFineRootLitter = interp1d(ageToBmNoLeaves(a_arr), K_fineroot_litter, fill_value=(K_fineroot_litter[0], K_fineroot_litter[-1]), bounds_error=False)
    
        bmToNWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), N_woody_litter, fill_value=(N_woody_litter[0], N_woody_litter[-1]), bounds_error=False)
        bmToPWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), P_woody_litter, fill_value=(P_woody_litter[0], P_woody_litter[-1]), bounds_error=False)
        bmToKWoodyLitter = interp1d(ageToBmNoLeaves(a_arr), K_woody_litter, fill_value=(K_woody_litter[0], K_woody_litter[-1]), bounds_error=False)
    
        bmToNMortalityLeaves =  interp1d(ageToBmNoLeaves(a_arr), N_mortality_leaves, fill_value=(N_mortality_leaves[0], N_mortality_leaves[-1]), bounds_error=False)
        bmToPMortalityLeaves =  interp1d(ageToBmNoLeaves(a_arr), P_mortality_leaves, fill_value=(P_mortality_leaves[0], P_mortality_leaves[-1]), bounds_error=False)
        bmToKMortalityLeaves =  interp1d(ageToBmNoLeaves(a_arr), K_mortality_leaves, fill_value=(K_mortality_leaves[0], K_mortality_leaves[-1]), bounds_error=False)
    
        bmToNMortalityFineRoot =  interp1d(ageToBmNoLeaves(a_arr), N_mortality_fineroot, fill_value=(N_mortality_fineroot[0], N_mortality_fineroot[-1]), bounds_error=False)
        bmToPMortalityFineRoot =  interp1d(ageToBmNoLeaves(a_arr), P_mortality_fineroot, fill_value=(P_mortality_fineroot[0], P_mortality_fineroot[-1]), bounds_error=False)
        bmToKMortalityFineRoot =  interp1d(ageToBmNoLeaves(a_arr), K_mortality_fineroot, fill_value=(K_mortality_fineroot[0], K_mortality_fineroot[-1]), bounds_error=False)
    
   
    
        allometry_f={}
        allometry_f['ageToHdom'] = ageToHdom 
        allometry_f['ageToBa'] = ageToBa
        allometry_f['ageToVol'] = ageToVol
        allometry_f['ageToYield'] = ageToYield
        allometry_f['ageToBm'] = ageToBm
        allometry_f['ageToBmNoLeaves'] = ageToBmNoLeaves
        allometry_f['ageToLeaves']= ageToLeaves
        
        allometry_f['bmToLeafMass'] = bmToLeafMass
        allometry_f['bmToLAI'] = bmToLAI
        allometry_f['bmToHdom'] = bmToHdom
        allometry_f['bmToYi'] = bmToYi
        allometry_f['bmToBa'] = bmToBa
        allometry_f['bmToDbm'] = bmToDbm
    
        allometry_f['bmToStems'] = bmToStems
        allometry_f['yiToVol'] = yiToVol
        allometry_f['yiToBm'] = yiToBm
        allometry_f['volToLogs'] = volToLogs
        allometry_f['volToPulp'] = volToPulp
        
        allometry_f['bmToDbm'] = bmToDbm
        allometry_f['bmToFinerootLitter'] = bmToFinerootLitter
        allometry_f['bmToWoodyLitter'] = bmToWoodyLitter
        allometry_f['bmToMortalityFineRoot'] =  bmToMortalityFineRoot
        allometry_f['bmToMortalityWoody'] = bmToMortalityWoody
        allometry_f['bmToMortalityLeaves'] = bmToMortalityLeaves
    
        allometry_f['bmToNdemand'] = bmToNdemand
        allometry_f['bmToPdemand'] = bmToPdemand
        allometry_f['bmToKdemand'] = bmToKdemand
        allometry_f['bmToNFineRootLitter'] = bmToNFineRootLitter
        allometry_f['bmToPFineRootLitter'] = bmToPFineRootLitter
        allometry_f['bmToKFineRootLitter'] = bmToKFineRootLitter
    
        allometry_f['bmToNWoodyLitter'] = bmToNWoodyLitter
        allometry_f['bmToPWoodyLitter'] = bmToPWoodyLitter
        allometry_f['bmToKWoodyLitter'] = bmToKWoodyLitter
    
        allometry_f['bmToNMortalityLeaves'] = bmToNMortalityLeaves
        allometry_f['bmToPMortalityLeaves'] = bmToPMortalityLeaves
        allometry_f['bmToKMortalityLeaves'] = bmToKMortalityLeaves
    
        allometry_f['bmToNMortalityFineRoot'] = bmToNMortalityFineRoot
        allometry_f['bmToPMortalityFineRoot'] = bmToPMortalityFineRoot
        allometry_f['bmToKMortalityFineRoot'] = bmToKMortalityFineRoot
        
        allometry_f['bmToWoodyLoggingResidues'] = bmToWoodyLoggingResidues
        allometry_f['bmToNWoodyLoggingResidues'] = bmToNWoodyLoggingResidues
        allometry_f['bmToPWoodyLoggingResidues'] = bmToPWoodyLoggingResidues
        allometry_f['bmToKWoodyLoggingResidues'] = bmToKWoodyLoggingResidues
        
        allometry_f['bmToFineRoots'] = bmToFineRoots
        allometry_f['bmToNFineRoots'] = bmToNFineRoots
        allometry_f['bmToPFineRoots'] = bmToPFineRoots
        allometry_f['bmToKFineRoots'] = bmToKFineRoots
        
        
        self.allometry_f = allometry_f
        self.sp=sp
        self.df = df
        #return allometry_f, sp, df 