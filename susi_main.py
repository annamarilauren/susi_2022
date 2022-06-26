# -*- coding: utf-8 -*-
"""
Created on Mon May 21 18:38:10 2018

@author: lauren
"""
import numpy as np
import pandas as pd
import datetime

from canopygrid import CanopyGrid
from mosslayer import MossLayer
from strip import StripHydrology, drain_depth_development
from temperature import PeatTemperature
#from docclass import DocModel
from gvegetation import Gvegetation
from esom import Esom
from stand import Stand
from methane import Methane
from fertilization import Fertilization
from outputs import Outputs

import susi_io
from susi_utils import  rew_drylimit
from susi_utils import  get_temp_sum

class Susi():
    def __init(self):
        pass

    def run_susi(self, forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation=None, mottifile=None, peat=None, 
                 photosite=None, folderName=None, hdomSim=None, volSim=None, ageSim=None, 
                 sarkaSim=None, sfc=None, susiPath=None, simLAI=None, kaista=None, sitename=None): 
        
        print ('******** Susi-peatland simulator v.10 (2022) c Annamari Laurén *********************')
        print ('           ')    
        print ('Initializing stand and site:') 
         
        dtc = cpara['dt']                                                          # canopy model timestep
    
        start_date = datetime.datetime(start_yr,1,1) 
        end_date=datetime.datetime(end_yr,12,31)
        length = (end_date - start_date).days +1                                   # simulation time in days
        yrs = end_yr - start_yr +1                                                 # simulation time in years
        ts = get_temp_sum(forc)                                                    # temperature sum degree days
        nscens = len(spara['ditch depth east'])                                    # number of scenarios          
        n = spara['n']                                                             # number of columns along the strip        
        
        outname = outpara['outfolder'] +  outpara['netcdf']                        # name and path for the netcdf4 file for output
        
        out = Outputs(nscens, n, length, yrs, spara['nLyrs'], outname)
        out.initialize_scens()
        out.initialize_paras()
        
        lat=forc['lat'][0]; lon=forc['lon'][0]                                     # location of weather file, determines the simulation location
        print ('      - Weather input:', wpara['description'], ', start:', start_yr, ', end:', end_yr) 
        print ('      - Latitude:', lat, ', Longitude:', lon )
    
    
        stand = Stand(nscens, yrs, spara['canopylayers'], spara['n'], sfc, ageSim, mottifile, photopara)   
        stand.update()
        #spara = stand.update_spara(spara)    
        
        out.initialize_stand()
        out.initialize_canopy_layer('dominant')
        out.initialize_canopy_layer('subdominant')
        out.initialize_canopy_layer('under')
        
        out.write_paras(spara['sfc'], stand.dominant.tree_species, stand.subdominant.tree_species, stand.under.tree_species)
    
        susi_io.print_site_description(spara)                                      # Describe site parameters for user
    
        groundvegetation = Gvegetation(spara['n'], lat, lon, sfc, stand.dominant.species)
        groundvegetation.run(stand.basalarea, stand.stems, stand.volume,\
                             stand.dominant.species, ts, ageSim['dominant'])
        out.initialize_gv()
    
        
        esmass = Esom(spara, sfc, 366*yrs, substance='Mass')                       # initializing organic matter decomposition instace for mass
        esN = Esom(spara, sfc, 366*yrs, substance='N')                             # initializing organic matter decomposition instace for N
        esP = Esom(spara, sfc, 366*yrs, substance='P')                             # initializing organic matter decomposition instace for P
        esK = Esom(spara, sfc, 366*yrs, substance='K')                             # initializing organic matter decomposition instace for K
        ferti = Fertilization(spara)
    
        out.initialize_esom('Mass')
        out.initialize_esom('N')
        out.initialize_esom('P')
        out.initialize_esom('K')
        out.initialize_fertilization()
        out.initialize_nutrient_balance('N')
        out.initialize_nutrient_balance('P')
        out.initialize_nutrient_balance('K')
        out.initialize_carbon_balance()
            
        #********* Above ground hydrology initialization ***************
        cmask = np.ones(spara['n'])                                                # compute canopy and moss for each soil column (0, and n-1 are ditches)
        cstate = cpara['state'].copy()
        for key in cstate.keys():
            cstate[key] *= cmask
        cpy = CanopyGrid(cpara, cstate, outputs=False)                             # initialize above ground vegetation hydrology model
        #cpy.update_amax(cpara['physpara'], stand.nut_stat)
        out.initialize_cpy()
        
        for key in org_para.keys():                                                 
            org_para[key] *= cmask
        moss = MossLayer(org_para, outputs=True)  
        print ('Canopy and moss layer hydrology initialized')
    
        #******** Soil and strip parameterization *************************
        stp = StripHydrology(spara)                                                # initialize soil hydrology model
        out.initialize_strip(stp)
        
        pt = PeatTemperature(spara, forc['T'].mean())                              # initialize peat temperature model       
        out.initialize_temperature()
            
        ch4s = Methane(n, yrs)
        out.initialize_methane()
    
        out.initialize_export()
        print ('Soil hydrology, temperature and DOC models initialized')
    
        ets = np.zeros((length, n))                                                # Evapotranspiration, mm/day
        
        #********initialize result arrays***************************
        scen = spara['scenario name']                                              # scenario name for outputs    
        rounds = len(spara['ditch depth east'])                                    # number of ditch depth scenarios (used in comparison of management)
        
        stpout = stp.create_outarrays(rounds, length, n)    
        peat_temperatures = pt.create_outarrays(rounds, length, spara['nLyrs'])     # daily peat temperature profiles
        intercs, evaps, ETs, transpis, efloors, swes = cpy.create_outarrays(rounds, length, n)
        
        #***********Scenario loop ********************************************************
        
        for r, dr in enumerate(zip(spara['ditch depth west'], spara['ditch depth 20y west'], \
                                   spara['ditch depth east'], spara['ditch depth 20y east'])):   
    
            dwt=spara['initial h']*np.ones(spara['n'])                             # set the initial WT for the scenario
            hdr_west, hdr20y_west,hdr_east, hdr20y_east = dr                       # drain depth [m] in the beginning and after 20 yrs
            h0ts_west = drain_depth_development(length, hdr_west, hdr20y_west)     # compute daily values for drain bottom boundary condition
            h0ts_east = drain_depth_development(length, hdr_east, hdr20y_east)     # compute daily values for drain bottom boundary condition
    
            
            # ---- Initialize integrative output arrays (outputs in nodewise sums) -------------------------------
    
            print ('***********************************')        
            print ('Computing canopy and soil hydrology ', length, ' days', 'scenario:', scen[r])
            
            stand.reset_domain(ageSim)  
            out.write_scen(r, hdr_west, hdr_east)
            
            out.write_stand(r, 0, stand)
            out.write_canopy_layer(r, 0, 'dominant', stand.dominant)
            out.write_canopy_layer(r, 0, 'subdominant', stand.subdominant)
            out.write_canopy_layer(r, 0, 'under', stand.under)
            
            groundvegetation.reset_domain()
            groundvegetation.run(stand.basalarea, stand.stems, stand.volume,\
                                 stand.dominant.species, ts, ageSim['dominant'])
            out.write_groundvegetation(r, 0, groundvegetation)
    
            
            esmass.reset_storages()
            esN.reset_storages()
            esP.reset_storages()
            esK.reset_storages()
    
            out.write_esom(r, 0, 'Mass', esmass, inivals=True)
            out.write_esom(r, 0, 'N', esN, inivals=True)
            out.write_esom(r, 0, 'P', esP, inivals=True)
            out.write_esom(r, 0, 'K', esK, inivals=True)
            
            stp.reset_domain()   
            pt.reset_domain()
            
            d = 0                                                                  # day index
            start = 0                                                              # day counter in annual loop
            year = 0                                                               # year counter in annual loop
            # *************** Annual loop *****************************************************************
            for yr in range(start_yr, end_yr + 1):                                 # year loop 
                days = (datetime.datetime(yr,12, 31) - datetime.datetime(yr,1, 1)).days + 1
                
                # CHECK THIS AND TEST
                #cpy.update_amax(cpara['physpara'], stand.nut_stat)
    
                #**********  Daily loop ************************************************************
                for dd in range(days):                                             # day loop   
                    #-------Canopy hydrology--------------------------            
                    reww = rew_drylimit(dwt)                                       # for each column: moisture limitation from ground water level (Feddes-function)            
                    doy = forc.iloc[d, 14]                                         # day of the year
                    ta =  forc.iloc[d, 4]                                          # air temperature deg C
                    vpd = forc.iloc[d, 13]                                         # vapor pressure deficit
                    rg = forc.iloc[d, 8]                                           # solar radiation
                    par = forc.iloc[d, 10]                                         # photosynthetically active radiation
                    prec=forc.iloc[d, 7]/86400.                                    # precipitation
        
                    potinf, trfall, interc, evap, ET, transpi, efloor, MBE, SWE = cpy.run_timestep(doy, dtc, ta, prec, rg, par, vpd, 
                                                                    hc=stand.hdom, LAIconif=stand.leafarea, Rew=reww, beta=moss.Ree)                       # canopy hydrology computation
                    
                    intercs, evaps, ETs, transpis, efloors, SWEs = cpy.update_outarrays(r, d, interc, evap, ET, transpi, efloor, SWE)
                    
                    potinf, efloor, MBE2 = moss.interception(potinf, efloor)       # ground vegetation and moss hydrology 
                    stpout['deltas'][r, d, :] = potinf - transpi                                   # water flux thru soil surface
                    ets[d] = efloor + transpi + interc                             # evapotranspiration components 
                    
                    
                    if d%365==0: print ('  - day #', d, ' hdom ', np.round(np.mean(stand.hdom),2), ' m, ',  
                                        'LAI ', np.round(np.mean(stand.leafarea),2), ' m2 m-2')
        
                    #--------Soil hydrology-----------------
                    stp.run_timestep(d, h0ts_west[d], h0ts_east[d], stpout['deltas'][r, d,:], moss)  # strip/peat hydrology
                    stpout = stp.update_outarrays(r, d, stpout) 
    
                    z, peat_temperature = pt.run_timestep(ta, np.mean(SWE), np.mean(efloor))   # peat temperature in different depths
                    peat_temperatures[r,d,:] = peat_temperature
                    
                    swes[r,d] = np.mean(SWE)                                       # snow water equivalent    
                    d += 1
               #******* End of daily loop***************************** 
            
            #----- Hydrology and temperature-related variables to time-indexed dataframes -----------------
                sday = datetime.datetime(yr, 1, 1)                                 # start day of the year 
                df_peat_temperatures = pd.DataFrame(peat_temperatures[r,start:start+days,:],index=pd.date_range(sday,periods=days)) # all peat temperatures 
                dfwt = pd.DataFrame(stpout['dwts'][r,start:start+days,:],index=pd.date_range(sday,periods=days))                              # daily water table data frame
                dfafp = pd.DataFrame(stpout['afps'][r,start:start+days,:],index=pd.date_range(sday,periods=days))                             # air filled porosity
                
                out.write_cpy(r, start, days, year+1, intercs[r,start:start+days,:], 
                                          evaps[r,start:start+days,:], ETs[r,start:start+days,:],
                                          transpis[r,start:start+days,:], efloors[r,start:start+days,:],
                                          SWEs[r,start:start+days,:])
                
                out.write_temperature(r, start, days, peat_temperatures[r,start:start+days,:])
    
                stp.update_residence_time(dfwt)            
                out.write_strip(r, start, days, yr, year+1, dfwt, stpout, outpara, stp)            
                
            # **************  Biogeochemistry ***********************************
                groundvegetation.run(stand.basalarea, stand.stems, stand.volume,\
                                     stand.dominant.species, ts, ageSim['dominant'])
    
                stand.assimilate(forc.loc[str(yr)], dfwt.loc[str(yr)], dfafp.loc[str(yr)])
                stand.update()
                
                # --------- Locate cuttings here--------------------
                # if yr == 2001: 
                #      stand.dominant.cutting(yr, nut_stat =  stand.nut_stat, to_ba = 8)            
                #      stand.update_lresid()
                #      """ ATTN distinct stem mortaility from the logging -> fate of stems different!!!"""            
                #---------- Organic matter decomposition and nutrient release-----------------
                    
                #---------------- Fertilization --------------------------------
                if yr >= spara['fertilization']['application year']:
                    pH_increment = ferti.ph_effect(yr)
                    esmass.update_soil_pH(pH_increment)
                    esN.update_soil_pH(pH_increment)
                    esP.update_soil_pH(pH_increment)
                    esK.update_soil_pH(pH_increment)
                ferti.nutrient_release(yr)    
                out.write_fertilization(r, year+1, ferti)
                
                """
                TODO:
                    logging resdues to output, to be joined wtih litter in figures
                    harvested volume and biomass to outputs
                    construct balcances at the end of simulation, join to outputs
                """
                nonwoodylitter = (stand.nonwoodylitter + stand.nonwoody_lresid + groundvegetation.nonwoodylitter)/10000.     #conversion kg/ha/yr -> kg/m2/yr
                woodylitter = (stand.woodylitter + stand.woody_lresid) /10000. 
                esmass.run_yr(forc.loc[str(yr)], df_peat_temperatures, dfwt, nonwoodylitter, woodylitter)
                esmass.compose_export(stp)            
                out.write_esom(r, year+1, 'Mass', esmass)
                
                n_nonwoodylitter = (stand.n_nonwoodylitter + stand.n_nonwoody_lresid + groundvegetation.n_litter)/10000. 
                n_woodylitter = (stand.n_woodylitter + stand.n_woody_lresid) /10000. 
                esN.run_yr(forc.loc[str(yr)], df_peat_temperatures, dfwt, n_nonwoodylitter, n_woodylitter)
                out.write_esom(r, year+1, 'N', esN)
                
                p_nonwoodylitter = (stand.p_nonwoodylitter + stand.p_nonwoody_lresid + groundvegetation.p_litter)/10000.
                p_woodylitter = (stand.p_woodylitter + stand.p_woody_lresid)/10000. 
                esP.run_yr(forc.loc[str(yr)], df_peat_temperatures, dfwt, p_nonwoodylitter, p_woodylitter)
                out.write_esom(r, year+1, 'P', esP)
    
                k_nonwoodylitter = (stand.k_nonwoodylitter + stand.k_nonwoody_lresid + groundvegetation.k_litter)/10000. 
                k_woodylitter = (stand.k_woodylitter + stand.k_woody_lresid)/10000. 
                esK.run_yr(forc.loc[str(yr)], df_peat_temperatures, dfwt, k_nonwoodylitter, k_woodylitter)
                out.write_esom(r, year+1, 'K', esK)
    
                stand.update_nutrient_status(groundvegetation, esN.out_root_lyr  + spara['depoN'] + ferti.release['N'] ,\
                                             esP.out_root_lyr + spara['depoP']+ ferti.release['P'] ,\
                                             esK.out_root_lyr + spara['depoK']+ ferti.release['K'] )
                
                # move stand.assimilate here, if first year, take foliage litter from 'table growth (interpolation functions)'    
                # stand.assimilate(forc.loc[str(yr)], dfwt.loc[str(yr)], dfafp.loc[str(yr)])
                # stand.update()
    
                    
                CH4, CH4mean, CH4asCO2eq = ch4s.run_ch4_yr(yr, dfwt)                   # annual ch4 nodewise (kg ha-1 yr-1), mean ch4, and mean ch4 as co2 equivalent
                out.write_methane(r, year+1, CH4)
    
                out.write_stand(r, year+1, stand)
                out.write_canopy_layer(r, year+1, 'dominant', stand.dominant)
                out.write_canopy_layer(r, year+1, 'subdominant', stand.subdominant)
                out.write_canopy_layer(r, year+1, 'under', stand.under)
                out.write_groundvegetation(r, year+1, groundvegetation)
                out.write_export(r, year + 1, esmass)
                
                out.write_nutrient_balance(r, year+1, 'N', esN, spara['depoN'], ferti.release['N'], 
                                           stand.n_demand + stand.n_leaf_demand, groundvegetation.nup)            
                out.write_nutrient_balance(r, year+1, 'P', esP, spara['depoP'], ferti.release['P'], 
                                           stand.p_demand + stand.p_leaf_demand, groundvegetation.pup)            
                out.write_nutrient_balance(r, year+1, 'K', esK, spara['depoK'], ferti.release['K'], 
                                           stand.k_demand + stand.k_leaf_demand, groundvegetation.kup)            
                
                out.write_carbon_balance(r, year+1, stand, groundvegetation, esmass, CH4)
                
                stand.reset_lresid()
                start = start + days                                               # starting point of the next year daily loop                                    
                year +=1                                                           # update the next year for the biogeochemistry loop 
    
        out.close()
        #del stand.dominant
        #del stand.subdominant
        #del stand.under
        #del stand, groundvegetation, esmass, esN, esP, esK, ferti, cpy, moss, stp, pt        
