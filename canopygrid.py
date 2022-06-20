# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:01:50 2017

@author: slauniai

******************************************************************************
CanopyGrid:

Gridded canopy and snow hydrology model for SpaFHy -integration
Based on simple schemes for computing water flows and storages within vegetation
canopy and snowpack at daily or sub-daily timesteps.

(C) Samuli Launiainen, 2016-
last edit: Oct 2018 / Samuli
******************************************************************************

"""
import numpy as np
import configparser
eps = np.finfo(float).eps


class CanopyGrid():
    def __init__(self, cpara, state, outputs=False):
        """
        initializes CanopyGrid -object

        Args:
            cpara - parameter dict:
            state - dict of initial state
            outputs - True saves output grids to list at each timestep

        Returns:
            self - object
            
        NOTE:
            Currently the initialization assumes simulation start 1st Jan, 
            and sets self._LAI_decid and self.X equal to minimum values.
            Also leaf-growth & senescence parameters are intialized to zero.
        """
        epsi = 0.01

        # physiology: transpi + floor evap
        self.physpara = cpara['physpara']

        # phenology
        self.phenopara = cpara['phenopara']
        
        # canopy parameters and state
        self.hc = state['hc'] + epsi
        self.cf = state['cf'] + epsi
        #self.cf = 0.1939 * ba / (0.1939 * ba + 1.69) + epsi
        # canopy closure [-] as function of basal area ba m2ha-1;
        # fitted to Korhonen et al. 2007 Silva Fennica Fig.2
            
        self._LAIconif = state['lai_conif'] + epsi # m2m-2
        self._LAIdecid = state['lai_decid_max']
        self.LAI = self._LAIconif + self._LAIdecid

        self._LAIdecid_max = state['lai_decid_max'] + epsi # m2m-2

        # self.cpara = cpara  # added new parameters self.cpara['kmt'],
        # self.cpara['kmr'] here for testing radiation-based snow melt model
        self.wmax = cpara['interc']['wmax']
        self.wmaxsnow = cpara['interc']['wmaxsnow']
        self.Kmelt = cpara['snow']['kmelt']
        self.Kfreeze = cpara['snow']['kfreeze']
        self.R = cpara['snow']['r']  # max fraction of liquid water in snow

        # --- for computing aerodynamic resistances
        self.zmeas = cpara['flow']['zmeas']
        self.zground =cpara['flow']['zground'] # reference height above ground [m]
        self.zo_ground = cpara['flow']['zo_ground'] # ground roughness length [m]
        self.gsoil = self.physpara['gsoil']
        
        # --- state variables
        self.W = np.minimum(state['w'], self.wmax*self.LAI)
        self.SWE = state['swe']
        self.SWEi = self.SWE
        self.SWEl = np.zeros(np.shape(self.SWE))

        # NOTE: this assumes simulations start 1st Jan each year !!!
        self.X = 0.0

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'PotInf': [], 'Trfall': [], 'Interc': [], 'Evap': [], 'ET': [],
                            'Transpi': [], 'Efloor': [], 'SWE': [], 'LAI': [],
                            'Mbe': [], 'LAIdecid': [], 'erate': [], 'Unload': [],
                            'fact': []}

    def update_amax(self, cpara, nutstat):
        self.physpara['amax'] = cpara['amax']*1.05 + (nutstat - 1) 
        
    def run_timestep(self, doy, dt, Ta, Prec, Rg, Par, VPD,  hc, LAIconif, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0):
        """
        Runs CanopyGrid instance for one timestep
        IN:
            doy - day of year
            dt - timestep [s]
            Ta - air temperature  [degC], scalar or (n x m) -matrix
            prec - precipitatation rate [mm/s]
            Rg - global radiation [Wm-2], scalar or matrix
            Par - photos. act. radiation [Wm-2], scalar or matrix
            VPD - vapor pressure deficit [kPa], scalar or matrix
            U - mean wind speed at ref. height above canopy top [ms-1], scalar or matrix
            CO2 - atm. CO2 mixing ratio [ppm]
            Rew - relative extractable water [-], scalar or matrix
            beta - term for soil evaporation resistance (Wliq/FC) [-]
            P - pressure [Pa], scalar or matrix
        OUT:
            updated CanopyGrid instance state variables
            flux grids PotInf, Trfall, Interc, Evap, ET, MBE [m]
        """
        # SL 18.9.2018: tässä oli kaiken pahan alku ja juuri: pitää päivittää myös tilanmuuttujat self._LAIconif ja self._LAIdecid!
        # näitä käytetään kun lasketaan gsref
        # muutin kutsua niin että siinä annetaan LAIconif ja päivitetään LAI sen kautta. Oletus että 
        self._LAIconif = LAIconif  # oletus on että _LAIdecid = vakio == 0.01
        self.LAI = self._LAIconif + self._LAIdecid
        self.hc = hc

        # Rn = 0.7 * Rg #net radiation
        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2,
                        0.55) * Rg  # Launiainen et al. 2016 GCB, fit to Fig 2a


        """ --- update phenology: self.ddsum & self.X ---"""
        fPheno = self._photoacclim(Ta)

        """ --- aerodynamic conductances --- """
        Ra, Rb, Ras, ustar, Uh, Ug = aerodynamics(self.LAI, self.hc, U, w=0.01, zm=self.zmeas,
                                                  zg=self.zground, zos=self.zo_ground)

        """ --- interception, evaporation and snowpack --- """
        PotInf, Trfall, Evap, Interc, MBE, erate, unload, fact = self.canopy_water_snow(dt, Ta, Prec, Rn, VPD, Ra=Ra)

        """--- dry-canopy evapotranspiration [mm s-1] --- """
        Transpi, Efloor, Gc = self.dry_canopy_et(VPD, Par, Rn, Ta, Ra=Ra, Ras=Ras, CO2=CO2, Rew=Rew, beta=beta, fPheno=fPheno)

        Transpi = Transpi * dt
        Efloor = Efloor * dt
        ET = Transpi + Efloor

        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['PotInf'].append(PotInf)
            self.results['Trfall'].append(Trfall)
            self.results['Interc'].append(Interc)
            self.results['Evap'].append(Evap)
            self.results['ET'].append(ET)
            self.results['Transpi'].append(Transpi)
            self.results['Efloor'].append(Efloor)
            self.results['SWE'].append(self.SWE)
            self.results['LAI'].append(self.LAI)
            self.results['Mbe'].append(np.nanmax(MBE))
            self.results['LAIdecid'].append(self._LAIdecid)
            self.results['erate'].append(erate)
            self.results['Unload'].append(unload)
            self.results['fact'].append(fact)

        return PotInf*1e-3, Trfall*1e-3, Interc*1e-3, Evap*1e-3, ET*1e-3, Transpi*1e-3, Efloor*1e-3, MBE*1e-3, self.SWE   # fQ,fD,fRew

    def _degreeDays(self, T, doy):
        """
        Calculates and updates degree-day sum from the current mean Tair.
        INPUT:
            T - daily mean temperature (degC)
            doy - day of year 1...366 (integer)
        """
        To = 5.0  # threshold temperature
        if doy == 1:  # reset in the beginning of the year
            self.DDsum = 0.
        else:
            self.DDsum += np.maximum(0.0, T - To)

    def _photoacclim(self, T):
        """
        computes new stage of temperature acclimation and phenology modifier.
        Peltoniemi et al. 2015 Bor.Env.Res.
        IN: object, T = daily mean air temperature
        OUT: fPheno - phenology modifier [0...1], updates object state
        """

        self.X = self.X + 1.0 / self.phenopara['tau'] * (T - self.X)  # degC
        S = np.maximum(self.X - self.phenopara['xo'], 0.0)
        fPheno = np.maximum(self.phenopara['fmin'],
                            np.minimum(S / self.phenopara['smax'], 1.0))
        return fPheno

    def _lai_dynamics(self, doy):
        """
        Seasonal cycle of deciduous leaf area

        Args:
            self - object
            doy - day of year

        Returns:
            none, updates state variables self.LAIdecid, self._growth_stage,
            self._senec_stage
        """
        lai_min = self.phenopara['lai_decid_min']
        ddo = self.phenopara['ddo']
        ddur = self.phenopara['ddur']
        sso = self.phenopara['sso']
        sdur = self.phenopara['sdur']

        # growth phase
        if self.DDsum <= ddo:
            f = lai_min
            self._growth_stage = 0.
            self._senesc_stage = 0.
        elif self.DDsum > ddo:
            self._growth_stage += 1.0 / ddur
            f = np. minimum(1.0, lai_min + (1.0 - lai_min) * self._growth_stage)

        # senescence phase
        if doy > sso:
            self._growth_stage = 0.
            self._senesc_stage += 1.0 / sdur
            f = 1.0 - (1.0 - lai_min) * np.minimum(1.0, self._senesc_stage)

        # update self.LAIdecid and total LAI
        self._LAIdecid = self._LAIdecid_max * f
        self.LAI = self._LAIconif + self._LAIdecid
        return f

    def dry_canopy_et(self, D, Qp, AE, Ta, Ra=25.0, Ras=250.0, CO2=380.0, Rew=1.0, beta=1.0, fPheno=1.0):
        """
        Computes ET from 2-layer canopy in absense of intercepted precipitiation,
        i.e. in dry-canopy conditions
        IN:
           self - object
           D - vpd in kPa
           Qp - PAR in Wm-2
           AE - available energy in Wm-2
           Ta - air temperature degC
           Ra - aerodynamic resistance (s/m)
           Ras - soil aerodynamic resistance (s/m)
           CO2 - atm. CO2 mixing ratio (ppm)
           Rew - relative extractable water [-]
           beta - relative soil conductance for evaporation [-]
           fPheno - phenology modifier [-]
        Args:
           Tr - transpiration rate (mm s-1)
           Efloor - forest floor evaporation rate (mm s-1)
           Gc - canopy conductance (integrated stomatal conductance)  (m s-1)
        SOURCES:
        Launiainen et al. (2016). Do the energy fluxes and surface conductance
        of boreal coniferous forests in Europe scale with leaf area?
        Global Change Biol.
        Modified from: Leuning et al. 2008. A Simple surface conductance model
        to estimate regional evaporation using MODIS leaf area index and the
        Penman-Montheith equation. Water. Resources. Res., 44, W10419
        Original idea Kelliher et al. (1995). Maximum conductances for
        evaporation from global vegetation types. Agric. For. Met 85, 135-147

        Samuli Launiainen, Luke
        Last edit: 13.6.2018: TESTING UPSCALING
        """

        # ---Amax and g1 as LAI -weighted average of conifers and decid.

        rhoa = 101300.0 / (8.31 * (Ta + 273.15)) # mol m-3
        
        Amax = 1./self.LAI * (self._LAIconif * self.physpara['amax']
                + self._LAIdecid *self.physpara['amax']) # umolm-2s-1

        g1 = 1./self.LAI * (self._LAIconif * self.physpara['g1_conif']
                + self._LAIdecid *self.physpara['g1_decid']) 

        kp = self.physpara['kp']  # (-) attenuation coefficient for PAR
        q50 = self.physpara['q50']  # Wm-2, half-sat. of leaf light response
        rw = self.physpara['rw']  # rew parameter
        rwmin = self.physpara['rwmin']  # rew parameter

        tau = np.exp(-kp * self.LAI)  # fraction of Qp at ground relative to canopy top

        """--- canopy conductance Gc (integrated stomatal conductance)----- """

        # fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
        fQ = 1./ kp * np.log((kp*Qp + q50) / (kp*Qp*np.exp(-kp * self.LAI) + q50 + eps) )

        # the next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to 
        # Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.        
        # fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

        # soil moisture response: Lagergren & Lindroth, xxxx"""
#        fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))
        fRew=Rew

        # CO2 -response of canopy conductance, derived from APES-simulations
        # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
        fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)
        
        # leaf level light-saturated gs (m/s)
        gs = 1.6*(1.0 + g1 / np.sqrt(D))*Amax / CO2 / rhoa
        
        # canopy conductance
        Gc = gs * fQ * fRew * fCO2 * fPheno
        Gc[np.isnan(Gc)] = eps

        """ --- transpiration rate --- """
        Tr = penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='mm')
        Tr[Tr < 0] = 0.0

        """--- forest floor evaporation rate--- """
        # soil conductance is function of relative water availability
        # gcs = 1. / self.soilrp * beta**2.0
        # beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES
        Gcs = self.gsoil
        
        Efloor = beta * penman_monteith(tau * AE, 1e3*D, Ta, Gcs, 1./Ras, units='mm')
        Efloor[self.SWE > 0] = 0.0  # no evaporation from floor if snow on ground or beta == 0

        return Tr, Efloor, Gc


    def canopy_water_snow(self, dt, T, Prec, AE, D, Ra=25.0, U=2.0):
        """
        Calculates canopy water interception and SWE during timestep dt
        Args: 
            self - object
            dt - timestep [s]
            T - air temperature (degC)
            Prec - precipitation rate during (mm d-1)
            AE - available energy (~net radiation) (Wm-2)
            D - vapor pressure deficit (kPa)
            Ra - canopy aerodynamic resistance (s m-1)
        Returns:
            self - updated state W, Wf, SWE, SWEi, SWEl
            Infil - potential infiltration to soil profile (mm)
            Evap - evaporation / sublimation from canopy store (mm)
            MBE - mass balance error (mm)
        Samuli Launiainen & Ari Laurén 2014 - 2017
        Last edit 12 / 2017
        """

        # quality of precipitation
        Tmin = 0.0  # 'C, below all is snow
        Tmax = 1.0  # 'C, above all is water
        Tmelt = 0.0  # 'C, T when melting starts

        # storage capacities mm
        Wmax = self.wmax * self.LAI
        Wmaxsnow = self.wmaxsnow * self.LAI

        # melting/freezing coefficients mm/s
        Kmelt = self.Kmelt - 1.64 * self.cf / dt  # Kuusisto E, 'Lumi Suomessa'
        Kfreeze = self.Kfreeze

        kp = self.physpara['kp']
        tau = np.exp(-kp*self.LAI)  # fraction of Rn at ground

        # inputs to arrays, needed for indexing later in the code
        gridshape = np.shape(self.LAI)  # rows, cols
    
        if np.shape(T) != gridshape:
            T = np.ones(gridshape) * T
            Prec = np.ones(gridshape) * Prec
            AE = np.ones(gridshape) * AE
            D = np.ones(gridshape) * D
            Ra = np.ones(gridshape) * Ra

        Prec = Prec * dt  # mm
        # latent heat of vaporization (Lv) and sublimation (Ls) J kg-1
        Lv = 1e3 * (3147.5 - 2.37 * (T + 273.15))
        Ls = Lv + 3.3e5

        # compute 'potential' evaporation / sublimation rates for each grid cell
        erate = np.zeros(gridshape)
        ixs = np.where((Prec == 0) & (T <= Tmin))
        ixr = np.where((Prec == 0) & (T > Tmin))
        Ga = 1. / Ra  # aerodynamic conductance

        # resistance for snow sublimation adopted from:
        # Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
        # Best et al. 2011 Geosci. Mod. Dev.
        # ri = (2/3*rhoi*r**2/Dw) / (Ce*Sh*W) == 7.68 / (Ce*Sh*W

        Ce = 0.01*((self.W + eps) / Wmaxsnow)**(-0.4)  # exposure coeff (-)
        Sh = (1.79 + 3.0*U**0.5)  # Sherwood numbner (-)
        gi = Sh*self.W*Ce / 7.68 + eps # m s-1
        # print ixs
        # print('ixs', np.shape(ixs), 'gi', np.shape(gi[ixs]), 'ga', np.shape(Ga[ixs]),
        #      'T', np.shape(T[ixs]), 'AE', np.shape(AE[ixs]), 'tau', np.shape(tau[ixs]))
        erate[ixs] = dt / Ls[ixs] * penman_monteith((1.0 - tau[ixs])*AE[ixs], 1e3*D[ixs], T[ixs], gi[ixs], Ga[ixs], units='W')
#        print('gi', gi, 'Ce', Ce, 'Sh', Sh)

        # evaporation of intercepted water, mm
        gs = 1e6
        erate[ixr] = dt / Lv[ixr] * penman_monteith((1.0 - tau[ixr])*AE[ixr], 1e3*D[ixr], T[ixr],  gs, Ga[ixr],units='W')
        
        # print('erate', erate)

        # ---state of precipitation [as water (fW) or as snow(fS)]
        fW = np.zeros(gridshape)
        fS = np.zeros(gridshape)

        fW[T >= Tmax] = 1.0
        fS[T <= Tmin] = 1.0

        ix = np.where((T > Tmin) & (T < Tmax))
        fW[ix] = (T[ix] - Tmin) / (Tmax - Tmin)
        fS[ix] = 1.0 - fW[ix]
        del ix

        # --- Local fluxes (mm)
        Unload = np.zeros(gridshape)  # snow unloading
        Interc = np.zeros(gridshape)  # interception
        Melt = np.zeros(gridshape)   # melting
        Freeze = np.zeros(gridshape)  # freezing
        Evap = np.zeros(gridshape)

        """ --- initial conditions for calculating mass balance error --"""
        Wo = self.W  # canopy storage
        SWEo = self.SWE  # Snow water equivalent mm

        """ --------- Canopy water storage change -----"""
        # snow unloading from canopy, ensures also that seasonal LAI development does
        # not mess up computations
        ix = (T >= Tmax)
        Unload[ix] = np.maximum(self.W[ix] - Wmax[ix], 0.0)
        self.W = self.W - Unload
        del ix
        # dW = self.W - Wo

        # Interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
        ix = (T < Tmin)
        Interc[ix] = (Wmaxsnow[ix] - self.W[ix]) \
                    * (1.0 - np.exp(-(self.cf[ix] / Wmaxsnow[ix]) * Prec[ix]))
        del ix
        
        # above Tmin, interception capacity equals that of liquid precip
        ix = (T >= Tmin)
        Interc[ix] = np.maximum(0.0, (Wmax[ix] - self.W[ix]))\
                    * (1.0 - np.exp(-(self.cf[ix] / Wmax[ix]) * Prec[ix]))
        del ix
        self.W = self.W + Interc  # new canopy storage, mm

        Trfall = Prec + Unload - Interc  # Throughfall to field layer or snowpack

        # evaporate from canopy and update storage
        Evap = np.minimum(erate, self.W)  # mm
        self.W = self.W - Evap

        """ Snowpack (in case no snow, all Trfall routed to floor) """
        ix = np.where(T >= Tmelt)
        Melt[ix] = np.minimum(self.SWEi[ix], Kmelt[ix] * dt * (T[ix] - Tmelt))  # mm
        del ix
        ix = np.where(T < Tmelt)
        Freeze[ix] = np.minimum(self.SWEl[ix], Kfreeze * dt * (Tmelt - T[ix]))  # mm
        del ix

        # amount of water as ice and liquid in snowpack
        Sice = np.maximum(0.0, self.SWEi + fS * Trfall + Freeze - Melt)
        Sliq = np.maximum(0.0, self.SWEl + fW * Trfall - Freeze + Melt)

        PotInf = np.maximum(0.0, Sliq - Sice * self.R)  # mm
        Sliq = np.maximum(0.0, Sliq - PotInf)  # mm, liquid water in snow

        # update Snowpack state variables
        self.SWEl = Sliq
        self.SWEi = Sice
        self.SWE = self.SWEl + self.SWEi
        
        # mass-balance error mm
        MBE = (self.W + self.SWE) - (Wo + SWEo) - (Prec - Evap - PotInf)

        # MBEsnow = self.SWE - SWEo - (Trfall - PotInf)
#        if np.nanmax(MBEsnow) > 0.1:
#            print('MBEsnow', MBEsnow)
#            print('Freeze', Freeze)
#            print('melt', Melt)
#            print('Potinf', PotInf)
#            print('Trfall', Trfall)
#            print('SWE', self.SWE)
#            print('SWo', SWEo)
#            print('dSWE', self.SWE-SWEo)
#            print('T', T)
#            print('fS', fS)
#            print('fW', fW)
#        MBEcan = (self.W - Wo) - (Interc - Evap - Unload)
#        if np.nanmax(MBEcan) > 0.1:
#            print('Mbecan', MBEcan)
#            print('Unload', Unload)
#            print('dW', dW)

        return PotInf, Trfall, Evap, Interc, MBE, erate, Unload, fS + fW



    def create_outarrays(self, nrounds, ndays, ncols):
        self.intercs = np.zeros((nrounds, ndays,ncols), dtype=float)           
        self.evaps = np.zeros((nrounds, ndays,ncols), dtype=float)                  
        self.ETs = np.zeros((nrounds, ndays,ncols), dtype=float)                   
        self.transpis =np.zeros((nrounds,ndays, ncols), dtype=float)                   
        self.efloors =np.zeros((nrounds,ndays, ncols), dtype=float) 
        self.SWEs =np.zeros((nrounds,ndays, ncols), dtype=float)
       
        return self.intercs, self.evaps, self.ETs, self.transpis, self.efloors, self.SWEs
    
    def update_outarrays(self, scen, d, interc, evap, ET, transpi, efloor, SWE):
        self.intercs[scen, d, :] = interc           
        self.evaps[scen, d, :] = evap                  
        self.ETs[scen, d, :] = ET                   
        self.transpis[scen, d, :] = transpi                   
        self.efloors[scen, d, :] = efloor 
        self.SWEs[scen, d, :] = SWE
        
        return self.intercs, self.evaps, self.ETs, self.transpis, self.efloors, self.SWEs
""" *********** utility functions ******** """

# @staticmethod
def degreeDays(dd0, T, Tbase, doy):
    """
    Calculates degree-day sum from the current mean Tair.
    INPUT:
        dd0 - previous degree-day sum (degC)
        T - daily mean temperature (degC)
        Tbase - base temperature at which accumulation starts (degC)
        doy - day of year 1...366 (integer)
    OUTPUT:
        x - degree-day sum (degC)
   """
    if doy == 1:  # reset in the beginning of the year
        dd0 = 0.
    return dd0 + max(0, T - Tbase)


# @staticmethod
def eq_evap(AE, T, P=101300.0, units='W'):
    """
    Calculates the equilibrium evaporation according to McNaughton & Spriggs,\
    1986.
    INPUT:
        AE - Available energy (Wm-2)
        T - air temperature (degC)
        P - pressure (Pa)
        units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT:
        equilibrium evaporation rate (Wm-2)
    """
    Mw = 18e-3  # kg mol-1
    # latent heat of vaporization of water [J/kg]
    L = 1e3 * (2500.8 - 2.36 * T + 1.6e-3 * T ** 2 - 6e-5 * T ** 3)
    # latent heat of sublimation [J/kg]
    if T < 0:
        L = 1e3 * (2834.1 - 0.29 * T - 0.004 * T ** 2)

    _, s, g = e_sat(T, P)

    x = np.divide((AE * s), (s + g))  # Wm-2 = Js-1m-2
    if units == 'mm':
        x = x / L  # kg m-2 s-1 = mm s-1
    elif units == 'mol':
        x = x / L / Mw  # mol m-2 s-1
    x = np.maximum(x, 0.0)
    return x


# @staticmethod
def e_sat(T, P=101300):
    """
    Computes saturation vapor pressure (Pa), slope of vapor pressure curve
    [Pa K-1]  and psychrometric constant [Pa K-1]
    IN:
        T - air temperature (degC)
        P - ambient pressure (Pa)
    OUT:
        esa - saturation vapor pressure in Pa
        s - slope of saturation vapor pressure curve (Pa K-1)
        g - psychrometric constant (Pa K-1)
    """
    NT = 273.15
    cp = 1004.67  # J/kg/K

    Lambda = 1e3 * (3147.5 - 2.37 * (T + NT))  # lat heat of vapor [J/kg]
    esa = 1e3 * (0.6112 * np.exp((17.67 * T) / (T + 273.16 - 29.66)))  # Pa

    s = 17.502 * 240.97 * esa / ((240.97 + T) ** 2)
    g = P * cp / (0.622 * Lambda)
    return esa, s, g


# @staticmethod
def penman_monteith(AE, D, T, Gs, Ga, P=101300.0, units='W'):
    """
    Computes latent heat flux LE (Wm-2) i.e evapotranspiration rate ET (mm/s)
    from Penman-Monteith equation
    INPUT:
       AE - available energy [Wm-2]
       VPD - vapor pressure deficit [Pa]
       T - ambient air temperature [degC]
       Gs - surface conductance [ms-1]
       Ga - aerodynamic conductance [ms-1]
       P - ambient pressure [Pa]
       units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT:
       x - evaporation rate in 'units'
    """
    # --- constants
    cp = 1004.67  # J kg-1 K-1
    rho = 1.25  # kg m-3
    Mw = 18e-3  # kg mol-1
    _, s, g = e_sat(T, P)  # slope of sat. vapor pressure, psycrom const
    L = 1e3 * (3147.5 - 2.37 * (T + 273.15))

    x = (s * AE + rho * cp * Ga * D) / (s + g * (1.0 + Ga / Gs))  # Wm-2

    if units is 'mm':
        x = x / L  # kgm-2s-1 = mms-1
    if units is 'mol':
        x = x / L / Mw  # mol m-2 s-1

    x = np.maximum(x, 0.0)
    return x


# @staticmethod
# def aerodynamic_conductance_from_ust(Ust, U, Stanton):
#    """
#    computes canopy aerodynamic conductance (ms-1) from frict. velocity
#    IN:
#       Ustar - friction velocity (ms-1)
#       U - mean wind speed at flux measurement heigth (ms-1)
#       Stanton - Stanton number (kB-1) for quasi-laminar boundary layer
#           resistance. Typically kB=1...12, use 2 for vegetation ecosystems
#           (Verma, 1989, Garratt and Hicks, 1973)
#    OUT:
#       Ga - aerodynamic conductance [ms-1]
#    """
#    kv = 0.4  # von Karman constant
#    ra = U / (Ust ** 2.0 + eps) + Stanton / (kv * (Ust + eps))  # sm-1
#    Ga = 1.0 / ra  # ms-1
#    return Ga

def aerodynamics(LAI, hc, Uo, w=0.01, zm=2.0, zg=0.5, zos=0.01):
    """
    computes wind speed at ground and canopy + boundary layer conductances
    Computes wind speed at ground height assuming logarithmic profile above and
    exponential within canopy
    Args:
        LAI - one-sided leaf-area /plant area index (m2m-2)
        hc - canopy height (m)
        Uo - mean wind speed at height zm (ms-1)
        w - leaf length scale (m)
        zm - wind speed measurement height above canopy (m)
        zg - height above ground where Ug is computed (m)
        zos - forest floor roughness length, ~ 0.1*roughness element height (m)
    Returns:
        ra - canopy aerodynamic resistance (sm-1)
        rb - canopy boundary layer resistance (sm-1)
        ras - forest floor aerod. resistance (sm-1)
        ustar - friction velocity (ms-1)
        Uh - wind speed at hc (ms-1)
        Ug - wind speed at zg (ms-1)
    SOURCE:
       Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
       Massman 1987, BLM 40, 179 - 197.
       Magnani et al. 1998 Plant Cell Env.
    """
    zm = hc + zm  # m
    kv = 0.4  # von Karman constant (-)
    beta = 285.0  # s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
    alpha = LAI / 2.0  # wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*hc  # m
    zom = 0.123*hc  # m
    zov = 0.1*zom
    zosv = 0.1*zos

    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / np.log((zm - d) / zom) 
    Uh = ustar / kv * np.log((hc - d) / zom)
    
    # U(zg) from exponential wind profile
    zn = np.minimum(zg / hc, 1.0)  # zground can't be above canopy top
    Ug = Uh * np.exp(alpha*(zn - 1.0))

    # canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    #ra = 1. / (kv*ustar) * np.log((zm - d) / zom)
    ra = 1./(kv**2.0 * Uo) * np.log((zm - d) / zom) * np.log((zm - d) / zov)    
    rb = 1. / LAI * beta * ((w / Uh)*(alpha / (1.0 - np.exp(-alpha / 2.0))))**0.5

    # soil aerodynamic resistance (sm-1)
    ras = 1. / (kv**2.0*Ug) * (np.log(zg / zos))*np.log(zg / (zosv))
    
    #print('ra', ra, 'rb', rb)
    ra = ra + rb
    return ra, rb, ras, ustar, Uh, Ug


#def aerodynamic_conductance(hc, U, zm=2.0):
#    """
#    Aerodynamic conductance from log-law neglecting stability effects.
#    Args:
#        hc - canopy height (m)
#        U - mean wind speed at zm
#        zm - measurement height above canopy (m)
#    Returns:
#        ga - aerodynamic conductance for water vapor (ms-1)
#
#    Source:
#        Leuning et al. 2008 WRR
#        Eq: ga = U*k**2 / [ln((zm - d) / zom) * ln((zm - d) / dov) ]
#    """
#
#    zm = zm + hc
#    d = 0.66*hc  # displacement heigth (m)
#    zom = 0.123*hc  # roughness lenght for momentum (m)
#    zov = 0.1*zom  # roughness length for H2O (m)
#    # k**2 = 0.16
#    ga = 0.16*U / (np.log((zm - d) / zom) * np.log((zm - d) / zov))
#
#    return ga


def wind_profile(LAI, hc, Uo, z, zm=2.0, zg=0.2):
    """
    Computes wind speed at ground height assuming logarithmic profile above and
    hyperbolic cosine profile within canopy
    INPUT:
        LAI - one-sided leaf-area /plant area index (m2m-2)
        hc - canopy height (m)
        Uo - mean wind speed at height zm (ms-1)
        zm - wind speed measurement height above canopy (m)
        zg - height above ground where U is computed
    OUTPUT:
        Uh - wind speed at hc (ms-1)
        Ug - wind speed at zg (ms-1)
    SOURCE:
       Cammalleri et al. 2010 Hydrol. Earth Syst. Sci
       Massman 1987, BLM 40, 179 - 197.
    """

    k = 0.4  # von Karman const
    Cd = 0.2  # drag coeff
    alpha = 1.5  # (-)

    zm = zm + hc
    d = 0.66*hc
    zom = 0.123*hc
    beta = 4.0 * Cd * LAI / (k**2.0*alpha**2.0)
    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * k / np.log((zm - d) / zom)  # m/s

    U = np.ones(len(z))*np.NaN

    # above canopy top wind profile is logarithmic
    U[z >= hc] = ustar / k * np.log((z[z >= hc] - d) / zom)

    # at canopy top, match log and exponential profiles
    Uh = ustar / k * np.log((hc - d) / zom)  # m/s

    # within canopy hyperbolic cosine profile
    U[z <= hc] = Uh * (np.cosh(beta * z[z <= hc] / hc) / np.cosh(beta))**0.5

    return U, ustar, Uh


def daylength(LAT, LON, DOY):
    """
    Computes daylength from location and day of year.

    Args:
        LAT, LON - in deg, float or arrays of floats
        doy - day of year, float or arrays of floats

    Returns:
        dl - daylength (hours), float or arrays of floats
    """
    CF = np.pi / 180.0  # conversion deg -->rad

    LAT = LAT*CF
    LON = LON*CF

    # ---> compute declination angle
    xx = 278.97 + 0.9856*DOY + 1.9165*np.sin((356.6 + 0.9856*DOY)*CF)
    DECL = np.arcsin(0.39785*np.sin(xx*CF))
    del xx

    # --- compute day length, the period when sun is above horizon
    # i.e. neglects civil twilight conditions
    cosZEN = 0.0
    dl = 2.0*np.arccos(cosZEN - np.sin(LAT)*np.sin(DECL) / (np.cos(LAT)*np.cos(DECL))) / CF / 15.0  # hours

    return dl


def read_ini(inifile):
    """read_ini(inifile): reads canopygrid.ini parameter file into pp dict"""

    cfg = configparser.ConfigParser()
    cfg.read(inifile)

    pp = {}
    for s in cfg.sections():
        section = s.encode('ascii', 'ignore')
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k.encode('ascii', 'ignore')
            val = v.encode('ascii', 'ignore')
            if section == 'General':  # 'general' section
                pp[section][key] = val
            else:
                pp[section][key] = float(val)

    pp['General']['dt'] = float(pp['General']['dt'])

    pgen = pp['General']
    cpara = pp['CanopyGrid']
    return pgen, cpara
