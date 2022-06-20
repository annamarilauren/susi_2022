    # -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 10:41:03 2022

@author: alauren
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.sparse import  diags
from susi_utils import peat_hydrol_properties, wrc

class Esom():
  def __init__(self, spara, sfc, days, substance = 'Mass'):
    """
    Input 
        shape_area shape of the computation domain, tuple (x, y)
        ash litter ash content %
        litter N content %
        sfc site fertility class, array of integers between 1 and 6, shape as shape_area
        days length of the simuilation in days to be used in the output array
        substance defines what is this instance for:
                                     can be 'Mass'              - organic matter decomposition
                                              'N'               - nitrogen dynamics
                                              'P'               - phosphorus dynamics
                                              'K'               - potassium dynamics
    TODO:
        change self.M update to only one place, so it can be called outside
        create an initialization loop where mor layer M[2,3,4,5,6] are iterated with a given wt, temperature and litterfall
        
    """
    
    # -------------- Parameters -----------------------------------------
    
    self.substance = substance
    # Mass instance accounts for combined CO2 and DOC release from organic matter 
    # Lauren et al 2012, Lappalainen et al. 2018 ja Laurén et al 2019
    # L: DOC / CO2-C = 0.1; F,H, peat DOC / CO2-C = 0.05
    # added here to Mass instance relase k1: L to out 0.1, k2 F to out 0.05 ; H to peat 0.05
                                                                            
    nutcpara = {'Mass': {'k1':1.05, 'k2': 1.05, 'k6':1.05},                     # modifiers for nutrient release in comparison to mass release
            #'N':{'k1':0.1, 'k2': 0.5, 'k6':0.5},
            'N':{'k1':0.1, 'k2': 0.1, 'k6':0.125},
            #'P':{'k1':1.1, 'k2': 1.1, 'k6':1.0},
            'P':{'k1':1.1, 'k2': 0.4, 'k6':0.3},
            'K':{'k1':1.5, 'k2': 1.5, 'k6':1.5}}

    self.nutc = nutcpara[substance]         

    self.contpara = {'Mass': {2: {1:100., 2:100}, 3: {1:100., 2:100}, 4: {1:100., 2:100}, 5: {1:100., 2:100}},        # Content of mass, N, P, K in peat, unit gravimettric %
                'N':{2: {1:1.9, 2:1.9}, 3: {1:1.6, 2:1.6}, 4: {1:1.4, 2:1.4}, 5: {1:1.2, 2:1.2}},                # Unit gravimetric % Mese study           
                'P':{2: {1:0.1,2:0.1}, 3: {1:0.08, 2:0.08}, 4: {1:0.06, 2:0.06}, 5: {1:0.05, 2:0.05}}, 
                'K':{2: {1:0.045, 2:0.045}, 3: {1:0.04, 2:0.038}, 4: {1:0.037, 2:0.034}, 5: {1:0.03, 2:0.03}}    
                }
    
    keys_in_spara ={'N':'peatN', 'P': 'peatP', 'K':'peatK'}
    if self.substance != 'Mass':
        sfcs = np.unique(sfc)
        if spara[keys_in_spara[self.substance]] is not None:
            for s in sfcs:
                self.contpara[self.substance][s][1] = spara[keys_in_spara[self.substance]]
                self.contpara[self.substance][s][2] = spara[keys_in_spara[self.substance]]
    self.enable_peattop = spara['enable_peattop']
    self.enable_peatmiddle = spara['enable_peatmiddle']
    self.enable_peatbottom = spara['enable_peatbottom']
    
    self.contpara_mor = {'Mass': 100.0, 'N': 1.5, 'P': 0.05, 'K': 0.03}         # Initial concentration of mass components in mor layer, gravimetric %

    self.dph = {1:3.5, 2:3.4, 3:3.3, 4:3.2, 5:3.1, 6:3.0}                      # pH according to site fertility class

    #------------These from spara dictionary
    x=1; y=spara['n']                          # shape of the domain    
    shape_area = (x,y)                         # input shape 

    self.nLyrs = spara['nLyrs']                   # number of soil layers
    self.dz = np.ones(self.nLyrs)*spara['dzLyr']            # thickness of layers, m
    self.z = np.cumsum(self.dz)-self.dz/2.                   # depth of the layer center point, m 
    if spara['vonP']:
        vpost = spara['vonP bottom']*np.ones(self.nLyrs)
        vpost[:len(spara['vonP top'])] = spara['vonP top']
        self.bd =  0.035 + 0.0159 * vpost                           # Päivänen 1973 page 36 Figure 9 unit g cm-3
                                                      
    else:
        self.bd = spara['bd bottom']*np.ones(self.nLyrs)
        self.bd[:len(spara['bd top'])] = spara['bd top']
    # ---------this to array------------
    self.sfc_specification = 1 #np.ones(shape_area, dtype=int)                        # MTkg1, MTkg2
    
    self.bound1 = 0.3                                                                # boundary between top and middle layer, m
    self.bound2 = 0.6                                                                # boundary between middle and bottom layers, m
    self.i = 0                                                                  # day counter
    self.x,self.y = shape_area                                                            # shape of the computation domain 
    self.mass = np.zeros((self.x,self.y,11,days))                                         # model output in four dimensions (area: x,y, storages 0:8, time)
    self.pH = np.zeros(shape_area)
    self.ash = np.ones(shape_area)* 5.0                                        # %
    self.litterN = np.ones(shape_area)* 1.2                                    # %    
    self.sfc = np.expand_dims(sfc, axis=0)
    self.reset_storages()

    # 
    gwl=np.linspace(0,-6,150)
    pF, _ = peat_hydrol_properties(self.bd[self.idtop],  var='bd', ptype='A') # peat hydraulic properties after Päivänen 1973    
    water_sto = [sum(wrc(pF, x = np.minimum(self.z[self.idtop]+g, 0.0))*self.dz[self.idtop]) for g in gwl]     #equilibrium head m
    volume_fraction_of_air = (water_sto[0] - water_sto)/water_sto[0]
    self.wtToVfAir_top = interp1d(gwl, volume_fraction_of_air, fill_value=(volume_fraction_of_air[0], volume_fraction_of_air[-1]), bounds_error=False)
    
    pF, _ = peat_hydrol_properties(self.bd[self.idmiddle],  var='bd', ptype='A') # peat hydraulic properties after Päivänen 1973    
    water_sto = [sum(wrc(pF, x = np.minimum(self.z[self.idmiddle]+g, 0.0))*self.dz[self.idmiddle]) for g in gwl]     #equilibrium head m
    water_sto
    volume_fraction_of_air = (water_sto[0] - water_sto)/water_sto[0]
    self.wtToVfAir_middle = interp1d(gwl, volume_fraction_of_air, fill_value=(volume_fraction_of_air[0], volume_fraction_of_air[-1]), bounds_error=False)
    
    pF, _ = peat_hydrol_properties(self.bd[self.idbottom],  var='bd', ptype='A') # peat hydraulic properties after Päivänen 1973    
    water_sto = [sum(wrc(pF, x = np.minimum(self.z[self.idbottom]+g, 0.0))*self.dz[self.idbottom]) for g in gwl]     #equilibrium head m
    water_sto
    volume_fraction_of_air = (water_sto[0] - water_sto)/water_sto[0]
    self.wtToVfAir_bottom = interp1d(gwl, volume_fraction_of_air, fill_value=(volume_fraction_of_air[0], volume_fraction_of_air[-1]), bounds_error=False)
    
    self.pF, _ = peat_hydrol_properties(self.bd,  var='bd', ptype='A') # peat hydraulic properties after Päivänen 1973    


    #temperature_functions:
    self.t2 = interp1d([-40.,-5., -1, 25., 35., 60.],         # effect of temperature on the decomposition rate
                [0.,   0.,  0.2, 1.53, 1.53, 0.])
    self.t3 = interp1d([-40., -3., 0., 7., 60.],
                   [0., 0., 1.3, 1.3, 0.])
    self.t4 = interp1d([-40., -5., 1., 20., 40., 80.],
                   [0., 0., 0.2, 1., 1., 0.])
    self.t5 = interp1d([-40., -5., 1., 13., 25., 50.],
                   [0., 0., 0.2, 1., 1., 0.])
    self.t6 = interp1d([-40., -5., 1., 27.5, 35., 60.],
                  [0., 0., 0.2, 1.95, 1.95, 0.])
    #self.t6 = interp1d([-40., -30., -20. ,-10.,   0.,  10.,  20.,  30.,  40.,  50.],             #Q10 = 2
    #               [ 0.03125,0.0625, 0.125, 0.25, 0.5, 1., 2., 4., 4., 4.]) 
    
    #self.t7 = interp1d([-40., -5., 1., 27.5, 35., 60.],
    #               [0., 0., 0.2, 1.95, 1.95, 0.])
    self.t7 = interp1d([-40., -30., -20. ,-10.,   0.,  10.,  20.,  30.,  40.,  50.],             #Q10 = 2
                   [ 0.03125,0.0625, 0.125, 0.25, 0.5, 1., 2., 4., 4., 4.]) 

    

    # moisture_functions:
    # Description of Romul model Table 2
    self.phi1236 = interp1d([0.02, 0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4, 0.417, 1.333, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8,4], 
                      [0.0, 0.004, 0.026, 0.074, 0.154, 0.271, 0.432, 0.64,  0.899, 1.0, 1.0, 0.844, 0.508, 0.305, 0.184, 0.111, 0.067, 0.04, 0.024, 0])
    self.phi4 = interp1d([0.0, 0.133, 1.333, 2.333,4.0],
                    [0.0, 1.0, 1.0 ,0.0, 0.0])
    self.phi5 = interp1d([0.0, 0.067, 0.5, 2.333, 4.0, 10.0],
                     [0.0, 0.0, 1.0, 1.0 ,0.0, 0.0])
    

  #def lignin_corrections(self, nitrogen = 0.7, lignin = 25.0):
      #N contents, branches Skonieczna,et al 2014 unit %
      #Lignin contents Kilpeläinen et al. 2003 unit %
    nitrogen = 0.7; lignin = 25.0
    adjust = 2. #4.
    self.mu_k1 = 0.092 * (lignin / nitrogen)**-0.7396 * adjust
    self.mu_k2 = 0.0027 * (lignin / nitrogen)**-0.3917 * adjust
    self.mu_k3 = 0.062 * (lignin / nitrogen)**-0.3972 * adjust
    #return  mu_k1, mu_k2, mu_k3
    self.out_root_lyr = np.zeros(self.y)
    self.out_below_root_lyr =  np.zeros(self.y)
    
  def reset_storages(self): 

      #initial values for the storages -> make these to dictionary and input values
      self.previous_mass = np.zeros(self.y)    
      
      self.i = 0
      h_humus = 0.01                                          # Mor humus thickness in (m)
      rho_humus = 120.0                                       # Bulk density of the mor (kg m3)
      frac_L = 0.1                                            # Share of undecomposed litter (L) from the mor thickness (fraction 0...1)
      frac_F = 0.2                                            # Share of partly decomposed F material from the mor thickness (fraction 0...1)
      frac_H = 0.7                                            # Share of humified H material from the mor thickness (fraction 0...1)
      frac_leaf = 0.5                                         # Share of non-woody material from L and F material (fraction 0...1)
      frac_woody = 0.5                                        # Share of non-woody material from L and F material (fraction 0...1)

      LL_mass = h_humus*frac_L*rho_humus * frac_leaf * self.contpara_mor[self.substance] / 100.
      LW_mass = h_humus*frac_L*rho_humus * frac_woody* self.contpara_mor[self.substance] / 100.

      FL_mass = h_humus*frac_F*rho_humus * frac_leaf* self.contpara_mor[self.substance] / 100.
      FW_mass = h_humus*frac_F*rho_humus * frac_woody* self.contpara_mor[self.substance] / 100.

      H_mass = h_humus*frac_H*rho_humus* self.contpara_mor[self.substance] / 100.

      # ****** Initialize the Mass matrix, that contains the current mass of matter in different stages of decomposition
      self.M_shape = (self.x,self.y,11)
      self.M = np.zeros(self.M_shape, dtype = float)
      self.M[:,:,0] = 0.                              # Fresh litter as input: needles and fine roots, kg m-2
      self.M[:,:,1] = 0.                              # Woody debris as input: branches and coarse roots, kg m-2
      self.M[:,:,2] = LL_mass                         # Leaf litter storage, kg m-2
      self.M[:,:,3] = LW_mass                         # Woody litter storage, kg m-2
      self.M[:,:,4] = FL_mass                         # F material storage (leaf & fine roots), kg m-2
      self.M[:,:,5] = FW_mass                         # Woody F material (braches and coarse roots), kg m-2
      self.M[:,:,6] = H_mass                          # Humus material storage, kg m-2
                                                        # Peat mass initialized in the loop below
      self.M[:,:,10] = 0.                              # Cumulative output in organic material kg m-2        

      for scode in np.unique(np.ravel(self.sfc)):
          ix = np.where(self.sfc==scode)
          cont_array = np.ones(self.nLyrs)*self.contpara[self.substance][scode][self.sfc_specification] / 100.  
          self.idtop = np.where(self.z < self.bound1)
          self.idmiddle = np.where((self.z < self.bound2) & ( self.z > self.bound1))
          self.idbottom = np.where(self.z > self.bound2)
          self.pH[ix[0],ix[1]] = self.dph[scode]
          self.M[ix[0],ix[1],7] = sum(self.bd[self.idtop]*self.dz[self.idtop]*1000.*cont_array[self.idtop])                     # Peat material storage, top, kg m-2
          self.M[ix[0],ix[1],8] = sum(self.bd[self.idmiddle]*self.dz[self.idmiddle]*1000*cont_array[self.idmiddle])     # Peat material storage, middle, kg m-2
          self.M[ix[0],ix[1],9] = sum(self.bd[self.idbottom]*self.dz[self.idbottom]*1000*cont_array[self.idbottom])   # Peat material storage, bottom, kg m-2

  def update_soil_pH(self, increment):
      for scode in np.unique(np.ravel(self.sfc)):
          ix = np.where(self.sfc==scode)
          self.pH[ix[0],ix[1]] = self.dph[scode] + increment

      
  def get_rates(self, tair, tp_top, tp_middle, tp_bottom, wn, peat_w1, peat_w2, peat_w3, H_w):
      """
      ash, ash content in gravimetric %
      nitrogen, N content in gravimetric %
      t temperature in deg C
      t2...t7 temperature functions
      wn normalaized water content w/wfc
      phi1236, phi4, phi5 moisture functions
      """
      nu = np.clip(0, 0.701*self.pH -1.6018 - 0.038*self.pH**2, 1)   # ph Romul documentation Table 1

      k1= (0.002 + 0.00009 * self.ash + 0.003 * self.litterN) * min(0.1754 * np.exp(0.0871 * tair), 1.)*self.phi1236(wn)*nu # adjusted decomposition rates
      k2= np.clip((0.00114 -0.00028 * self.litterN)*self.t2(tair)*self.phi1236(wn)*nu, 0., 1.)     #
      k3= np.clip((0.04 - 0.003*self.litterN) * self.t3(tair) * self.phi1236(wn), 0., 1.) 
      k4= 0.005 * self.litterN * self.t4(tair) * self.phi4(wn) 
      k5= 0.007 * self.t5(tair) * self.phi5(wn)
      k6= 0.0006 * self.t6(tp_top) * self.phi1236(wn) #* 0.5
      #k6= 0.0006*self.t6(tp_top)*H_w 
      
      #THESE can be modified by you
      k7= 0.0001 * self.t7(tp_top) * peat_w1 * self.enable_peattop * 1.0 #2.5 #4.0  #5.0       #Change this                                # Lappalainen et al 2018, gamma/VfAir slightly decomposed peat
      k8 = 0.0001 * self.t7(tp_middle) * peat_w2 * self.enable_peatmiddle * 1.0 #2.0 #1.0                                               # Lappalainen et al. 2018 gamma/VfAir highly decomposed
      k9 = 0.0001 *self.t7(tp_bottom) * peat_w3 *self.enable_peatbottom * 0.33  #0.5

      # -> temperature separately for top peat 30 cm (take from 15 cm)
      # -> air filled porosity separately for the top and bottom ()

      return (k1, k2, k3, k4, k5, k6, k7, k8, k9)


  def decompose(self, k1, k2,  k3,  k4, k5, k6, k7, k8, k9, M):
      """
      Main matrix contains 11 storages 
      0 - L0L input of leaf and fine root litter
      1 - L0W input of branch and coarse root, i.e. woody litter
      2 - LL storage of leaf litter
      3 - LW storage of woody litter
      4 - FL storage of leaf F material
      5 - FW storage of woody F material 
      6 - H storage of humus material from 4 and 5
      7 - P1 storage of peat; depth 0-30 cm
      8 - P2 storage of peat; depth 30-60 cm
      9 - P3 storage of peat; depth 60 - bottom cm
      10 - Out cumulative output of mass
      The diagonals and subdiagonals are indexed using the above codes, arrangement to the sparse matrix takes into account 
      the different length of subdiagonals 
      """
      
      diagonal_shape = (self.x,self.y,11)
      #Staying fraction
      k_diag = np.zeros(diagonal_shape, dtype = float)
      k_diag[:,:,2] = 1-(self.nutc['k1'] * k1 +  k3)                  
      k_diag[:,:,3] = 1-(self.nutc['k1'] * k1 * self.mu_k1 +  k3 * self.mu_k3)                 
      k_diag[:,:,4] = 1-(self.nutc['k2'] * k2 + k4 +  k5)                          
      k_diag[:,:,5] = 1-(self.nutc['k2'] * k2 * self.mu_k2 + k4 +  k5)                           
      k_diag[:,:,6] = 1-(self.nutc['k6'] * k6)                          
      k_diag[:,:,7] = 1-(self.nutc['k6'] * k7)                        
      k_diag[:,:,8] = 1-(self.nutc['k6'] * k8)                         
      k_diag[:,:,9] = 1-(self.nutc['k6'] * k9)                        
      k_diag[:,:,10] = 1                          
      k_diag = np.ravel(k_diag)
      
      #1st subdiagonal
      k_low0 = np.zeros(diagonal_shape, dtype = float)
      k_low0[:,:,6] = k4 + k5                                # to H from FW 
      k_low0[:,:,10] = (self.nutc['k6'] * k9)                # to Out from P3 
      k_low0 = np.ravel(k_low0)

      #2nd subdiagonal
      k_low1 = np.zeros(diagonal_shape, dtype = float)
      k_low1[:,:,2] = 1                                      # to LL from L0L 
      k_low1[:,:,3] = 1                                      # to LW from L0W
      k_low1[:,:,4] = k3                                     # to FL from LL
      k_low1[:,:,5] = k3 * self.mu_k3                        # to FW from LW
      k_low1[:,:,6] = k4 + k5                                # to H from FL
      k_low1[:,:,10] = self.nutc['k6'] * k8                  # to Out from P2
      k_low1 = np.ravel(k_low1)

      #3rd subdiagonal, only to Out
      k_low2 = np.zeros(diagonal_shape, dtype = float)
      k_low2[:,:,10] = self.nutc['k6'] * k7                  # to Out from P1
      k_low2 = np.ravel(k_low2)
      
      k_low3 = np.zeros(diagonal_shape, dtype = float)
      k_low3[:,:,10] = self.nutc['k6'] * k6                  # to Out from H
      k_low3 = np.ravel(k_low3)

      k_low4 = np.zeros(diagonal_shape, dtype = float)
      k_low4[:,:,10] = self.nutc['k2'] * k2 * self.mu_k2     # to Out from FW
      k_low4 = np.ravel(k_low4)

      k_low5 = np.zeros(diagonal_shape, dtype = float)
      k_low5[:,:,10] = self.nutc['k2'] * k2                  # to Out from FL
      k_low5 = np.ravel(k_low5)

      k_low6 = np.zeros(diagonal_shape, dtype = float)
      k_low6[:,:,10] = self.nutc['k1'] * k1 * self.mu_k1     # to Out from LW
      k_low6 = np.ravel(k_low6)

      k_low7 = np.zeros(diagonal_shape, dtype = float)
      k_low7[:,:,10] = self.nutc['k1'] * k1                  # to Out from LL
      k_low7 = np.ravel(k_low7)

      # Now we locate the elements to a matrix
      length = len(k_diag)
      kmat =  diags(diagonals=[k_diag, k_low0[1:], k_low1[2:],
                              k_low2[3:], k_low3[4:], k_low4[5:], 
                              k_low5[6:], k_low6[7:], k_low7[8:],],       # Create the main matrix
          offsets=[0, -1, -2, -3, -4, -5, -6, -7, -8], shape=(length,length),
          format='csr')

      #with np.printoptions(edgeitems=30, linewidth=100000,formatter=dict(float=lambda x: "%.3g" % x)):
      #    kdense = kmat.todense()
      #    print (kdense)
      
      M_tmp = np.ravel(self.M)
      M_tmp = kmat@M_tmp
      self.M = np.reshape(M_tmp,(self.x, self.y,11)) 
      return self.M

  def run_yr(self, weather, df_peat_temperatures, water_tables, nonwoodylitter, woodylitter):
    
    ini_i = self.i                                                             # Day calculator, set the first day of the year
    P1_ini = self.M[:,:, 7]*10000.                                             # top layer peat mass in the in the beginning of yr 
    P2_ini = self.M[:,:, 8]*10000.                                             # middle layer peat mass in the in the beginning of yr
    P3_ini = self.M[:,:, 9]*10000.                                             # bottom layer peat mass in the in the beginning of yr 
    
    L0L = np.zeros((self.x, self.y))                                           # Litterfall time series (leaf & fineroots), kg m-2, along the strip 
    L0L[0,:] = nonwoodylitter                                                  # Leaf and fine root litter, kg m-2, locate the litter to array
    self.nonwoodylitter = nonwoodylitter                                       # litter to instance variables       
    self.woodylitter = woodylitter                                             # litter to instance variables
    
    L0W = np.zeros((self.x, self.y))                                           # Litterfall time series (woody litter), kg m-2, empty strip shape array
    L0W[0,:] = woodylitter                                                     # Leaf and fine root litter, kg m-2, locate the intput to arrays
    air_ts = weather['T'].values                                               # Daily air temperatures, deg C
    tp_top_ts = df_peat_temperatures.iloc[:,2].values                          # Peat temperature -0.125 m depth
    tp_middle_ts = df_peat_temperatures.iloc[:,8].values                       # Peat temperature -0.4 m depth
    tp_bottom_ts = df_peat_temperatures.iloc[:,15].values                      # Peat temperature -0.75 m depth
    
    for n, (tair, tp_top, tp_middle, tp_bottom, wts) in enumerate(zip(air_ts, tp_top_ts, tp_middle_ts, tp_bottom_ts, water_tables.values)):    
        # Physical conditions in the peat profile
        wn =   wrc(self.pF[0], wts) / wrc(self.pF[0], -0.3)                    # Relative water content with respect to field capacity, pF[0] refers to water retention characterisitcs in the topmost layer
        H_w = wrc(self.pF[0], 0.0) - wrc(self.pF[0], wts)                      # Air filled pore space  
        peat_w1 = self.wtToVfAir_top(wts)                                      # Call interpolation function WT -> volume fraction of air 
        peat_w2 = self.wtToVfAir_middle(wts)                                   # Call interpolation function WT -> volume fraction of air
        peat_w3 = self.wtToVfAir_bottom(wts)                                   # Call interpolation function WT -> volume fraction of air   

        try:                                                
            k1, k2, k3, k4, k5, k6, k7, k8, k9 = self.get_rates(tair, tp_top, tp_middle,tp_bottom, wn,  peat_w1, peat_w2, peat_w3, H_w)
        except:
          print ('fail in rates, esom run_yr')

        if n==243:                                                             # n is day of the year
            self.M[:,:,0] = L0L                                                 # fresh litter leaves and fine roots, kg m-2, locate end of August
            self.M[:,:,1] = L0W                                                 # woody litter branches and coarse roots, kg m-2, locate end of August

        self.M = self.decompose(k1,  k2, k3, k4, k5, k6, k7, k8, k9, self.M)
        self.mass[:,:,:,self.i] = self.M                                       # locate mass to output array  
        self.i +=1                                                             # day counter
    end_i = self.i
    

    self.out = self.M[:,:, 10]*10000. - self.previous_mass
    self.previous_mass = self.M[:,:, 10]*10000.
    self.P1_out =  P1_ini - self.M[:,:, 7]*10000.
    self.P2_out =  P2_ini - self.M[:,:, 8]*10000.  
    self.P3_out =  P3_ini - self.M[:,:, 9]*10000.  
    self.out_root_lyr = self.out - self.P2_out - self.P3_out
    self.out_below_root_lyr =  self.P2_out + self.P3_out
    
    #return self.out
  
  def compose_export(self, stp):
        # To get total export, sum the left and right ditches 
        docshare = 0.05
        lmwtohmwshare = 0.04
        mass_to_c = 0.5
        hmw = (1-lmwtohmwshare)*self.out * 1/1.05 * docshare * mass_to_c  
        lmw = self.out * 1/1.05 * docshare * mass_to_c * lmwtohmwshare  
        self.hmw = hmw 
        self.lmw = lmw
        self.hmwtoditch = hmw*np.exp(-0.0004*stp.residence_time)                    # biodegradation parameters from Kalbiz et al 2003
        self.lmwtoditch= lmw*np.exp(-0.15*stp.residence_time)
        self.hmw_to_west = len(np.ravel(stp.ixwest))/stp.n * np.mean(self.hmwtoditch[0, np.ravel(stp.ixwest)])
        self.hmw_to_east = len(np.ravel(stp.ixeast))/stp.n * np.mean(self.hmwtoditch[0, np.ravel(stp.ixeast)])
        self.lmw_to_west = len(np.ravel(stp.ixwest))/stp.n * np.mean(self.lmwtoditch[0, np.ravel(stp.ixwest)])
        self.lmw_to_east = len(np.ravel(stp.ixeast))/stp.n * np.mean(self.lmwtoditch[0, np.ravel(stp.ixeast)])
        
        
        