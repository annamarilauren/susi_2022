# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 08:59:36 2019

@author: kersti haahti
"""
import numpy as np
eps = np.finfo(float).eps

class MossLayer(object):
    """
    Moss/organic layer on top of soil.
    """
    def __init__(self, org_para, outputs=False):
        """
        Initializes MossLayer:
        Args:
            org_para (dict):
                # parameters
                'org_depth': depth of organic top layer (m)
                'org_poros': porosity (-)
                'org_fc': field capacity (-)
                'org_rw': critical vol. moisture content (-) for decreasing phase in Ef
                'pond_storage_max': maximum pond depth [m]
                # initial states
                'org_sat': organic top layer saturation ratio (-)
                'pond_storage': initial pond depth at surface [m]
        """
        # top layer is interception storage, which capacity is depends on its depth [m]
        # and field capacity
        self.dz_top = org_para['org_depth']  # depth, m3 m-3
        self.poros_top = org_para['org_poros']  # porosity, m3 m-3
        self.fc_top = org_para['org_fc']  # field capacity m3 m-3
        self.rw_top = org_para['org_rw']  # ree parameter m3 m-3
        self.Wsto_top_max = self.fc_top * self.dz_top  # maximum storage m

        # pond storage
        self.h_pond_max = org_para['pond_storage_max']

        # initialize state
        # toplayer storage and relative conductance for evaporation
        self.Wsto_top = self.Wsto_top_max * org_para['org_sat']
        self.Wliq_top = self.poros_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(
                0.98*self.Wliq_top / self.rw_top, 1.0)) # relative evaporation rate (-)
        # pond storage
        self.h_pond = org_para['pond_storage']

        # create dictionary of empty lists for saving results
        if outputs:
            self.results = {'Mbe': [], 'PondSto': [], 'Wliq_top': [], 'Ree': [], 'SRunoff': []}

    def interception(self, potinf=0.0, evap=0.0):
        """ 
        Solves top layer interception
        Args:
            potinf (float): throughfall [m]
            evap (float): evaporation from top layer [m]
        Returns:
            potinf (float): potential infiltration [m]
            evap (float): evaporation from top layer [m]
        """
        potinf0 = potinf.copy()

        # initial state
        Wsto_top_ini = self.Wsto_top.copy()
        pond_ini = self.h_pond.copy()

        # add current pond storage to rr & update storage
        potinf += pond_ini
        self.h_pond -= pond_ini

        # top layer interception & water balance
        interc = np.maximum(0.0, (self.Wsto_top_max - self.Wsto_top))\
                    * (1.0 - np.exp(-(potinf / self.Wsto_top_max)))
        potinf -= interc  # to soil profile
        self.Wsto_top += interc
        evap = np.minimum(evap, self.Wsto_top)
        self.Wsto_top -= evap

        self.mbe = ((pond_ini - self.h_pond) + (Wsto_top_ini - self.Wsto_top) +
                    potinf0 - potinf - evap)

        return potinf, evap, self.mbe

    def returnflow(self, rflow):
        """ 
        Water not fitting in soil first fills water storage in moss layer, then pond, and
        rest is routed to surface runoff. States updated.
        Args:
            rflow (float): water that cannot fit in soil [m]
        Returns:
            surface_runoff (float): surface runoff [m]
        """
        # first fill top layer
        # water that can fit in top layer
        to_top_layer = np.minimum(rflow, self.Wsto_top_max - self.Wsto_top)
        self.Wsto_top += to_top_layer
        # then pond storage
        to_pond = np.minimum(rflow - to_top_layer, self.h_pond_max - self.h_pond)
        self.h_pond += to_pond

        #if max(to_pond) > 0.0:
        #    print('pond', max(self.h_pond), max(to_pond))
        
        # and route remaining to surface runoff
        surface_runoff = rflow - to_top_layer - to_pond

        self.mbe += rflow - to_top_layer - to_pond - surface_runoff

        # update moss layer state
        self.Wliq_top = self.fc_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(0.98*self.Wliq_top / self.rw_top, 1.0))

        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['Mbe'].append(self.mbe)
            self.results['PondSto'].append(self.h_pond.copy())
            self.results['Wliq_top'].append(self.Wliq_top)
            self.results['Ree'].append(self.Ree)
            self.results['SRunoff'].append(surface_runoff)

        return surface_runoff
