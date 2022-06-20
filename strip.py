# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 19:59:31 2019

@author: lauren
"""
import numpy as np
from susi_utils import peat_hydrol_properties, CWTr

#from susi_utils import gmeanTr, Hadjacent, Amatrix, boundConst, rightSide


class StripHydrology():
    def __init__(self, spara):
        self.nLyrs = spara['nLyrs']                                                 # number of soil layers
        dz = np.ones(self.nLyrs)*spara['dzLyr']                                     # thickness of layers, m
        z = np.cumsum(dz)-dz/2.                                                # depth of the layer center point, m 
        self.spara = spara
        if spara['vonP']:
            lenvp=len(spara['vonP top'])    
            vonP = np.ones(self.nLyrs)*spara['vonP bottom'] 
            vonP[0:lenvp] = spara['vonP top']                                      # degree of  decomposition, von Post scale
            ptype = spara['peat type bottom']*spara['nLyrs']
            lenpt = len(spara['peat type']); ptype[0:lenpt] = spara['peat type']    
            self.pF, self.Ksat = peat_hydrol_properties(vonP, var='H', ptype=ptype) # peat hydraulic properties after Päivänen 1973    
        else:
            lenbd=len(spara['bd top'])    
            bd = np.ones(self.nLyrs)*spara['bd bottom'] 
            bd[0:lenbd] = spara['bd top']                                      # degree of  decomposition, von Post scale
            ptype = spara['peat type bottom']*spara['nLyrs']
            lenpt = len(spara['peat type']); ptype[0:lenpt] = spara['peat type']    
            self.pF, self.Ksat = peat_hydrol_properties(bd, var='bd', ptype=ptype) # peat hydraulic properties after Päivänen 1973    
        
        for n in range(self.nLyrs):
            if z[n] < 0.41: 
                self.Ksat[n]= self.Ksat[n]*spara['anisotropy']
            else:
                self.Ksat[n]= self.Ksat[n]*1.
                
        self.dwtToSto, self.stoToGwl, self.dwtToTra, self.C, self.dwtToRat, self.dwtToAfp = CWTr(self.nLyrs, z, dz, self.pF, 
                                                               self.Ksat, direction='negative') # interpolated storage, transmissivity, diff water capacity, and ratio between aifilled porosoty in rooting zone to total airf porosity  functions
        
    
        self.L= spara['L']                                                     # compartemnt width, m
        self.n= spara['n']                                                     # number of computation nodes
        self.dy = float(self.L/self.n)                                         # node width m
        sl= spara['slope']                                                     # slope %
        lev = 1.                                                               # basic level of soil surface
        
        self.ele = np.linspace(0,self.L*sl/100., self.n) + lev                 # surface rise in y direction, m
        self.dt = 1                                                             # time step, days
        self.implic = 1.# 0.5                                                  # 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson
        self.DrIrr = False
        self.dwt = spara['initial h']                                          # h in the compartment
        self.H = (self.ele + self.dwt)  #np.empty(self.n)                      # depth to water table, negative down m

        self.Kmap = np.tile(self.Ksat, (self.n,1))                             # Ksat map (col, lyr)
        self.residence_time = np.zeros(self.n)                                 # residence time from column to the ditch, days
        
        print ('Peat strip initialized')
        
    def reset_domain(self):
        self.A = np.zeros((self.n,self.n))                                     # computation matrix 
        self.dwt = np.ones(self.n)*self.spara['initial h']                                     # right hand side vector
        self.H = self.ele + self.dwt                                             # head with respect to absolute reference level, m
        #self.sruno = 0.
        self.roff = 0.
        print ('Resetting strip scenario')


    def run_timestep(self,d, h0ts_west, h0ts_east, p, moss):
        """
        IN: 
            d day number
            h0ts boudary (ditch depth, m) in time series
            p rainfall-et m, arrayn n length
            moss as object
        """
        n = self.n        
        Htmp = self.H.copy(); Htmp1=self.H.copy()
        self.dwt = (Htmp - self.ele)
        #S = p/1000.*np.ones(n)                                                # source/sink, in m
        S = p.copy() #*np.ones(n)                                              # source/sink, in m
        self.dwt[0] = h0ts_west
        self.dwt[n-1] = h0ts_east                                              # symmetrical boundaries, set water level
        
        #TESTING here
        #airv = self.hToSto(self.ele)-self.hToSto(Htmp-self.ele)                # air volume, in m
        airv = np.maximum(self.dwtToSto(np.zeros(n)) - self.dwtToSto(self.dwt), np.zeros(n))
        
        S = np.where(S > airv, airv, S)                        
        exfil = p - S
        
        self.surface_runoff = moss.returnflow(exfil)

        #self.sruno += self.surface_runoff
        #self.sruno += np.where(np.ones(len(airv))*(p)/1000. > airv, np.ones(len(airv))*(p)/1000.-airv, 0.0)  #cut the surface water above to runoff
        Tr0 = self.dwtToTra(self.dwt)                                     # Transmissivity from the previous time step 
        
        
        Trminus0, Trplus0 = self.gmeanTr(Tr0)                                  # geometric mean of adjacent node transmissivities
        Hminus, Hplus = self.Hadjacent(self.H)                                 # vector of adjacent node H

        for it in range(100):                                                  # iteration loop for implicit solution
            Tr1 = self.dwtToTra(Htmp1-self.ele)                                  # transmissivity in new iteration
            Tr0 = np.maximum(self.dwtToTra(self.H-self.ele),0.0)               
            Tr1 = np.maximum(self.dwtToTra(Htmp1-self.ele),0.0)  
            CC = self.C(Htmp1-self.ele)                                        # storage coefficient in new iteration
            Trminus1, Trplus1 = self.gmeanTr(Tr1)                              # geometric mean of adjacent node transmissivity                
            alfa = CC*self.dy**2/self.dt
            self.A = self.Amatrix(self.A, n, self.implic, Trminus1, Trplus1, alfa)              # construct tridiaginal matrix
            self.A = self.boundConst(self.A, n)                                # constant head boundaries to A matrix
            hs=self.rightSide(S,self.dt,self.dy,self.implic,alfa, self.H,Trminus0, Hminus, Trplus0, Hplus,\
                self.DrIrr, Htmp1, self.ele, h0ts_west, h0ts_east)                             # right hand side of the equation
            Htmp1 = np.linalg.multi_dot([np.linalg.inv(self.A),hs])            # solve equation                    
            #Htmp1 = np.matmul(np.linalg.inv(self.A),hs)            # solve equation                    

            Htmp1=np.where(Htmp1>self.ele, self.ele,Htmp1)                     # cut the surface water 
            conv = max(np.abs(Htmp1-Htmp))                                     # define convergence
            Htmp=Htmp1.copy()                                                  # new wt to old for new iteration
            if conv < 1.e-7: 
                if d%365==0: print ('  - day #',d, 'iterations', it)                    
                break
        self.H=Htmp1.copy()                     
        
        #**********************construction*****************
        self.roffwest, self.roffeast = self.runoff(self.H, Trminus1, Trplus1,\
                        self.dt, self.dy, self.L)        
        self.surface_runoff
        self.roff = self.roffwest + self.roffeast + np.mean(self.surface_runoff)    
        
        self.dwt= self.H-self.ele
        self.air_ratio = self.dwtToRat(self.dwt)
        self.afp = self.dwtToAfp(self.dwt)
        #**************************************************
        
        return self.dwt, self.H, self.roff, self.air_ratio, self.afp

    def Hadjacent(self,H):
        """
        Input:
            H vector, H in each node
        Output:
            Hwest H(i-1), Heast H(i+1)
        """
        n=len(H)
        Hwest = H[0:n-1]; Hwest=np.append(Hwest, 0.0)
        Heast = H[1:]; Heast=np.insert(Heast, 0, 0.0)
        return Hwest, Heast  

    def Amatrix(self, A, n, implic, Trwest, Treast, alfa):
        """
        Construction of tridiagonal matrix
        """     
        i,j = np.indices(A.shape)
        A[i==j]= implic*(Trwest+Treast)+alfa                                   # diagonal element
        A[i==j+1]=-implic*Trwest[:n-1]                                         # West element
        A[i==j-1]=-implic*Treast[1:]                                           # East element     
    
        return A
    
    def boundConst(self, A, n):
        """
        Diriclet (constant head boundary conditions)
        """    
        A[0,0]=1; A[0,1]=0.                                                    # Dirichlet, west boundary
        A[n-1,n-1]=1.; A[n-1, n-2]=0.                                          # Dirichlet, east boundary
        return A
    
    def boundNoFlow(A, n, implic, Trwest, Treast, alfa):
        """
        Diriclet (constant head boundary conditions)
        """    
                                                                                
        A[0,0]= 2.*implic*(Treast[0])+alfa[0]                                  # Diagonal element
        A[0,1]=-2*implic*Trwest[0]                                             # East element     
        A[n-1,n-1]=2.*implic*(Trwest[n-1])+alfa[0]
        A[n-1, n-2]=-2*implic*Treast[n-1]                                      # West element
    
        return A
    
    
    def rightSide(self, S,dt,dy, implic,alfa, H, Trminus0,Hminus,Trplus0,Hplus, DrIrr, Htmp1, ele, h0_west, h0_east):
        
        hs = S*dt*dy**2 + alfa*H + (1-implic)*(Trminus0*Hminus) - (1-implic)*(Trminus0 + Trplus0)*H  + (1-implic)*(Trplus0*Hplus)
        n=len(Htmp1)
                
        if DrIrr==False:                
            hs[0]=Htmp1[1] if Htmp1[0]>Htmp1[1] else min(ele[0]+h0_west, Htmp1[1])
            hs[n-1]=Htmp1[n-2] if Htmp1[n-1]>Htmp1[n-2] else min(ele[n-1]+h0_east, Htmp1[n-2])    #if wt below canal water level, lower the canal wl to prevent water inflow to compartment
        else:
            hs[0]=ele[0]+h0_west
            hs[n-1]=ele[n-1]+h0_east 
        return hs
    
    
    def gmeanTr(self, Tr):
        """
        Input: 
            Transmissivity vector, tr in node center point
        Output:
            Transmissivity, tr in west surface sqrt(Tr(i-1)*Tr(i)) and east sqrt(Tr(i)*Tr(i+1)) 
        """
        n = len(Tr)
        trwest = np.maximum(Tr[:n-1]*Tr[1:], 0.0)
        #Trwest = np.sqrt(Tr[:n-1]*Tr[1:])
        Trwest = np.sqrt(trwest); Trwest=np.append(Trwest, 0.0)
        treast = np.maximum(Tr[1:]*Tr[:n-1], 0.0)
        #Treast = np.sqrt(Tr[1:]*Tr[:n-1])
        Treast = np.sqrt(treast)   
        Treast=np.insert(Treast, 0, 0.0)                
        return Trwest, Treast

    def runoff(self, H, Trminus, Trplus, dt, dy,L):
        roffwest = ((H[1]-H[0])/dy*Trminus[0]*dt)/L
        roffeast = ((H[-2]-H[-1])/dy*Trplus[-1]*dt/L)        
        return roffwest, roffeast
    

    def create_outarrays(self, nrounds, ndays, ncols):
        stpout ={}
        stpout['dwts'] = np.zeros((nrounds, ndays,ncols), dtype=float)                   # water table depths, m,  ndarray(scenarios, days, number of nodes)
        stpout['afps'] = np.zeros((nrounds, ndays,ncols), dtype=float)                   # air-filled porosity (m3 m-3),  ndarray(scenarios, days, number of nodes)
        stpout['deltas'] = np.zeros((nrounds, ndays,ncols), dtype=float)
        stpout['hts'] = np.zeros((nrounds, ndays,ncols), dtype=float)                    # water table depths, m,  ndarray(scenarios, days, number of nodes)
        stpout['runoff'] = np.zeros((nrounds,ndays), dtype=float)                        # daily total runoff, here in m, sum of west, east, surface runoff
        stpout['runoffwest'] = np.zeros((nrounds,ndays), dtype=float)                    # daily runoff, here in m, from west ditch
        stpout['runoffeast'] = np.zeros((nrounds,ndays), dtype=float)                    # daily runoff, here in m, from east ditch
        stpout['surfacerunoff'] = np.zeros((nrounds,ndays, ncols), dtype=float)                 # daily surfacerunoff, here in m, from each column

        return stpout

    def update_outarrays(self, r, d, stpout):
        stpout['dwts'][r,d,:] = self.dwt                                                 # daily water tables
        stpout['hts'][r,d,:] = self.H                                                    # water table height in comparison to stabile datum (elevation)
        stpout['afps'][r,d,:] = self.afp                                                 # air filled porosity
        stpout['runoff'][r,d] = self.roff                                                # daily runoff 
        stpout['runoffwest'][r,d] = self.roffwest
        stpout['runoffeast'][r,d] = self.roffeast
        stpout['surfacerunoff'][r,d,:]  =  self.surface_runoff                                                      # daily surfacerunoff, here in m, from each column

        return stpout
    
    def update_residence_time(self, dfwt):
       
       timetoditch = np.zeros(self.n)                                          # residence time from column to ditch, days
       porosity = 0.9
       K = 10**(-4)*86400                                                      # generic Koivusalo et al, 2008
       
       dist = np.arange(0,self.n*self.dy,self.dy)                              # distance array 
       H = self.ele + dfwt.mean(axis = 0)                                                       # water table height to common datum 
       rtime =self.dy / (K * np.gradient(H, dist)/porosity)                    # residence time within a column
       self.ixwest = np.where(rtime > 0)                                       # separate with directions, west, east
       self.ixeast = np.where(rtime < 0)                                 
    
       timetoditch[self.ixwest] = np.cumsum(rtime[self.ixwest])                       
       timetoditch[self.ixeast] = np.flip(np.cumsum(np.flip(rtime[self.ixeast]*-1)))
   
       self.residence_time = timetoditch
       
def drain_depth_development(length, hdr, hdr20y):
    """
    Computes daily level of drain bottom thru the time of the simulation. Model adjusted from Hannu Hökkä drain model.
    Input:
        - drain depth in the beginning of the simulation (m, negative down)
        - drain depth after 20 yrs (m, negative down)
        - length of simulation in days
    Output:
        - daily drain bottom level (m, negative down) 
    """    
    timeyrs = np.linspace(0,length/365., length)                                     # time vector telling the time elapsed from the simulation beginning time, yrs        
    h0ts = ((-100*hdr20y +100.*hdr) /np.log(20.0)*np.log(timeyrs+1.)-100*hdr)/-100.  # Ditch model Hökkä
    return h0ts
  