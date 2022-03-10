# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 19:31:48 2020

@author: alauren
"""
import numpy as np
import matplotlib.pylab as plt
import pandas as pd

class DocModel():
    def __init__(self, lyrs, cols, bd_top, bd_bottom, dz, length):
        """
        construct cross section of a strip (cols, lyrs)
        """

        self.wet_limit = 0.3                                                   # m
        self.gamma_dry = 1.331                                                 # micrograms / gram / day = mg/kg/day
        self.gamma_wet = 2.814                                                 # DOC

        self.gamma_dry_hmw = 1.331-0.561                                       # micrograms / gram / day = mg/kg/day
        self.gamma_wet_hmw = 2.814-0.272                                       # DOC
        self.DOC_Q10 = 1.98                                                    # 2.52

        self.length=length
        self.cols = cols                                                       # columns in strip
        self.lyrs = lyrs
        self.z = -np.cumsum(np.ones(self.lyrs)*dz)                             # depth in column, one column
        self.zmap = np.tile(self.z, (self.cols,1))                             # depth map (col, lyr)

        rho = np.ones(self.lyrs)*bd_bottom*1000.                               # bulk density convert from g/cm3 to kg m-3, one column
        if bd_top is not None:
            rho[:len(bd_top)] = np.array(bd_top) * 1000.                           # separate top layer bulk densities convert from g/cm3 to kg m-3
            
        dry_mass = rho*dz                                                      # organic matter kg in layer
        self.mass_map = np.tile(dry_mass, (self.cols,1))                       # layer dry mass as map (cols, lyrs) 
        self.releasedDOC = np.zeros((cols, lyrs))                              # cumulative release of total DOC
        self.releasedHMW = np.zeros((cols, lyrs))                              # cumulative release of HMW DOC

    def reset_domain(self):
        self.releasedDOC = np.zeros((self.cols, self.lyrs))                              # cumulative release of total DOC
        self.releasedHMW = np.zeros((self.cols, self.lyrs))                              # cumulative release of HMW DOC
        self.HMW_ts = np.zeros(self.length)    
        self.LMW_ts = np.zeros(self.length)
        self.i = 0
        
    def doc_release(self,temp, dfwt):
        days, _ = np.shape(temp)                                               # dimensions days and layers in column

        q10 = lambda t: self.DOC_Q10**((t-10.)/10.)
    
        #changing values wt and temp
        wt = dfwt.to_numpy()
        temp[temp<0.0]=-40.   #restrict when frozen
        t = temp.to_numpy()
        
        doc = np.zeros(self.cols)
        hmw = np.zeros(self.cols)
        
        for d in range(days):
            wt_map = np.tile(wt[d],(self.lyrs,1)).T                            # water table map in shape (cols, lyrs)
            t_map = np.tile(t[d], (self.cols,1))                               # temperature map in shape (cols, lyrs)

            above_mask = np.ones((self.cols, self.lyrs))
            above_mask = np.where(self.zmap > wt_map, above_mask, 0.0)    
            
            dry_mask = np.ones((self.cols, self.lyrs))
            dry_mask = np.where(self.zmap  > wt_map + self.wet_limit, dry_mask, 0.0)    
            
            wet_mask = above_mask-dry_mask
            
            q10_map = q10(t_map)
            doc_dry = self.gamma_dry * q10_map * self.mass_map * dry_mask * 10000. / 1000000.        #kg/ha/day 
            doc_wet = self.gamma_wet * q10_map * self.mass_map * wet_mask * 10000. / 1000000.        #kg/ha/day 
            hmw_dry = self.gamma_dry_hmw * q10_map * self.mass_map * dry_mask * 10000. / 1000000. 
            hmw_wet = self.gamma_wet_hmw * q10_map * self.mass_map * wet_mask * 10000. / 1000000. 
            
            self.releasedDOC = self.releasedDOC + doc_dry + doc_wet
            self.releasedHMW = self.releasedHMW + hmw_dry + hmw_wet
            
            doc = doc + np.sum(doc_dry, axis=1)
            doc = doc + np.sum(doc_wet, axis=1)
            hmw = hmw + np.sum(hmw_dry, axis=1)
            hmw = hmw + np.sum(hmw_wet, axis=1)
            
            self.HMW_ts[self.i] = np.sum(hmw_dry, axis=(0,1)) + np.sum(hmw_wet, axis=(0,1)) 
            self.LMW_ts[self.i] = np.sum(doc_dry, axis=(0,1)) + np.sum(doc_wet, axis=(0,1)) - self.HMW_ts[self.i]
            self.i += 1
            
        return doc, hmw
 
    def draw_doc(self,sitename, yrs):
        fig = plt.figure(num = 'doc release '+sitename, figsize =(15,9) )
        plt.subplot(311)
        plt.imshow((self.releasedHMW.T)/yrs, cmap='Reds', aspect='auto'); plt.colorbar()  
        plt.title('HMW DOC, $ kg ha^{-1} yr^{-1}$ ')
        plt.subplot(312)
        plt.imshow((self.releasedDOC.T - self.releasedHMW.T)/yrs, cmap='Purples', aspect='auto'); plt.colorbar()     
        plt.title('HMW DOC, $ kg ha^{-1} yr^{-1}$ ')
        plt.subplot(313)
        plt.plot(self.HMW_ts, 'r-', label='HMW')
        plt.plot(self.LMW_ts, 'b-', label='LMW')
        plt.legend(loc='upper right')
        outOpt=False
        if outOpt:
            folder = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Susi_doc/'
            plt.savefig(folder + sitename +'.png')
            dout = {'hmw_ts':self.HMW_ts, 'lmw_ts': self.LMW_ts}
            dfout = pd.DataFrame(dout)
            dfout.to_excel(folder + 'doc_ts.xlsx')
        
