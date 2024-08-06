# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import datetime
from susi.susi_utils import read_FMI_weather
from inputs.susi_para import get_susi_para
from susi.susi_main import Susi


#***************** local call for SUSI*****************************************************
folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/' #'sensitivity/'
wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/inputs/'


mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/codes/Susi_10/inputs/',
              'dominant':{1: 'motti viitasaari_mtkg.xls'},
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 
"""
mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/krycklan/',
              'dominant':{1: 'susi_krycklan_input_lyr_0.xlsx'},
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 
"""
wdata='parkano_weather.csv'

start_date = datetime.datetime(1995,1,1)
end_date=datetime.datetime(2014,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 40.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 40.*np.ones(n),
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                                         # age of the stand in each node

sfc =  np.ones(n, dtype=int)*4                                                                        # site fertility class

#ageSim['dominant'][int(n/2):] = 2.
#ageSim[4:-4] = 2.

site = 'develop_scens'

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          n=n)
#spara['canopylayers']['dominant'][int(n/2):] = 2                                                                        
#spara['canopylayers']['subdominant'][:int(n/2)] = 1                                                                        

spara['ditch depth west'] = [-0.5]   #nLyrs kerrosten lkm, dzLyr kerroksen paksuus m, saran levys m, n laskentasolmulen lukumäärä, ditch depth pjan syvyys simuloinnin alussa m  
spara['ditch depth east'] = [-0.5]
spara['ditch depth 20y west'] = [-0.5]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
spara['ditch depth 20y east'] = [-0.5]                                            #ojan syvyys 20 vuotta simuloinnin aloituksesta
spara['scenario name'] =  ['D60']                                #kasvunlisaykset


susi = Susi()
 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc)
    
          
             
#%%
from susi.figures import *
ff = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/susi.nc'
#ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/vesiensuojelu_2023/susi_base.nc'
#ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/vesiensuojelu_2023/susi_ash.nc'
#ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/vesiensuojelu_2023/susi_partial_d.nc'
#ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/vesiensuojelu_2023/susi_partial_d_ash.nc'


scen = 2
# hydrology(ff, scen)
#stand(ff, scen)
# mass(ff, scen)
# carbon(ff, scen)
# nutrient_balance(ff, 'N', scen)
# nutrient_balance(ff, 'P', scen)
# nutrient_balance(ff, 'K', scen)
# compare_1(ff, [0,3])
#compare_scens(ff)
#%%
from netCDF4 import Dataset  
scen = 0
ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    
#-------------- Water tables cross section-------------------------------------
facecolor = '#f2f5eb'
fs = 15
fig = plt.figure(num='susi', figsize=(15,18), facecolor=facecolor)   #width, height
gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=2.2, hspace=0.6)


out = np.array(ncf['stand']['volume'][0,:,1:-1]) 
ax = fig.add_subplot(gs[:4, :4])
ax.plot(out)
ax.set_ylabel('volume', fontsize=fs)

out = np.mean(np.array(ncf['strip']['dwtyr'][0,1:,:]), axis=1)
ax = fig.add_subplot(gs[:4, 4:8])
ax.plot(out)
ax.set_ylabel('WT', fontsize=fs)

out = np.mean(np.array(ncf['strip']['dwtyr_latesummer'][0,1:,:]), axis=1)
ax = fig.add_subplot(gs[:4, 8:])
ax.plot(out)
ax.set_ylabel('WT latesummer', fontsize=fs)

out = np.array(ncf['balance']['C']['stand_c_balance_c'][0,1:,:])
ax = fig.add_subplot(gs[4:8, :4])
ax.plot(out)
#ax.plot(np.cumsum(out))
ax.set_ylabel('stand c bal', fontsize=fs)

out = np.array(ncf['balance']['C']['soil_c_balance_c'][0,1:,:])
ax = fig.add_subplot(gs[4:8, 4:8])
ax.plot(out)
ax.set_ylabel('stand c bal', fontsize=fs)

out = np.array(ncf['balance']['C']['stand_litter_in'][0,1:,:] + ncf['balance']['C']['gv_litter_in'][0,1:,:])
ax = fig.add_subplot(gs[4:8, 8:])
ax.plot(out)
ax.set_ylabel('litter in', fontsize=fs)


out = np.array(ncf['balance']['C']['co2c_release'][0,1:,:])
ax = fig.add_subplot(gs[8:12, :4])
ax.plot(out)
ax.set_ylabel('co2 release', fontsize=fs)


substance = 'Mass'
mor =  ncf['esom'][substance]['LL'][scen, :, :] \
    + ncf['esom'][substance]['LW'][scen, :, :] \
    + ncf['esom'][substance]['FL'][scen, :, :]\
    + ncf['esom'][substance]['FW'][scen, :, :]\
    + ncf['esom'][substance]['H'][scen, :, :]

#mor = np.sum(mor, axis =1)
ax = fig.add_subplot(gs[8:12, 4:8])
ax.plot(mor)
ax.set_ylabel('mor', fontsize=fs)


peat = ncf['esom'][substance]['P1'][scen, :, :] \
        + ncf['esom'][substance]['P2'][scen, :, :] \
        + ncf['esom'][substance]['P3'][scen, :, :]

#peat = np.sum(peat, axis=1)
ax = fig.add_subplot(gs[8:12, 8:])
ax.plot(peat+mor)
ax.set_ylabel('peat+mor', fontsize=fs)

out = np.array(ncf['balance']['C']['co2c_release'][0,1:,:])
ax = fig.add_subplot(gs[8:12, :4])
ax.plot(out)
ax.set_ylabel('co2 release', fontsize=fs)

#plt.plot(mor+peat)
#plt.plot(peat)
#print(np.shape(out))
#out = np.mean(out, axis=1)

#cols =  (np.shape(out)[-1])
#scenario = 0            # 0 = 30 cm syvä oja, 1 = 60 cm syvä oja ja 2 = 90 cm syvä oja
#df = pd.DataFrame(data=out[:,:], columns=list(range(cols)))
#df = pd.DataFrame(data=out)

ncf.close()
#plt.plot(np.cumsum(out))
#plt.plot(out)
#print (df)
#df.plot()
#%%
"""
from netCDF4 import Dataset  
ff = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/susi.nc'

ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode

wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)
ncf.close()                                 
"""
