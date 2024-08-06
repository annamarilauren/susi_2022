# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 09:49:39 2022

@author: alauren
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:52:35 2022

@author: alauren
"""

import numpy as np
import pandas as pd
import datetime
from susi.susi_utils import  get_motti, read_FMI_weather, get_mese_input
from inputs.susi_para import get_susi_para
from susi.susi_main import Susi
import susi.susi_io
from netCDF4 import Dataset 

#***************** local call for SUSI*****************************************************
folderName=r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'

wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'
wdata='katila_weather.csv'


mottifile = {'path':r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/',
              'dominant':{1: 'ullika_1.xlsx', 2: 'ullika_2.xlsx', 3: 'ullika_3.xlsx' },
              'subdominant':{0:'susi_motti_input_lyr_1.xlsx'},
              'under':{0:'susi_motti_input_lyr_2.xlsx'}} 


start_date = datetime.datetime(2005,1,1)
end_date=datetime.datetime(2014,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25

sarkaSim = 120.                                                                  # strip width, ie distance between ditches, m
n = int(sarkaSim / 2)                                                           # number of computation nodes in the strip

ageSim = {'dominant': 120.*np.ones(n),
          'subdominant': 0*np.ones(n),
          'under': 0*np.ones(n)}                                                         # age of the stand in each node

sfc =  np.ones(n, dtype=int)*5                                                                        # site fertility class

#ageSim['dominant'][int(n/2):] = 2.

#ageSim[4:-4] = 2.
#site = 'develop_scens'
site = 'ullika'


forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=None,
                                                                          n=n)
spara['canopylayers']['dominant'][int(n/3):int(n/3*2)] = 2
spara['canopylayers']['dominant'][int(n/3*2):] = 3



spara['bd top'] = [0.1, 0.1, 0.1, 0.12, 0.12, 0.12, 0.12, 0.14]
spara['bd bottom'] = [0.14]

outpara['netcdf'] = 'susi_ullikka_120m.nc'


susi = Susi() 
susi.run_susi(forc, wpara, cpara, org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=None)
    
          
#%%
from susi.figures import *
folderName =r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'
ff = folderName + 'susi_ullikka_100m.nc'
scen = 1
#hydrology(ff, scen)
#stand(ff, scen)
#mass(ff, scen)
carbon(ff, scen)
    
#scens = [0,1]
#compare_1(ff, scens)

#%%

folderName =r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'

ncf=Dataset(folderName+'susi_ullikka_120m.nc', mode='r')
for scen in [0,1,2]:

    roff = ncf['strip']['roff'][scen, :]*1000
    sflow = np.mean(ncf['strip']['surfacerunoff'][scen, :, :], axis=1)*1000
    totr = roff + sflow
    nofdays = len(roff)
    docexport = ncf['export']['hmwtoditch'][scen,:, :] + ncf['export']['lmwtoditch'][scen,:, :]
    
    
    import pandas as pd
    dates = pd.date_range(str(start_date.date()),periods=nofdays)
    dfrunoff = pd.DataFrame(data={'date':dates, 'runoff':totr})
    dfrunoff = dfrunoff.set_index('date')
    water = dfrunoff['runoff'].resample('y').sum()
    docs = np.mean(docexport, axis=1)[1:]
    conc = (docs*1e6 )/( water.values*1e4)
    print ('scenario',scen)
    print ('exp',np.mean(docs), np.std(docs))
    print ('conc', np.mean(conc), np.std(conc))
ncf.close()

#%%

folderName =r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'

ncf=Dataset(folderName+'susi_ullikka_80m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print ('-----80m----------')
print ('exports',hmwtoditch0, hmwtoditch1,hmwtoditch2)
wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)

ncf.close()                                 

ncf=Dataset(folderName+'susi_ullikka_100m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print ('-----100m----------')
print ('exports', hmwtoditch0, hmwtoditch1,hmwtoditch2)
wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)
ncf.close()                                 


ncf=Dataset(folderName+'susi_ullikka_120m.nc', mode='r')
hmwtoditch0 = np.mean(ncf['export']['hmwtoditch'][0,:, :] + ncf['export']['lmwtoditch'][0,:, :])
hmwtoditch1 = np.mean(ncf['export']['hmwtoditch'][1,:, :] + ncf['export']['lmwtoditch'][1,:, :])
hmwtoditch2 = np.mean(ncf['export']['hmwtoditch'][2,:, :] + ncf['export']['lmwtoditch'][2,:, :])
print ('-----120m----------')
print ('exports',hmwtoditch0, hmwtoditch1,hmwtoditch2)
wt0 =  np.mean(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1 =  np.mean(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2 =  np.mean(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

wt0sd =  np.std(ncf['strip']['dwtyr_growingseason'][0,:, :]) 
wt1sd =  np.std(ncf['strip']['dwtyr_growingseason'][1,:, :]) 
wt2sd =  np.std(ncf['strip']['dwtyr_growingseason'][2,:, :]) 

print ('wts', wt0,wt1,wt2)
print ('wtssd', wt0sd,wt1sd,wt2sd)


ncf.close()                                 

#%%
wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ullikka/'
wdata='katila_weather.csv'
start_date = datetime.datetime(2005,1,1)
end_date=datetime.datetime(2014,12,31)

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)    
print (np.mean(forc['T'].resample('Y').mean()))
print (np.std(forc['T'].resample('Y').mean()))

#print(forc['Prec'].resample('Y').sum())

print (np.mean(forc['Prec'].resample('Y').sum()))
print (np.std(forc['Prec'].resample('Y').sum()))
