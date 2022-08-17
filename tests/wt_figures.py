# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 18:04:05 2020

@author: alauren
"""
import datetime
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from scipy import stats
from itertools import product
from dwts_para import para
from netCDF4 import Dataset  

#from sklearn.metrics import r2_score
import seaborn as sns
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
sns.set()


folder_meas = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Stenberg/puusto_ym_data_Arille/Pohjavesiaineistot/'
folder_sim = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/'


params = para(period='start-end')
print (params.keys())
s_ansa = 40.
s_koira = 37.
s_neva = 30.
s_jaakkoin = 40.
wt_meas ={
          'ansa21':{'file':'muhos_2_pohjavesi_koottu.xlsx',
                    'tubes':[ 1, 2, 3, 4, 5, 6],
                    'dist':[5.,s_ansa/2.,5.,5.,s_ansa/2.,5.]},
          'ansa26':{'file':'muhos_2_pohjavesi_koottu.xlsx',
                    'tubes':[28, 29, 30],
                    'dist':[5., s_ansa/2., 5.]},
                    
          'koira11':{'file':'koiraoja_pohjavesi_koottu.xlsx',
                    'tubes':[ 13, 14, 15],
                    'dist':[5., s_koira/2., 5.]},
          'koira12':{'file':'koiraoja_pohjavesi_koottu.xlsx',
                    'tubes':[1, 2, 3],
                    'dist':[5., s_koira/2., 5.]},
          'koira21':{'file':'koiraoja_pohjavesi_koottu.xlsx',
                    'tubes':[19, 20, 21],
                    'dist':[5., s_koira/2., 5.]},
          'koira22':{'file':'koiraoja_pohjavesi_koottu.xlsx',
                    'tubes':[31, 32, 33],
                    'dist':[5., s_koira/2., 5.]},
                    
          'neva11':{'file':'nevajarvi_1_pohjavesi_koottu.xlsx',
                    'tubes':[ 5, 6, 7, 8, 9, 10],
                    'dist':[5.,5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.,]},
          'neva14':{'file':'nevajarvi_1_pohjavesi_koottu.xlsx',
                    'tubes':[36, 41, 42, 43, 44, 45],
                    'dist':[5.,5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.,]},
          'neva21':{'file':'nevajarvi_2_pohjavesi_koottu.xlsx',
                    'tubes':[5, 6, 7, 8, 9, 10, 11],
                    'dist':[5.,5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.,]},
          'neva24':{'file':'nevajarvi_2_pohjavesi_koottu.xlsx',
                    'tubes':[36, 41, 42, 43, 44, 45],
                    'dist':[5.,5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.,]},
          'neva31':{'file':'nevajarvi_3_pohjavesi_koottu.xlsx',
                    'tubes':[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    'dist':[5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.,5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.]},
          'neva34':{'file':'nevajarvi_3_pohjavesi_koottu.xlsx',
                    'tubes':[ 41, 42, 43, 44, 45, 46, 47, 48, 49, 50],
                    'dist':[5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.,5.,s_neva*0.25, s_neva*0.5, s_neva*0.25, 5.]},

          'jaakkoin61':{'file':'jaakkoinsuo_pohjavesi_koottu.xlsx',
                    'tubes':[ 2, 3, 12, 13, 14, 15],
                    'dist':[5., 5., 16.,16., 17.5,27.]},
          'jaakkoin62':{'file':'jaakkoinsuo_pohjavesi_koottu.xlsx',
                    'tubes':[ 28, 30, 36, 39, 42, 43, 44, 45],
                    'dist':[32., 26., 21., 13.,15., 23., 37.5, 46. ]},
            
          'parkano11':{'file':'parkano_1_pohjavesi_koottu.xlsx',
                    'tubes':[ 22,23,24,27,28,29,32,33,34,37,38,39]},
          'parkano12':{'file':'parkano_3_pohjavesi_koottu.xlsx',
                    'tubes':[1,2,9,10,11,12,19,20,21,22, 29,30,31,32]},
        }


#%%
# Measured and modelled WT as time series
""" WT time series figure  """
cols = 2
rows = 6
#coordinates = list(product(range(rows), range(cols)))


#fig, axs = plt.subplots(rows, cols, figsize = (10.,15.))
fig = plt.figure(constrained_layout=True, figsize = (10.,15.))
gs = GridSpec(rows, cols, figure=fig)
coordinates = [gs[0,0], gs[0,1],
               gs[1,0], gs[1,1],
               gs[2,0], gs[2,1],
               gs[3,0], gs[3,1],
               gs[4,0], gs[4,1],
               gs[5,0]
               ]

sites =['koira11', 'koira12',
        'ansa21', 'ansa26', 
        'neva11', 'neva14', 
        'neva31', 'neva34',
        'jaakkoin61', 'jaakkoin62', 
        'parkano11']

names = ['Koirasuo11','Koirasuo12',
         'Ansasaari21','Ansasaari26', 
         'Nevajärvi11','Nevajärvi14', 
         'Nevajärvi31','Nevajärvi34',
         'Jaakkoinsuo61','Jaakkoinsuo62', 
         'Parkano11']

txt=['a', 'b', 'c', 'd', 'e', 'f', 
    'g', 'h', 'i', 'j', 'k']

#sites =['neva11']
folder_meas = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Henkilöt/Stenberg/puusto_ym_data_Arille/Pohjavesiaineistot/'
i=0
for crd, k, tx, n in zip(coordinates, sites, txt, names):
    print (crd, k) 

    file_meas = wt_meas[k]['file']
    tubes = wt_meas[k]['tubes']
    dfmeas = pd.read_excel(folder_meas+file_meas, sheet_name ='CSV')
    dfmeas['date']=pd.to_datetime(dict(year=dfmeas.vuosi, month=dfmeas.kk, day=dfmeas.pv))
    dfmeas = dfmeas.set_index('date')
    dfmeas = dfmeas.drop(['vuosi','kk', 'pv'], axis=1)
    dfmeas = dfmeas.drop(columns=[col for col in dfmeas if col not in tubes])
    dfmeas = dfmeas/100.*-1
    
    ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/' + k + '.nc'
    params = para(period='start-end')
    sday = params[k]['start_date']
    end_date =params[k]['end_date']
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    dwt = ncf['strip']['dwt'][0,:,1:-1] 
    days, cols = np.shape(dwt)
    dfsim = pd.DataFrame(dwt, columns=range(cols), index=pd.date_range(sday,periods=days)) 
    ncf.close()
    

    fs=8        
    wt_max = np.max(dfsim, axis=1)
    wt_min = np.min(dfsim, axis=1)
    ax = fig.add_subplot(crd)
    ax.fill_between(dfsim.index, wt_max, wt_min, color='grey', alpha = 0.3)
    ax.set_ylim([-1.0,0.0])
    ax.plot(dfmeas, 'bo', markersize=2)
    ax.set_xlim([sday, end_date])
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    
    ax.text(0.02, 0.95, tx + '  ' + n, horizontalalignment='left',
               verticalalignment='top', fontsize=14, transform = ax.transAxes)

    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(11.0)
        xlabel_i.set_y(-0.01)

        
    if i in [1,3,5,7,9]:
        for ylabel_i in ax.get_yticklabels():
            ylabel_i.set_fontsize(0.0)
            ylabel_i.set_visible(False)
    else:
        for ylabel_i in ax.get_yticklabels():
            ylabel_i.set_fontsize(11.0)
            ylabel_i.set_x(-0.025)
            
        ax.set_ylabel('WT, m', fontsize = 14)


    if i in [9,10]:
        ax.set_xlabel('Time', fontsize = 14)
    i +=1
plt.tight_layout(h_pad=1.5)

#%%
# calculate measured and modelled late-summer WT and save to out{} dictionary
sites =['ansa21', 'ansa26', 'jaakkoin61', 'jaakkoin62', 'koira11', 'koira12',
        'neva11', 'neva14', 'neva31', 'neva34','parkano11']
out={}
for k in sites:
    file_meas = wt_meas[k]['file']
    tubes = wt_meas[k]['tubes']
    dfmeas = pd.read_excel(folder_meas+file_meas, sheet_name ='CSV')
    dfmeas['date']=pd.to_datetime(dict(year=dfmeas.vuosi, month=dfmeas.kk, day=dfmeas.pv))
    dfmeas = dfmeas.set_index('date')
    dfmeas = dfmeas.drop(['vuosi','kk', 'pv'], axis=1)
    dfmeas = dfmeas.drop(columns=[col for col in dfmeas if col not in tubes])
    dfmeas = dfmeas/100.*-1

    ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/' + k + '.nc'
    params = para(period='start-end')
    start_date = params[k]['start_date']
    end_date =params[k]['end_date']

    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    dwt = ncf['strip']['dwt'][0,:,1:-1] 
    days, cols = np.shape(dwt)
    dfsim = pd.DataFrame(dwt, columns=range(cols), index=pd.date_range(start_date,periods=days)) 
    ncf.close()
    
    
    print (k, start_date, end_date)
    dfmeas = dfmeas.sort_index()[str(dfsim.index[0]): str(dfsim.index[-1])]
    
    out[k]={}
    meansim = []
    stdsim =[]
    meanmeas=[]
    stdmeas=[]
    for yr in range(start_date.year, end_date.year+1):   
        wtsim = dfsim[str(yr)+'-07-01':str(yr)+'-08-31'].mean().values     #.values[:-1])
        meansim.append(wtsim.mean())
        stdsim.append(wtsim.std())
        wtmeas = dfmeas[str(yr)+'-07-01':str(yr)+'-08-31'].mean().values     #.values[:-1])
        meanmeas.append(wtmeas.mean())
        stdmeas.append(wtmeas.std())
    out[k]['meansim']=meansim
    out[k]['stdsim']=stdsim
    out[k]['meanmeas']=meanmeas
    out[k]['stdmeas']=stdmeas
print (out)



#%%   
# Collect biomass and volume growth from netcdf files
gro = np.zeros(len(sites))
grosd = np.zeros(len(sites))
bm_gr = np.zeros(len(sites))
bm_gr_sd = np.zeros(len(sites))

for n, k in enumerate(sites):
    ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/' + k + '.nc'
    params = para(period='start-end')
    sday = params[k]['start_date']
    end_date =params[k]['end_date']
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    vol = np.array(ncf['stand']['volume'][0,:,1:-1]) 
    growth = np.diff(vol, axis=0) 
    yrs, cols = np.shape(vol)
    dfgrowth = pd.DataFrame(growth, columns=range(cols)) 
    
    bm_growth = np.array(ncf['stand']['dominant']['NPP'][0,:,1:-1]*ncf['stand']['stems'][0,:,1:-1])
    dfbm = pd.DataFrame(bm_growth, columns=range(cols)) 
    
    ncf.close()
    gro[n] = np.mean(dfgrowth.mean(axis=0).values)
    grosd[n] = np.std(dfgrowth.mean(axis=0).values)
    bm_gr[n] = np.mean(dfbm.mean(axis=0).values)
    bm_gr_sd[n] = np.std(dfbm.mean(axis=0).values)


dvol = {'sites': sites, 'grsim':gro, 'grsd': grosd, 'bmgr':bm_gr, 'bmgrsd': bm_gr_sd}
dfvols = pd.DataFrame(data=dvol)
dfvols = dfvols.set_index('sites')

print (dfvols)

#%%     
def ojanen_2010(sfc, stand_v, t_peat, gs_wt):

      bd_d = {2: 0.14, 3: 0.11, 4: 0.10, 5: 0.08}                                      # Mese study: bulk densities in different fertility classes     g/cm3                                                            # peat layer thickness, cm            
      bd = bd_d[sfc] *1000.                                                            # set the bulk density according to site fertility class        kg/m3
      B = 350. ; T5zero = -46.02; T5ref = 10.
      Rref = (0.0695 + 3.7*10**(-4)*stand_v + 5.4*10**(-4)*bd + 1.2*10**(-3)*gs_wt) * 24.    #unit g/m2/h CO2-> g/m2/day
      Rhet = [rref*np.exp(B*((1./(T5ref-T5zero))-(1./(t_peat-T5zero)))) for rref in Rref]          
      Rhet = np.array(Rhet).T
      return np.sum(Rhet, axis=0)

#%%
"""  KUVA N: wt, biomass, vol  """
"""
needed: 
    grsim = vend - vini
    yrs
"""
import matplotlib.gridspec as gridspec

sites =['koira11', 'koira12',
        'ansa21', 'ansa26', 
        'neva11', 'neva14', 
        'neva31', 'neva34',
        'jaakkoin61', 'jaakkoin62', 
        'parkano11']

names = ['Koirasuo11','Koirasuo12',
         'Ansasaari21','Ansasaari26', 
         'Nevajärvi11','Nevajärvi14', 
         'Nevajärvi31','Nevajärvi34',
         'Jaakkoinsuo61','Jaakkoinsuo62', 
         'Parkano11']



fbio = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Susi combined/gr_bio.xlsx'
dfbio = pd.read_excel(fbio, index_col='site')

obs=[]
pred =[]
dfvols['grobs']=dfvols['grsim']*0.0
dfvols['yrs']=dfvols['grsim']*0.0
dfvols['bioobs']=dfvols['grsim']*0.0

for s,n in zip(sites, names):    
    print (s, np.round(params[s]['vol'][1] - params[s]['vol'][0]), dfvols.loc[s]['grsim'])
    dfvols.at[s, 'grobs']= params[s]['vol'][1] - params[s]['vol'][0]
    dfvols.at[s, 'yrs']= params[s]['end_date'].year -params[s]['start_date'].year +1.
    dfvols.at[s,'bioobs']=dfbio.at[s,'gr_bio']

fs=16
nsites = 11
figsi = (10,10)
fig= plt.figure(num = 'growth', figsize=figsi)   #Figsize(w,h), tuple inches 
col1 = 'yellow'
col2 = 'grey'
gs = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, wspace=0.25, hspace=0.25)


colors = plt.cm.jet(np.linspace(0,1,nsites))
#-----------WT figure ------------------------------------------------
#ax0 = fig.add_axes([0.08, 0.15, 0.25, 0.75]) #left, bottom, width, height)
ax0 = fig.add_subplot(gs[0,0])

mval=0.
col1 = 'yellow'
col2 = 'grey'
plt.fill_between([-1., mval],[mval,mval],[-1.,mval], color=col1, alpha=0.3)
plt.fill_between([-1., mval],[-1.,0.], [-1.,-1.], color=col2, alpha=0.3)

colors = plt.cm.jet(np.linspace(0,1,len(sites)))
si=[]
ob =[]
c = 0
for k, n in zip(sites, names):
    data = out[k]
    obs = np.array(data['meanmeas'])
    ob.extend(obs)
    obserr=np.array(data['stdmeas'])
    pre = np.array(data['meansim'])
    si.extend(pre)
    preerr=np.array(data['stdsim'])
    plt.plot(obs,pre,'o', markersize=10, label=n, color=colors[c])
    plt.errorbar(obs,pre,preerr*2,obserr*2, 'none', color=colors[c], capsize=4)
    c+=1
#plt.legend(loc='lower right', ncol=3)
si = np.array(si)
ob = (np.array(ob))

diff = ob - si
diff = diff[~np.isnan(diff)]
rmse = np.round(np.sqrt(np.square(diff).mean()),3)

data = {'meas':ob,'sim':si }
dfdata = pd.DataFrame.from_dict(data)
dfdata=dfdata.dropna()
#slope, intercept, r_value, p_value, std_err = stats.linregress(dfdata['meas'].values, dfdata['sim'].values)

fs = 16

obs_wt = dfdata['meas'].values[:,np.newaxis]
sim_wt =  dfdata['sim'].values
a, _, _, _ = np.linalg.lstsq(obs_wt,sim_wt)
eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(-0.9 , -0.03, 0.05)
plt.plot(x, a*x, 'k--', linewidth=2)


plt.xlim([-1.,0.0])
plt.ylim([-1.,0.0])
plt.text(-0.9, -0.1, eq, fontsize = fs-1)
plt.text(-0.9, -0.15, 'RMSE ' + str(rmse), fontsize=fs-1)
plt.xlabel('Observed $\it{WT}$, m', fontsize =fs )
plt.ylabel('Predicted $\it{WT}$, m', fontsize =fs)

for tick in ax0.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax0.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

ax0.text(0.04, 0.95, 'a', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax0.transAxes,  fontweight='bold')

#----------BM figure ---------------------------------------------------
ax1 = fig.add_subplot(gs[0,1])

mval=10000.
plt.fill_between([0., mval],[mval,mval],[0.,mval], color=col1, alpha=0.3)
plt.fill_between([0., mval],[0.,mval], [0.,0.], color=col2, alpha=0.3)

#dfvols.to_excel(folder_out+'dfvols.xlsx')
c=0
for s,n in zip(sites, names):
    obs = dfvols.loc[s]['bioobs']
    pre = dfvols.loc[s]['bmgr']
    plt.plot(obs,pre,'o', markersize=10, label=n, color=colors[c])
    preerr=dfvols.loc[s]['bmgrsd']
    plt.errorbar(obs,pre,preerr*2, 0, 'none', color=colors[c], capsize=4)
    

    c+=1
plt.xlabel('Observed bm growth, $kg \ ha^{-1} yr^{-1}$ ', fontsize =fs )
plt.ylabel('Predicted bm growth, $kg \ ha^{-1} yr^{-1}$ ', fontsize =fs, labelpad=-7.5)
plt.xlim([0,mval])
plt.ylim([0,mval])
obs_growths = dfvols['bioobs'].values[:,np.newaxis]
sim_growths = dfvols['bmgr'].values

rmse = np.round(np.sqrt(np.square(obs_growths-sim_growths).mean()),0)
a, _, _, _ = np.linalg.lstsq(obs_growths, sim_growths)
eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(1000.0 , 9000.0, 100)
plt.plot(x, a*x, 'k--', linewidth=2)
plt.text(5000, 8300, 'RMSE ' + str(rmse), fontsize=fs-1)


plt.text(5000, 8900, eq , fontsize=fs-1)


ax1.text(0.04, 0.95, 'b', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax1.transAxes, fontweight='bold')

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

#-------------------Vol figure ------------------
ax2 = fig.add_subplot(gs[1,0])
        
mval=12.
plt.fill_between([0., mval],[mval,mval],[0.,mval], color=col1, alpha=0.3)
plt.fill_between([0., mval],[0.,mval], [0.,0.], color=col2, alpha=0.3)
c=0
obsvols = []
prevols = []

for s,n in zip(sites, names):
    yrs = params[s]['end_date'].year -params[s]['start_date'].year +1.
    obs = dfvols.loc[s]['grobs']/yrs
    pre = dfvols.loc[s]['grsim']
    plt.plot(obs,pre,'o', markersize=10, label=n, color=colors[c])
    preerr=dfvols.loc[s]['grsd']
    plt.errorbar(obs,pre,preerr*2, 0, 'none', color=colors[c], capsize=4)
    obsvols.append(obs)
    prevols.append(pre)
    c+=1

ax1.legend(loc='upper center', bbox_to_anchor=(-0.15, 1.3),
          ncol=4,  fontsize=fs-3)
#ax1.legend(loc='lower left', ncol=1, fontsize=fs-3)
plt.xlim([0,mval])
plt.ylim([0,mval])

plt.xlabel('Observed $\it{i_V}$, $m^{3} ha^{-1} yr^{-1}$ ', fontsize =fs )
plt.ylabel('Predicted $\it{i_V}$, $m^{3} ha^{-1} yr^{-1}$ ', fontsize =fs)

obs_growths = dfvols['grobs'].values/dfvols['yrs'].values
sim_growths =dfvols['grsim'].values
slope, intercept, r_value, p_value, std_err = stats.linregress(obs_growths, 
                                                               sim_growths)

diffv = dfvols['grobs'].values/dfvols['yrs'].values - dfvols['grsim'].values
rmse = np.round(np.sqrt(np.square(diffv).mean()),3)
obs_growths = obs_growths[:,np.newaxis]
a, _, _, _ = np.linalg.lstsq(obs_growths, sim_growths)


eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(1.0 , 8.0, 0.5)
plt.text(6, 8.3, 'RMSE ' + str(rmse), fontsize=fs-1)

plt.plot(x, a*x, 'k--', linewidth=2)
predicted = np.ravel(a[0]*obs_growths)

plt.text(6, 9, eq , fontsize=fs-1)

ax2.text(0.04, 0.95, 'c', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax2.transAxes,  fontweight='bold')

for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

#-------------------CO2 figure ------------------
ax3 = fig.add_subplot(gs[1,1])
        
mval=3000.
plt.fill_between([0., mval],[mval,mval],[0.,mval], color=col1, alpha=0.3)
plt.fill_between([0., mval],[0.,mval], [0.,0.], color=col2, alpha=0.3)
esarr = np.empty(0)
emps = np.empty(0)
c=0
for k in (sites):
    ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/' + k + '.nc'
    params = para(period='start-end')
    sday = params[k]['start_date']
    end_date =params[k]['end_date']
    sfc = params[k]['sfc']
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    vol = np.array(ncf['stand']['volume'][0,:,1:-1]) 
    yrs, cols = np.shape(vol)
    dfvol = pd.DataFrame(vol, columns=range(cols)) 
    
    dwt = ncf['strip']['dwtyr_growingseason'][0,:,1:-1] 
    yrs, cols = np.shape(dwt)
    dfwt = pd.DataFrame(dwt, columns=range(cols)) 
    
    peat_t = ncf['temperature']['T'][0,:,3]
    days = np.shape(peat_t)[0]
    dft = pd.DataFrame(peat_t, columns=['T'], index=pd.date_range(sday,periods=days))    
    
    esom_co2 = ncf['esom']['Mass']['co2'][0,:, 1:-1]/10. #kg ha-1 yr-1 ->g m-2 yr-1
    
    ncf.close()
    empirical = np.zeros(np.shape(dwt))
    for n, yr in enumerate(range(sday.year, end_date.year+1)):
        t_peat_yr = np.ravel(dft.loc[str(yr)].values)
        rhet = ojanen_2010(sfc, vol[n+1], t_peat_yr, dfwt.loc[n+1].values*-100.) 
        empirical[n+1,:] = rhet       
    for m in range(1,yrs):
        plt.plot(empirical[m, :], esom_co2[m, :], color=colors[c])
    
    c+=1
    
    esarr = np.append(esarr,np.ravel(esom_co2[1:,:]) )
    emps = np.append(emps, np.ravel(empirical[1:,:]))

plt.xlim([0,mval])
plt.ylim([0,mval])

plt.xlabel('Empirical $CO_2$, $g m^{-2} yr^{-1}$ ', fontsize =fs )
plt.ylabel('Esom $CO_2$, $g m^{-3} yr^{-2}$ ', fontsize =fs, labelpad=-5.0)

obs_growths = emps
sim_growths =esarr
slope, intercept, r_value, p_value, std_err = stats.linregress(obs_growths, 
                                                               sim_growths)

#diffv = dfvols['grobs'].values/dfvols['yrs'].values - dfvols['grsim'].values
#rmse = np.round(np.sqrt(np.square(diffv).mean()),3)
obs_growths = obs_growths[:,np.newaxis]
a, _, _, _ = np.linalg.lstsq(obs_growths, sim_growths)


eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(200.0 , 2000.0, 0.5)
#plt.text(5, 8.5, 'RMSE ' + str(rmse), fontsize=fs-1)

plt.plot(x, a*x, 'k--', linewidth=2)
predicted = np.ravel(a[0]*obs_growths)

plt.text(2000, 2000, eq , fontsize=fs-1)

ax3.text(1.3, 0.95, 'd', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax2.transAxes,  fontweight='bold')

for tick in ax3.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax3.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 


#%%


for n, k in enumerate(sites[:1]):
    ff = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/susi_22_out/' + k + '.nc'
    params = para(period='start-end')
    sday = params[k]['start_date']
    end_date =params[k]['end_date']
    sfc = params[k]['sfc']
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    vol = np.array(ncf['stand']['volume'][0,:,1:-1]) 
    yrs, cols = np.shape(vol)
    dfvol = pd.DataFrame(vol, columns=range(cols)) 
    
    dwt = ncf['strip']['dwtyr_growingseason'][0,:,1:-1] 
    yrs, cols = np.shape(dwt)
    dfwt = pd.DataFrame(dwt, columns=range(cols)) 
    
    peat_t = ncf['temperature']['T'][0,:,3]
    days = np.shape(peat_t)[0]
    dft = pd.DataFrame(peat_t, columns=['T'], index=pd.date_range(sday,periods=days))    
    
    esom_co2 = ncf['esom']['Mass']['co2'][0,:, 1:-1]/10. #kg ha-1 yr-1 ->g m-2 yr-1
    
    ncf.close()
    empirical = np.zeros(np.shape(dwt))
    for n, yr in enumerate(range(sday.year, end_date.year+1)):
        t_peat_yr = np.ravel(dft.loc[str(yr)].values)
        rhet = ojanen_2010(sfc, vol[n+1], t_peat_yr, dfwt.loc[n+1].values*-100.) 
        empirical[n+1,:] = rhet       
   
dfesom = pd.DataFrame(data=esom_co2)
dfempirical = pd.DataFrame(data=empirical)
#%%
sfcs = [2,2,3,3,3,3,3,3,5,5,5]  
print ('***********************')
dfresid = pd.DataFrame(list(zip(names,sfcs,obsvols, prevols )),
                       columns=['name', 'sfc', 'obs', 'pre'])
dfresid['residual'] = dfresid['pre'] - dfresid['obs']
print (dfresid.groupby(['sfc']).mean())
print (dfresid)