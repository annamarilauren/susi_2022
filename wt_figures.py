# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 18:04:05 2020

@author: alauren
"""
#%%
import datetime
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from scipy import stats
from itertools import product
from dwts_para import para
#from sklearn.metrics import r2_score
import seaborn as sns
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
sns.set()
#%%
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

"""
          'ansa11':{'file':'muhos_1_pohjavesi_koottu.xlsx',
                    'tubes':[ 1, 2, 3, 4, 5, 6],
                    'dist':[5.,s_ansa/2.,5.,5.,s_ansa/2.,5.]},
          'ansa16':{'file':'muhos_1_pohjavesi_koottu.xlsx',
                    'tubes':[28, 29, 30],
                    'dist':[5., s_ansa/2., 5.]},
"""
#%%
folder_meas = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Stenberg/puusto_ym_data_Arille/Pohjavesiaineistot/'
folder_sim = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/'

#%%
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


#frame1 = plt.gca()

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
folder_meas = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Stenberg/puusto_ym_data_Arille/Pohjavesiaineistot/'
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
    
    folder_sim = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/'
    file_sim = 'gwl_'+ k +'.xlsx'
    dfsim = pd.read_excel(folder_sim+file_sim)
    start_date = params[k]['start_date']
    end_date = params[k]['end_date']
    dfsim = dfsim.set_index(pd.date_range(start = start_date, end= end_date, freq='D'))
    _,cols = np.shape(dfsim)
    dropcols = [0,1, cols-1]
    dfsim = dfsim.drop(dfsim.columns[dropcols], axis = 1)

    #dfmeas = dfmeas[str(dfsim.index[0]): str(dfsim.index[-1])]

    fs=8        
    wt_max = np.max(dfsim, axis=1)
    wt_min = np.min(dfsim, axis=1)
    ax = fig.add_subplot(crd)
    ax.fill_between(dfsim.index, wt_max, wt_min, color='grey', alpha = 0.3)
    ax.set_ylim([-1.0,0.0])
    ax.plot(dfmeas, 'bo', markersize=2)
    ax.set_xlim([start_date, end_date])
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
sites =['ansa21', 'ansa26', 'jaakkoin61', 'jaakkoin62', 'koira11', 'koira12',
        'neva11', 'neva14', 'neva31', 'neva34','parkano11', 'parkano12']
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
    
    folder_sim = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/'
    file_sim = 'gwl_'+ k +'.xlsx'
    dfsim = pd.read_excel(folder_sim+file_sim)
    start_date = params[k]['start_date']
    end_date = params[k]['end_date']
    dfsim = dfsim.set_index(pd.date_range(start = start_date, end= end_date, freq='D'))
    _,cols = np.shape(dfsim)
    dropcols = [0,1, cols-1]
    dfsim = dfsim.drop(dfsim.columns[dropcols], axis = 1)

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
folder_meas = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Stenberg/puusto_ym_data_Arille/Pohjavesiaineistot/'
for k in params.keys():
#for k in ['neva11']:
    #if k != 'jaakkoin62':
    #if not params[k]['thinning']:
    #if params[k]['thinning']:
        file_meas = wt_meas[k]['file']
        tubes = wt_meas[k]['tubes']
        dfmeas = pd.read_excel(folder_meas+file_meas, sheet_name ='CSV')
        dfmeas['date']=pd.to_datetime(dict(year=dfmeas.vuosi, month=dfmeas.kk, day=dfmeas.pv))
        dfmeas = dfmeas.set_index('date')
        dfmeas = dfmeas.drop(['vuosi','kk', 'pv'], axis=1)
        dfmeas = dfmeas.drop(columns=[col for col in dfmeas if col not in tubes])
        dfmeas = dfmeas/100.*-1
        
        folder_sim = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/'
        file_sim = 'gwl_'+ k +'.xlsx'
        dfsim = pd.read_excel(folder_sim+file_sim)
        start_date = params[k]['start_date']
        end_date = params[k]['end_date']
        dfsim = dfsim.set_index(pd.date_range(start = start_date, end= end_date, freq='D'))
        _,cols = np.shape(dfsim)
        dropcols = [0,1, cols-1]
        dfsim = dfsim.drop(dfsim.columns[dropcols], axis = 1)
 
        
        print (start_date.year, end_date.year)
        print (dfmeas)
        #import sys; sys.exit()
        try:
            dfmeas = dfmeas[start_date.year:end_date.year]
        except:
            pass
        fs=15        
        fig = plt.figure(num=k, figsize=(15,8))
        plt.subplot(111)
        wt_mean = np.mean(dfsim, axis=1)
        wt_max = np.max(dfsim, axis=1)
        wt_min = np.min(dfsim, axis=1)
        
        #plt.plot(np.mean(dfsim, axis=1), 'b-')
        plt.fill_between(dfsim.index, wt_max, wt_min, color='green', alpha = 0.3)
        plt.xlim([start_date, end_date])
        plt.ylim([-1.0,0.0])
        plt.plot(dfmeas, 'ro')

        plt.xlabel('Time', fontsize = fs)
        plt.ylabel('WT, m', fontsize = fs)
        
        



#%%
        print (k)
        
        cols = np.around(np.array(wt_meas[k]['dist'])/2-1)
        cols = cols.astype(int)
        
        fig = plt.figure(num=k + 'tubes', figsize=(16,8))
        
        colors = plt.cm.jet(np.linspace(0,1,len(cols)))
        fs=15
        
        
        ax0 = fig.add_axes([0.05, 0.15, 0.55, 0.75]) #left, bottom, width, height)
        for n,(s,m) in enumerate(zip(cols,tubes)):
            plt.plot(dfsim[s], color = colors[n], label= m)
            plt.plot(dfmeas[m], 'o', color=colors[n], label=m)
        wt_max = np.max(dfsim, axis=1)
        wt_min = np.min(dfsim, axis=1)
        plt.fill_between(dfsim.index, wt_max, wt_min, color='green', alpha = 0.3)
        plt.legend(loc='lower left', ncol=2)
        plt.xlim([start_date, end_date])
        plt.xlabel('Time', fontsize = fs)
        plt.ylabel('WT, m', fontsize = fs)
        
        dfsim2 = dfsim.drop(columns=[col for col in dfsim if col not in cols])
        meas_days = list(dfmeas.index)
        dfsim2 = dfsim2.loc[meas_days]
        
        
        #fig = plt.figure(num=k+'rmse', figsize=(8,8))
        ax1 = fig.add_axes([0.65, 0.15, 0.3, 0.75]) #left, bottom, width, height)
        
        mval=0.
        col1 = 'yellow'
        col2 = 'grey'
        plt.fill_between([-1., mval],[mval,mval],[-1.,mval], color=col1, alpha=0.3)
        plt.fill_between([-1., mval],[-1.,0.], [-1.,-1.], color=col2, alpha=0.3)
        
        meas =[]; sim=[]
        for n, (s,t) in enumerate(zip(cols, tubes)):
            plt.plot(dfmeas[t], dfsim2[s], 'o', color=colors[n], label=t)
            meas.append(dfmeas[t].values)
            sim.append(dfsim2[s].values)
        plt.legend(loc='lower left')
        
        meas=np.ravel((np.array(meas)))
        sim =np.ravel((np.array(sim)))
        diff = meas - sim
        diff = diff[~np.isnan(diff)]
        rmse = np.round(np.sqrt(np.square(diff).mean()),3)
        
        
        data = {'meas':meas,'sim':sim }
        dfdata = pd.DataFrame.from_dict(data)
        dfdata=dfdata.dropna()
        slope, intercept, r_value, p_value, std_err = stats.linregress(dfdata['meas'].values, dfdata['sim'].values)
        
        eq = 'y = ' + str(np.round(intercept,3)) + ' +  ' + str(np.round(slope,3)) + 'x'
        x =np.arange(-0.9 , -0.03, 0.05)
        plt.plot(x, intercept + slope*x, 'k--', linewidth=2)
        plt.xlim([-1.,0.0])
        plt.ylim([-1.,0.0])
        plt.text(-0.8, -0.05, eq, fontsize = fs-1)
        plt.text(-0.8, -0.1, 'RMSE ' + str(rmse), fontsize=fs-1)
        plt.text(-0.8, -0.15, '$R^{2}$  ' + str(np.round(r_value**2,3)), fontsize = fs-1)
        plt.xlabel('Observed WT, m', fontsize =fs )
        plt.ylabel('Predicted WT, m', fontsize =fs)

#%%        
"""  KUVA N: wt, biomass, vol  """


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


folder_out = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/'
file_vol = 'vols.xlsx'
dfvols = pd.read_excel(folder_out+file_vol, index_col=[0])
dfvols['grsim']=dfvols['v']-dfvols['v_ini']
fbio = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Susi combined/gr_bio.xlsx'
dfbio = pd.read_excel(fbio, index_col='site')

obs=[]
pred =[]
dfvols['grobs']=dfvols['grsim']*0.0
dfvols['yrs']=dfvols['grsim']*0.0
dfvols['bioobs']=dfvols['grsim']*0.0

#for s in params.keys():
for s,n in zip(sites, names):    
    print (s, np.round(params[s]['vol'][1] - params[s]['vol'][0]), dfvols.loc[s]['grsim'])
    dfvols.at[s, 'grobs']= params[s]['vol'][1] - params[s]['vol'][0]
    dfvols.at[s, 'yrs']= params[s]['end_date'].year -params[s]['start_date'].year +1.
    dfvols.at[s,'bioobs']=dfbio.at[s,'gr_bio']

fs=16
nsites = 11
figsi = (18,6)
fig= plt.figure(num = 'growth', figsize=figsi)   #Figsize(w,h), tuple inches 
col1 = 'yellow'
col2 = 'grey'


colors = plt.cm.jet(np.linspace(0,1,nsites))
#-----------WT figure ------------------------------------------------
ax0 = fig.add_axes([0.08, 0.15, 0.25, 0.75]) #left, bottom, width, height)

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
#ax0 = fig.add_axes([0.08, 0.15, 0.375, 0.75]) #left, bottom, width, height)
ax1 = fig.add_axes([0.40, 0.15, 0.25, 0.75]) #left, bottom, width, height)
mval=10000.
plt.fill_between([0., mval],[mval,mval],[0.,mval], color=col1, alpha=0.3)
plt.fill_between([0., mval],[0.,mval], [0.,0.], color=col2, alpha=0.3)

dfvols.to_excel(folder_out+'dfvols.xlsx')
c=0
for s,n in zip(sites, names):
    obs = dfvols.loc[s]['bioobs']
    pre = dfvols.loc[s]['bmgr']
    plt.plot(obs,pre,'o', markersize=10, label=n, color=colors[c])
    c+=1
plt.xlabel('Observed biomass growth, $kg \ ha^{-1} yr^{-1}$ ', fontsize =fs )
plt.ylabel('Predicted biomass growth, $kg \ ha^{-1} yr^{-1}$ ', fontsize =fs)
plt.xlim([0,mval])
plt.ylim([0,mval])
obs_growths = dfvols['bioobs'].values[:,np.newaxis]
sim_growths = dfvols['bmgr'].values

rmse = np.round(np.sqrt(np.square(obs_growths-sim_growths).mean()),0)
a, _, _, _ = np.linalg.lstsq(obs_growths, sim_growths)
eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(1000.0 , 9000.0, 100)
plt.plot(x, a*x, 'k--', linewidth=2)
plt.text(5000, 8400, 'RMSE ' + str(rmse), fontsize=fs-1)


#plt.text(5000, 8500, eq + '  $R^{2}$ ' + str(np.round(r2,2)), fontsize=fs-1)
plt.text(5000, 8900, eq , fontsize=fs-1)


ax1.text(0.04, 0.95, 'b', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax1.transAxes, fontweight='bold')

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

#-------------------Vol figure ------------------
#ax1 = fig.add_axes([0.55, 0.15, 0.375, 0.75]) #left, bottom, width, height)
ax2 = fig.add_axes([0.72, 0.15, 0.25, 0.75]) #left, bottom, width, height)
        
mval=10.
plt.fill_between([0., mval],[mval,mval],[0.,mval], color=col1, alpha=0.3)
plt.fill_between([0., mval],[0.,mval], [0.,0.], color=col2, alpha=0.3)
c=0
for s,n in zip(sites, names):
    yrs = params[s]['end_date'].year -params[s]['start_date'].year +1.
    obs = dfvols.loc[s]['grobs']/yrs
    pre = dfvols.loc[s]['grsim']/yrs
    plt.plot(obs,pre,'o', markersize=10, label=n, color=colors[c])
    c+=1
ax1.legend(loc='lower left', ncol=1, fontsize=fs-3)
plt.xlim([0,mval])
plt.ylim([0,mval])

plt.xlabel('Observed $\it{i_V}$, $m^{3} ha^{-1} yr^{-1}$ ', fontsize =fs )
plt.ylabel('Predicted $\it{i_V}$, $m^{3} ha^{-1} yr^{-1}$ ', fontsize =fs)

obs_growths = dfvols['grobs'].values/dfvols['yrs'].values
sim_growths =dfvols['grsim'].values/dfvols['yrs'].values
slope, intercept, r_value, p_value, std_err = stats.linregress(obs_growths, 
                                                               sim_growths)

diffv = dfvols['grobs'].values/dfvols['yrs'].values - dfvols['grsim'].values/dfvols['yrs'].values
rmse = np.round(np.sqrt(np.square(diffv).mean()),3)
obs_growths = obs_growths[:,np.newaxis]
a, _, _, _ = np.linalg.lstsq(obs_growths, sim_growths)


eq = 'y = ' +  str(np.round(a[0],3)) + 'x'
x =np.arange(1.0 , 8.0, 0.5)
#plt.plot(x, intercept + slope*x, 'k--', linewidth=2)
plt.text(5, 8.5, 'RMSE ' + str(rmse), fontsize=fs-1)

plt.plot(x, a*x, 'k--', linewidth=2)
predicted = np.ravel(a[0]*obs_growths)

plt.text(5, 9, eq , fontsize=fs-1)

ax2.text(0.04, 0.95, 'c', horizontalalignment='left',
               verticalalignment='top', fontsize=18, transform = ax2.transAxes,  fontweight='bold')

for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(14) 
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(14) 

#print (np.corrcoef(obs_growths, predicted)[0, 1]**2)

#print (np.mean(predicted))

#%%
