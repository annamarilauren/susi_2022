# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 17:37:58 2022

@author: alauren
"""
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset  
import numpy as np
import pandas as pd

def create_profile_line(ax, wt, wtmin, sd, cols, ylabel, label, fs, facecolor, colorin, hidex=False, hidey=False, elevation=None):
    if elevation is not None: ax.plot(elevation, color='brown', label = 'soil surface')
    ax.plot(wt, color=colorin, label = label)
    ax.fill_between(range(cols), wt+sd*2, wt-sd, color=colorin, alpha=0.075)
    if elevation is not None: 
        drainnorm = elevation - 0.35
    else:
        drainnorm = -0.35
    ax.hlines(y= drainnorm, xmin=0, xmax = cols, color='red',linestyles='--')
    ax.get_xaxis().set_visible(False) 
    ax.tick_params(axis='y', labelsize=fs)
    if elevation is None: ax.set_ylim([wtmin,0])
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.legend()
    ax.grid(visible=False)
    ax.set_facecolor(facecolor)
    if hidex: 
        ax.get_xaxis().set_visible(False) 
    else:
        ax.tick_params(axis='x', labelsize=fs)
    if hidey: 
        ax.get_yaxis().set_visible(False)
    else:
        ax.get_yaxis().set_visible(True)
        

    return ax

def create_profile_boxplot(ax, datain, cols, colorin, title, label, fs, facecolor, zero=False, hidex=True, hidey=False):
    
    df = pd.DataFrame(data=datain, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes=colorin, whiskers=colorin, medians=colorin, caps=colorin),
                  boxprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color=colorin, alpha=0.6),
                  showfliers=False, grid=False, rot=0)

    if zero: ax.hlines(y= 0, xmin=0, xmax = cols, color='red',linestyles='--')
    meanval = df.mean(axis=0)
    meanval = np.round(meanval.mean(),2)
    title = title + ': mean '+ str(meanval)
    ax.set_title(title, fontsize = fs)
    ax.set_ylabel(label, fontsize=fs)
    ax.set_facecolor(facecolor)
    #ax.get_xaxis().set_visible(False) 
    if hidex: 
        ax.get_xaxis().set_visible(False) 
    else:
        ax.tick_params(axis='x', labelsize=fs)
    if hidey: 
        ax.get_yaxis().set_visible(False)
    else:
        ax.get_yaxis().set_visible(True)
        ax.set_ylabel(label, fontsize=fs)
        ax.tick_params(axis='y', labelsize=fs)

    return ax

def hydrology(ff, scen):
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    
    #-------------- Water tables cross section-------------------------------------
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='hydro', figsize=(15,18))   #width, height
    gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.25, hspace=0.25)
    
    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    wtgs = np.mean(ncf['strip']['dwtyr_growingseason'][scen,:, :], axis = 0)
    sdgs = np.std(ncf['strip']['dwtyr_growingseason'][scen,:, :], axis = 0)
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    wtmin = min(wtls) -0.2
    
    ax = fig.add_subplot(gs[10:, :4])
    ax = create_profile_line(ax, wt, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, 'blue', hidex=True, hidey=False)

    ax = fig.add_subplot(gs[10:, 4:8])
    ax = create_profile_line(ax, wtgs, wtmin, sdgs, cols, None, 'annual', fs, facecolor, 'green', hidex=True, hidey=True)

    ax = fig.add_subplot(gs[10:, 8:])
    ax = create_profile_line(ax, wtls, wtmin, sdls, cols, None, 'annual', fs, facecolor, 'orange', hidex=True, hidey=True)

    #----------Water tables time series
    
    wt = np.mean(ncf['strip']['dwt'][scen,:, :], axis = 1)
    days = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwt'][scen,:, :], axis = 1)# sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    
    axwtts = fig.add_subplot(gs[8:10, :])
    axwtts.plot(wt, color='green', label = 'WT')
    #axwtts.fill_between(range(days), wt+sd*2, wt-sd*2, color='green', alpha=0.25)
    axwtts.hlines(y= -0.35, xmin=0, xmax = days, color='red',linestyles='--')
    for c in range(1,cols-1):
        axwtts.plot(range(days), ncf['strip']['dwt'][scen,:, c], alpha=0.2)
        
    #axwtts.get_xaxis().set_visible(False) 
    axwtts.tick_params(axis='y', labelsize=fs)
    axwtts.set_ylim([wtmin,0])
    axwtts.set_ylabel('WT m', fontsize=fs)
    axwtts.legend(loc= 'upper left')
    axwtts.grid(visible=False)
    axwtts.set_facecolor(facecolor)
    
    runoff = np.cumsum(ncf['strip']['roff'][scen,:]) 
    ulimruno = max(runoff)*1.1*1000.
    axruno = fig.add_subplot(gs[7, :]) #axwtts.twinx()
    axruno.plot(range(len(runoff)), runoff*1000., color='blue', label='total runoff')
    axruno.set_ylim([0.,ulimruno ])
    axruno.fill_between(range(len(runoff)), 0.0, runoff*1000., color='blue', alpha=0.3)
    axruno.grid(visible=False)
    axruno.set_ylabel('mm', fontsize=fs)
    axruno.tick_params(axis='y', labelsize=fs)
    axruno.legend(loc='upper left')
    axruno.set_facecolor(facecolor)
    axruno.get_xaxis().set_visible(False) 

    runoff = np.cumsum(ncf['strip']['roffwest'][scen,:])        
    axruno = fig.add_subplot(gs[6, :]) #axwtts.twinx()
    axruno.plot(range(len(runoff)), runoff*1000., color='green', label='west runoff')
    axruno.set_ylim([0., ulimruno])
    axruno.fill_between(range(len(runoff)), 0.0, runoff*1000., color='green', alpha=0.3)
    axruno.grid(visible=False)
    axruno.set_ylabel('mm', fontsize=fs)
    axruno.tick_params(axis='y', labelsize=fs)
    axruno.legend(loc='upper left')
    axruno.set_facecolor(facecolor)
    axruno.get_xaxis().set_visible(False) 
    
    runoff = np.cumsum(ncf['strip']['roffeast'][scen,:])            
    axruno = fig.add_subplot(gs[5, :]) #axwtts.twinx()
    axruno.plot(range(len(runoff)), runoff*1000., color='red', label='east runoff')
    axruno.set_ylim([0., ulimruno])
    axruno.fill_between(range(len(runoff)), 0.0, runoff*1000., color='red', alpha=0.3)
    axruno.grid(visible=False)
    axruno.set_ylabel('mm', fontsize=fs)
    axruno.tick_params(axis='y', labelsize=fs)
    axruno.legend(loc='upper left')
    axruno.set_facecolor(facecolor)
    axruno.get_xaxis().set_visible(False) 
    
    runoff = np.cumsum(np.mean(ncf['strip']['surfacerunoff'][scen,:,:],axis=1))               
    axruno = fig.add_subplot(gs[4, :]) #axwtts.twinx()
    axruno.plot(range(len(runoff)), runoff*1000., color='orange', label='surface runoff')
    axruno.set_ylim([0., ulimruno])
    axruno.fill_between(range(len(runoff)), 0.0, runoff*1000., color='orange', alpha=0.3)
    axruno.grid(visible=False)
    axruno.set_ylabel('mm ', fontsize=fs)
    axruno.tick_params(axis='y', labelsize=fs)
    axruno.legend(loc='upper left')
    axruno.set_facecolor(facecolor)
    axruno.get_xaxis().set_visible(False) 
    
    
    #------deltas-----------------
    deltas = ncf['strip']['deltas'][scen,1:, :]*1000.
    dfdeltas = pd.DataFrame(data=deltas, columns=list(range(cols)))
 
    ax = fig.add_subplot(gs[2:4, :4])
    ax = create_profile_boxplot(ax, dfdeltas, cols, 'blue', 'Through soil surface', 'Water flux mm $yr^{-1}$', fs, facecolor, zero=False)

    #------ETs-----------------
    ET = ncf['cpy']['ET_yr'][scen,1:, :]*1000.
    dfET = pd.DataFrame(data=ET, columns=list(range(cols)))
    
    ax = fig.add_subplot(gs[2:4, 4:8])
    ax = create_profile_boxplot(ax, dfET, cols, 'green', 'ET', '', fs, facecolor, zero=False, hidex=True)
    
    
    #------transpi-----------------
    transpi = ncf['cpy']['transpi_yr'][scen,1:, :]*1000.
    dftranspi = pd.DataFrame(data=transpi, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, 8:])
    ax = create_profile_boxplot(ax, dftranspi, cols, 'orange', 'Transpiration', '', fs, facecolor, zero=False, hidex=True)
    
    
    #------efloor-----------------
    efloor = ncf['cpy']['efloor_yr'][scen,1:, :]*1000.
    dfefloor = pd.DataFrame(data=efloor, columns=list(range(cols)))
    
    ax = fig.add_subplot(gs[:2, :4])    
    ax = create_profile_boxplot(ax, dfefloor, cols, 'blue', 'Soil evaporation', 'Water flux mm $yr^{-1}$', fs, facecolor, zero=False, hidex=True)
    
    #------SWE max-----------------
    swe = ncf['cpy']['SWEmax'][scen,1:, :]
    dfswe = pd.DataFrame(data=swe, columns=list(range(cols)))

    ax= fig.add_subplot(gs[:2, 4:8])
    ax = create_profile_boxplot(ax, dfswe, cols, 'green', 'Max snow water equivalent', '', fs, facecolor, zero=False, hidex=True)
    
    #------Interc max-----------------
    interc = ncf['cpy']['interc_yr'][scen,1:, :]*1000.
    dfinterc = pd.DataFrame(data=interc, columns=list(range(cols)))
   
    ax= fig.add_subplot(gs[:2, 8:])
    ax = create_profile_boxplot(ax, dfinterc, cols, 'orange', 'Mean interception storage', '', fs, facecolor, zero=False, hidex=True)
    

    ncf.close()



def stand(ff, scen):
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode

    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='stand', figsize=(15,18))   #width, height
    gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.25, hspace=0.25)
  
    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    wtgs = np.mean(ncf['strip']['dwtyr_growingseason'][scen,:, :], axis = 0)
    sdgs = np.std(ncf['strip']['dwtyr_growingseason'][scen,:, :], axis = 0)
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    wtmin = min(wtls) -0.2
  
    #-------water table as a reference-----------------------------
    ax = fig.add_subplot(gs[10:, :4])
    ax.plot(wtls, color='orange', label = 'late summer')
    ax.fill_between(range(cols), wtls+sdls*2, wtls-sd*2, color='orange', alpha=0.075)
    ax.hlines(y= -0.35, xmin=0, xmax = cols, color='red',linestyles='--')
    ax.get_xaxis().set_visible(False) 
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylim([wtmin,0])
    ax.set_ylabel('WT m', fontsize=fs)
    ax.legend()
    ax.grid(visible=False)
    ax.set_facecolor(facecolor)
    
    
    #------------stand growth--------------------
    vol = ncf['stand']['volume'][scen,:, :]
    growth = np.diff(vol, axis=0)
    dfgrowth = pd.DataFrame(data=growth, columns=list(range(cols)))
    axgrowth = fig.add_subplot(gs[8:10, :4])
    dfgrowth.boxplot(ax = axgrowth,
                  color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    axgrowth.set_title('Stand growth')
    axgrowth.set_ylabel('$m^3 ha^{-1} yr^{-1}$', fontsize=fs)
    
    axgrowth.get_xaxis().set_visible(False)
    axgrowth.tick_params(axis='y', labelsize=fs)
    axgrowth.set_facecolor(facecolor)
    
    #-------------end volume-------------------
    ax = fig.add_subplot(gs[10:, 4:8])
    totvol = vol[-1,:]
    domvol = ncf['stand']['dominant']['volume'][scen,-1, :]
    subdomvol = ncf['stand']['subdominant']['volume'][scen,-1, :]
    undervol= ncf['stand']['under']['volume'][scen,-1, :]
    df = pd.DataFrame({'total': totvol, 'dominant': domvol,
                       'subdominant': subdomvol, 
                       'under':undervol}, index=range(cols))
    df.plot(kind='bar', width=1.1, ax=ax)
    ax.set_title('Stand volume')
    #ax.set_ylabel('$m^3 ha^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    ax.legend(loc='upper right', fontsize=7)
    
    #----- volume growth------------------------
    ax = fig.add_subplot(gs[10:, 8:])
    totvol = vol[:,:]
    domvol = ncf['stand']['dominant']['volume'][scen,:, :]
    subdomvol = ncf['stand']['subdominant']['volume'][scen,:, :]
    undervol= ncf['stand']['under']['volume'][scen,:, :]
    for c in range(cols):
        ax.plot(totvol[:,c], alpha=0.2)
    for c in range(cols):
        ax.plot(domvol[:,c], alpha=0.2)
    for c in range(cols):
        ax.plot(subdomvol[:,c], alpha=0.2)
    for c in range(cols):
        ax.plot(undervol[:,c], alpha=0.2)
        
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title('Stand volume increment')
        
  
    #----- volume growth ------------------------
    ax = fig.add_subplot(gs[8:10, 4:8])
    vol = ncf['stand']['volume'][scen,0:, :]
    for c in range(cols):
        ax.plot(np.diff(vol[:,c]), alpha=0.2)
        
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title('Volume')
    ax.get_xaxis().set_visible(False) 
    
    #----- log volume growth------------------------
    ax = fig.add_subplot(gs[8:10, 8:])
    logvol = ncf['stand']['logvolume'][scen,:, :]
    pulpvol = ncf['stand']['pulpvolume'][scen,:, :]
    for c in range(cols):
        yrs = len(logvol[:,c])
        ax.plot(range(1,yrs), logvol[1:,c], alpha=0.2)
    for c in range(cols):
        ax.plot(range(1,yrs), pulpvol[1:,c], alpha=0.2)
        
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title('Log and pulp volume')
    ax.get_xaxis().set_visible(False) 
    
    #-----------------leaf mass------------------
    lmass = ncf['stand']['leafmass'][scen,:, :]
    df = pd.DataFrame(data=lmass, columns=list(range(cols)))
    ax = fig.add_subplot(gs[6:8, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Leaf mass')
    ax.set_ylabel('$kg \ ha^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    
    #----- leaf mass time series------------------------
    ax = fig.add_subplot(gs[6:8, 4:8])
    for c in range(cols):
        ax.plot(lmass[:,c], alpha=0.2)
        
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title('Leaf mass')
    ax.get_xaxis().set_visible(False) 
 
    #----- leaf mass time series------------------------
    dlmass = ncf['stand']['dominant']['leafmass'][scen,:, :]
    upperlim = ncf['stand']['dominant']['leafmax'][scen,:, :]
    lowerlim = ncf['stand']['dominant']['leafmin'][scen,:, :]
    ax = fig.add_subplot(gs[6:8, 8:])
    for c in range(cols):    
        yrs = len(dlmass[:,c])
        ax.fill_between(range(1,yrs), upperlim[1:,c], lowerlim[1:,c], color='green', alpha=0.01)
        ax.plot(dlmass[:,c], alpha=0.5, color='green')
    ax.plot(dlmass[:,c], alpha=0.5, color='green', label= 'dominant')

        #ax.plot(range(1,yrs), upperlim[1:,c],  color='blue', alpha=0.01)
    dlmass = ncf['stand']['subdominant']['leafmass'][scen,:, :]
    upperlim = ncf['stand']['subdominant']['leafmax'][scen,:, :]
    lowerlim = ncf['stand']['subdominant']['leafmin'][scen,:, :]
    for c in range(cols):    
        yrs = len(dlmass[:,c])
        ax.fill_between(range(1,yrs), upperlim[1:,c], lowerlim[1:,c], color='cyan', alpha=0.01)
        ax.plot(dlmass[:,c], alpha=0.5, color='cyan')
    ax.plot(dlmass[:,c], alpha=0.5, color='cyan', label = 'subdominant')

    dlmass = ncf['stand']['under']['leafmass'][scen,:, :]
    upperlim = ncf['stand']['under']['leafmax'][scen,:, :]
    lowerlim = ncf['stand']['under']['leafmin'][scen,:, :]
    for c in range(cols):    
        yrs = len(dlmass[:,c])
        ax.fill_between(range(1,yrs), upperlim[1:,c], lowerlim[1:,c], color='orange', alpha=0.01)
        ax.plot(dlmass[:,c], alpha=0.5, color='orange')
    ax.plot(dlmass[:,c], alpha=0.5, color='orange', label = 'under')
    
    ax.legend(loc='upper left')
    
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_title('Leaf mass in canopy layers')
    ax.get_xaxis().set_visible(False) 
    

    #-----------nutrient status-------------------------
    ax = fig.add_subplot(gs[4:6, :4])
    
    ns = ncf['stand']['nut_stat'][scen,:, :]
    for c in range(cols):    
        yrs = len(ns[:,c])
        ax.plot(ns[:,c], alpha=0.5, color='orange')
    ax.plot(ns[:,c], alpha=0.5, color='orange', label = 'nutrient status')
    
    ax.legend(loc='upper left')
    
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    #ax.set_title('Nutrient status')
    ax.get_xaxis().set_visible(False) 
    


    #-----------physcal growth restrictions-------------------------
    ax = fig.add_subplot(gs[4:6, 4:8])
    dom_phys_r = (ncf['stand']['dominant']['NPP'][scen, :, :]/ncf['stand']['dominant']['NPP_pot'][scen, :, :])
    df = pd.DataFrame(data=dom_phys_r, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Physical restrictions')
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    

    #-----------nutrient status-------------------------
    ax = fig.add_subplot(gs[4:6, 8:])
    
 
    for c in range(cols):    
        yrs = len(ns[:,c])
        ax.plot(ns[:-1,c], growth[:,c], 'go', alpha=0.5)
    
    ax.legend(loc='upper left')
    
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    #ax.set_title('Nutrient status')
    ax.set_ylabel('volume growth')
    ax.set_xlabel('nutrient status')


    ndemand = ncf['stand']['n_demand'][scen,:, :]
    df = pd.DataFrame(data=ndemand, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('N demand')
    ax.set_ylabel('$kg \ ha^{-1} \ yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    pdemand = ncf['stand']['p_demand'][scen,:, :]
    df = pd.DataFrame(data=pdemand, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, 4:8])
    df.boxplot(ax = ax,
                  color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('P demand')
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    
    kdemand = ncf['stand']['k_demand'][scen,:, :]
    df = pd.DataFrame(data=kdemand, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, 8:])
    df.boxplot(ax = ax,
                  color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('K demand')
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
   
    
    ax = fig.add_subplot(gs[:2, :4])
 
    
    
    
    for c in range(cols):    
        yrs = len(ns[:,c])
        ax.plot(range(1,yrs), ndemand[1:,c],  alpha=0.5, color='blue')
    ax.plot(range(1,yrs), ndemand[1:,c],  alpha=0.5, color='blue', label='N demand')    
    ax.legend(loc='upper left')
    
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
    #ax.set_title('Nutrient status')
    ax.set_ylabel('kg $ha^{-1} yr^{-1}$', fontsize = fs)

#-----------------------------------    
    ax = fig.add_subplot(gs[:2, 4:8])
 
    for c in range(cols):    
        yrs = len(ns[:,c])
        ax.plot(range(1,yrs), pdemand[1:,c],  alpha=0.5, color='green')
    ax.plot(range(1,yrs), pdemand[1:,c],  alpha=0.5, color='green', label='P demand')
    ax.legend(loc='upper left')
    
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)
#-----------------------------------
    ax = fig.add_subplot(gs[:2, 8:])
 
    for c in range(cols):    
        yrs = len(ns[:,c])
        ax.plot(range(1,yrs), kdemand[1:,c],  alpha=0.5, color='orange')
    ax.plot(range(1,yrs), kdemand[1:,c],  alpha=0.5, color='orange', label = 'K demand')
    ax.legend(loc='upper left')
    
    ax.set_facecolor(facecolor)
    ax.tick_params(axis='y', labelsize=fs)

    
    ncf.close()


def mass(ff, scen):
    
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='peat', figsize=(15,18))   #width, height
    gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.25, hspace=0.25)
  
    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    wtmin = min(wtls) -0.2
    
    #-------water table as a reference-----------------------------
    ax = fig.add_subplot(gs[10:, :4])
    ax.plot(wtls, color='orange', label = 'late summer')
    ax.fill_between(range(cols), wtls+sdls*2, wtls-sd*2, color='orange', alpha=0.075)
    ax.hlines(y= -0.35, xmin=0, xmax = cols, color='red',linestyles='--')
    ax.get_xaxis().set_visible(False) 
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylim([wtmin,0])
    ax.set_ylabel('WT m', fontsize=fs)
    ax.legend()
    ax.grid(visible=False)
    ax.set_facecolor(facecolor)
    
    #-------mass loss in kg/ha/yr---------------------------
    litter = ncf['groundvegetation']['ds_litterfall'][scen,:, :]/10000.\
        + ncf['groundvegetation']['h_litterfall'][scen,:, :]/10000.\
        + ncf['groundvegetation']['s_litterfall'][scen,:, :]/10000.\
        + ncf['stand']['nonwoodylitter'][scen, :, :]/10000.\
        + ncf['stand']['woodylitter'][scen, :, :]/10000.

    soil = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1\
        + litter

    soilout = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1

    df = pd.DataFrame(data=soil, columns=list(range(cols)))
    ax = fig.add_subplot(gs[8:10, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Soil mass change')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
 
    #-------------end volume-------------------
    ax = fig.add_subplot(gs[10:, 4:8])
    esoms =['L0L', 'L0W', 'LL', 'LW', 'FL', 'FW', 'H', 'P1', 'P2', 'P3']
    inipeat =  np.zeros(cols)
    for sto in esoms[7:]:
        inipeat += ncf['esom']['Mass'][sto][scen,0, :]/10000.
    endpeat =  np.zeros(cols)
    for sto in esoms[7:]:
        endpeat += ncf['esom']['Mass'][sto][scen,-1, :]/10000.
    inimor =  np.zeros(cols)
    for sto in esoms[:7]:
        inimor += ncf['esom']['Mass'][sto][scen,0, :]/10000.
    endmor =  np.zeros(cols)
    for sto in esoms[:7]:
        endmor += ncf['esom']['Mass'][sto][scen,-1, :]/10000.

    
    maxval = np.max(np.vstack((inipeat+inimor, endpeat+endmor))) 
    minval = np.min(np.vstack((inipeat+inimor, endpeat+endmor)))
               
    #df = pd.DataFrame({'peat initial': inipeat, 'peat end': endpeat}, index=range(cols))
    df = pd.DataFrame({'peat initial': inipeat, 'mor initial': inimor, 'peat end': endpeat, 'mor end': endmor}, index=range(cols))
    df[['peat initial','mor initial']].plot.bar(stacked = True, position=-0.35, width=0.25, ax=ax,  colormap='copper', edgecolor='k', alpha=0.6)
    df[['peat end','mor end']].plot.bar(stacked = True, position= -1.95, width=0.25, ax=ax, colormap='copper', edgecolor='k', alpha=0.6)
    
    
    ax.set_title('Organic soil mass, kg $m^{-2}$')
    ax.set_ylim([minval*0.95, maxval*1.025])    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    ax.legend(loc='lower right', fontsize=7)

    #------Litter input --------------------------------------------
    
    df = pd.DataFrame(data=litter, columns=list(range(cols)))
    ax = fig.add_subplot(gs[6:8, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Litter input')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #-------ground vegetation biomass change  kg/ha/m2---------------------------
    gv = ncf['groundvegetation']['gv_tot'][scen,:, :]/10000.
    grgv = np.diff(gv, axis = 0)
    df = pd.DataFrame(data=grgv, columns=list(range(cols)))
    ax = fig.add_subplot(gs[4:6, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Ground vegetation biomass change')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

        

    #-------stand biomass change  kg/m2/yr---------------------------
    stand = ncf['stand']['biomass'][scen,:, :]/10000.
    gr = np.diff(stand, axis = 0)
    df = pd.DataFrame(data=gr, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Stand biomass change')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)


    #-------Site mass balance kg/m2/yr---------------------------
    site = gr + grgv + soilout[1:, :] + litter[1:,:]
    df = pd.DataFrame(data=site, columns=list(range(cols)))
    ax = fig.add_subplot(gs[:2, :4])
    df.boxplot(ax = ax,
                  color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Site mass change')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)


    #------- leaf litter kg/m2/yr---------------------------
    leaflitter =  ncf['stand']['nonwoodylitter'][scen, :, :]\
            - ncf['stand']['finerootlitter'][scen, :, :]

    df = pd.DataFrame(data=leaflitter, columns=list(range(cols)))
    ax = fig.add_subplot(gs[:2, 4:8])
    df.boxplot(ax = ax,
                  color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Leaf litter')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #------- leaf litter kg/m2/yr---------------------------
    finerootlitter =  ncf['stand']['finerootlitter'][scen, :, :]

    df = pd.DataFrame(data=finerootlitter, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, 4:8])
    df.boxplot(ax = ax,
                  color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('fineroot litter')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)


    #------- woody litter kg/m2/yr---------------------------
    woodylitter =  ncf['stand']['woodylitter'][scen, :, :]

    df = pd.DataFrame(data=woodylitter, columns=list(range(cols)))
    ax = fig.add_subplot(gs[4:6, 4:8])
    df.boxplot(ax = ax,
                  color=dict(boxes='brown', whiskers='brown', medians='brown', caps='brown'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Woody litter')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #------- ground vegetation litter kg/m2/yr---------------------------
    gvlitter =  ncf['groundvegetation']['ds_litterfall'][scen, :, :]\
        +ncf['groundvegetation']['h_litterfall'][scen, :, :]\
        +ncf['groundvegetation']['s_litterfall'][scen, :, :]

    df = pd.DataFrame(data=gvlitter, columns=list(range(cols)))
    ax = fig.add_subplot(gs[6:8, 4:8])
    df.boxplot(ax = ax,
                  color=dict(boxes='red', whiskers='red', medians='red', caps='red'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='red', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='red', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='red', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Ground vegetation litter')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #---soil mass outflux----------
    out = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1
       
    df = pd.DataFrame(data=out, columns=list(range(cols)))
    ax = fig.add_subplot(gs[8:10, 4:8])
    df.boxplot(ax = ax,
                  color=dict(boxes='grey', whiskers='grey', medians='grey', caps='grey'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='grey', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='grey', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='grey', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Soil mass loss')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)


    #----- soil mass change timeseries------------------------
    
    yrs = np.shape(ncf['esom']['Mass'][sto][scen,:, :])[0] 
    ax = fig.add_subplot(gs[10:, 8:])
    esoms =['L0L', 'L0W', 'LL', 'LW', 'FL', 'FW', 'H', 'P1', 'P2', 'P3']
    for c in range(cols):
        peat =  np.zeros(yrs)
        for sto in esoms[7:]:
            peat += ncf['esom']['Mass'][sto][scen,:, c]/10000.
        mor =  np.zeros(yrs)
        for sto in esoms[:7]:
            mor += ncf['esom']['Mass'][sto][scen,:, c]/10000.

        ax.plot(np.diff(peat), color='brown')
        ax.plot(np.diff(mor), color='orange')
    ax.plot(np.diff(peat), color='brown', label = 'peat')
    ax.plot(np.diff(mor), color='orange', label= 'mor')
    ax.legend()
    ax.set_title('Organic soil mass chagne, kg $m^{-2} yr^{-1}$')
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    ax.legend(loc='lower right', fontsize=7)



    #---soil mass outflux----------
    out = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1
       
    ax = fig.add_subplot(gs[8:10, 8:])
    for c in range(cols):

        ax.plot(out[1:,c], color='orange')
    
    ax.set_title('Soil mass loss')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)


    #---soil mass change----------
    out = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1 + litter
       
    ax = fig.add_subplot(gs[6:8, 8:])
    for c in range(cols):
        ax.plot(out[1:,c], color='red')
    
    ax.set_title('Soil mass change')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)



    #---Litter----------
       
    ax = fig.add_subplot(gs[4:6, 8:])
    for c in range(cols):
        ax.plot(litter[1:,c], color='orange')
    
    ax.set_title('Litterfall')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #---Fineroot litter----------
   
    ax = fig.add_subplot(gs[2:4, 8:])
    for c in range(cols):
        ax.plot(finerootlitter[1:,c]/10000., color='brown')
    
    ax.set_title('FinerootLitter')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #---Leaf litter----------
   
    ax = fig.add_subplot(gs[:2, 8:])
    for c in range(cols):
        ax.plot(leaflitter[1:,c]/10000., color='green')
    
    ax.set_title('LeafLitter')
    #ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    ncf.close()
   



def carbon(ff, scen):

    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='carbon', figsize=(15,18))   #width, height
    fig.suptitle('Carbon balance components', fontsize = fs+2)
    gs = gridspec.GridSpec(ncols=12, nrows=14, figure=fig, wspace=0.5, hspace=0.5)
    mass_to_c = 0.5
    
    #-------soil C balance kg/ha/yr---------------------------
    litter = (ncf['groundvegetation']['ds_litterfall'][scen,:, :]/10000.\
        + ncf['groundvegetation']['h_litterfall'][scen,:, :]/10000.\
        + ncf['groundvegetation']['s_litterfall'][scen,:, :]/10000.\
        + ncf['stand']['nonwoodylitter'][scen, :, :]/10000.\
        + ncf['stand']['woodylitter'][scen, :, :]/10000.)*mass_to_c

    soil = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1 * mass_to_c + litter
    
    out = ncf['esom']['Mass']['out'][scen,:, :]/10000.*-1 * mass_to_c

    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    wtmin = min(wtls) -0.2
    
    #-------water table as a reference-----------------------------
    ax = fig.add_subplot(gs[12:, :6])
    ax = create_profile_line(ax, wt, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, 'blue', hidex=True, hidey=False)

    #----------- water table in H------------------------------   
    elevation = np.array(ncf['strip']['elevation'])
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sd = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    h = elevation + wtls
    
    ax = fig.add_subplot(gs[12:, 6:])
    ax = create_profile_line(ax, h, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, \
                             'blue', hidex=True, hidey=False, elevation = elevation)
    
    
    #-------------LMW to Ditch----------------------------
    lmwtoditch = ncf['balance']['C']['LMWdoc_to_water'][scen,:, :] *-1
    ax = fig.add_subplot(gs[10:12, :6])
    df = pd.DataFrame(data=lmwtoditch, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'brown', 'LMW to ditch', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------HMW to Ditch----------------------------
    hmwtoditch = ncf['balance']['C']['HMW_to_water'][scen,:, :]*-1
    ax = fig.add_subplot(gs[10:12, 6:])
    df = pd.DataFrame(data=hmwtoditch, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'brown', 'HMW to ditch', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-----------LMW to atmosphere--------------------
    lmwtoatm = ncf['balance']['C']['LMWdoc_to_atm'][scen,:, :]*-1
    ax = fig.add_subplot(gs[8:10, :6])
    df = pd.DataFrame(data=lmwtoatm, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'orange', 'LMW to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------HMW to atmosphere--------------------
    hmwtoatm = ncf['balance']['C']['HMW_to_atm'][scen,:, :]*-1
    ax = fig.add_subplot(gs[8:10, 6:])
    df = pd.DataFrame(data=hmwtoatm, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'orange', 'HMW to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------CO2C to atmosphere--------------------
    co2 = ncf['balance']['C']['co2c_release'][scen,:, :]*-1
    ax = fig.add_subplot(gs[6:8, :6])
    df = pd.DataFrame(data=co2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'grey', '$CO_2C$ to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------CH4C to atmosphere--------------------
    co2 = ncf['balance']['C']['ch4c_release'][scen,:, :]*-1
    ax = fig.add_subplot(gs[6:8, 6:])
    df = pd.DataFrame(data=co2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'grey', '$CH_4C$ to atmosphere', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-----------stand litter in--------------------
    standl = ncf['balance']['C']['stand_litter_in'][scen,:, :]
    ax = fig.add_subplot(gs[4:6, :6])
    df = pd.DataFrame(data=standl, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Stand litter', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
 
    #-----------ground vegetation litter in--------------------
    gvl = ncf['balance']['C']['gv_litter_in'][scen,:, :]
    ax = fig.add_subplot(gs[4:6, 6:])
    df = pd.DataFrame(data=gvl, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Groundvegetation litter', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-----------soil balance c--------------------
    soilc = ncf['balance']['C']['soil_c_balance_c'][scen,:, :]
    ax = fig.add_subplot(gs[2:4, :6])
    df = pd.DataFrame(data=soilc, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', 'Soil C balance in C', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-----------soil balance co2 equivalents--------------------
    soilco2 = ncf['balance']['C']['soil_c_balance_co2eq'][scen,:, :]
    ax = fig.add_subplot(gs[2:4, 6:])
    df = pd.DataFrame(data=soilco2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', 'Soil C balance in $CO_2$ eq', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-----------stand balance c--------------------
    standc = ncf['balance']['C']['stand_c_balance_c'][scen,:, :]
    ax = fig.add_subplot(gs[:2, :6])
    df = pd.DataFrame(data=standc, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Stand C balance in C', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-----------stand balance co2 equivalents--------------------
    standco2 = ncf['balance']['C']['stand_c_balance_co2eq'][scen,:, :]
    ax = fig.add_subplot(gs[:2, 6:])
    df = pd.DataFrame(data=standco2, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'Stand C balance in $CO_2$ eq', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    ncf.close()
    
def nutrient_balance(ff, substance, scen):

    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num=substance, figsize=(15,18))   #width, height
    tx = substance + ' balance components'
    fig.suptitle(tx, fontsize = fs+2)
    gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.5, hspace=0.5)


    wt = np.mean(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    cols = np.shape(wt)[0]
    sd = np.std(ncf['strip']['dwtyr'][scen,:, :], axis = 0)
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sdls = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    wtmin = min(wtls) -0.2
    
    #-------water table as a reference-----------------------------
    ax = fig.add_subplot(gs[10:, :6])
    ax = create_profile_line(ax, wt, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, 'blue', hidex=True, hidey=False)

    #----------- water table in H------------------------------   
    elevation = np.array(ncf['strip']['elevation'])
    wtls = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sd = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    h = elevation + wtls
 
    ax = fig.add_subplot(gs[10:, 6:])
    ax = create_profile_line(ax, h, wtmin, sd, cols, 'WT m', 'annual', fs, facecolor, \
                             'blue', hidex=True, hidey=False, elevation = elevation)

    #-------------to Ditch----------------------------
    towater = ncf['balance'][substance]['to_water'][scen,:, :] 
    ax = fig.add_subplot(gs[8:10, :6])
    df = pd.DataFrame(data=towater, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'brown', substance + ' to water', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-------------below root layer----------------------------
    brl = ncf['balance'][substance]['decomposition_below_root_lyr'][scen,:, :] 
    ax = fig.add_subplot(gs[8:10, 6:])
    df = pd.DataFrame(data=brl, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'brown', substance + ' below root layer', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------Release in decomposition----------------------------
    de = ncf['balance'][substance]['decomposition_tot'][scen,:, :] 
    ax = fig.add_subplot(gs[6:8, :6])
    df = pd.DataFrame(data=de, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', substance + ' release in decomposition', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------Release in decomposition----------------------------
    dert = ncf['balance'][substance]['decomposition_root_lyr'][scen,:, :] 
    ax = fig.add_subplot(gs[6:8, 6:])
    df = pd.DataFrame(data=dert, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', substance + ' release in root layer', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------- Total supply ----------------------------------------    
    supply = ncf['balance'][substance]['decomposition_root_lyr'][scen,:, :] \
        + ncf['balance'][substance]['deposition'][scen,:, :]  \
        + ncf['balance'][substance]['fertilization_release'][scen,:, :] 

    ax = fig.add_subplot(gs[4:6, :6])
    df = pd.DataFrame(data=supply, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'blue', substance + ' supply', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------Release in decomposition----------------------------
    fert = ncf['balance'][substance]['fertilization_release'][scen,:, :] 
    ax = fig.add_subplot(gs[4:6, 6:])
    df = pd.DataFrame(data=fert, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'black', substance + ' release in fertilizers', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-------------Stand uptake----------------------------
    dem = ncf['balance'][substance]['stand_demand'][scen,:, :] 
    ax = fig.add_subplot(gs[2:4, :6])
    df = pd.DataFrame(data=dem, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', substance + ' stand uptake', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------ground vegetation uptake----------------------------
    dem = ncf['balance'][substance]['gv_demand'][scen,:, :] 
    ax = fig.add_subplot(gs[2:4, 6:])
    df = pd.DataFrame(data=dem, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', substance + ' groundvegetation uptake', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)

    #-------------Stand nutrient balance----------------------------
    dem = ncf['balance'][substance]['balance_root_lyr'][scen,:, :] 
    ax = fig.add_subplot(gs[:2, :6])
    df = pd.DataFrame(data=dem, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', substance + ' balance', 'kg $ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    #-------------Stand volume growth----------------------------
    vg = ncf['stand']['volumegrowth'][scen,:, :] 
    ax = fig.add_subplot(gs[:2, 6:])
    df = pd.DataFrame(data=vg, columns=list(range(cols)))
    ax = create_profile_boxplot(ax, df, cols, 'green', 'volume growth', '$m^{3} ha^{-1} yr^{-1}$', fs, facecolor, zero=False)
    
    
    ncf.close()
#***************************************************************
#***************************************************************
#              COMPARISON FIGURES
#***************************************************************

def compare_1(ff, scens):
    
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='comparison', figsize=(15,18))   #width, height
    gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.25, hspace=0.25)
    mass_to_c = 0.5
    
    wt0 = np.mean(ncf['strip']['dwtyr'][scens[0],:, :], axis = 0)
    sd0 = np.std(ncf['strip']['dwtyr'][scens[0],:, :], axis = 0)
    wtmin = min(wt0) - 0.7
    cols = np.shape(wt0)[0]
    #------------------------------
    ax = fig.add_subplot(gs[10:, :4])
    ax = create_profile_line(ax, wt0, wtmin, sd0, cols, 'WT m', 'annual', fs, facecolor, 'blue')

    wt1 = np.mean(ncf['strip']['dwtyr'][scens[1],:, :], axis = 0)
    sd1 = np.std(ncf['strip']['dwtyr'][scens[1],:, :], axis = 0)
 
    ax = fig.add_subplot(gs[10:, 4:8])
    ax = create_profile_line(ax, wt1, wtmin, sd1, cols, '','annual', fs, facecolor, 'orange')

    deltawt = ncf['strip']['dwtyr'][scens[1],:, :] - ncf['strip']['dwtyr'][scens[0],:, :]
    ax = fig.add_subplot(gs[10:, 8:])
    ax = create_profile_boxplot(ax, deltawt,cols,'green', 'WT difference', 'WT m', fs, facecolor, zero=True)
    #--------------------------------
    vol = ncf['stand']['volume'][scens[0],:, :]
    growth0 = np.diff(vol, axis=0)
    ax = fig.add_subplot(gs[8:10, :4])
    ax = create_profile_boxplot(ax, growth0, cols,'blue', 'Stand growth', 'm3ha-1yr-1', fs, facecolor)

    vol = ncf['stand']['volume'][scens[1],:, :]
    growth1 = np.diff(vol, axis=0)
    ax = fig.add_subplot(gs[8:10, 4:8])
    ax = create_profile_boxplot(ax, growth1, cols,'orange', 'Stand growth', 'm3ha-1yr-1', fs, facecolor)

    deltagr = growth1 - growth0
    ax = fig.add_subplot(gs[8:10, 8:])
    ax = create_profile_boxplot(ax, deltagr,cols,'green', 'Growth difference', 'WT m', fs, facecolor, zero=True)
    #-----------------------------------------
    hmwtoditch0 = ncf['export']['hmwtoditch'][scens[0],:, :]
    ax = fig.add_subplot(gs[6:8, :4])
    ax = create_profile_boxplot(ax, hmwtoditch0, cols,'blue', 'HMW to ditch', '', fs, facecolor)

    hmwtoditch1 = ncf['export']['hmwtoditch'][scens[1],:, :]
    ax = fig.add_subplot(gs[6:8, 4:8])
    ax = create_profile_boxplot(ax, hmwtoditch1, cols,'orange', 'HMW to ditch', '', fs, facecolor)

    deltahmw = hmwtoditch1 - hmwtoditch0
    ax = fig.add_subplot(gs[6:8, 8:])
    ax = create_profile_boxplot(ax, deltahmw,cols,'green', 'HMW difference', '', fs, facecolor, zero=True)

    #-----------------------------------------
    lmwtoditch0 = ncf['export']['lmwtoditch'][scens[0],:, :]
    ax = fig.add_subplot(gs[4:6, :4])
    ax = create_profile_boxplot(ax, lmwtoditch0, cols,'blue', 'LMW to ditch', '', fs, facecolor)

    lmwtoditch1 = ncf['export']['lmwtoditch'][scens[1],:, :]
    ax = fig.add_subplot(gs[4:6, 4:8])
    ax = create_profile_boxplot(ax, lmwtoditch1, cols,'orange', 'LMW to ditch', '', fs, facecolor)

    deltalmw = lmwtoditch1 - lmwtoditch0
    ax = fig.add_subplot(gs[4:6, 8:])
    ax = create_profile_boxplot(ax, deltalmw,cols,'green', 'LMW difference', '', fs, facecolor, zero=True)

    #-------soil C balance kg/ha/yr---------------------------
    litter0 = (ncf['groundvegetation']['ds_litterfall'][scens[0],:, :]/10000.\
        + ncf['groundvegetation']['h_litterfall'][scens[0],:, :]/10000.\
        + ncf['groundvegetation']['s_litterfall'][scens[0],:, :]/10000.\
        + ncf['stand']['nonwoodylitter'][scens[0], :, :]/10000.\
        + ncf['stand']['woodylitter'][scens[0], :, :]/10000.)*mass_to_c

    soil0 = ncf['esom']['Mass']['out'][scens[0],:, :]/10000.*-1 * mass_to_c + litter0

    ax = fig.add_subplot(gs[2:4, :4])
    ax = create_profile_boxplot(ax, soil0, cols,'blue', 'Soil C balance', '', fs, facecolor)

    litter1 = (ncf['groundvegetation']['ds_litterfall'][scens[1],:, :]/10000.\
        + ncf['groundvegetation']['h_litterfall'][scens[1],:, :]/10000.\
        + ncf['groundvegetation']['s_litterfall'][scens[1],:, :]/10000.\
        + ncf['stand']['nonwoodylitter'][scens[1], :, :]/10000.\
        + ncf['stand']['woodylitter'][scens[1], :, :]/10000.)*mass_to_c

    soil1 = ncf['esom']['Mass']['out'][scens[1],:, :]/10000.*-1 * mass_to_c + litter1

    ax = fig.add_subplot(gs[2:4, 4:8])
    ax = create_profile_boxplot(ax, soil1, cols,'orange', 'Soil C balance', '', fs, facecolor)

    deltasoil = soil1 - soil0
    ax = fig.add_subplot(gs[2:4, 8:])
    ax = create_profile_boxplot(ax, deltasoil,cols,'green', 'Soil C difference', '', fs, facecolor, zero=True)

    #------------Site C balance----------------------
    gv = ncf['groundvegetation']['gv_tot'][scens[0],:, :]/10000. * mass_to_c
    grgv0 = np.diff(gv, axis = 0)
    stand = ncf['stand']['biomass'][scens[0],:, :]/10000.* mass_to_c
    gr0 = np.diff(stand, axis = 0)
    out0 = ncf['esom']['Mass']['out'][scens[0],:, :]/10000.*-1 * mass_to_c

    site0 = gr0 + grgv0 + out0[1:, :] + litter0[1:,:]

    ax = fig.add_subplot(gs[:2, :4])
    ax = create_profile_boxplot(ax, site0, cols,'blue', 'Site C balance', '', fs, facecolor)

    gv = ncf['groundvegetation']['gv_tot'][scens[1],:, :]/10000. * mass_to_c
    grgv1 = np.diff(gv, axis = 0)
    stand = ncf['stand']['biomass'][scens[1],:, :]/10000.* mass_to_c
    gr1 = np.diff(stand, axis = 0)
    out1 = ncf['esom']['Mass']['out'][scens[1],:, :]/10000.*-1 * mass_to_c

    site1 = gr1 + grgv1 + out1[1:, :] + litter0[1:,:]

    ax = fig.add_subplot(gs[:2, 4:8])
    ax = create_profile_boxplot(ax, site1, cols,'orange', 'Site C balance', '', fs, facecolor)

    deltasoil = site1 - site0
    ax = fig.add_subplot(gs[:2, 8:])
    ax = create_profile_boxplot(ax, deltasoil,cols,'green', 'Site C difference', '', fs, facecolor, zero=True)
   

    ncf.close()


def compare_scens(ff):

    def draw_comparison(ax, x, y, sd, label, ylabel, title, color, facecolor, fs, xlabel=False):
        ax.plot(x,y, 'o-', color=color, label=label)
        ax.set_ylabel(ylabel, fontsize=fs)
        if xlabel: ax.set_xlabel('Ditch depth, m', fontsize = fs)
        ax.fill_between(x, y+sd, y-sd, color = color, alpha = 0.2)
        ax.set_facecolor(facecolor)
        ax.set_title(title, fontsize = fs)
        #ax.legend()
        return ax
    
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num='compare_scens', figsize=(12,12))   #width, height
    gs = gridspec.GridSpec(ncols=3, nrows=4, figure=fig, wspace=0.35, hspace=0.35)


    ditch_depths = ncf['scen']['ditch_depths_mean'][:] * -1
    
    grresponse = ncf['stand']['volumegrowth'][:,:, :] - ncf['stand']['volumegrowth'][0,:, :] 
    grr = np.mean(grresponse, axis=(1,2))
    grrsd = np.std(grresponse, axis=(1,2))
    ax = fig.add_subplot(gs[0,0])
    ax = draw_comparison(ax, ditch_depths, grr, grrsd, 'm3 yr-1', '$m^3 ha^{-3} yr^{-1}$', 'Growth response', 'green', facecolor, fs)
    
    wt = np.mean(ncf['strip']['dwtyr_latesummer'][:,:, :], axis = (1,2))   # mean annual wr dim: nscens
    sd = np.std(ncf['strip']['dwtyr_latesummer'][:,:, :], axis = (1,2))     #
    ax = fig.add_subplot(gs[0,1])
    ax = draw_comparison(ax, ditch_depths, wt, sd, 'water table', 'WT, m', 'Water table', 'blue', facecolor, fs)
    
    standco2bal = np.mean(ncf['balance']['C']['stand_c_balance_co2eq'][:, :, :], axis=(1,2))
    standco2balsd = np.std(ncf['balance']['C']['stand_c_balance_co2eq'][:, :, :], axis=(1,2))
    ax = fig.add_subplot(gs[0,2])
    ax = draw_comparison(ax, ditch_depths, standco2bal, standco2balsd, '', '$kg \ ha^{-1} yr^{-1}$', 'Stand $CO_2$ balance', 'grey', facecolor, fs)
    
    
    soilco2bal = np.mean(ncf['balance']['C']['soil_c_balance_co2eq'][:, :, :], axis=(1,2))
    soilco2balsd = np.std(ncf['balance']['C']['soil_c_balance_co2eq'][:, :, :], axis=(1,2))
    ax = fig.add_subplot(gs[1,0])
    ax = draw_comparison(ax, ditch_depths, soilco2bal, soilco2balsd, '', 
                         '$kg \ ha^{-1} yr^{-1}$', 'Soil $CO_2$ balance', 'grey', 
                         facecolor, fs, xlabel=True)
    
    
    ntowater = np.mean(ncf['balance']['N']['to_water'][:,:,:], axis = (1,2))
    ntowater = np.maximum(ntowater, 0.0)
    ntowatersd = np.std(ncf['balance']['N']['to_water'][:,:,:], axis = (1,2))
    ax = fig.add_subplot(gs[1,1])
    ax = draw_comparison(ax, ditch_depths, ntowater, ntowatersd, '',
                         '$kg \ ha^{-1} yr^{-1}$', 'N to water', 'red', 
                         facecolor, fs, xlabel= True)
    
    ptowater = np.mean(ncf['balance']['P']['to_water'][:,:,:], axis = (1,2))
    ptowater = np.maximum(ptowater, 0.0)
    ptowatersd = np.std(ncf['balance']['P']['to_water'][:,:,:], axis = (1,2))
    ax = fig.add_subplot(gs[1,2])
    ax = draw_comparison(ax, ditch_depths, ptowater, ptowatersd, '', 
                         '$kg \ ha^{-1} yr^{-1}$', 'P to water', 'orange', 
                         facecolor, fs,  xlabel= True)
    """
    #ax = fig.add_subplot(gs[2,:2])
    data = {'Name': ['Growth response', 'Water table', 'Stand CO2 balance', 'Soil CO2 balence', 'N export' , 'P export'], 
            '0.3 -> 0.7': [grr[2]-grr[0], wt[2]-wt[0], standco2bal[2]-standco2bal[0], soilco2bal[2]-soilco2bal[0], 
                           ntowater[2]-ntowater[0], ptowater[2]-ptowater[0] ], 
            '0.5 -> 0.7': [grr[2]-grr[1], wt[2]-wt[1], standco2bal[2]-standco2bal[1], soilco2bal[2]-soilco2bal[1], 
                           ntowater[2]-ntowater[1], ptowater[2]-ptowater[1] ], 
            '0.3 -> 0.9': [grr[3]-grr[0], wt[3]-wt[0], standco2bal[3]-standco2bal[0], soilco2bal[3]-soilco2bal[0], 
                           ntowater[3]-ntowater[0], ptowater[3]-ptowater[0] ],
            '0.5 -> 0.9': [grr[3]-grr[1], wt[3]-wt[1], standco2bal[3]-standco2bal[1], soilco2bal[3]-soilco2bal[1], 
                           ntowater[3]-ntowater[1], ptowater[3]-ptowater[1] ]}    
    df = pd.DataFrame.from_dict(data = data)
    df = df.round(decimals=2)
    df.set_index('Name', inplace=True)
    import dataframe_image as dfi
    dfi.export(df, 'table.png')
    im = plt.imread('table.png')
    ax = fig.add_subplot(gs[2:,:])
    ax.imshow(im)
    ax.axis('off')
    """
    ncf.close()


#ff = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_3_py37/outputs/susi.nc'
# ff = 'D:/Immala_simulations/Metsakeskus/immala__38304191_50vuotta_susi_lyr_0.nc'
#compare_scens(ff)
# hydrology(ff, 0)
# stand(ff, 0)
# mass(ff, 0)
# nutrient_balance(ff, 'N', 0)