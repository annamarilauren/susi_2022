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
    
    axwt0 = fig.add_subplot(gs[10:, :4])
    axwt0.plot(wt, color='blue', label = 'annual')
    axwt0.fill_between(range(cols), wt+sd*2, wt-sd*2, color='blue', alpha=0.075)
    axwt0.hlines(y= -0.35, xmin=0, xmax =cols, color='red',linestyles='--')
    axwt0.get_xaxis().set_visible(False) 
    axwt0.tick_params(axis='y', labelsize=fs)
    axwt0.set_ylim([wtmin,0])
    axwt0.set_ylabel('WT m', fontsize=fs)
    axwt0.legend()
    axwt0.grid(visible=False)
    axwt0.set_facecolor(facecolor)
    
    axwt1 = fig.add_subplot(gs[10:, 4:8])
    axwt1.plot(wtgs, color='green', label = 'growing season')
    axwt1.fill_between(range(cols), wtgs+sdgs*2, wtgs-sdgs*2, color='green', alpha=0.075)
    axwt1.hlines(y= -0.35, xmin=0, xmax = cols, color='red',linestyles='--')
    axwt1.get_xaxis().set_visible(False) 
    axwt1.get_yaxis().set_visible(False) 
    axwt1.set_ylim([wtmin,0])
    axwt1.legend()
    axwt1.grid(visible=False)
    axwt1.set_facecolor(facecolor)
    
    axwt2 = fig.add_subplot(gs[10:, 8:])
    axwt2.plot(wtls, color='orange', label = 'late summer')
    axwt2.fill_between(range(cols), wtls+sdls*2, wtls-sdls*2, color='orange', alpha=0.075)
    axwt2.hlines(y= -0.35, xmin=0, xmax = cols, color='red',linestyles='--')
    axwt2.get_xaxis().set_visible(False) 
    axwt2.get_yaxis().set_visible(False) 
    axwt2.set_ylim([wtmin,0])
    axwt2.legend()
    axwt2.grid(visible=False)
    axwt2.set_facecolor(facecolor)
    
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
    axdeltas = fig.add_subplot(gs[2:4, :4])
    dfdeltas.boxplot(ax = axdeltas,
                 color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                 boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                 flierprops=dict(linestyle='-', linewidth=1.5),
                 medianprops=dict(linestyle='-', linewidth=1.5),
                 whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                 capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                 showfliers=False, grid=False, rot=0)
    
    axdeltas.set_title('Through soil surface')
    axdeltas.set_ylabel('Water flux mm $yr^{-1}$', fontsize=fs)
    
    axdeltas.get_xaxis().set_visible(False)
    axdeltas.tick_params(axis='y', labelsize=fs)
    axdeltas.set_facecolor(facecolor)
    
    #------ETs-----------------
    ET = ncf['cpy']['ET_yr'][scen,1:, :]*1000.
    dfET = pd.DataFrame(data=ET, columns=list(range(cols)))
    axET = fig.add_subplot(gs[2:4, 4:8])
    dfET.boxplot(ax = axET,
                 color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                 boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                 flierprops=dict(linestyle='-', linewidth=1.5),
                 medianprops=dict(linestyle='-', linewidth=1.5),
                 whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                 capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                 showfliers=False, grid=False, rot=0)
    
    axET.set_title('ET')
    
    axET.get_xaxis().set_visible(False)
    axET.tick_params(axis='y', labelsize=fs)
    axET.set_facecolor(facecolor)
    
    
    #------transpi-----------------
    transpi = ncf['cpy']['transpi_yr'][scen,1:, :]*1000.
    dftranspi = pd.DataFrame(data=transpi, columns=list(range(cols)))
    axtranspi = fig.add_subplot(gs[2:4, 8:])
    dftranspi.boxplot(ax = axtranspi,
                 color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                 boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                 flierprops=dict(linestyle='-', linewidth=1.5),
                 medianprops=dict(linestyle='-', linewidth=1.5),
                 whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                 capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                 showfliers=False, grid=False, rot=0)
    
    axtranspi.set_title('Transpiration')
    
    axtranspi.get_xaxis().set_visible(False)
    axtranspi.tick_params(axis='y', labelsize=fs)
    axtranspi.set_facecolor(facecolor)
    
    
    #------efloor-----------------
    efloor = ncf['cpy']['efloor_yr'][scen,1:, :]*1000.
    dfefloor = pd.DataFrame(data=efloor, columns=list(range(cols)))
    axefloor = fig.add_subplot(gs[:2, :4])
    dfefloor.boxplot(ax = axefloor,
                 color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                 boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                 flierprops=dict(linestyle='-', linewidth=1.5),
                 medianprops=dict(linestyle='-', linewidth=1.5),
                 whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                 capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                 showfliers=False, grid=False, rot=0)
    
    axefloor.set_title('Soil evaporation')
    axefloor.set_ylabel('Water flux mm $yr^{-1}$', fontsize=fs)
    
    axefloor.get_xaxis().set_visible(False)
    axefloor.tick_params(axis='y', labelsize=fs)
    axefloor.set_facecolor(facecolor)
    
    #------SWE max-----------------
    swe = ncf['cpy']['SWEmax'][scen,1:, :]
    dfswe = pd.DataFrame(data=swe, columns=list(range(cols)))
    axswe= fig.add_subplot(gs[:2, 4:8])
    dfswe.boxplot(ax = axswe,
                 color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                 boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                 flierprops=dict(linestyle='-', linewidth=1.5),
                 medianprops=dict(linestyle='-', linewidth=1.5),
                 whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                 capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                 showfliers=False, grid=False, rot=0)
    
    axswe.set_title('Max snow water equivalent')
    
    axswe.get_xaxis().set_visible(False)
    axswe.tick_params(axis='y', labelsize=fs)
    axswe.set_facecolor(facecolor)
    
    #------Interc max-----------------
    interc = ncf['cpy']['interc_yr'][scen,1:, :]*1000.
    dfinterc = pd.DataFrame(data=interc, columns=list(range(cols)))
    axinterc= fig.add_subplot(gs[:2, 8:])
    dfinterc.boxplot(ax = axinterc,
                 color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                 boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                 flierprops=dict(linestyle='-', linewidth=1.5),
                 medianprops=dict(linestyle='-', linewidth=1.5),
                 whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                 capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                 showfliers=False, grid=False, rot=0)
    
    axinterc.set_title('Mean interception storage')
    
    axinterc.get_xaxis().set_visible(False)
    axinterc.tick_params(axis='y', labelsize=fs)
    axinterc.set_facecolor(facecolor)
    
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
    gs = gridspec.GridSpec(ncols=12, nrows=12, figure=fig, wspace=0.25, hspace=0.25)
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
    
    df = pd.DataFrame(data=soil, columns=list(range(cols)))
    ax = fig.add_subplot(gs[2:4, :6])
    df.boxplot(ax = ax,
                  color=dict(boxes='blue', whiskers='blue', medians='blue', caps='blue'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='blue', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Soil carbon balance')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
    #-------------------------
    
    #-------Site carbon balance kg/m2/yr---------------------------
    gv = ncf['groundvegetation']['gv_tot'][scen,:, :]/10000. * mass_to_c
    grgv = np.diff(gv, axis = 0)
    stand = ncf['stand']['biomass'][scen,:, :]/10000.* mass_to_c
    gr = np.diff(stand, axis = 0)


    site = gr + grgv + out[1:, :] + litter[1:,:]
    df = pd.DataFrame(data=site, columns=list(range(cols)))
    ax = fig.add_subplot(gs[:2, :6])
    df.boxplot(ax = ax,
                  color=dict(boxes='green', whiskers='green', medians='green', caps='green'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='green', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Site carbon balance')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    #-------Soil carbon emissions kg/m2/yr---------------------------
    cemiss = ncf['esom']['Mass']['c_to_atm'][scen,:, :]/10000.*-1 
    df = pd.DataFrame(data=cemiss, columns=list(range(cols)))
    ax = fig.add_subplot(gs[4:6, :6])
    df.boxplot(ax = ax,
                  color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Soil carbon emissions')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
  
    #-------Soil methane emissions kg/m2/yr---------------------------    
    
    ch4 = ncf['methane']['ch4'][scen,:, :]/10000.
        
    
    df = pd.DataFrame(data=ch4, columns=list(range(cols)))
    ax = fig.add_subplot(gs[6:8, :6])
    df.boxplot(ax = ax,
                  color=dict(boxes='magenta', whiskers='magenta', medians='magenta', caps='magenta'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='magenta', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='magenta', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='magenta', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Methane')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    ax.hlines(y= 0, xmin=0, xmax = cols, color='red',linestyles='--')
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
  
    
    #-------peat DOC relase kg/m2/yr---------------------------
    doc = ncf['esom']['Mass']['doc'][scen,:, :]/10000. *-1
    df = pd.DataFrame(data=doc, columns=list(range(cols)))
    ax = fig.add_subplot(gs[8:10, :6])
    df.boxplot(ax = ax,
                  color=dict(boxes='brown', whiskers='brown', medians='brown', caps='brown'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('DOC release')
    ax.set_ylabel('$kg m^{-2} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    
    #-------water table as a reference-----------------------------
    ax = fig.add_subplot(gs[10:, :6])
    ax.plot(wt, color='blue', label = 'late summer')
    ax.fill_between(range(cols), wt+sd*2, wt-sd, color='blue', alpha=0.075)
    ax.hlines(y= -0.35, xmin=0, xmax = cols, color='red',linestyles='--')
    ax.get_xaxis().set_visible(False) 
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylim([wtmin,0])
    ax.set_ylabel('WT m', fontsize=fs)
    ax.legend()
    ax.grid(visible=False)
    ax.set_facecolor(facecolor)
    
    #----------- water table in H------------------------------   
    wt = np.mean(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    sd = np.std(ncf['strip']['dwtyr_latesummer'][scen,:, :], axis = 0)
    ax = fig.add_subplot(gs[10:, 6:])

    elevation = np.array(ncf['strip']['elevation'])
    ax.plot(elevation, color='brown', label = 'soil surface')
    ax.plot(elevation + wt, color='blue', label = 'late summer')
    ax.fill_between(range(cols), elevation + wt+sd*2, elevation + wt-sd, color='blue', alpha=0.075)
    ax.plot(range(cols), elevation-0.35, color='red')
    ax.get_xaxis().set_visible(False) 
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylabel('WT m', fontsize=fs)
    ax.legend()
    ax.grid(visible=False)
    ax.set_facecolor(facecolor)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
 
    
    #-------residence time days-----------------------------
    timetoditch = ncf['strip']['residencetime'][scen,:, :]
    ax = fig.add_subplot(gs[8:10, 6:])
    df = pd.DataFrame(data=timetoditch, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes='brown', whiskers='brown', medians='brown', caps='brown'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='brown', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('Residence time')
    ax.set_ylabel('days', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)
 

    hmw = ncf['esom']['Mass']['doc'][scen,:, :] - ncf['esom']['Mass']['lmwdoc'][scen,:, :]
    ax = fig.add_subplot(gs[6:8, 6:])
    df = pd.DataFrame(data=hmw, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes='black', whiskers='black', medians='black', caps='black'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='black', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='black', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='black', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('HMW DOC')
    ax.set_ylabel('kg $ha^{-1} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    hmwtoditch = ncf['export']['hmwtoditch'][scen,:, :]
    yrs, cols = np.shape(hmw)
    
    ax = fig.add_subplot(gs[4:6, 6:])
    df = pd.DataFrame(data=hmwtoditch, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes='red', whiskers='red', medians='red', caps='red'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='red', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='red', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='red', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('HMW to ditch')
    ax.set_ylabel('kg $ha^{-1} yr^{-1}$', fontsize=fs)
    ax.set_facecolor(facecolor)
    ax.get_xaxis().set_visible(False)
    
    
    lmw = ncf['esom']['Mass']['lmwdoc'][scen,:, :]
    ax = fig.add_subplot(gs[2:4, 6:])
    df = pd.DataFrame(data=lmw, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes='cyan', whiskers='cyan', medians='cyan', caps='cyan'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='cyan', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='cyan', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='cyan', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('LMW DOC')
    ax.set_ylabel('kg $ha^{-1} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    lmwtoditch = ncf['export']['lmwtoditch'][scen,:, :]

    ax = fig.add_subplot(gs[:2, 6:])
    df = pd.DataFrame(data=lmwtoditch, columns=list(range(cols)))
    df.boxplot(ax = ax,
                  color=dict(boxes='orange', whiskers='orange', medians='orange', caps='orange'),
                  boxprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  flierprops=dict(linestyle='-', linewidth=1.5),
                  medianprops=dict(linestyle='-', linewidth=1.5),
                  whiskerprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  capprops=dict(linestyle='-', linewidth=1.5, color='orange', alpha=0.6),
                  showfliers=False, grid=False, rot=0)
    
    ax.set_title('LMW to ditch')
    ax.set_ylabel('kg $ha^{-1} yr^{-1}$', fontsize=fs)
    
    ax.get_xaxis().set_visible(False)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_facecolor(facecolor)

    ncf.close()

