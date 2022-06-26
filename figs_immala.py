# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 18:56:24 2022

@author: alauren
"""


import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset  
import numpy as np
import pandas as pd
import dataframe_image as dfi
import glob
import os
from os import listdir

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
    
    fname = ff.split('/')[-1].split('.')[0].split('_')[2]
    
    ncf=Dataset(ff, mode='r')                                        # water netCDF, open in reading mode
    facecolor = '#f2f5eb'
    fs = 15
    fig = plt.figure(num=fname, figsize=(12,12))   #width, height
    plt.suptitle('Kuvio ' + str(fname), fontsize=fs+2)
    
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

    dfi.export(df, 'table.png')
    im = plt.imread('table.png')
    ax = fig.add_subplot(gs[2:,:])
    ax.imshow(im)
    ax.axis('off')

    ncf.close()
    outfolder =r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Immala/figs/'
    plt.savefig(outfolder + fname )

#ifiles = glob.glob(r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Immala/susi_22_out/susi_out_1/' +"*.nc")
ifolder = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/Immala/susi_22_out/susi_out_1/'

for ifile in listdir(ifolder):
    ff = os.path.join(ifolder, ifile)
    compare_scens(ff)
    print (ifile)
