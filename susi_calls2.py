# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:10:42 2020

@author: alauren
"""
import numpy as np
import pandas as pd
import datetime
from susi_utils import  get_motti, read_FMI_weather, get_mese_input
from susi_para import get_susi_para
from susi83 import run_susi
import susi_io
#%%
def call_local_susi():
    #from dwts_para import para
    """
    single run, all input here
    switches in susi83: Opt_strip= True, no output, only figs
    """    
    #***************** local call for SUSI*****************************************************
    folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/' #'sensitivity/'
    susiPath = r'C:/Users/alauren/Documents/Susi_9/'
    wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/input_aikakauskirja/'
    mottipath =  r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/input_aikakauskirja/'
   
    mf='Rovaniemi_Mtkg_manty.xls'
    wdata='Rovaniemi.csv'

    start_date = datetime.datetime(2005,1,1)
    end_date=datetime.datetime(2006,12,31)
    start_yr = start_date.year 
    end_yr = end_date.year
    yrs = (end_date - start_date).days/365.25
    mottifile = mottipath + mf
    df = get_motti(mottifile)
    sfc =  3                                                                         #soil fertility class
    ageSim = 48.       #Lohja 31; 48 Rovaniemi 46 100  
    sarkaSim = 40. 
    n = int(sarkaSim / 2)
            
    site = 'develop_scens'
    
    forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
    
            
    wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                              folderName=folderName, hdomSim=None,  
                                                                              ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                              susiPath=susiPath,
                                                                              n=n)
                                                                              
    v_ini, v, iv5,  cbt, dcbt, cb, dcb, w,dw,logs,pulp, dv,dlogs,dpulp,yrs, bmgr,  \
                                    Nleach, Pleach, Kleach, DOCleach, runoff, \
                                    nrelease, prelease,krelease, ch4release  = run_susi(forc, wpara, cpara, 
                                    org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'other', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath)
    
    print ( '************************************************')
    
    outf = r'C:/Users/alauren/OneDrive - University of Eastern Finland\Susi/input_aikakauskirja/'
    viniarr = np.ones(len(v))*v_ini
    dout = {'ditch depth':spara['ditch depth east'],'v_ini':viniarr,'end volume': v, 'annual growrth': iv5,
        'carbon balance': cbt, 'water table': w, 'log volume': logs, 'pulp vol': pulp, 'Nleach': Nleach,
        'Pleach':Pleach, 'Kleach': Kleach, 'DOCleach': DOCleach, 'runoff':runoff}
    dfout = pd.DataFrame(data=dout)
    dfout = dfout.set_index('ditch depth')
    dfout.to_excel(outf + 'standout.xlsx')
    
            
             
call_local_susi()
#%%


#%%
def call_vesitase_susi():
    from dwts_para import para
    """
    signle runs, input from water balance trials
    switches in susi83: Opt_strip= True, vesitase=True
    """    
    #***************** local call for SUSI*****************************************************
    folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/'#'sensitivity/'
    susiPath ='C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/'
    wpath = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/susi_experim/weather/'
    mottipath =  r'C:/Users/alauren/OneDrive - University of Eastern Finland/Stenberg/motti_Arille/'

    sites =['ansa21', 'ansa26', 'jaakkoin61', 'jaakkoin62', 'koira11', 'koira12',
        'neva11', 'neva14', 'neva31', 'neva34','parkano11']
 
    out = {}    
    params = para(period='start-end')
    #for s in params.keys():
    #for s in ['koira11']:
    for s in sites:
        if not params[s]['thinning']: 
            print (s)
            start_date = params[s]['start_date']
            end_date=params[s]['end_date']
            #end_date = datetime.datetime(2009,12,31)
            start_yr = start_date.year; end_yr = end_date.year
            yrs = (end_date - start_date).days/365.25
            length = (end_date - start_date).days +1
            mottifile = mottipath + params[s]['mottifile']
            df = get_motti(mottifile)
            wdata  =  params[s]['wfile']
            sfc =  params[s]['sfc']                                                                         #soil fertility class
            ageSim=  params[s]['Aini'] #age #         
        
            sarkaSim = params[s]['Swidth']
            n = int(sarkaSim / 2)
            ddwest = params[s]['ddepth']
            ddeast = params[s]['ddepth'] 
            bd = params[s]['bulk dens']/1000. #0.14 #idata.T[nro]['bd']
            
            #site = 'develop_scens'
            site = 'develop'               
            
            forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
            peatN = None #idata.T[nro]['n_mg/g']/10.  
            peatP = None #idata.T[nro]['p_mg/g']/10.
            peatK = None# idata.T[nro]['k_mg/g']/10.
            kaista = None #idata.T[nro]['kaista']
            
            #susi_io.weather_fig(forc)
            #susi_io.motti_fig(df, ageSim, yrs, sfc)
            wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                              folderName=folderName, hdomSim=None,  
                                                                              ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                              susiPath=susiPath,
                                                                              ddwest=ddwest, ddeast=ddeast, n=n, bd=bd,
                                                                              peatN=peatN, peatP=peatP, peatK=peatK)
            spara['scenario name']=[s]
            spara['vonP top'] = params[s]['vonP']
            
            spara['peat type'] = params[s]['ptype']
            #spara['peat type'] = ['S','S','S','S','S','S','S','S']
            #spara['peat type bottom']=['S']

            #spara['peat type'] = ['A','A','A','A','A','A','A','A']
            #spara['peat type bottom']=['A']

            spara['depoN'] = params[s]['depoN']       #
            spara['depoP'] = params[s]['depoP']       #+dry deposition Ruoho-Airola et al 2015
            spara['depoK'] = params[s]['depoK']*1.1   #+dry deposition  Ruoho-Airola et al 2003 
            
            # v_ini, v, iv5,  cbt, dcbt, cb, dcb, w,dw,logs,pulp, dv,dlogs,dpulp,yrs, bmgr, pot_gr, phys_r, chem_r, \
            #                         Nleach, Pleach, Kleach, DOCleach, HMWleach, runoff, N_gr, P_gr, K_gr 

            v_ini, v, iv5,  cbt, dcbt, cb, dcb, w,dw,logs,pulp, dv,dlogs,dpulp,yrs, bmgr,  \
                Nleach, Pleach, Kleach, DOCleach,  runoff, nrelease, prelease, krelease, ch4release  = run_susi(forc, wpara, cpara, 
                                    org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'other', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath, 
                                    kaista=kaista, sitename = s)
            
                                                        
            
            out[s]={} 
            out[s]['v_ini']=v_ini
            out[s]['v']=v[0]
            out[s]['iv5']=iv5[0]            
            out[s]['bmgr']=bmgr[0]
            #out[s]['pot_gr'] = pot_gr
            #out[s]['phys_r'] = phys_r
            #out[s]['chem_r'] = chem_r
            out[s]['wt'] = np.mean(w)
            out[s]['Nleach'] = Nleach[0]
            out[s]['Pleach'] = Pleach[0]
            out[s]['Kleach'] = Kleach[0]
            out[s]['DOCleach'] = DOCleach[0]
            #out[s]['HMWleach'] = HMWleach
            out[s]['runoff'] = runoff[0]
            out[s]['cbt'] = cbt[0]
            out[s]['cb'] = cb[0]
            
                       
            
            print ( '************************************************')
            print (v_ini, v, iv5, Kleach)
    dfOut = pd.DataFrame.from_dict(out)
    dfOut=dfOut.T
    print (dfOut)
    files =['vols_swidth+20.xlsx', 'vols_swidth-20.xlsx', 'vols_ddepth+20.xlsx', 'vols_ddepth-20.xlsx', 
        'vols_nroot_30.xlsx', 'vols_nroot_50.xlsx',
        'vols_afp_25.xlsx', 'vols_afp_5.xlsx', 'vols_nut+20.xlsx', 'vols_nut-20.xlsx', 'vols_anisotropy+20.xlsx',
        'vols_anisotropy-20.xlsx', 'vols_sphagnum.xlsx', 'vols_carex.xlsx']

    dfOut.to_excel(folderName+'vols.xlsx')
    #dfOut.to_excel(folderName+files[13])
    
call_vesitase_susi()

#%%
def susi_mese():
    folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/outputs/'
    susiPath = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/'
    wpath = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_7_0_py27/wfilesmese/'
    mottipath = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_8_4_py3/MESE-motti/'
    #demfile= r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/Susi_7_0_py27/nutdem.xls'
    idata = get_mese_input(susiPath + 'mese_input.xlsx')
    
    start_date = datetime.datetime(1980,1,1); end_date=datetime.datetime(1984,12,31)
    start_yr = start_date.year; end_yr = end_date.year
    yrs = (end_date - start_date).days/365.25
    length = (end_date - start_date).days +1
    
    #site = 'develop_scens'
    site = 'develop'
    
    for nro in range(1,208):  #208
            
        mottifile = mottipath + 'motti_'+str(nro)+   '.xls'
        df = get_motti(mottifile)
        
        wdata = idata.T[nro]['wfile']
        forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
        
        sfc =  idata.T[nro]['sfc']                                                                         #soil fertility class
        ageSim=  idata.T[nro]['age_ini']                                                                     #90
        sarkaSim = idata.T[nro]['stripw']  
        n = int(sarkaSim / 2)
        ddwest = -idata.T[nro]['dd_west']/100.
        ddeast = -idata.T[nro]['dd_east']/100.
        bd = idata.T[nro]['bd']
        peatN = idata.T[nro]['n_mg/g']/10.  
        peatP = idata.T[nro]['p_mg/g']/10.
        peatK = idata.T[nro]['k_mg/g']/10.
        kaista = idata.T[nro]['kaista']
        
        
        #susi_io.weather_fig(forc)
        #susi_io.motti_fig(df, ageSim, yrs, sfc)
        wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=susiPath,
                                                                          ddwest=ddwest, ddeast=ddeast, n=n, bd=bd,
                                                                          peatN=peatN, peatP=peatP, peatK=peatK)

        spara['depoN'] = idata.T[nro]['depoN']       #
        spara['depoP'] = idata.T[nro]['depoP']       #+dry deposition Ruoho-Airola et al 2015
        spara['depoK'] = idata.T[nro]['depoK']*1.1   #+dry deposition  Ruoho-Airola et al 2003 
    
        runopt = True  
        if runopt:
            if len(spara['scenario name']) == 1:
                
                v_ini, v, iv5, Nrel, Prel, Krel, Crel, dwt_loc,cb, cbt, runoff, drunoff, w, dw, dv, \
                                        gr_pot, N_gr, P_gr, K_gr, phys_gr = run_susi(forc, wpara, cpara, 
                                        org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                        mottifile=mottifile, peat= 'other', photosite='All data', 
                                        folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath, kaista=kaista)

                fout = susiPath + 'meseout_2021.xls'

                print (nro, '************************************************')
                
                
                susi_io.write_mese(fout, nro, v_ini, v, iv5, np.mean(Nrel), \
                                   np.mean(Prel), np.mean(Krel), np.mean(Crel), \
                                   np.mean(dwt_loc), np.mean(cb), np.mean(cbt), spara['sfc'], gr_pot, N_gr, P_gr, K_gr, phys_gr)
                
        
    
            else:    #
            
                v_ini, v_end, gr, w, dw, dv = run_susi(forc, wpara, cpara, 
                                        org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                        mottifile=mottifile, peat= 'other', photosite='All data', 
                                        folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath, kaista=kaista)



                fout = susiPath + 'mese_scen_out_2021.xls'            
                #insert runoof and drunoff to call
                susi_io.write_mese_scen(fout, nro, v_ini, v_end, gr,  w, dw)
    
                print (nro, '************************************************')
                print ('         + Vini', v_ini)
                print ('         + V_end',v_end) 
                print ('         + annual gr',gr)
                #print ('         + annual carbon b',cb) 
                #print ('         + change in cb', dcb)
                #print ('         + annual carbon b trees',cbt) 
                #print ('         + change in cbt', dcbt)
                print ('         + grow s wt',w)
                print ('         + change in wt',dw)
    

#susi_mese()
#%%                
def susi_jaali():

    #***************** local call for SUSI*****************************************************
    folderName=r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/prepare_susi/outs/' #sensitivity/'
    susiPath ='C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/prepare_susi/'
    wpath = r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/prepare_susi//weather/'
    mottipath =  r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/prepare_susi/mottiouts/'

    infold =r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/prepare_susi/'
    f = r'test_vals.csv'
    dfin = pd.read_csv(infold+f, sep=';', usecols =['ID','h', 'Ahat', 'BA','N','d', 'V', 'xlname', 'spe', 'sfc', 'swidth'])
    rows, cols = np.shape(dfin)
    print (dfin.loc[0:2])

    out = {}    
    for r in range(569, 949): #range(rows):    
        print ('++++++++++++++++++++++++++')
        print ('ROW', r)
        print ('++++++++++++++++++++++++++')
        
        if dfin['V'][r]<1.0:
            continue
        start_date = datetime.datetime(2000,1,1) 
        end_date=datetime.datetime(2019,12,31)
        start_yr = start_date.year; end_yr = end_date.year
        yrs = (end_date - start_date).days/365.25
        length = (end_date - start_date).days +1
        mottifile = mottipath + dfin['xlname'][r]
        df = get_motti(mottifile)
        wdata  =  'weather_oulu_2493_filled.csv'
        sfc = dfin['sfc'][r]                                                                          #soil fertility class
        ageSim=  dfin['Ahat'][r] #age #         
        sarkaSim = dfin['swidth'][r] 
        n = int(sarkaSim / 2)
        site = 'jaali'               
        
        forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
        
        peatN = None #idata.T[nro]['n_mg/g']/10.  
        peatP = None #idata.T[nro]['p_mg/g']/10.
        peatK = None# idata.T[nro]['k_mg/g']/10.
        kaista = 1 #idata.T[nro]['kaista']
        
        #susi_io.weather_fig(forc)
        #susi_io.motti_fig(df, ageSim, yrs)
        wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=susiPath, n=n, 
                                                                          peatN=peatN, peatP=peatP, peatK=peatK)
        #spara['scenario name']=[dfin['ID'][r]]
        
        single = False
        if single:        
            v_ini, v, iv5,  cbt, dcbt, cb, dcb, w,dw,logs,pulp, dv,dlogs,dpulp,yrs, bmgr, pot_gr, phys_r, chem_r, \
                                    Nleach, Pleach, Kleach, DOCleach, runoff = run_susi(forc, wpara, cpara, 
                                    org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'other', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath, kaista=kaista)
        else:
            v_ini, v, iv5, Nrel, Prel, Krel, Crel, dwt_loc,cb, cbt, runoff, drunoff, wt, dwt, dv = run_susi(forc, wpara, cpara, 
                                    org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                    mottifile=mottifile, peat= 'other', photosite='All data', 
                                    folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath, kaista=kaista)


        print ( '************************************************')
        print (v_ini, v, iv5, runoff, drunoff, wt, dwt, dv)

        fout = susiPath + 'jaali_scen_out.xls'            
        #insert runoof and drunoff to call
        susi_io.write_jaali_scen(fout, r, dfin['ID'][r], v_ini, v, iv5, cb, cbt, wt, dwt, runoff, drunoff)


    """
        out[r]={} 
        out[r]['v_ini']=v_ini
        out[r]['v']=v[0]
        out[r]['iv5']=iv5            
        out[r]['bmgr']=bmgr[0]
        out[r]['pot_gr'] = pot_gr
        out[r]['phys_r'] = phys_r
        out[r]['chem_r'] = chem_r
        out[r]['wt'] = np.mean(w)
        out[r]['Nleach'] = Nleach
        out[r]['Pleach'] = Pleach
        out[r]['Kleach'] = Kleach
        out[r]['DOCleach'] = DOCleach
        out[r]['runoff'] = runoff
                   
        
    dfOut = pd.DataFrame.from_dict(out)
    dfOut=dfOut.T
    print (dfOut)
    #dfOut.to_excel(folderName+'vols.xlsx')    
    """
#susi_jaali() 
#%%     
 
"""
single run, all input here
switches in susi83: Opt_strip= True, no output, only figs
"""    
#***************** local call for SUSI*****************************************************
folderName= r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/output/' #'sensitivity/'
susiPath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/output/' #r'C:/Users/alauren/Documents/Susi_9/'
wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/wfiles/'
mottipath =  r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/mottifiles/'

mf='susi_stand_koirasuo11_nofert.xlsx'

wdata='koirasuo_weather.csv'

start_date = datetime.datetime(2014,1,1)
end_date=datetime.datetime(2020,12,31)
start_yr = start_date.year 
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25
mottifile = mottipath + mf
df = get_motti(mottifile)
sfc =  2                                                                         #soil fertility class
ageSim = 94.         
sarkaSim = 37. 
n = int(sarkaSim / 2)
        
site = 'develop_scens'

forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
            
wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=susiPath,
                                                                          n=n)

outpara['hydfig']=True                  # produces a figure for the each scenario of the hydrology
outpara['to_file']=True                 # saves csv and excel outputs to susi-folder           
outpara['figs']=True                    # figure for comparing growth and water tables

spara['L']= 37.0                        # strip width, m 25...55 
 
spara['depoK']= 0.67                     # atmospheric deposition kg / ha / yr
spara['depoN']= 4.0                       
spara['depoP']= 0.1

#++++++ Ash fertilization +++++++++++++++++++++++++++++
spara['fertilization']['application year'] = 2014         # year of application, check consistency with simulation time
spara['fertilization']['P']['dose'] = 0.0 #45.0                 # kg / ha of elemental P 
spara['fertilization']['K']['dose'] = 0.0 #100.0               # kg / ha of elemental K

#+++++++ These go together, same amount of entries in the list
spara['ditch depth east']= [-0.85]       # ditch depths, m. The first one is the reference depth
spara['ditch depth west']= [-0.85]
spara['ditch depth 20y east']= [-0.85]     # ditch depth after 20 yrs
spara['ditch depth 20y west']= [-0.85]
spara['scenario name']= ['koira']    # give the scenario names    
v_ini, v, iv5,  cbt, dcbt, cb, dcb, w,dw,logs,pulp, dv,dlogs,dpulp,yrs, bmgr,  \
                                Nleach, Pleach, Kleach, DOCleach, runoff, \
                                nrelease, prelease,krelease, ch4release = run_susi(forc, wpara, cpara, 
                                org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath)
    
viniarr = np.ones(len(v))*v_ini
dout = {'ditch depth':spara['ditch depth east'],'initial volume': [v_ini]*len(spara['ditch depth east']), 'end volume': v, 'annual growrth': iv5,
        'stand C bal': cbt, 'peat C bal': cb, 'ch4 emiss': ch4release,
        'water table': w, 'log volume': logs, 'pulp vol': pulp, 'Nleach': Nleach }
dfout = pd.DataFrame(data=dout)
dfout = dfout.set_index('ditch depth')
dfout.to_excel('standout.xlsx')
dfout
#%%
def call_ash_susi():
    from dwts_para import para_ash
    """
    signle runs, input from water balance trials
    switches in susi83: Opt_strip= True, vesitase=True
    """    
    #***************** local call for SUSI*****************************************************
    folderName= r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/output/' #'sensitivity/'
    susiPath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/output/' #r'C:/Users/alauren/Documents/Susi_9/'
    wpath = r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/wfiles/'
    mottipath =  r'C:/Users/alauren/OneDrive - University of Eastern Finland/Susi/ash_article/mottifiles/'

    sites =['ansa_cont', 'ansa_fert_0','ansa_fert_1', 'ansa_fert_2', 'ansa_fert_3',
            'koira_cont', 'koira_fert_0', 'koira_fert_1', 'koira_fert_2', 'koira_fert_3',
            'neva_cont', 'neva_fert_0', 'neva_fert_1', 'neva_fert_2', 'neva_fert_3']
 
    out = {}    
    params = para_ash()
    for s in sites:
        print (s)
        start_date = params[s]['start_date']
        end_date=params[s]['end_date']
        start_yr = start_date.year; end_yr = end_date.year
        yrs = (end_date - start_date).days/365.25
        length = (end_date - start_date).days +1
        mottifile = mottipath + params[s]['mottifile']
        df = get_motti(mottifile)
        wdata  =  params[s]['wfile']
        sfc =  params[s]['sfc']                                                                         #soil fertility class
        ageSim=  params[s]['Aini'] #age #         
    
        sarkaSim = params[s]['Swidth']
        n = int(sarkaSim / 2)
        ddwest = params[s]['ddepth']
        ddeast = params[s]['ddepth'] 
        bd = params[s]['bulk dens']/1000. #0.14 #idata.T[nro]['bd']
        
        #site = 'develop_scens'
        site = 'develop'               
        
        forc=read_FMI_weather(0, start_date, end_date, sourcefile=wpath+wdata)           # read weather input
        
        peatN = None #idata.T[nro]['n_mg/g']/10.  
        peatP = None #idata.T[nro]['p_mg/g']/10.
        peatK = None# idata.T[nro]['k_mg/g']/10.
        kaista = None #idata.T[nro]['kaista']
        
        #susi_io.weather_fig(forc)
        #susi_io.motti_fig(df, ageSim, yrs, sfc)
        wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation='undefined', peat=site, 
                                                                          folderName=folderName, hdomSim=None,  
                                                                          ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, 
                                                                          susiPath=susiPath,
                                                                          ddwest=ddwest, ddeast=ddeast, n=n, bd=bd,
                                                                          peatN=peatN, peatP=peatP, peatK=peatK)
        spara['scenario name']=[s]
        spara['vonP top'] = params[s]['vonP']
        spara['peat type'] = params[s]['ptype']
        #spara['peat type'] = ['S','S','S','S','S','S','S','S']
        #spara['peat type bottom']=['S']

        #spara['peat type'] = ['A','A','A','A','A','A','A','A']
        #spara['peat type bottom']=['A']

        spara['depoN'] = params[s]['depoN']       #
        spara['depoP'] = params[s]['depoP']       #+dry deposition Ruoho-Airola et al 2015
        spara['depoK'] = params[s]['depoK']*1.1   #+dry deposition  Ruoho-Airola et al 2003 
        
        spara['fertilization']['application year'] = params[s]['fertilization']['application year']         # year of application, check consistency with simulation time
        spara['fertilization']['P']['dose'] =  params[s]['fertilization']['P']['dose']               # kg / ha of elemental P 
        spara['fertilization']['K']['dose'] = params[s]['fertilization']['K']['dose']               # kg / ha of elemental K

        photopara['beta'] = photopara['beta'] * params[s]['LUE_multiplier']

        
        outpara['hydfig']=False                  # produces a figure for the each scenario of the hydrology
        outpara['to_file']=True                  # saves csv and excel outputs to susi-folder           
        outpara['figs']=False                    # figure for comparing growth and water tables
        
        
        
        v_ini, v, iv5,  cbt, dcbt, cb, dcb, w,dw,logs,pulp, dv,dlogs,dpulp,yrs, bmgr,  \
                                Nleach, Pleach, Kleach, DOCleach, runoff, \
                                nrelease, prelease,krelease, ch4release = run_susi(forc, wpara, cpara, 
                                org_para, spara, outpara, photopara, start_yr, end_yr, wlocation = 'undefined', 
                                mottifile=mottifile, peat= 'other', photosite='All data', 
                                folderName=folderName,ageSim=ageSim, sarkaSim=sarkaSim, sfc=sfc, susiPath=susiPath)
                                            
        
        out[s]={} 
        out[s]['ditch_depth']=spara['ditch depth east']
        out[s]['v_ini']=v_ini
        out[s]['v']=v[0]
        out[s]['iv5']=iv5[0]            
        out[s]['wt'] = w[0]
        out[s]['stand C bal'] = cbt[0]
        out[s]['peat C bal'] = cb[0]
        out[s]['ch4 emiss'] = ch4release[0]
        out[s]['Nleach'] = Nleach[0]
        out[s]['Pleach'] = Pleach[0]
        out[s]['Kleach'] = Kleach[0]
        out[s]['DOCleach'] = DOCleach[0]
        out[s]['runoff'] = runoff[0]
        out[s]['nrelease'] = nrelease[0]
        out[s]['prelease'] = prelease[0]
        out[s]['krelease'] = krelease[0]



        # N,P,K release           
        
        print ( '************************************************')
        print (v_ini, v, iv5, Kleach)
    dfOut = pd.DataFrame.from_dict(out)
    dfOut=dfOut.T
    print (dfOut)
    files =['vols_swidth+20.xlsx', 'vols_swidth-20.xlsx', 'vols_ddepth+20.xlsx', 'vols_ddepth-20.xlsx', 
        'vols_nroot_30.xlsx', 'vols_nroot_50.xlsx',
        'vols_afp_25.xlsx', 'vols_afp_5.xlsx', 'vols_nut+20.xlsx', 'vols_nut-20.xlsx', 'vols_anisotropy+20.xlsx',
        'vols_anisotropy-20.xlsx', 'vols_sphagnum.xlsx', 'vols_carex.xlsx']

    dfOut.to_excel(folderName+'susi_ash.xlsx')
    #dfOut.to_excel(folderName+files[13])
    
call_ash_susi()
#%%
from dwts_para import para
params = para(period='start-end')
print (params)