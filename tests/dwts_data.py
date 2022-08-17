# -*- coding: utf-8 -*-
"""
Created on Thu Aug 08 13:36:07 2019

@author: alauren
"""
import pandas as pd
import numpy as np
import susi_utils as su
import matplotlib.pylab as plt
import datetime
import dwts_para
from scipy.interpolate import interp1d


class Dwt_site():
    def __init__(self, site, name):
        self.name = name
        self.folder= r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/susi_experim/wt/'
        self.wfolder =r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/susi_experim/weather/'
        self.site = site
        self.heads =['time', 'tube1', 'tube2','tube3','tube4','tube5','tube6','tube7','tube8', 
                     'tube9', 'tube10','tube11','tube12','tube13','tube14','tube15','tube16',
                     'tube17', 'tube18','tube19','tube20','tube21','tube22','tube23','tube24',
                     'tube25', 'tube26','tube27','tube28','tube29','tube30','tube31','tube32',
                     'tube33', 'tube34','tube35','tube36','mean']
                     
        self.forc = su.read_FMI_weather(0, self.site['start_date'], self.site['end_date'], sourcefile=self.wfolder+self.site['wfile'])           # read weather input0

        self.length = len(self.forc)        
        self.vol = self.site['vol']                                                                               # stand volume in the beginning and end of period, list []
        self.bd = self.site['bulk dens']                                                                          # bulk density of peat kg m-3
        
        self.z = [-0.05,-0.15,-0.25,-0.35]
        self.ppara = dwts_para.get_photo_para(site='All data')
        #if 'koira' in self.name:        
        #    self.ppara = dwts_para.get_photo_para(site='Sodankyla')
        #    print 'sodankyla'
    def prepare_data(self, start_date, end_date, wtfile, weatherfile, sheads, gseason=True ):
        dfw = pd.read_csv(wtfile,sep=';', skiprows=1, names=sheads)
        try:    
            dfw['date']= pd.to_datetime(dfw['time'], format='%Y-%m-%d')
        except:
            dfw['date']= pd.to_datetime(dfw['time'], format='%d.%m.%Y')       
        dfw = dfw.set_index(dfw['date'])
        dfw = dfw.drop(['time', 'date'], axis = 1)
        dfw = dfw[dfw.index < end_date]
        dfw = dfw[dfw.index > start_date]
        forc=su.read_FMI_weather(0, start_date, end_date, sourcefile=weatherfile)           # read weather input
        dfw['airT']=forc['T']
        dfw['Par']=forc['Par']                                                    #W/m2
        
        if gseason: dfw =dfw[dfw['airT']>5.0]
        return dfw
    
    def het_respiration(self, vol, bd, start_date, end_date, wtfile, weatherfile, sheads):
        dfwt = pd.read_csv(wtfile,sep=';', skiprows=1, names=sheads)
        try:    
            dfwt['date']= pd.to_datetime(dfwt['time'], format='%Y-%m-%d')
        except:
            dfwt['date']= pd.to_datetime(dfwt['time'], format='%d.%m.%Y')       
        dfwt = dfwt.set_index(dfwt['date'])
        dfwt = dfwt.drop(['time', 'date'], axis = 1)
        dfwt = dfwt[dfwt.index < end_date]
        dfwt = dfwt[dfwt.index > start_date]
    
        forc=su.read_FMI_weather(0, start_date, end_date, sourcefile=weatherfile)           # read weather input
        dfwt['airT']=forc['T']
        B = 350. ; T5zero = -46.02; T5ref = 10.
        years = float(len(range(start_date.year, end_date.year+1)))
        vol_inter = interp1d([start_date.year, end_date.year], vol)    
        out = 0.    
        for yr in  range(start_date.year, end_date.year+1):
            wt = np.mean(dfwt[str(yr)+'-05-01':str(yr)+'-10-31'].mean().values[:-1])
            T5 = dfwt[str(yr)+'-05-01':str(yr)+'-10-31']['airT'].values
            T5 = np.minimum(T5, 16.)
            v = vol_inter(yr)        
            Rref = (0.0695+3.7*10**(-4)*v +5.4*10**(-4)*bd + 1.2*10**(-3)*wt)*24.    #unit g/m2/h CO2-> g/m2/day
            Rhet = Rref*np.exp(B*((1./(T5ref-T5zero))-(1./(T5-T5zero))))        
            out+= sum(Rhet)
            #print yr, sum(Rhet)
        return out*10./years, years   #unit kg/ha in the time interval, in CO2 -> /year
    
    
    def water_table(self, tube, dfwa):
        return dfwa[tube].mean()
    
    def get_rew(self, ptype, vonP, df):
        pF, _ = su.peat_hydrol_properties(vonP, var='H', ptype=ptype) # peat hydraulic properties after Päivänen 1973    
        df = df.drop(['airT'], axis=1) 
        h0 =  self.z[0]- df.mean(axis=1)*-1e-2
        h1 =  self.z[1]- df.mean(axis=1)*-1e-2
        wcont0 = su.wrc(pF[0],x=-h0,var='Psii')
        wcont1=su.wrc(pF[1],x=-h1,var='Psii')
        fc0 = su.wrc(pF[0],x=-0.3+h0 ,var='Psii')
        wp0 =su.wrc(pF[0],x=-150.+h0 ,var='Psii')
        fc1 = su.wrc(pF[1],x=-0.3+h1 ,var='Psii')
        wp1 =su.wrc(pF[1],x=-150.+h1 ,var='Psii')
        df['rew0'] = (wcont0-wp0)/(fc0-wp0)
        df['rew1'] = (wcont1-wp1)/(fc1-wp1)
        df['rew'] = df[['rew0','rew1']].mean(axis=1)
        df['ones'] = 1.
        df['zeros'] = 0.
        df['cut1'] = df[['rew','ones']].min(axis=1)
        df['cut2'] = df[['cut1','zeros']].max(axis =1)
        return  df['cut2']
    
    def drought(self, tube, df):
        limits = [70., 80., 90., 100.]
        gsdays = float(len(df))
        fs=[]    
        for limit in limits:
            crit = df[tube].values > limit
            fs.append(crit.sum()/gsdays)
        return fs[0],fs[1],fs[2],fs[3]
        
    def air_filled(self, ptype, vonP, tube, df):
        pF, _ = su.peat_hydrol_properties(vonP, var='H', ptype=ptype) # peat hydraulic properties after Päivänen 1973    
        h0 =  self.z[0]- df[tube]*-1e-2
        h1 =  self.z[1]- df[tube]*-1e-2
        wcont0 = su.wrc(pF[0],x=-h0,var='Psii')
        airf0 = pF[0][0]-wcont0
        wcont1=su.wrc(pF[1],x=-h1,var='Psii')
        airf1 = pF[1][0]-wcont1
        return (airf0+airf1)/2.
    
    def frequency(self, afp):
        limits = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.15, 0.18, 0.2]
        gsdays = float(len(afp))
        fs=[]    
        for limit in limits:
            crit = afp.values < limit
            fs.append(crit.sum()/gsdays)
        return fs[2], fs[4], fs[7]        
        
    def nutrient_demand(self, vol, leaf_mass, yrs):
        #vol as list: start, end, leaf_mass kg/ha at the end, yrs is length of the study period
        MarjoNut = lambda vol, lna, b, k: np.exp(lna + b*np.log(vol) + k)
        from_abovegr_to_total = 1.2    
    
        #litter_nut={'N':12.0, 'P':1.0, 'K':3.5}      #mg/g Laiho 1997 page 49
        litter_nut={'N':6.9, 'P':0.81, 'K':2.81}      #mg/g Palviainen & Finer 2012
        retrans = {'N': 0.69, 'P': 0.73, 'K':0.8}     #Nieminen Helmisaari 1996 Tree Phys
        longevityLeaves = {'pine':4.5, 'spruce':5.} #yrs, Lamppu & Huttunen 2001 CJFR life span of leaves and fine roots
    
        par = {                                                                 # Palviainen & Finer, 2012 Eur J For Res 131:945-964, eq 2, Table 7
        'N':{'pine':[1.856,0.631,0.050]},
        'P':{'pine':[-2.387,0.754,0.158]},
        'K':{'pine':[0.839,0.650,0.070]}        
        }
        lna,b,k=par['N']['pine']                                                # pine
        N = (MarjoNut(vol[-1], lna, b, k) - MarjoNut(vol[0], lna, b, k))* from_abovegr_to_total + \
            leaf_mass*litter_nut['N']/1000.*(1.-retrans['N'])*yrs/longevityLeaves['pine']
    
        lna,b,k=par['P']['pine']                                                # pine
        P = (MarjoNut(vol[-1], lna, b, k) - MarjoNut(vol[0], lna, b, k))* from_abovegr_to_total + \
            leaf_mass*litter_nut['P']/1000.*(1.-retrans['P'])*yrs/longevityLeaves['pine']
    
        lna,b,k=par['K']['pine']                                                # pine
        K = (MarjoNut(vol[-1], lna, b, k) - MarjoNut(vol[0], lna, b, k))* from_abovegr_to_total + \
            leaf_mass*litter_nut['K']/1000.*(1.-retrans['K'])*yrs/longevityLeaves['pine']
    
    
        return N, P, K
    
    def nutrient_supply(self, hr, sfc):
        #hr as kg CO2/ha/yr
        #sfc is soil fertility class
        Ncont =  {2: 1.9, 3: 1.6, 4: 1.4, 5: 1.2}                                       # Mese study: N cont in OM % dm
        Pcont =  {2: 0.1, 3: 0.08, 4: 0.06, 5: 0.05}                                    # Mese study: P cont in OM % dm
        Kcont=  {2: 0.041, 3: 0.04, 4: 0.035, 5: 0.03}                                 # Mese study: P cont in OM % dm (alkuperäinen); 2: 0.041, 3: 0.04, 4: 0.035, 5: 0.03
        C_in_OM = 0.55                                                              # C content in OM kg kg-1
        CO2_to_C = 12./44.
        Nt = hr * CO2_to_C / C_in_OM * Ncont[sfc] / 100.    # Nrelease kg ha-1 yr-1
        Pt = hr * CO2_to_C / C_in_OM * Pcont[sfc] / 100.    # Prelease kg ha-1 yr-1
        Kt = hr * CO2_to_C / C_in_OM * Kcont[sfc] / 100.    # Prelease kg ha-1 yr-1
        return Nt, Pt, Kt
    
    #def get_weather_data(self, start_date, end_date, weatherfile):
    #    forc=su.read_FMI_weather(0, start_date, end_date, sourcefile=weatherfile)           # read weather input0
    #    return forc
    
    def assimilation(self, ppara, rg, vpd, Ta_minus1,Ta, rew, LAI, Xk, hdom):
        """
        Computes photosynthesis and respiration of the stand in daily time step
        Mäkelä et al. 2008. Empirical model of stand GPP LUE approach. Global Change Biology 14: 92-108
        IN:
            (develpment of all sided leaf area index (LAI) m2 m-2)
            meteorological input data:
            rg in W/m2
            vpd in kPa
            Ta previous day tamperature
            Ta  in deg C
            Xk acclimatisation temperature, Ta in the first time step
            Ns,Ps,Ks nutrient release in time step kg/ha/timestep
            Nconc, Pconc, Kconc nutrient concentration in stand 
            hdom dominant height of the stand (used in computation of respiration)
        OUT:
            Net primary production NPP in kg/ha/yr        
            Xa
        """
        par = rg*4.6*0.5 /1000000.*86400.                                          # Unit conversion to mol/m2/day, 0.5 is the share of par from rg 
    
        """ Eq 2 """
        fL = 1./(ppara['gamma']*par+1.)                                             #
    
        """Eq. 3a and 3b"""
        Xk = Xk+(Ta-Xk)/ppara['tau']    
        Sk = np.maximum(Xk-ppara['X0'], 0.0)                           
    
        """Eq 4"""
        fs = np.minimum(Sk/ppara['Smax'], 1.0)                         # Final temperature adustment
    
        """Eq 5"""                                                                  # Vapor pressure deficit function: limits photos when stomata are closed
        fd = np.exp(ppara['kappa']*vpd)                                             # Eq 5
    
        """Eq 6"""                                                                  # Soil water content modifyer 
        fREW=((1.+ ((1.-rew)/ppara['alfa'])**ppara['nu']))**(-1.)
        #fREW = rew
        """ Beer-Lambert function: LAI function"""                                  # Fabrika 2013
        kext = 0.19
        
        """Eq 1"""
        #LAI = LAI*11./6.4     #to all sided LAI Mäkelä et al. 2001 (LAI 11)
        Pk=ppara['beta'] *(1.- np.exp(-kext*LAI)) * par * fL * fs * fd * fREW   #canopy GPP gC m-2 day-1
        """Respi"""
        #Conversion micro mol/m2/s ->g/m2/day
        #Rtrees= np.maximum((r0*q10**((Ta-10.)/10.)-cr)*rlai, np.zeros(sim_len))
        #Rtrees = 0.47*Pk  #Waring et al. 1998, Tree Physiology 18, 129-132
        fr=interp1d(np.array([0.0,10.0,20.0,30.0,]), 
                    np.array([0.6, 0.50,0.36,0.20]), fill_value='extrapolate')  #Mäkelä & Valentine 2001 Tree Physiology 21, 1015-1030, Härkönen et al. 2010 ForEco 259, 524-533    
    
        NPP = Pk*fr(hdom)* 2.      #Change from gC m-2 to g organic matter
        #NPP = (Pk-Rtrees)*2. 
        
        return NPP*10., Xk    #kg/ha organic matter
    
    
    def compute(self):
        #afpfrs6=[]; afpfrs10=[];afpfrs15=[];wts=[]; droughts70=[];droughts80=[];droughts90=[];droughts100=[]
        sheads = self.heads[0:self.site['ntubes']+2]
        dfwater = self.prepare_data(self.site['start_date'], self.site['end_date'], self.folder+self.site['file'], \
            self.wfolder+self.site['wfile'], sheads)
        dfall = self.prepare_data(self.site['start_date'], self.site['end_date'], self.folder+self.site['file'], \
            self.wfolder+self.site['wfile'], sheads, gseason=False)
        rews = self.get_rew(self.site['ptype'],self.site['vonP'], dfall)
        afpfrs6=[]; afpfrs10=[];afpfrs15=[];wts=[]; droughts70=[];droughts80=[];droughts90=[];droughts100=[]
        for tube in sheads[1:]:                                            # loop through wt tubes in the data
            afp = self.air_filled(self.site['ptype'],self.site['vonP'], tube, dfwater)                                  # convert water table to airfilled porosity in 20 cm layer
            afpfr = self.frequency(afp)                                                                          # compute frequeny of air filled porosity against afp criteria
            afpfrs6.append(afpfr[0])    
            afpfrs10.append(afpfr[1])    
            afpfrs15.append(afpfr[2])    
            wt = self.water_table(tube, dfwater)
            wts.append(wt)
            dro =self.drought(tube, dfwater)                                                                          # frequency of drought in growing season according to different criteria
            droughts70.append(dro[0])
            droughts80.append(dro[1])
            droughts90.append(dro[2])
            droughts100.append(dro[3])

        Ta_minus1 = self.forc['T'][0]; Xk =0.0; Xk2=0.0; rew=1.0; year = self.forc.index.year[0]
        fhdom = interp1d(np.array([float(self.site['start_date'].year), float(self.site['end_date'].year)]), \
                                    np.array(self.site['hdom']))     
        fvol = interp1d(np.array([float(self.site['start_date'].year), float(self.site['end_date'].year)]), np.array(self.vol))     
        fleaf = interp1d(np.array([0.,50.,100.,150.,200.,250.,300.]), np.array([0., 3., 4.3, 5.4, 6.2, 6.5, 6.7]), fill_value='extrapolate')        
        fbm = interp1d(np.array([0.,50.,100.,150.,200.,250.,300.]), np.array([0., 33.7, 64.2, 96., 125., 153., 186.])*1000., fill_value='extrapolate') #motti biomass units from tn/ha to kg /ha   
        fbmtovol = interp1d(np.array([0., 33.7, 64.2, 96., 125., 153., 186.])*1000., np.array([0.,50.,100.,150.,200.,250.,300.]), fill_value = 'extrapolate') 
        bm_growth = fbm(self.vol[-1])-fbm(self.vol[0])    

        NPP=0.; NPPwater =0.
        for n,row in enumerate(self.forc.iterrows()):
            if year != float(row[0].year): 
                Xk = 0.
                Xk2 = 0.
            year =float(row[0].year)        
            V = fvol(year)
            bm = fbm(V)
            leaf_mass = fleaf(V)
            sla=12. #(12.5+15.)/2.             #Härkönen et al. 2015 BER 20, 181-195, Table 3, All sided LAI      
            LAI = leaf_mass/10.*sla       
            hdom = fhdom(year)
            rg, vpd, Ta = row[1][8], row[1][13], row[1][4]
            npp, Xk = self.assimilation(self.ppara, rg, vpd, Ta_minus1, Ta, rew, LAI, Xk, hdom)
            NPP = NPP + npp
            rew2 = rews.values[n] if n < len(rews.values) else 1.0
            nppwater, Xk2 = self.assimilation(self.ppara, rg, vpd, Ta_minus1, Ta, rew2, LAI, Xk2, hdom)       
            NPPwater = NPPwater + nppwater
            Ta_minus1 = Ta
            #print s, V, LAI
        hr, yrs = self.het_respiration(self.vol, self.bd, self.site['start_date'], self.site['end_date'], self.folder+self.site['file'], \
            self.wfolder+self.site['wfile'], sheads)                                                             # heterotrophic respiration
        Nsup, Psup, Ksup = self.nutrient_supply(hr, self.site['sfc'])    
        Ndem, Pdem, Kdem = self.nutrient_demand(self.vol,leaf_mass*1000.,yrs)
        vol_computed_end = fbmtovol(fbm(self.vol[0]) + NPPwater)    
        growth = self.vol[-1]-self.vol[0]
        _,_,Kdem_pot = self.nutrient_demand([self.vol[0], vol_computed_end], leaf_mass*1000., yrs)
    
        print self.name, NPP/yrs, dfall['airT'].mean(), dfall['Par'].mean() 
        reslist = [self.name, self.site['status'],np.mean(afpfrs6), np.mean(afpfrs10),np.mean(afpfrs15), \
                np.mean(wts)/100., \
                np.mean(droughts70),np.mean(droughts80),np.mean(droughts90),np.mean(droughts100), \
                hr, growth/yrs, Ndem/yrs, Pdem/yrs, Kdem/yrs, Nsup, Psup, Ksup, bm_growth/yrs, NPP/yrs, vol_computed_end, \
                min(1.0,Ksup/(Kdem_pot/yrs)), min(1.0,Ksup/(Kdem_pot/yrs))*NPP/yrs, self.site['bulk dens'], self.site['vol'][0], \
                dfall['Par'].mean()]
    
        lab = self.site['mark'] 
        plt.subplot(331) 
        plt.title('resp')
        plt.plot(hr, growth/yrs, lab)
        plt.subplot(332)
        plt.title('dwt')
        plt.plot(np.mean(wts)/100.,growth/yrs, lab)
        plt.subplot(333)
        plt.title('afp')
        plt.plot(np.mean(afpfrs10),growth/yrs, lab)
        plt.subplot(334)
        plt.title('N demanand')
        plt.plot(Nsup, Ndem/yrs, lab)
        plt.subplot(335)
        plt.title('P demanand')
        plt.plot(Psup, Pdem/yrs, lab)
        plt.subplot(336)
        plt.title('K demanand')
        plt.plot(Ksup, Kdem/yrs, lab)
        plt.plot(Ksup,Ksup,'ko')
        plt.xlabel('supply')
        plt.ylabel('demand')        
        plt.subplot(337)
        
        plt.title('NPP')
        plt.plot(bm_growth/yrs, NPPwater/yrs,  lab)
        plt.plot(NPPwater/yrs,NPPwater/yrs, 'ko')
        plt.xlabel('observed')
        plt.ylabel('predicted')
        #plt.plot(bm_growth/yrs, NPPwater/yrs,  'co')

        plt.subplot(338)
        plt.title('K demand pot')
        plt.plot(Ksup, Kdem_pot/yrs, lab)
        plt.plot(Ksup,Ksup, 'ko')
        
        plt.subplot(339)
        plt.title('K restricted growth')
        kmodif = min(1.0,Ksup/(Kdem_pot/yrs))
        #plt.scatter(bm_growth/yrs, kmodif*NPPwater/yrs, marker = r"$ {} $".format(self.name[:3]), s=350.)
        plt.plot(bm_growth/yrs, NPPwater/yrs*kmodif, lab)        
        plt.plot(bm_growth/yrs,bm_growth/yrs, 'ko')
        plt.xlabel('observed bm growth kg/yr')
        plt.ylabel('computed growth')
        #plt.plot(NPPwater/yrs,NPPwater/yrs,'co')

        return reslist

sites = dwts_para.para('start-end')  #'start-end'  'start-thin'  'thin-end'

sitelist=['koira11', 'koira12','koira21','koira22','ansa11','ansa16','ansa21','ansa26','neva11',
          'neva14','neva21','neva24','neva31','neva34','jaakkoin61','jaakkoin62','parkano11','parkano12']  #,'katila'
results =[]
koira11 = Dwt_site(sites['koira11'], 'koira11')
koira12 = Dwt_site(sites['koira12'], 'koira12')
koira21 = Dwt_site(sites['koira21'], 'koira21')
koira22 = Dwt_site(sites['koira22'], 'koira22')
ansa11 = Dwt_site(sites['ansa11'], 'ansa11')
ansa16 = Dwt_site(sites['ansa16'], 'ansa16')
ansa21 = Dwt_site(sites['ansa21'], 'ansa21')
ansa26 = Dwt_site(sites['ansa26'], 'ansa26')

neva11 = Dwt_site(sites['neva11'], 'neva11')
neva14 = Dwt_site(sites['neva14'], 'neva14')
neva21 = Dwt_site(sites['neva21'], 'neva21')
neva24 = Dwt_site(sites['neva24'], 'neva24')
neva31 = Dwt_site(sites['neva31'], 'neva31')
neva34 = Dwt_site(sites['neva34'], 'neva34')
jaakkoin61 = Dwt_site(sites['jaakkoin61'], 'jaakkoin61')
jaakkoin62 = Dwt_site(sites['jaakkoin62'], 'jaakkoin62')

parkano11 = Dwt_site(sites['parkano11'], 'parkano11')
parkano12 = Dwt_site(sites['parkano12'], 'parkano12')
katila = Dwt_site(sites['katila'], 'katila')

siteins =[koira11, koira12, koira21, koira22, ansa11, ansa16, ansa21, ansa26,
          neva11, neva14, neva21, neva24, neva31, neva34, jaakkoin61, jaakkoin62, parkano11, parkano12]
for s in siteins:
    reslist=s.compute()    
    results.append(reslist)
    
plt.show()
f =  r'C:/Users/alauren/Documents/WinPython-64bit-2.7.10.3/susi_experim/out.csv'

with open(f,'wb') as fi:
    titles=  ['site', 'status', 'afp6', 'afp10', 'afp15', 'meanwts','drought70', 'drought80', 'drought90', \
                    'drought100' , 'hrespi', 'volgr', 'Ndem', 'Pdem', 'Kdem', 'Nsup','Psup','Ksup','bm_growth', 'NPP', \
                    'vol_computed_end', 'kmodif', 'krestrictedgr', 'bulk density', 'inivol', 'meanparWm-2' ]  
    for entries in titles:
        fi.write(str(entries))
        fi.write(';')
    fi.write('\n')
    for line in results:
        for entries in line:
            fi.write(str(entries))
            fi.write(';')
        fi.write('\n')
                                                          # number of headings    
#fig = plt.subplot(331)

    