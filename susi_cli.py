## 2020-11-23

## simple cli interface with cli arguments for stand (motti) and weather data

import os
##import pandas as pd
import datetime
from susi_utils import get_motti, read_FMI_weather, get_mese_input
from susi_para import get_susi_para
from susi83 import run_susi
import susi_io
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='testi')
    parser.add_argument('--stand', type = str, required = True)
    parser.add_argument('--weather', type = str, required = True)
    return parser.parse_args()

args = parse_arguments()

## read in stand data, not used currently
## df = get_motti(args.stand)

## read in weather data
start_date = datetime.datetime(2008,1,1) ## ?
end_date = datetime.datetime(2009,12,31) ## ?
forc = read_FMI_weather(0, start_date, end_date, sourcefile = args.weather)

start_yr = start_date.year
end_yr = end_date.year
yrs = (end_date - start_date).days/365.25
sfc: int = 3                               ## soil fertility class
ageSim: float = 40.0
sarkaSim: float = 40.0
n: int = int(sarkaSim / 2)
site: str = 'develop_scens'
folderName = os.getcwd()
susiPath = os.getcwd()

wpara, cpara, org_para, spara, outpara, photopara = get_susi_para(wlocation = 'undefined',
                                                                  peat = site,
                                                                  folderName = folderName,
                                                                  hdomSim = None,
                                                                  ageSim = ageSim,
                                                                  sarkaSim = sarkaSim,
                                                                  sfc = sfc,
                                                                  susiPath = susiPath,
                                                                  n = n)

v_ini, v, iv5, cbt, dcbt, cb, dcb, w, dw, logs, pulp, dv, dlogs, dpulp, yrs, bmgr, \
    Nleach, Pleach, Kleach, DOCleach, runoff = run_susi(forc, wpara, cpara,
                                                        org_para, spara, outpara,
                                                        photopara, start_yr, end_yr,
                                                        wlocation = 'undefined',
                                                        mottifile = args.stand,
                                                        peat = 'other',
                                                        photosite = 'All data',
                                                        folderName = os.getcwd(),  ## ?
                                                        ageSim = ageSim,
                                                        sarkaSim = sarkaSim,
                                                        sfc = sfc,
                                                        susiPath = os.getcwd())  ## ?
