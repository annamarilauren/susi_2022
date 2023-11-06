# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 15:15:55 2023

@author: alauren
"""
import numpy as np
from scipy.integrate import quad
# Copied from Vauhkonen (2019):

''' Class StemCurve is structured to include stem taper models and stem bucking based on analyzing the stem tapering. '''
''' The attributes of the class include both allowable log dimensions for stem bucking and required parameters for Laaasasenaho's models. '''
''' The methods of the class are both for 1) determining the stem curve and extracting values from it; 2) virtually bucking the stem '''
''' based on the tapering and allowable log dimensions, with an aim to produce as much saw wood logs as possible. ''' #(Vauhkonen 2019)

class StemCurve:
   def __init__(self):
      self.min_diams = {'log':(0,15,16,18), 'pulp':(0,7,7,7)}
      self.lengths =  { 'log':{ 1:(0,0,0,0,0,4.3,4.6,4.9,5.2,5.5,5.8,6.1), \
                              2:(0,0,0,0,0,4.3,4.6,4.9,5.2,5.5,5.8,6.1), \
                              3:(0,3.1,3.4,3.7,4.0,4.3,4.6,4.9,5.2,5.5,5.8,6.1)},
                       'pulp':{ 1:(0,3.0,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,6.0), \
                                2:(0,3.0,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,6.0), \
                                3:(0,3.0,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,6.0) } 
                      }
      self.stemcurve_coefs = { 1:(2.1288,-0.63157,-1.6082,2.4886,-2.4147,2.3619,-1.7539,1.0817), \
                               2:(2.3366,-3.2684,3.6513,-2.2608,0,2.1501,-2.7412,1.8876), \
                               3:(0.93838,4.106,-7.8517,7.8993,-7.5018,6.3863,-4.3918,2.1604) }
      self.stumpheight_coefs = { 1:(0.09522,0.4456), 2:(0.56,0.5089), 3:(0.497936,0.4862) }

   def stemCurve(self, hx, h, species):
      x = 1 - hx / h
      b = self.stemcurve_coefs[species]
      return b[0]*x + b[1]*x**2 + b[2]*x**3 + b[3]*x**5 + b[4]*x**8 + b[5]*x**13 + b[6]*x**21 + b[7]*x**34

   def integrand(self,hx,h,d20,sp):
      return np.pi / 4 * (self.stemCurve(hx,h,sp) * d20 / 100) **2 * 1000

   def calculateVolume(self,d13,h,sp,lower,upper):
      d20 = d13 / self.stemCurve(1.3, h, sp)
      return quad(self.integrand, lower, upper, args=(h,d20,sp))[0]

   def detectCutHeight(self, dcut, d13, h, sp):
      d20 = d13 / self.stemCurve(1.3, h, sp)
      htest = 0
      diff = h / float(2)
      hsol = -1
      while hsol != htest:
         if htest > h:
            htest = h
         dtest = self.stemCurve(htest, h, sp) * d20
         if dtest < dcut + 0.001 and dtest > dcut - 0.001:
            hsol = htest
            break
         if dtest > dcut:
            htest = htest + diff
         else:
            htest = htest - diff
         diff = diff * 0.5
      return hsol

   def predictStumpHeight(self,d13,h,sp):
      stumph = self.stumpheight_coefs[sp][0] * h + self.stumpheight_coefs[sp][1] * d13
      stumph = stumph / 100
      if stumph < 0.1:
         return 0.10
      else: return stumph

   def buckStem(self,d13,h,sp):
      cutPoints = { 'stump':[0,0,0,0,0], 'pulp':[0,0,0,0,0], 'log':[0,0,0,0,0] }
      cutPoints['stump'][0] = self.predictStumpHeight(d13,h,sp)

      log_flag, pulp_flag = 1,1
      if d13 < self.min_diams['log'][sp]: log_flag = 0
      if d13 < self.min_diams['pulp'][sp]: pulp_flag=0

      if log_flag==1:
         hlog = self.detectCutHeight(self.min_diams['log'][sp], d13, h, sp)
         logprop = hlog - self.predictStumpHeight(d13,h,sp)
         for l in self.lengths['log'][sp]:
            for k in self.lengths['log'][sp]:
               for j in self.lengths['log'][sp]:
                  for i in self.lengths['log'][sp]:
                     apu = i+j+k+l
                     if apu > logprop:
                        break
                     if apu > cutPoints['log'][0]:
                        cutPoints['log'][0] = apu
                        cutPoints['log'][1] = i
                        cutPoints['log'][2] = j
                        cutPoints['log'][3] = k
                        cutPoints['log'][4] = l
      cutPoints['log'][0] = cutPoints['log'][0] + self.predictStumpHeight(d13,h,sp)

      if pulp_flag==1:
         hpulp = self.detectCutHeight(self.min_diams['pulp'][sp], d13, h, sp)
         pulpprop = hpulp - cutPoints['log'][0]
         for l in self.lengths['pulp'][sp]:
            for k in self.lengths['pulp'][sp]:
               for j in self.lengths['pulp'][sp]:
                  for i in self.lengths['pulp'][sp]:
                     apu = i+j+k+l
                     if apu > pulpprop:
                        break
                     if apu > cutPoints['pulp'][0]:
                        cutPoints['pulp'][0] = apu
                        cutPoints['pulp'][1] = i
                        cutPoints['pulp'][2] = j
                        cutPoints['pulp'][3] = k
                        cutPoints['pulp'][4] = l
      cutPoints['pulp'][0] = cutPoints['pulp'][0] + cutPoints['log'][0]

      return cutPoints

   def predictAssortmentVolumes(self,d13,h,sp):
      cutPoints = self.buckStem(d13,h,sp)
      volume = {'log':0, 'pulp':0, 'residual':0, 'total':0 }
      volume['log'] = self.calculateVolume(d13, h, sp, cutPoints['stump'][0], cutPoints['log'][0])
      volume['pulp'] = self.calculateVolume(d13, h, sp, cutPoints['log'][0], cutPoints['pulp'][0])
      volume['residual'] = self.calculateVolume(d13, h, sp, 0, cutPoints['stump'][0]) + self.calculateVolume(d13, h, sp, cutPoints['pulp'][0], h)
      volume['total'] = volume['log'] + volume['pulp'] + volume['residual']
      volume['residual'] = volume['residual'] - self.calculateVolume(d13, h, sp, 0, cutPoints['stump'][0])
      return volume