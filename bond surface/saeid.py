#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 18:15:18 2022

@author: dicksonnkwantabisa
"""

import os
import sys
import numpy as np
import pandas as pd
import math
# from scipy.stats import truncnorm
# from scipy.stats import weibull_min
# from scipy.stats import invgamma
# from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import random
import seaborn as sns
import math
# from scipy.optimize import minimize
from math import inf 
# from numba import jit
# import time as time
# from IPython import get_ipython
# from datetime import datetime as datetime
# import multiprocessing


def inst_to_ann(r):
    """
    Converts short rate to annualized rate
    """
    return np.expm1(r)

def ann_to_inst(r):
    """
    Converts annualized to short rate
    """
    return np.log1p(r)



# Simulation of Bond prices zero-coupon
def ann_to_inst(r):
    """
    Converts annualized to short rate
    """
    return np.log1p(r)


def cir_price(n_years, n_scen, θ, m, σ, θ_star, 
              m_star, λ_r, steps_per_year,
              r_0=None, risk_n=True):
   """
   Generate random interest rate and bond price evolution over 
   time using the CIR model.
   m and r_0 (which are annualized) are changed to the instantaneous short
   rate and the returned values are short rates as well.
   n_years : number of years on the horizon
   n_scen : number of simulations or scenarios
   θ : mean reversion rate
   m : long-run mean of interest rate
   σ : diffusion parameter
   θ_star : mean reversion parameter under the risk neutral measure
   m_star : long-run mean of interest rate under the risk neutral measure
   λ_r : market price of risk claiberated in R
   r_0 : initial (start) interest rate annualized
   """
   if r_0 is None: r_0 = m 
   r_0 = ann_to_inst(r_0)
   dt = 1/steps_per_year
   num_steps = int(n_years*steps_per_year) + 1 # int because n_years might be a float 
      
   shock = np.random.normal(0, scale=np.sqrt(dt), size=(num_steps, n_scen))
   rates = np.empty_like(shock)
   rates[0] = r_0
   
   ## For Price Generation
   h = math.sqrt(θ**2 + 2*σ**2)
   h_r = math.sqrt((θ+λ_r)**2 + 2*σ**2) # this incorporate the market price of risk (λ_r)
   prices = np.empty_like(shock)
   ####
   def price(ttm, r):
       if risk_n==True:           
           _A = ((2*h_r*math.exp((θ+λ_r+h_r)*ttm/2))/(2*h+(θ+λ_r+h_r)*(math.exp(h_r*ttm)-1)))**(2*θ_star*m_star/σ**2)
           _B = (2*(math.exp(h_r*ttm)-1))/(2*h_r + (θ+λ_r+h_r)*(math.exp(h_r*ttm)-1))
           _P = _A*np.exp(-_B*r)
       else:
           _A = ((2*h*math.exp((h+θ)*ttm/2))/(2*h+(h+θ)*(math.exp(h*ttm)-1)))**(2*θ*m/σ**2)
           _B = (2*(math.exp(h*ttm)-1))/(2*h + (h+θ)*(math.exp(h*ttm)-1))
           _P = _A*np.exp(-_B*r)
       return _P
   
   prices[0] = price(n_years, r_0)
   ####
   for step in range(1, num_steps):
       r_t = rates[step-1]
       if risk_n==True:           
           d_r_t = θ_star*(m_star-r_t)*dt + σ*np.sqrt(r_t)*shock[step]
           rates[step] = abs(r_t + d_r_t)
           # generate prices at time t as well ...
           prices[step] = price(n_years-step*dt, rates[step])
       else:
           d_r_t = θ*(m-r_t)*dt + σ*np.sqrt(r_t)*shock[step]
           rates[step] = abs(r_t + d_r_t)
           prices[step] = price(n_years-step*dt, rates[step])
       
   rates = pd.DataFrame(data=(rates), index=range(num_steps))
   ### for prices
   prices = pd.DataFrame(data=prices, index=range(num_steps))
 
   return rates, prices


    
#Example below

θ = 0.9552238
m = 0.7100216 
σ = 1.1974436
λ_r = -0.5547602

θ_star = θ + λ_r
m_star = (θ * m) / (θ + λ_r)
riskn_rates, riskn_prices = cir_price(n_years=2.,θ=θ, 
                              m=m, σ=σ, θ_star=θ_star,\
                              m_star=m_star, λ_r=λ_r,n_scen=100,\
                              steps_per_year=50,risk_n=True)


    
phy_rates, phy_prices = cir_price(n_years=2.,θ=θ, 
                              m=m, σ=σ, θ_star=θ_star,\
                              m_star=m_star, λ_r=λ_r,n_scen=100,\
                              steps_per_year=50,risk_n=False)


# I am changing the matrix into a vector here
zc_price = np.array(riskn_prices.iloc[1:101,:])
cir_ir = np.array(riskn_rates.iloc[1:101, :])

ir_vec = []
for i in range(0,len(cir_ir)):
    ir_vec.append(cir_ir[:, i])

pd.DataFrame(np.array(ir_vec).reshape(10000,1)).to_csv('interest_rate_cir_vec_saeid.csv', index=False)


# read in a catastrophe risk file

