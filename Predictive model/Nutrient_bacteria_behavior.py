# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:36:25 2023

@author: Usuario
"""

import numpy as np
import matplotlib.pyplot as plt

def monod_growth_law(c, kn, gnmax=0.034):
    return gnmax * c / (kn + c)


        
        
gmax = 0.079
B = 1.0e6
n0 = 1.0e7
delta_t = 0.01
tmax = 70

t,nutrients = nutrient_integral(gmax,B,n0,delta_t,tmax)       
plt.plot(t,nutrients)
