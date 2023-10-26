# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 14:31:01 2023

@author: Usuario
"""
import numpy as np
import matplotlib.pyplot as plt

def get_sigmoidal_burst_size(t, c=0.5, p=0.5, t0 = 30.0, a=1.0):
    tlim = 38.0
    if t<=tlim:
        return float(c)/(1.0+np.exp(-(t-t0)*p))
    elif t>tlim:
        return float(c)/(1.0+np.exp(-(tlim-t0)*p))+ a*(t-tlim)

p = 0.4
c = 20.0
t0 = 27.0
a = 0.1

t = np.linspace(0.0,55.0,1000)
offspring = [get_sigmoidal_burst_size(i,c,p,t0,a) for i in t]

plt.figure()
plt.plot(t,offspring)
plt.ylim(0.0,30.0)

