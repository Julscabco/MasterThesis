# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:15:24 2023

@author: Usuario
"""

import matplotlib.pyplot as plt
import numpy as np


def get_burst_size(t, tlim, a):
    if t < tlim:
        return 0.0
    elif t >= tlim:
        return a * t - a*tlim
    
    
time = np.linspace(0.0,40.0,1000)
burst_size = [get_burst_size(t,15.0,10.0) for t in time]

plt.plot(time,burst_size)
plt.xlabel(r'${\tau}_{l} (min)$',fontsize=14)
plt.ylabel('Burst size',fontsize=14)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)