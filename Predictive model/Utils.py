# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:23:34 2023

@author: Usuario
"""


import numba

import numpy as np


@numba.jit(nopython=True)
def get_linear_burst(t, tlim, a):
    if t < tlim:
        return 0.0
    elif t >= tlim:
        return a * t - a*tlim

@numba.jit(nopython=True)
def taul_to_minutes(t, n, k):
    return t * (float(k) / float(n)) * 23.0


def minutes_to_taul(t, n, k):
    return t * (float(n) / float(k)) / 23.0


@numba.jit(nopython=True)
def monod_growth_law(c, kn, gnmax=0.034):
    return gnmax * c / (kn + c)

def nutrient_integral(gmax,B,n0,delta_t,tmax):

    niter = int(tmax/delta_t)
    time = np.zeros(niter)
    n = np.zeros(niter)
    b = np.zeros(niter)
    
    n[0] = n0
    time[0] = 0.0
    b[0] = B

    for ite in range(1,niter):
        n[ite] =  n[ite -1] - (1.0)*b[ite-1]*monod_growth_law(n[ite-1],n0/5.0,gmax)*delta_t
        b[ite] = b[ite-1] + monod_growth_law(n[ite-1],n0/5.0,gmax)*b[ite-1]*delta_t
        time[ite] = time[ite-1] + delta_t
        
    return time,n,b






