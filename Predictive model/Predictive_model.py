# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:35:15 2023

@author: Usuario
"""

import numba

import numpy as np

from Utils import get_linear_burst, monod_growth_law, minutes_to_taul, taul_to_minutes



@numba.njit
def system_evolution_nlg(delta_t, Nb0, B, P, N, k, nd, nu, n0, doubling_time, niter):
    # Volume of the system
    V = np.float64(Nb0 / B)

    # How many phages do we have in this voalume according to the initial concentration?
    Np0 = np.float64(P * V)

    # how many units of nutrients do we have?
    nutrients = n0 * V
    gnmax = np.log(2.0) / doubling_time
    kn = n0 / 5.0
    
    # What is the maximum size that cells can achieve?
    max_size = int(10)

    # Initialization of vectors that contain information about bacteria
    states = np.zeros(Nb0, dtype=np.int64)
    first_infection_time = np.zeros(Nb0, dtype=np.float64)
    Nt = np.zeros(Nb0, dtype=np.int64)
    infection_times = np.zeros(1, dtype=np.float64)
    infected_cell = np.zeros(1, dtype=np.int64)
    lysis_times = np.zeros(Nb0, dtype=np.float64)
    
    sizes = np.zeros(Nb0,dtype=np.int64)
    for j in range(0,Nb0):
        sizes[j] = np.random.randint(0,max_size)

    # Initialization of vectors that grow with time steps
    time = []
    phages = []
    bacteria = []
    nut = []


    # We store the first values of the concentratioin of bacteria and phage
    bacteria.append(float(Nb0))
    phages.append(float(Np0))
    nut.append(nutrients)

    

    # We initialize time
    t = 0.0
    time.append(t)

    # We initialize the number of phages
    Np = Np0

    # probability of new timer molecules for infected bacteria
    prob_new_timer = k * delta_t

    
    n = 0

    while n<niter:

        # We iterate over the already infected cells
        infected_indices = np.where(states == 1)[0]
        
        prob_infection = nu * delta_t * Np / V
        
        for ii in infected_indices:
            r = np.random.rand()
            if r < prob_new_timer:
                Nt[ii] += np.int8(1)
                if Nt[ii] >= N:
                    tau = t - first_infection_time[ii]
                    tau_mins = taul_to_minutes(tau,N,k)
                    Np = Np + int(get_linear_burst(tau_mins, 15.0, 10.0))
                    states[ii] = 2
                    
                    lysis_times[ii] = tau_mins
                    

            elif prob_new_timer < r < prob_new_timer + prob_infection:
                
                infection_times = np.append(infection_times,t)
                infected_cell = np.append(infected_cell,ii)
                
                if Np > 0:
                    Np = Np - 1
                else:
                    Np = 0

                if Nt[ii] < nd:
                    Nt[ii] = 0
                elif Nt[ii] >= nd:
                    Nt[ii] = Nt[ii] - nd



        # We calculate concentration of nutrients and growth rate
        nn = float(nutrients)/V
        alpha = monod_growth_law(c=nn, kn=kn, gnmax=gnmax)
        if alpha <= 0.0:
            alpha = 0.0
        
        Y = 1.35
        # probability of growing
        prob_growth = alpha*delta_t*max_size*Y
    

        # We iterate over the non infected cells to see if they get infected
        non_infected_indices = np.where(states == 0)[0]

        for i in non_infected_indices:
            r = np.random.rand()
            if r < prob_infection:
                # They are infected now
                states[i] = 1
                Np = Np - 1
                # We store their infection time
                first_infection_time[i] = t
                infection_times = np.append(infection_times,t)
                infected_cell = np.append(infected_cell,i)

            elif prob_infection< r < (prob_infection + prob_growth):
                
                if sizes[i] < max_size: 
                    if nutrients>0.0:
                        sizes[i] = sizes[i] + 1
                        nutrients = nutrients - 1.0/(max_size*Y)
                
                if sizes[i] >= max_size:
                    # We re-initialize all vectors that depend on Nb
                    states = np.append(states, 0)
                    first_infection_time = np.append(first_infection_time, 0.0)
                    Nt = np.append(Nt, 0)
                    sizes = np.append(sizes, 0)
                    lysis_times = np.append(lysis_times,0.0)
                    
                    sizes[i] = 0



        t = t + delta_t
        phages.append((Np))
        time.append(t)
        bacteria.append((len(np.where(states == 0)[0]) + len(np.where(states == 1)[0])))
        nut.append((nutrients))

        n = n + 1

    return time, phages, bacteria, nut, lysis_times, infection_times, infected_cell


if __name__ == '__main__':
    
    delta_t = 0.01
    Nb0 = 1000
    B = 1.0e6
    P = 1.0e8
    N = 60
    k = 6.0
    nd = 30
    nu = 5.0e-10
    n0 = 1.0e7

    # The doubling time in minutes
    doubling_time_mins = 20.0
    doubling_time = minutes_to_taul(doubling_time_mins, N, k)
    
    niter = 80000

    time_vector, nphages, nbacteria, nnutrients, lysist, ninfect = system_evolution_nlg(delta_t, Nb0, B, P, N, k, nd, nu, n0, doubling_time, niter)
