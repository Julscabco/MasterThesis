# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:35:15 2023

@author: Usuario
"""

import numba

import numpy as np

from Utils import get_linear_burst, monod_growth_law, minutes_to_taul, taul_to_minutes



@numba.njit
def system_evolution_nlg(delta_t, Nb0, B, P, N, k, nd, nu, n0, doubling_time, twindow, Nimax, Nimax_std, niter, time_window_track=True, constant_burst = False, ctburstsize=50.0):
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
    ninfect = np.zeros(Nb0, dtype = np.int64)
    lysis_times = np.zeros(Nb0, dtype=np.float64)

    
    # calculate number of iterations corresponding to time window
    iter_window = int(twindow/delta_t)
    std_nimax = Nimax_std
    
    # Create a matrix with size (number of bacteria, number of iterations)
    last_infections = np.zeros((Nb0,iter_window), dtype=np.int64)
    
    
    sizes = np.zeros(Nb0,dtype=np.int64)
    for j in range(0,Nb0):
        sizes[j] = np.random.randint(0,max_size)

    # Initialization of vectors that grow with time steps
    time = []
    phages = []
    healthy_bacteria = []
    infected_bacteria = []
    dead_bacteria = []
    nut = []
    average_mosi = []
    window_mosi = []
    burst_events_md = []
    burstsize_events_md = []


    # We store the first values of the concentratioin of bacteria and phage
    healthy_bacteria.append(float(Nb0))
    infected_bacteria.append(0.0)
    dead_bacteria.append(0.0)
    phages.append(float(Np0))
    nut.append(nutrients)
    burst_events_md.append(0)
    burstsize_events_md.append(0.0)
    
    # We initialize time
    t = 0.0
    time.append(t)

    # We initialize the number of phages
    Np = Np0

    # probability of new timer molecules for infected bacteria
    prob_new_timer = k * delta_t

    n = 0
    burst_md_counter = 0


    while n<niter:
        
        burstsize_md_counter = 0

        # We iterate over the already infected cells
        infected_indices = np.where(states == 1)[0]
        
        prob_infection = nu * delta_t * Np / V
        if n==0: 
            print(prob_infection)
        
        new_infections = np.zeros(len(ninfect),dtype=np.int64)
        newborn_bacteria = 0

        
        for ii in infected_indices:
            r = np.random.rand()
            if r < prob_new_timer:
                Nt[ii] += np.int8(1)
                if Nt[ii] >= N:
                    tau_mins = taul_to_minutes(t,N,k) - first_infection_time[ii]
                    if constant_burst == False:
                        burst = int(get_linear_burst(tau_mins, 15.0, 10.0))
                    elif constant_burst == True:
                        burst = ctburstsize
                    Np = Np + burst
                    states[ii] = 2
                    lysis_times[ii] = tau_mins
                    

            elif prob_new_timer < r < prob_new_timer + prob_infection:
                
                if time_window_track == True:
                    relevant_infections = np.sum(last_infections[ii,:])
                elif time_window_track == False:
                    relevant_infections = ninfect[ii]
                
                if  relevant_infections>= np.random.uniform(low=Nimax-(Nimax_std/2.0),high=Nimax+(Nimax_std/2.0)):
                    tau_mins = taul_to_minutes(t,N,k) - first_infection_time[ii]
                    if constant_burst==False:
                        burst = int(get_linear_burst(tau_mins, 15.0, 10.0))
                    elif constant_burst == True:
                        burst = ctburstsize
                    Np = Np + burst
                    burstsize_md_counter += burst
                    states[ii] = 2
                    lysis_times[ii] = tau_mins
                    burst_md_counter += 1
                   
                   
   
                else:   
                    new_infections[ii] = 1
                    ninfect[ii] += 1
                
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
        
        Y = 1.00
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
                first_infection_time[i] = taul_to_minutes(t,N,k)
                new_infections[i] = 1
                ninfect[i] += 1


            elif prob_infection< r < (prob_infection + prob_growth):
                
                if sizes[i] < int(7): 
                    if nutrients>0.0:
                        sizes[i] = sizes[i] + 1
                        nutrients = nutrients - 1.0/(max_size*Y)
                
                if sizes[i] >= int(7):
                    newborn_bacteria += 1
                    # We re-initialize all vectors that depend on Nb
                    states = np.append(states, 0)
                    first_infection_time = np.append(first_infection_time, 0.0)
                    Nt = np.append(Nt, 0)
                    sizes = np.append(sizes, 0)
                    lysis_times = np.append(lysis_times,0.0)
                    ninfect = np.append(ninfect, 0)
                    
                    sizes[i] = 0

        

        last_infections = last_infections[:,1:]
        # add the new infections array to the right of the matrix
        last_infections = np.column_stack((last_infections,new_infections))
        alive_indices = np.hstack((infected_indices, non_infected_indices))
        if len(alive_indices)>0:
            window_mosi.append(np.mean(np.sum(last_infections[alive_indices,:],axis=1)))
        else:
            window_mosi.append(0.0)
        # add as many rows as new bacteria there are
        last_infections = np.vstack((last_infections,np.zeros((newborn_bacteria,iter_window),dtype=np.int64)))

        

        t = t + delta_t
        phages.append((Np))
        time.append(t)
        healthy_bacteria.append(len(np.where(states == 0)[0]))
        infected_bacteria.append(len(np.where(states == 1)[0]))
        dead_bacteria.append(len(np.where(states == 2)[0]))
        nut.append((nutrients))
        burst_events_md.append(burst_md_counter)
        burstsize_events_md.append(burstsize_md_counter)
        
        alive_cells = np.where(states == 1)[0]
        if len(alive_cells)>0:
            average_mosi.append(np.mean(ninfect[alive_cells]))
        else:
            average_mosi.append(0.0)

        n = n + 1

    return time, phages, healthy_bacteria, infected_bacteria, dead_bacteria, nut, average_mosi, window_mosi, burst_events_md, burstsize_events_md, first_infection_time, lysis_times, ninfect


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
    Nimax = 10
    
    # time parameters in minutes
    doubling_time_mins = 20.0
    time_window_mins = 2.0
    
    # Collapse time in hours
    collapse_time_hours = 2.5
    
    # Time parameters in tau_l units
    doubling_time = minutes_to_taul(doubling_time_mins, N, k)
    time_window = minutes_to_taul(time_window_mins, N, k)
    collapse_time = minutes_to_taul((collapse_time_hours*60.0),N,k)
    
    niter = 20000

    time_vector, nphages, nhbacteria, nibacteria, ndbacteria, nnutrients, avg_ninf, avg_burst, burstmd, burstsizemd, firstinf, lysist, ninfect = system_evolution_nlg(delta_t, Nb0, B, P, N, k, nd, nu, n0, doubling_time, time_window, Nimax, Nimax/4.0, niter)
