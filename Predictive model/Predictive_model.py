# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:55:23 2023

@author: Usuario
"""
import os 

import numba
import numpy as np
import matplotlib.pyplot as plt


def get_linear_burst(t,tlim,a):
    if t<tlim:
        return 1.0
    elif t>=tlim:
        return a*t

@numba.jit    
def system_evolution(delta_t,Nb0,B,Np0,alpha,mmax,N,k,nd,nu,niter):
    
    # Volume of the system
    V = float(Nb0)/float(B)
    
    # Initialization of vectors that contain information about bacteria
    states = np.zeros(Nb0,dtype=np.int8)
    infection_time = np.zeros(Nb0, dtype=np.float16)
    sizes = np.ones(Nb0, dtype=np.float16)
    Nt = np.zeros(Nb0, dtype=np.int8)
    
    
    # Initialization of vectors that grow with time steps
    time = []
    phages = []
    bacteria = []
    
    #We store the first values of the concentratioin of bacteria and phage
    bacteria.append(float(Nb0)/float(V))
    phages.append(float(Np0)/float(V))
    
    # We initialize time
    t = 0.0
    time.append(t)
    
    # We initialize the number of phages
    Np = Np0
    
    # probability of new timer molecules for infected bacteria
    prob_new_timer = k*delta_t
    
    iteration = 0
    
    for n in range(0,niter):
        
        # we count the infected and uninfected bacteria
        n_infected = len(np.where(states==1)[0])
        n_uninfected = len(np.where(states==0)[0])
        Nb = n_infected + n_uninfected
        
        # Probability of infection at this time step
        prob_infection = (nu*delta_t*(float(Np)/V))
        prob_superinfection = (nu*delta_t*(float(Np)/V))
        sum_probs = prob_superinfection + prob_new_timer
        
        # We iterate over the already infected cells
        for i in np.where(states==1)[0]:
            r = np.random.uniform()
            if r < prob_new_timer:
                Nt[i] = Nt[i] + 1
                #print('New timer')
                if Nt[i] >= N:
                    tau = t - infection_time[i]
                    Np = Np + int(get_linear_burst(tau,15.0,1.0))
                    states[i] = 2
                    print('death',iteration,Nt[i],sizes[i],tau,int(get_linear_burst(tau,15.0,1.0)))
                    
                    # We delete the dead cells from all the variables dependent on the number of bacteria
                    np.delete(states,i)
                    np.delete(infection_time,i)
                    np.delete(sizes,i)
                    np.delete(Nt,i)
                    
            elif prob_new_timer < r < sum_probs:
                if Np > 0:
                    Np = Np - 1
                else:
                    Np = 0
                    
                if Nt[i] < nd:
                    Nt[i] = 0 
                elif Nt[i] >= nd:
                    Nt[i] = Nt[i] - nd
            
        
        # We iterate over the non infected cells to see if they get infected
        for i in np.where(states==0)[0]:
            r = np.random.uniform()
            if r<prob_infection:
                #print('first infection')
                # They are infected now
                states[i] = 1
                Np = Np-1
                # Wr store their infection time
                infection_time[i] = t
            else:
                # If the cell is too big it divides otherwise just grows
                if sizes[i]>=mmax:
                    #print('division',iteration)
                    # The cell that has divided has a new size
                    sizes[i] = 1.0
                    # We initialize all vectors that depend on Nb
                    sizes = np.append(sizes,1.0)
                    states = np.append(states,0)
                    infection_time = np.append(infection_time,0.0)
                    Nt = np.append(Nt,0)
                    
                elif sizes[i]<mmax:
                    # We make the cell grow on size
                    sizes[i] = sizes[i] + alpha*delta_t
                    
        
        t = t + delta_t
        phages.append(float(Np)/V)
        time.append(t)
        bacteria.append(float(len(np.where(states==0)[0]) + len(np.where(states==1)[0]))/V)

        
        iteration = iteration + 1

                    
    return time,phages,bacteria


delta_t = 0.01
initial_nbacteria = 100
bacteria_concentration = 2.0e7
initial_nphages = 10
alpha = 100.0
mmax = 2000.0
N = 60
k = 6.0
nd = 30
nu = 5.0e-9

niter=500000

time_vector,nphages,nbacteria = system_evolution(delta_t, initial_nbacteria, bacteria_concentration,
                                                 initial_nphages, alpha, mmax, N, k, nd, nu, niter)

save_plot = False

fig, ax = plt.subplots(2,1,figsize=(10,6))

ax[0].plot(time_vector,nphages,label='Phage')
ax[0].set_title(r'$\Delta t=$' + str(delta_t) + ', ' + r'${N}_{iter}$=' + str(niter)
          +', ' + r'$\alpha=$' + str(alpha)+ ', ' + r'${N}_{b0}=$' +str(initial_nbacteria)+
          ', ' + r'${N}_{p0}=$' + str(initial_nphages))
ax[0].set_ylabel('Number of phage')

ax[1].plot(time_vector,nbacteria,label='Bacteria',color='orange')

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Number of bacteria')


if save_plot == True:
    figure_path = os.path.join(os.getcwd(),'FIGURES','First figures')
    figure_name = str(delta_t)+'_'+str(niter)+'_'+str(initial_nbacteria)+'_'+str(initial_nphages)+'.png'
    fig.savefig(os.path.join(figure_path,figure_name))


