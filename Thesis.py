# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math

""" Functions to use in Gillespie algorithm """


def get_time_step(rates):
    total_rate = np.sum(rates)
    tau = (-1.0/float(total_rate))*np.log(np.random.uniform())
    return tau


def get_next_reaction(rates):
    normalized_rates = rates/np.sum(rates)
    random_number = np.random.uniform()
    for i in range(0,len(normalized_rates)):
        if random_number < np.sum(normalized_rates[:i+1]):
            return i
        
        

""" Parameters """

def one_bacteria_lysis(k=3.0,nu=1.0,n=3,initial_Nt=0.0,initial_Np=1.0,plot=False):
    
    initial_time = 0.0

    # Vectors to store evolution

    time_vector = []
    phage = []
    timers = []

    # We want simulation to start after the first infection
    first_infection = False

    Nt = initial_Nt
    Np = initial_Np
    time = initial_time

    # Vector with the rates of each process

    rates_vect = np.array([k,nu*Np])


    # Loop until Nt = 3
    while Nt<n or first_infection==False:
    
        time_interval = get_time_step(rates_vect)
        reaction = get_next_reaction(rates_vect)
    
        if reaction==0:
            # A new timer is born
            Nt = Nt + 1.0
            Np = Np
            if first_infection==True:
                time = time + time_interval
        elif reaction==1:
            # New infection: free phages is reduced by 1 and Nt = 0
            #Np = Np - 1.0
            Nt = 0
            if first_infection==True:
                time = time + time_interval
            first_infection=True
    
        rates_vect = np.array([k,nu*Np])
    
        if first_infection == True and plot==True:
            time_vector.append(time)
            phage.append(Np)
            timers.append(Nt)
    
    if plot==True:
        phage = np.array(phage)
        time_vector = np.array(time_vector)
        timers = np.array(timers)
    
        """ Plottings """

        plt.figure()
        plt.plot(time_vector,phage,label='Phage',color='red',marker='o')
        plt.plot(time_vector,timers,label='Timer molecules',color='green',marker='o')
        plt.legend(loc='best')
        plt.ylabel('Number')
        plt.xlabel('Time (a.u)')
    
    return time

def erlang(t,mu,k):
    
    return mu*((mu*t)**(k-1))*np.exp(-mu*t)/math.factorial(k-1)
        


def get_tau_pd(ndata,k=3.0,nu=1.0,n=3,initial_Nt=0.0,initial_Np=1.0,plot_histogram=False):
    
    tau_vector = []
    
    for i in range(0,ndata):
        
        tau_vector.append(one_bacteria_lysis())
        
        
    if plot_histogram == True:
        
        tau_vector = np.array(tau_vector)
        
        # Calculation and normamlization of the histogram
        hist,bins=np.histogram(tau_vector,bins=30)
        bins = (bins[1:] + bins[:-1])/2.0
        binwidth = bins[2]-bins[1]
        hist = hist/(len(tau_vector)*binwidth)
        
        # Erlang distribution 
        erlang_dist = [erlang(i,3.0,3.0) for i in bins]
    
    
        # The plot 
        fig_hist,ax_hist = plt.subplots(figsize=(10,6))
        
        ax_hist.plot(bins,hist,label= r'${\tau}_{l}$ histogram')
        ax_hist.plot(bins,erlang_dist,label='Erlang(3,3)')
        
        ax_hist.set(xlim=(0.0,3.5))
        ax_hist.set_xlabel(r'${\tau}_{l} (a.u)$',fontsize=18,labelpad=8)
        ax_hist.set_ylabel('Normalized Frequency',fontsize=18,labelpad=8)
        ax_hist.tick_params(axis='both', which='major', labelsize=20);
        ax_hist.legend(loc='best',fontsize=18)
        
    return None


get_tau_pd(5000,plot_histogram=True)


""" --------------------- <Tau> vs nu plot ---------------------------------"""
    
    
def get_tau_vs_nu(ndata,nuvect,histogram=False):
    
    taumean = []
    taustd = []
    
    if histogram:
        fig2,ax2 = plt.subplots(nrows=int(len(nuvect)/2.0),ncols=2,figsize=(13,18))
        col = 0
    
    for u in range(0,len(nuvect)):
        
        tau = []
        for i in range(0,ndata):
            tau.append(one_bacteria_lysis(nu=nuvect[u],initial_Np=4.0))
        
        tau = np.array(tau)
        taumean.append(np.mean(tau))
        taustd.append(np.var(tau))
        
        if histogram:
            hist,bins=np.histogram(tau,bins=20)
            bins = (bins[1:] + bins[:-1])/2.0
            binwidth = bins[2]-bins[1]
            hist = hist/(len(tau)*binwidth)

            ax2[math.floor(u/2.0),col].plot(bins,hist)
            ax2[math.floor(u/2.0),col].set_title(r'$\eta$ = ' + f'{nuvect[u]:.3}',fontsize=15)
            ax2[math.floor(u/2.0),col].tick_params(axis='both', which='major', labelsize=12)
            ax2[math.floor(u/2.0),col].set_ylabel('Frequency',fontsize=14,labelpad=8)
            ax2[math.floor(u/2.0),col].set_xlabel(r'$\tau$ (${\tau}_{l}$)',fontsize=14,labelpad=8)
            ax2[math.floor(u/2.0),col].set_yscale('log')
            if col==0:
                col=1
            elif col==1:
                col=0
            
    fig, ax = plt.subplots()
    ax.errorbar(nuvect,taumean, yerr=(taustd/np.sqrt(ndata)), fmt='o', ecolor = 'black',color='red')
    #ax.scatter(nuvect,taumean)
    ax.set_xlabel(r'$\eta$ ( $\frac{ml}{{\tau}_{l}}$ )',fontsize=14,labelpad=8)
    ax.set_ylabel(r'$\langle \tau \rangle$ (${\tau}_{l}$)',fontsize=14,labelpad=8)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_yscale('log')
    


            
        

    
    return taumean,taustd

nu_vect = np.logspace(0.0,0.9,6)



tau_mean,tau_std = get_tau_vs_nu(1000,nu_vect,histogram=True)




 
        
    
    
    
    
    









