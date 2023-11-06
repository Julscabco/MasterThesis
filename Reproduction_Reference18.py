# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 09:06:06 2023

@author: Usuario
"""

""" THIS FILE CONTAINS THE CODE TO REPRODUCE THE RESULTS IN REFERENCE 18 """

# For the moment I have a function that takes the initial number of phages and 
# the initial number of timer molecules and the values of the parameters and
# calculates the evolution of the system, only from when the first infection takes 
# place

# We now control when the first infection takes place, and we want to control when
# the superinfection takes place: at some point of the time, more free phages are
# added to the system

# We leave the phages for some finite aount of time and then we just filtrate 
# the free phages in the environment of the bacteria. We have to be able to control
# the time at which we filtrate the free phages

# The lysis of the bacteria should also be included, and the phages inside the 
# bacteria should burst

# GOAL: We want to obtain the number of free fages as a function of time after
# the primary infection

# STEPS

# 1. Simulation starts recording data from the first moment, because simulation
# starts when we do the primary infection. OK

# 2. Bacteria have to lyse when timer molecules have reached the maximum timer
# molecules (that should be 10). OK

# 3. At some point t_superinfection, the number of phages increases OK

# 4. After some time of superinfection, then the number of free phages goes to 0 OK

# 5. Now we need to compare the "time" we have in the graphs with real minutes

# 5. We do a function that for a certain combination of parameters, it 
# calculates the evolution and stores the time of lysis. OK

# 6. If we repeat step 5 but with, for example different values of the time of
# superinfection, we have figure 2 from the paper. falta: Fer vector de diferents 
# MOSI, vector de diferents temps inicials i vector de diferents duracions
# Provar com queda amb diferents valors de nd

# Com fer aixo? necessito numero de bacteries mortes en funcio del temps: he
# de calcular tau per a moltes bacteries i llavors fer l'histograma.FET



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

def one_bacteria_lysis(nd=8,MOSI=0.0,ti_superinfection = 0.6, len_superinfection = 5.0/23.0,plot=False):
    
    # parameters
    k = 10.0
    nu = 5.0
    n = 10
    initial_Nt = 0.0
    initial_Np = 0.0
    initial_time = 0.0
    
    
    # To control if there is superinfection or not
    superinfection = False

    # Vectors to store evolution

    time_vector = []
    phage = []
    timers = []

    Nt = initial_Nt
    Np = initial_Np
    time = initial_time

    # Vector with the rates of each process

    rates_vect = np.array([k,nu*Np])
    
    MOSI = int(MOSI)


    # Loop until Nt = 10
    while Nt<=n:
        
        # If the time arrives to the initial superinfection time, we add phage
        if ti_superinfection<time and superinfection == False:
            Np = Np + MOSI
            superinfection = True
            rates_vect = np.array([k,nu*Np])
            
        if ti_superinfection+len_superinfection<time and superinfection == True:
            superinfection = False
            Np = 0.0
            rates_vect = np.array([k,nu*Np])
            
        time_interval = get_time_step(rates_vect)
        reaction = get_next_reaction(rates_vect)
    
        if reaction==0:
            # A new timer is born
            Nt = Nt + 1.0
            Np = Np
            time = time + time_interval
        elif reaction==1 and superinfection==True:
            # New infection: free phages is reduced by 1 and Nt by nd
            Np = Np - 1.0
            if Nt < nd:
                Nt = 0.0
            elif Nt >= nd:
                Nt = Nt - nd
            time = time + time_interval
            
        rates_vect = np.array([k,nu*Np])
    
        if plot==True:
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


""" MODIFICAR AQUESTA FUNCIO PER PODER PASSAR-LI ELS DIFERENTS PARAMETRES AMB
ELS QUALS VOLEM CRIDAR LA FUNCIO ONE_BACTERIA"""


def get_tau_pd(ndata,mosi=7.0,nd=2.0,ti_superinfection=15.0/23.0,len_superinfectio=5.0/23.0,plot_histogram=False):
    
    tau_vector = []
    
    for i in range(0,ndata):
        
        tau = 23.0*one_bacteria_lysis(nd=nd,ti_superinfection=15.0/23.0,len_superinfection=5.0/23.0,MOSI=mosi)
        tau_vector.append(tau)
        f.write(str(tau))
        
    if plot_histogram == True:
        
        tau_vector = np.array(tau_vector)
        
        # Calculation and normamlization of the histogram
        hist,bins=np.histogram(tau_vector,bins=60)
        bins = (bins[1:] + bins[:-1])/2.0
        binwidth = bins[2]-bins[1]
        hist = hist/(len(tau_vector)*binwidth)
        
        # The plot 
        fig_hist,ax_hist = plt.subplots(figsize=(10,6))
        
        ax_hist.scatter(bins,hist,label= 'ndata='+ str(ndata))
        
        ax_hist.set_title(r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=\frac{15}{23}$' + ' ,' + r'${\tau}_{s}=\frac{5}{23}$' + ' ,' + 'n=' + str(nd),fontsize=20)
        ax_hist.set(xlim=(20.0,70.0))
        ax_hist.set_xlabel(r'${\tau}_{l} (a.u)$',fontsize=18,labelpad=8)
        ax_hist.set_ylabel('Normalized Frequency',fontsize=18,labelpad=8)
        ax_hist.tick_params(axis='both', which='major', labelsize=20);
        ax_hist.legend(loc='best',fontsize=18)
        
    return None



mosi = 7.0
nd = 2.0
ti_superinfection = 15.0/23.0
len_superinfection = 5.0/23.0

f = open("file1.txt","a")
f.write("MOSI="+str(mosi)+"n="+str(nd)+"ti_sup="+str(ti_superinfection)+"len_sup="+str(len_superinfection)+'\n')
get_tau_pd(20000,mosi,nd,ti_superinfection,len_superinfection,plot_histogram=False)
f.close()

""" FER UNA FUNCIO QUE CALCULI L'HISTOGRAMA PER A DIFERENTS VALORS DE ND, MOSI,
TI_SUPERINFECTION I DURACIO DE LA SUPERINFECTION """