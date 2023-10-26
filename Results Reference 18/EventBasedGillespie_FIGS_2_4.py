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


""" Functions to use in Gillespie algorithm """



import os

import numpy as np
import matplotlib.pyplot as plt


def get_time_step(rates):
    total_rate = np.sum(rates)
    tau = (-1.0/float(total_rate))*np.log(np.random.uniform())
    return tau


def get_next_reaction(rates):
    normalized_rates = rates/np.sum(rates)
    random_number = np.random.uniform()
    for i in range(0, len(normalized_rates)):
        if random_number < np.sum(normalized_rates[:i+1]):
            return i
        
def integrate_exponential(nu,B,initial_time,end_time):
    a = np.exp(-nu*B*initial_time)
    b = np.exp(-nu*B*end_time)
    return (a-b)*(1.0/(nu*B))

def get_initial_number_of_phage(mosi,nu,B,ti_superinfection,len_superinfection):
    P_0 = mosi*B/integrate_exponential(nu,B,ti_superinfection,ti_superinfection + len_superinfection)
    return P_0

""" Parameters """


def one_bacteria_lysis(n=10,nd=8,k=10.0, MOSI=0.0, ti_superinfection=0.6, len_superinfection=5.0/23.0,adsorption_rate=5.0e-9,B=2.0e7,plot=False):

    # parameters
    #nu = 0.1
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
    
    P0 = get_initial_number_of_phage(mosi,adsorption_rate,B,ti_superinfection,ti_superinfection + len_superinfection)

    #nu = adsorption_rate*P0*integrate_exponential(adsorption_rate,B,ti_superinfection,ti_superinfection+len_superinfection)
    nu = adsorption_rate*P0
    # Vector with the rates of each process

    rates_vect = np.array([k, nu*Np])

    phages_superinfecting = int(P0/B)
    
    while Nt <= n:

        # If the time arrives to the initial superinfection time, we add phage
        if ti_superinfection < time <(ti_superinfection+len_superinfection) and superinfection == False:
            Np = Np + phages_superinfecting
            superinfection = True
            rates_vect = np.array([k, nu*Np])

        if ti_superinfection+len_superinfection < time and superinfection == True:
            superinfection = False
            Np = 0.0
            rates_vect = np.array([k, nu*Np])

        time_interval = get_time_step(rates_vect)
        reaction = get_next_reaction(rates_vect)

        if reaction == 0:
            # A new timer is born
            Nt = Nt + 1.0
            Np = Np
            time = time + time_interval
        elif reaction == 1 and superinfection == True:
            # New infection: free phages is reduced by 1 and Nt by nd
            Np = Np - 1.0
            if Nt < nd:
                Nt = 0.0
            elif Nt >= nd:
                Nt = Nt - nd
            time = time + time_interval

        rates_vect = np.array([k, nu*Np])

        if plot == True:
            time_vector.append(time)
            phage.append(Np)
            timers.append(Nt)

    if plot == True:
        phage = np.array(phage)
        time_vector = np.array(time_vector)
        timers = np.array(timers)

        """ Plottings """

        plt.figure()
        plt.plot(time_vector, phage, label='Phage', color='red', marker='o')
        plt.plot(time_vector, timers, label='Timer molecules',
                 color='green', marker='o')
        plt.legend(loc='best')
        plt.ylabel('Number')
        plt.xlabel('Time (a.u)')

    return time


""" THE FOLLOWING FUNCTION CALLS ONE_BACTERIA_LYSIS FUNCTION NDATA TIMES AND
STORES THE DATA OF THE LYSIS TIME IN A CERTAIN LOCATION """


def simulation(ndata,n,nd,kk,mosi,ti_superinfection,len_superinfection,ff,folder):
    name_of_file = str(n)+'_'+str(nd)+'_'+str(kk)+'_'+ str(mosi)+'_'+str(ti_superinfection)+'_'+str(len_superinfection)+'_'+str(ndata)+'.txt'
    file_path=os.path.join(os.getcwd(),'output_files',folder)
    
    if not os.path.exists(file_path):
        os.mkdir(file_path,0o666)
    
    file_path = os.path.join(file_path,name_of_file)
        
    f = open(file_path, "a")

    for i in range(0, ndata):
        
        tau = (1.0/float(ff))*one_bacteria_lysis(n=n,nd=nd,k=kk,MOSI=mosi,ti_superinfection=ti_superinfection*ff, len_superinfection=len_superinfection*ff)
        f.write(str(tau)+'\n')

    f.close()
    return


""" GRAPH 1: REPRESENTATION OF HISTOGRAM OF TAU AS A FUNCTION OF TIME AFTER 
FRIST INFECTION FOR DIFFERENT VALUES OF THE PARAMETERS """

ndata = 20001
mosi = 1.6
ti_superinfection = 15.0
len_superinfection = 3.0
k = 6.0
n = 60
nd_vect = np.linspace(0,n,10,dtype=np.int8)

conversion_factor = (float(n)/float(k))/23.0

output_folder = str(n)+'_'+str(int(k))+'_'+ str(mosi)+'_'+str(int(ti_superinfection))+'_'+str(int(len_superinfection))+'_'+str(ndata)
#output_folder='test'
for nd in nd_vect:
    simulation(ndata,n,nd,k,mosi,ti_superinfection,len_superinfection,conversion_factor,output_folder)
    

""" GRAPH 2 (FIGURE 4 REF 18): REPRESENTATION OF HISTOGRAM OF TAU AS A FUNCTION OF TIME AFTER 
PRIMARY INFECTION FOR DIFFERENT LENGTHS AND START TIMES OF SUPERINFECTION """
"""
ndata = 20000
mosi = 10.0
ti_superinfection_vect = np.array([15.0,14.0,13.0])
len_superinfection_vect = np.array([1.5,3.0,4.5])
k = 6.0
n = 60
nd = 33

conversion_factor = (float(n)/float(k))/23.0
output_folder = str(n)+'_'+str(nd)+'_'+str(int(k))+'_'+ str(mosi)
print(output_folder)
#for tis,lens in zip(ti_superinfection_vect,len_superinfection_vect):
    
    #simulation(ndata,n,nd,k,mosi,tis,lens,conversion_factor,output_folder)
    
"""