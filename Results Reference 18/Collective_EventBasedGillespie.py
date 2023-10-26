# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:13:22 2023

@author: Usuario
"""

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


def collective_bacteria_lysis(Nb=100, n=10, nd=8, k=10.0, MOSI=0.0, ti_superinfection=0.6, len_superinfection=5.0/23.0, nu=0.7, burst_size=180):

    # This function will return the lysis time of the bacteria (Nb) of the system.
    # Inputs: Number of bacteria, MOSI, nu, n,k, nd, ti_superinfection,len_superinfection

    Np = 0.0
    Nt = np.zeros(Nb, dtype=np.int8)
    time = 0.0

    # To control if there is superinfection or not
    superinfection = False

    # We initialize arrays that we will fill
    death_times = []

    # We calculate the initial number of bacteria we have to put in the system
    phages_superinfecting = int(
        (((MOSI*2.0e7)/np.exp(-0.1*len_superinfection))/2.0e7))*Nb
    print(phages_superinfecting)
    # The rates for the Gillespie algorithm. They are the same for every bacteria
    rates_vect = np.array([k, nu*Np])

    while any(Nt <= n):

        # If the time arrives to the initial superinfection time, we add phage
        if ti_superinfection < time < (ti_superinfection+len_superinfection) and superinfection == False:
            Np = Np + phages_superinfecting
            superinfection = True
            rates_vect = np.array([k, nu*Np])

        if ti_superinfection+len_superinfection < time and superinfection == True:
            superinfection = False
            Np = 0.0
            rates_vect = np.array([k, nu*Np])

        time_interval = get_time_step(rates_vect)
        # We calculate which is the next reaction occuring
        reaction = get_next_reaction(rates_vect)

        bacteria = np.random.randint(0, Nb)

        if reaction == 0:
            # A new timer is born
            Nt[bacteria] = Nt[bacteria] + 1
            Np = Np

        elif reaction == 1 and superinfection == True:
            # New infection: free phages is reduced by 1 and Nt by nd
            Np = Np - 1.0
            if Nt[bacteria] < nd:
                Nt[bacteria] = 0.0
            elif Nt[bacteria] >= nd:
                Nt[bacteria] = Nt[bacteria] - nd

        # We update the time:
        time = time + time_interval

        # We check if any bacteria has died
        if any(Nt >= n):
            index = np.where(Nt >= n)[0]
            Nt = np.delete(Nt, index)
            Nb = Nb - 1
            Np = Np + burst_size
            death_times.append(time)

        # We update the rates
        rates_vect = np.array([k, nu*Np])

    death_times = np.array(death_times)

    return death_times


mosi = 1.6
ti = 15.0
length = 3.0
k = 6.0
n = 60
nd = 33
Nb = 70
nu = 5e-9*2.0e8*np.exp(-0.1*length)



burst_size = 0

ff = ((float(n)/float(k)))/(23.0)

t = (1/ff)*collective_bacteria_lysis(Nb=Nb, n=n, nd=nd, k=k, MOSI=mosi,
                                     ti_superinfection=ff*ti, len_superinfection=ff*length, nu=nu, burst_size=burst_size)


output_folder = 'CollectiveGillespie_data'
file_name = str(Nb)+'_'+str(n)+'_'+str(int(k))+'_'+str(nd)+'_' + \
    str(mosi)+'_'+str(int(ti))+'_'+str(int(length))+'_'+str(burst_size)
print(file_name)
file_path = os.path.join(os.getcwd(), output_folder, file_name)

f = open(file_path, "a")

for i in range(0, len(t)):

    f.write(str(t[i])+'\n')

f.close()

# TIME NORMALIZATION: DIVIDE BY Nb, DIVIDE BY K/N AND MULTIPLY BY 23
