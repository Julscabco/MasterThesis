# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:05:56 2023

@author: Usuario
"""

import os
 
import numpy as np

def get_burst_size(t, tlim, a):
    if t < tlim:
        return 0.0
    elif t >= tlim:
        return a * t - a*tlim


def Nb_bacteria_lysis_2(deltat, Nb, n=60, nd=33, k=6.0, MOSI=7.6, ti_superinfection=15.0,len_superinfection=3.0):
    
    nu = 5e-9
    B = 1.0e7
    
    # Volume of the system
    V = Nb/B

    
    # we calculate the initial concentration of phages according to the desired MOSI
    P0 = (MOSI*2.0e7)/(np.exp(-0.1*len_superinfection))
    print(P0)
    
    # We calculate the corresponding number of phages
    Np = int(P0*V)
    print(Np)
    
    
    prob_new_timer = deltat*k
    prob_superinfection = 0.0
    
    Nt = np.zeros(Nb, dtype = np.int8)
    time = 0.0
    superinfection = False
    
    death_times = []

    
    while any(Nt <= n):
        if ti_superinfection<time<(ti_superinfection+len_superinfection) and superinfection == False:
            superinfection=True

            
        if ti_superinfection + len_superinfection < time and superinfection==True:
            superinfection=False
            prob_superinfection = 0.0
            
            
        
        if superinfection==True:
            prob_superinfection = deltat*nu*P0*np.exp(-0.1*((time-ti_superinfection)+(deltat/2.0)))

        
        for b in range(0,Nb):
            r = np.random.uniform()
        
            sum_probs = prob_new_timer + prob_superinfection
            if sum_probs>1.0:
                print('Probabilities are higher than one')
                #print(sum_probs)
        
            if r < prob_new_timer:
                Nt[b] = Nt[b] + 1.0
            elif prob_new_timer < r < sum_probs:
                if Nt[b] < nd:
                    Nt[b] = 0 
                elif Nt[b] >= nd:
                    Nt[b] = Nt[b] - float(nd)
        
        if any(Nt>=n):
            index = np.where(Nt>=n)[0]
            Nt = np.delete(Nt, index)
            Nb = Nb - len(index)
            for i in range(0,len(index)):
                death_times.append(time)
       
        
        time = time + delta_t
    
    
    death_times = np.array(death_times)  

                      
    return death_times
            
        
 
delta_t = 0.01
ndata = 20000
mosi = 0.4
ti= 15.0
length = 3.0
k = [2.0,4.0,6.0,8.0]
N = 60
nd =0
Nb = 10000



for i in range(0,len(k)):
    
    ff = ((float(N)/float(k[i])))/23.0
    t =  Nb_bacteria_lysis_2(deltat = delta_t, Nb = Nb, n = N, nd=nd, k = k[i], MOSI=mosi, ti_superinfection=ff*ti,len_superinfection=ff*length )


    output_file= str(delta_t)+'_'+str(Nb)+'_'+str(N)+'_'+str(int(k[i]))+'_'+ str(nd)+'_'+str(mosi)+'_'+str(int(ti))+'_'+str(int(length))
    file_path = os.path.join(os.getcwd(),os.path.join('Figure1ForOralPresentation','Prova3'),output_file)
    print(output_file)

    f = open(file_path, "w+")

    for i in range(0, len(t)):
    
        f.write(str((1/ff)*t[i])+'\n')

    f.close()



def integrate_exponential(nu,B,initial_time,end_time):
    a = np.exp(-nu*B*initial_time)
    b = np.exp(-nu*B*end_time)
    return (a-b)*(1.0/(nu*B))

def get_initial_number_of_phage(mosi,nu,B,ti_superinfection,len_superinfection):
    P_0 = mosi*B/integrate_exponential(nu,B,ti_superinfection,ti_superinfection + len_superinfection)
    return P_0






"""________________________ WITH BURSTING OF BACTERIA _______________________"""


def Nb_bacteria_lysis_3(deltat, Nb, N = 60,nd = 33,k = 6.0,MOSI=7.6, ti_superinfection=15.0,len_superinfection=3.0,burst_size=10):
    
    nu = 5.0e-9
    B = 1.0e6
    
    # Volume of the system
    V = Nb/B
    
    # we calculate the initial concentration of phages according to the desired MOSI
    P0 = get_initial_number_of_phage(MOSI,nu,B,ti_superinfection,len_superinfection)
    
    # We calculate the corresponding number of phages
    Np = int(P0*V)
    
    
    prob_new_timer = deltat*k
    prob_superinfection = 0.0
    
    Nt = np.zeros(Nb, dtype = np.int8)
    time = 0.0
    superinfection = False
    
    death_times = []
    
    while any(Nt <= N):
        if ti_superinfection<time<(ti_superinfection+len_superinfection) and superinfection == False:
            superinfection=True

          
        if ti_superinfection + len_superinfection < time and superinfection==True:
            superinfection=False
            prob_superinfection = 0.0
           
            
        
        if superinfection==True:
            prob_superinfection = deltat*nu*(float(Np)/V)

        
        for b in range(0,Nb):
            r = np.random.uniform()
        
            sum_probs = prob_new_timer + prob_superinfection
            if sum_probs>1.0:
                print('Probabilities are higher than one')
                #print(sum_probs)
        
            if r < prob_new_timer:
                Nt[b] = Nt[b] + 1.0
            elif prob_new_timer < r < sum_probs:
                Np = Np - 1
                if Nt[b] < nd:
                    Nt[b] = 0 
                elif Nt[b] >= nd:
                    Nt[b] = Nt[b] - float(nd)
        
        if any(Nt>=N):
            index = np.where(Nt>=N)[0]
            Nt = np.delete(Nt, index)
            Nb = Nb - len(index)
            for i in range(0,len(index)):
                death_times.append(time)
                Np = Np + burst_size
        
        
        time = time + deltat
    
    
    death_times = np.array(death_times)   
                      
    return death_times
"""

delta_t = 0.01
ndata = 20000
mosi = 7.6
ti= 15.0
length =3.0
k = 6.0
n = 60
nd = 33 
Nb = 10000
burst = 10

ff = ((float(n)/float(k)))/23.0



t =  Nb_bacteria_lysis_3(deltat = delta_t, Nb = Nb, N = n, nd=nd, k = k, MOSI=mosi, ti_superinfection=ff*ti,len_superinfection=ff*length, burst_size=burst )


output_file= str(delta_t)+'_'+str(Nb)+'_'+str(n)+'_'+str(int(k))+'_'+ str(nd)+'_'+str(mosi)+'_'+str(int(ti))+'_'+str(int(length))+'_'+str(burst)
file_path = os.path.join(os.getcwd(),'Nu_exploring_CollectiveTimeInterval',output_file)
print(output_file)

f = open(file_path, "a")

for i in range(0, len(t)):
    
    f.write(str((1/(ff)*t[i]))+'\n')

f.close()

"""
    