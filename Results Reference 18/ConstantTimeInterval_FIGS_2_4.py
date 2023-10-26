# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 16:48:24 2023

@author: Usuario
"""
import os

import numpy as np

""" REPRODUCTION OF RESULTS IN REFERENCE 18. CONSTANT TIME INTERVAL (NO GILLESPIE) """

# Main function:
# Simulates the evolution: for each time interval calculates the probabilities of events
# Chooses which event will happen
# If there is a superinfection we have one less phage and nd less timer molecules
# If there is the birht of a timer molecule then we have another timer molecule
# All this starts to happen after ti_superinfection and stops after len_superinfection

# INPUTS: N, nd, k, mosi, ti_superinfection, len_superinfection,
# OUTPUT: time at which the cell dies

def integrate_exponential(nu,B,initial_time,end_time):
    a = np.exp(-nu*B*initial_time)
    b = np.exp(-nu*B*end_time)
    return (a-b)*(1.0/(nu*B))

def get_initial_number_of_phage(mosi,nu,B,ti_superinfection,len_superinfection):
    P_0 = mosi*B/integrate_exponential(nu,B,ti_superinfection,ti_superinfection + len_superinfection)
    return P_0


    


def one_bacteria_lysis_2(deltat, N = 60,nd = 33,k = 6.0,MOSI=7.6, ti_superinfection=15.0,len_superinfection=3.0):
    
    nu = 5e-9
    B = 2.0e7
    
    # we calculate the initial concentration of phages according to the desired MOSI
    P0 = get_initial_number_of_phage(mosi,nu,B,ti_superinfection,len_superinfection)
    
    prob_new_timer = deltat*k
    prob_superinfection = 0.0
    
    Nt = 0.0
    time = 0.0
    superinfection = False
    
    while Nt <= N:
        if ti_superinfection<time<(ti_superinfection+len_superinfection) and superinfection == False:
            superinfection=True

            
        if ti_superinfection + len_superinfection < time and superinfection==True:
            superinfection=False
            prob_superinfection = 0.0
            
            
        
        if superinfection==True:
            prob_superinfection = deltat*nu*P0*np.exp(-0.1*((time-ti_superinfection)+(deltat/2.0)))

        
        r = np.random.uniform()
        sum_probs = prob_new_timer + prob_superinfection

        if sum_probs>1.0:
            print('Probabilities are higher than one')
        #print(sum_probs)
        if r < prob_new_timer:
            Nt = Nt + 1.0
        elif prob_new_timer < r < sum_probs:
            if Nt < float(nd):
                Nt = 0.0
            elif Nt >= float(nd):
                Nt = Nt - float(nd)
        time = time + delta_t
                      
    return time
            
        
def simulation(ndata,deltat,n,nd,kk,mosi,ti_superinfection,len_superinfection,ff,folder):
    name_of_file = str(delta_t)+'_'+str(n)+'_'+str(nd)+'_'+str(kk)+'_'+ str(mosi)+'_'+str(ti_superinfection)+'_'+str(len_superinfection)+'_'+str(ndata)+'.txt'
    file_path=os.path.join(os.getcwd(),'ConstantTimeInterval',folder)
    
    if not os.path.exists(file_path):
        os.mkdir(file_path,0o666)
    
    file_path = os.path.join(file_path,name_of_file)
        
    f = open(file_path, "a")

    for i in range(0, ndata):
        
        tau = (1.0/float(ff))*one_bacteria_lysis_2(deltat,N=n,nd=nd,k=kk,MOSI=mosi,ti_superinfection=ti_superinfection*ff, len_superinfection=len_superinfection*ff)
        f.write(str(tau)+'\n')

    f.close()
    return  

delta_t = 0.01
ndata = 10000
mosi = 0.4
ti_superinfection = 15.0
len_superinfection = 3.0
k = 6.0
n = 60
nd = 26

conversion_factor = (float(n)/float(k))/23.0


output_folder = str(delta_t)+'_'+str(n)+'_'+str(int(k))+'_'+ str(mosi)+'_'+str(int(ti_superinfection))+'_'+str(int(len_superinfection))+'_'+str(ndata)

simulation(ndata,delta_t,n,nd,k,mosi,ti_superinfection,len_superinfection,conversion_factor,output_folder)     
    
    