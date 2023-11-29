# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 11:56:16 2023

@author: Usuario
"""
import numpy as np

def get_collapse_time(time,nbacteria):
    index = np.argmax(nbacteria[10:]==0)
    collapse_time = time[index]
    return collapse_time

def get_init_collapse_time(df,time,avg_infections):
    Nimax = np.array(df['Value'][df['Name in code']=='Nimax'])[0]
    index = np.argmax(avg_infections>float(Nimax))
    tindex = time[index]
    return tindex