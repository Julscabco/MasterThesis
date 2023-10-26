# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:37:23 2023

@author: Usuario
"""

import numpy as np

MOSI = 7.6
len_superinfection = 3.0
phages_superinfecting = int(((MOSI*2.0e7)/(np.exp(-0.1*len_superinfection))/2.0e7))
print(phages_superinfecting)

n= 30
print(np.linspace(0,n,10,dtype=np.int8))

Nt = np.array([1,2,3,4,5,6,7,8,9])

print(np.delete(Nt,np.where(Nt==9)))


    

    