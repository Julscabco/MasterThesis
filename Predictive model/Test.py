import matplotlib.pyplot as plt
import numpy as np

vector = np.zeros((1,2))
print(vector)

vector = np.append(vector,[[1.0,1.0]],axis = 0 )
print(vector)

array = np.array([1,2,3,3,3,3])
print(np.where(array==3))