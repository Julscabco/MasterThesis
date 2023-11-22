import matplotlib.pyplot as plt
import numpy as np

# Shape has to be (number of cells, number of iterations in two minutes)
matrix = np.zeros((4,2),dtype=np.int64)
print(matrix)

matrix = matrix[:,1:]
print(matrix)

new_arra = np.zeros(4,dtype=np.int64)
matrix = np.column_stack((matrix,new_arra))
print(matrix)
