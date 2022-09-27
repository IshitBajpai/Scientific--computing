import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla


partial_sum = 0
for k in range(1,5001):
    partial_sum += (1/k) 
    if(k%100 == 0):
        eulerConstant = partial_sum - np.log(k)
        print("Euler's Constant at ",k," iteration : ",eulerConstant)