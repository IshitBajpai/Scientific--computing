from turtle import forward
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla

# condition number 

input  =  []
output =  []
condition_number = np.zeros(21,dtype=float)
forward_error = np.zeros(21,dtype=float)
backward_error = np.zeros(21,dtype=float)
absolute_input_error = np.zeros(21,dtype=float)

for i in range(0,21):
    x = ( np.pi / 4.) + (2.*np.pi * (10.**i) )

    # condition_number[i] = abs(1. - np.tan(x))/abs((np.pi / 4.) - np.arctan(np.tan(np.pi/4.)))
    condition_number[i] = abs(x/(np.sin(x)*np.cos(x)))
    forward_error[i] = abs(1 - np.tan(x))
    backward_error[i] = abs(forward_error[i]/condition_number[i])

    absolute_input_error[i] = abs(backward_error[i] * x)
    print("J = ",i)
    print('(x, tan(x)) = (%1.16f, %1.16f)'%(x, np.tan(x)))
    print("Condition number",condition_number[i])
    print("Relative Forward Error",forward_error[i])
    print("Relative Backward Error",backward_error[i])
    print("Absolute Input Error",absolute_input_error[i])


idx = np.argmax(absolute_input_error)
print("---------------------------------------------")
print("Maximum Input Error:",absolute_input_error[idx])

