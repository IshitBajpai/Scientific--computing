import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla
import math

def f(x):
    return (100/x)*(np.sin(10./x))
    # return x

def Quad(n):
    a = 1
    b = 3
    h = ((b-a)/n)
    # xi = np.linspace(a,b,n)

    result = 0.

    for i in range(1,int(n/2)+1):
        lower = a + 2*(i-1)*h
        upper = a + 2*i*h
        
        result+= GaussianQuad(lower,upper)

    true_value = -18.79829683678703
    rel_error = abs(true_value - result)/abs(true_value)
    print("n = ",n," Integeration Estimate = ",result, " Relative Error = ",rel_error)



    
# will find wi and xi and plug it in Q(f) += wi*f(xi)
def GaussianQuad(a,b):
    # changing the variable as done in Q2.
    w1 = 1.
    w2 = 1.
    x1 = -1/math.sqrt(3.)
    x2 = 1/math.sqrt(3.)
    t1 = ((b-a)*x1 + a + b )/(2.) 
    t2 = ((b-a)*x2 + a + b )/(2.) 


    Q_f = ((b-a)/2) * (w1*f(t1)+w2*f(t2))
   
    return Q_f

n_value = [2.,4.,8.,16.,32.,64.]
for i in n_value:
    print("---------------------------------------------------------------------------")
    Quad(i)