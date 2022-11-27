import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla
import math

def checkConvergence(A,A_prev,tol):
    if((np.linalg.norm(A - A_prev)/np.linalg.norm(A_prev)) < tol):
        return True
    return False      

def shifted_inverse(A,x0,mu):
    x_k = x0
    i = 0
    tol = 1e-9
    I = np.identity(len(A))
    convergence = True
    while(convergence):
        y_k = npla.solve((A - mu*I),x_k)
        x_prev = x_k
        x_k = y_k / np.linalg.norm(y_k,math.inf)
        if(checkConvergence(x_k,x_prev,tol)):
            break
        i+=1

    eigen_value = 1/np.linalg.norm(y_k,math.inf)
    print("Eigen Value = ",mu + eigen_value)
    print("Eigen Vector(Normalized) = ",y_k/np.linalg.norm(y_k),'Eigen Vector(Without Normalization) = ',y_k)
        
def ComputeAll(A):
    eigen_values,eigen_vectors = npla.eig(A)
    for i in range(len(eigen_values)):
        print("eigen value = ",eigen_values[i]," eigen_vector = ", eigen_vectors[:,i])


A =  np.array( [[6,2,1],[2,3,1],[1,1,1]] )
x0 = np.array([0,0,1])

mu = 2.

print('------------------------')
print('Shifted Inverse Iteration')
shifted_inverse(A,x0,mu)
print('------------------------')
print('Using Libraries')
ComputeAll(A)
print('------------------------')