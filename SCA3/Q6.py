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

# for largest magnitude eigenvalue
def normal_power_iteration(A,x0):
    
    x_k = x0
    i = 0
    tol = 1e-9
    convergence = True
    while(convergence):
        y_k = np.matmul(A,x_k)
        x_prev = x_k
        x_k = y_k / np.linalg.norm(y_k,math.inf)
        if(checkConvergence(x_k,x_prev,tol)):
            break
        i+=1

    print("Eigen Value = ",np.linalg.norm(y_k,math.inf))
    print("Eigen Vector(Normalized) = ",y_k/np.linalg.norm(y_k),'Eigen Vector(Without Normalization) = ',y_k)

def inverse_iteration(A,x0):
    x_k = x0
    i = 0
    while(i<1000):
        y_k = npla.solve(A,x_k)
        x_k = y_k / np.linalg.norm(y_k,math.inf)
        i+=1
    # print(y_k)
    eigen_value = 1/np.linalg.norm(y_k,math.inf)
    print("Eigen Value = ",eigen_value)
    print("Eigen Vector(Normalized) = ",y_k/np.linalg.norm(y_k),'Eigen Vector(Without Normalization) = ',y_k)
    
def ComputeAll(A):
    eigen_values,eigen_vectors = npla.eig(A)
    for i in range(len(eigen_values)):
        print("eigen value = ",eigen_values[i]," eigen_vector = ", eigen_vectors[:,i])


A =  np.array( [[2,3,2],[10,3,4],[3,6,1]] )
x0 = np.array([0,0,1])
print('------------------------')

print('Normal Power Iteration')
normal_power_iteration(A,x0)
print('------------------------')
print('Inverse Iteration')
inverse_iteration(A,x0)
print('------------------------')
print('Using Libraries')
ComputeAll(A)
print('------------------------')