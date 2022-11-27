import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla
import seaborn as sns
import math

def checkConvergence(A,A_prev,tol):
    # print(np.linalg.norm(A - A_prev)/np.linalg.norm(A_prev))
    if((np.linalg.norm(A - A_prev)/np.linalg.norm(A_prev)) < tol):
        return True
    return False       

def RayleighQuotientIteration(A,x):
    n = len(A)
    convergence =   True
    x_k = x
    i = 0
    I = np.identity(n)
    tol = 1e-9
    convergence = True
    while(convergence):
        sigma = (np.dot(x_k.T,np.dot(A,x_k)))/(np.dot(x_k.T,x_k))
        y_k = npla.solve(A-sigma*I,x_k)
        x_prev = x_k
        x_k = y_k/np.linalg.norm(y_k,math.inf)
        if(checkConvergence(x_k,x_prev,tol)):
            break
        i+=1
    eigen_value = 1/np.linalg.norm(y_k,math.inf)
    print("Eigen Value = ",eigen_value+sigma)
    print("Eigen Vector(Normalized) = ",y_k/np.linalg.norm(y_k),'Eigen Vector(Without Normalization) = ',y_k)

def ComputeAll(A):
    eigen_values,eigen_vectors = npla.eig(A)
    for i in range(len(eigen_values)):
        print("eigen value = ",eigen_values[i]," eigen_vector = ", eigen_vectors[:,i])

def convergence_rate(A,x,lamda):
    n = len(A)
    error_k = 1.
    x_k = x
    i = 0
    I = np.identity(n)
    error = []
    iteration = []
    roc = []
    while(error_k!=0):
        rayleighQuotient =  (np.dot(x_k.T,np.dot(A,x_k)))/(np.dot(x_k.T,x_k))
        # print(rayleighQuotient)
        y_k = npla.solve(A-rayleighQuotient*I,x_k)
        x_k = y_k/np.linalg.norm(y_k,math.inf)   

        error_k = abs(lamda - rayleighQuotient)     
        
        error.append(error_k)
        iteration.append(i)

        if(i>1 and error_k!=0):
            # print(len(error))
            rate_of_convergence = np.log(error[-1]/error[-2])/np.log(error[-2]/error[-3])
            roc.append(rate_of_convergence)
        i+=1
        # print(error)
    print('------------------------------------------------')
    for j in range(len(roc)):
        print("Iteration :",j+1," Convergence rate:",roc[j])
    print("--------------------------------------------------")
    # print(roc)
    return error,iteration

A =  np.array( [[2,3,2],[10,3,4],[3,6,1]] )
x0 = np.array([1,3,2])

print("----------------------")
print("Rayeleigh Quotient Iteration")

RayleighQuotientIteration(A,x0)
print("----------------------")
print('Using Libraries')
ComputeAll(A)
print("----------------------")
largest_eigenvalue = 11 # using libraries
y,x = convergence_rate(A,x0,largest_eigenvalue)
# print(len(x))
plt.scatter(x=x,y=y,marker='X')
plt.xlabel('Iterations')
plt.ylabel('Error')
plt.show()


