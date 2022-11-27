import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla
import math
def checkConvergence(A,A_prev,tol):
    # print(np.linalg.norm(A - A_prev)/np.linalg.norm(A_prev))
    if((np.linalg.norm(A - A_prev)/np.linalg.norm(A_prev)) < tol):
        return True
    return False       
def QR_iteration(A):
    n = len(A)
    convergence =   True
    i = 0
    I = np.identity(n)
    tol = 1e-9
    while(convergence):
        sigma = A[n-1][n-1]
        Q,R = spla.qr(A - sigma*I)
        A_prev = A
        A = np.dot(R,Q)+sigma*I
        i+=1
        if(checkConvergence(A,A_prev,tol)):
            break
    eigenvalues = A.diagonal()
    print("Eigenvalues using QR iteration  = ",eigenvalues)
    # print(i)

def ComputeEigenvalues(A):
    eigen_values,eigen_vectors = npla.eig(A)
    print("Eigenvalues using libraries = ",eigen_values)

print("-------------------")
print("Testing on Q6")
A =  np.array( [[2,3,2],[10,3,4],[3,6,1]] )
QR_iteration(A)
ComputeEigenvalues(A)
print("-------------------")
print("Testing on Q7")
A =  np.array( [[6,2,1],[2,3,1],[1,1,1]] )
QR_iteration(A)
ComputeEigenvalues(A)
print("-------------------")


