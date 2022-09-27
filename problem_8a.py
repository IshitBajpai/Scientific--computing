import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla

def matrix_generator(n,type): # returns a matrix of n*n and type = random,hilbert
    matrix = np.zeros((n,n),dtype = float) # initializing matrx
    if(type == 1): # random
        matrix = np.random.random_sample((n,n))
    elif(type == 2): # Hilbert matrix
        for i in range(n):
            for j in range(n):
                matrix[i][j] = (1./(float(i)+float(j)+1.))
    else:
        for i in range(n):
            for j in range(n):
                if(i>j):
                    matrix[i][j] = -1
                else:
                    matrix[i][j] = 1
    return matrix

def solve(A,b):
    n = len(A)
    x = np.zeros(n,dtype=float)

    for k in range(1,n):
        for i in range(k+1,n+1):
            a_temp = A[i-1][k-1]/A[k-1][k-1]
            A[i-1][k-1] = a_temp
            for j in range(k+1,n+1):
                A[i-1][j-1] = A[i-1][j-1] - ( a_temp * A[k-1][j-1] )
            b[i-1] = b[i-1]-(a_temp*b[k-1])
    x[n-1] = b[n-1] / A[n-1][n-1]
    for i in range(n-1,0,-1):
        sum = b[i-1]
        for j in range(i+1,n+1):
            sum = sum -  (A[i-1][j-1] * x[j-1])
        x[i-1] = sum / A[i-1][i-1]

    return x


size = [10,20,30,40]
for i in range(len(size)):
    x = np.ones(size[i],dtype=float)
    print("n =",size[i])
    for j in range(1,4):
        A = matrix_generator(size[i],j)
        b = np.matmul(A,x)
        
        A_copy = np.copy(A)
        if(j == 1):
            print("Solving Random Matrix\n")
        if(j == 2):
            print("Solving Hilbert's Matrix\n")
        if(j == 3):
            print("Solving Part 3\n")

    
        x_solve = solve(A,b)
        print("Condition number: (%1.16f)"%(np.linalg.cond(A_copy,2)))
        print("Error from un-pivoted solve is : %1.16f"%(np.linalg.norm(x-x_solve,2)/np.linalg.norm(x,2)))
        print("Residual from un-pivoted solve is : %1.16f"%(np.linalg.norm(np.matmul(A_copy,x_solve) - b,2)/np.linalg.norm(b,2)))
        print("Error from np.linalg.solve : %1.16f"%(np.linalg.norm(x-np.linalg.solve(A_copy,b),2)/np.linalg.norm(x,2)))
        print("Residual from np.linalg.solve : %1.16f"%(np.linalg.norm(np.matmul(A_copy,np.linalg.solve(A_copy,b))-b,2)/np.linalg.norm(b,2)))
        print()
    print("------------------------------------------------------------------------------------------------------------------")


# A = np.array([[1,2,-1,1],[-1,1,2,-1],[2,-1,2,2],[1,1,-1,2]],dtype=float) # n*n
# b = np.array([6,3,14,8],dtype=float) # n*1
# print(solve(A,b))
# n = len(A)
# print(x)


