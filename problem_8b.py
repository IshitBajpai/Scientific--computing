from re import L
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
                    matrix[i][j] =  1
    return matrix

def solve(A,b):
    n = len(A)
    x = [0.]*n

    l = [0.]*n
    s = [0.]*n

    for i in range(1,n+1):
        l[i-1] = i-1;smax = 0.
        for j in range(1,n+1):
            smax = max(smax,abs(A[i-1][j-1]))
        s[i-1] = smax

    for k in range(1,n):
        rmax = 0.
        for i in range(k,n+1):
            r = abs(A[l[i-1]][k-1]/s[l[i-1]])
            if(r > rmax):
                rmax =r
                j=i
        l_temp = l[k-1]
        l[k-1] = l[j-1]
        l[j-1] = l_temp
    
        for i in range(k+1,n+1):
            a_mult = A[l[i-1]][k-1] / A[l[k-1]][k-1]
            A[l[i-1]][k-1] = a_mult

            for j in range(k+1,n+1):
                A[l[i-1]][j-1] = A[l[i-1]][j-1] - (a_mult*A[l[k-1]][j-1])

    for k in range(1,n):
        for i in range(k+1,n+1):
            b[l[i-1]] = b[l[i-1]] - (A[l[i-1]][k-1] * b[l[k-1]])

    x[n-1] = b[l[n-1]]/A[l[n-1]][n-1]

    for i in range(n-1,0,-1):
        sum = b[l[i-1]]
        for j in range(i+1,n+1):
            sum = sum - (A[l[i-1]][j-1] * x[j-1])
        x[i-1] = sum/A[l[i-1]][i-1]
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
            print("Solving Random Matrix")
        if(j == 2):
            print("Solving Hilbert's Matrix")
        if(j == 3):
            print("Solving Part 3")


        try:
            x_solve = solve(A,b)
        except:
            print("Error...Solving Matrix again")
            continue

        
        print("Condition number: (%1.16f)"%(np.linalg.cond(A_copy,2)))
        print("Error from partially-pivoted solve is : %1.16f"%(np.linalg.norm(x-x_solve,2)/np.linalg.norm(x,2)))
        print("Residual from partially-pivoted solve is : %1.16f"%(np.linalg.norm(np.matmul(A_copy,x_solve) - b,2)/np.linalg.norm(b,2)))
        print("Error from np.linalg.solve : %1.16f"%(np.linalg.norm(x-np.linalg.solve(A_copy,b),2)/np.linalg.norm(x,2)))
        print("Residual from np.linalg.solve : %1.16f"%(np.linalg.norm(np.matmul(A_copy,np.linalg.solve(A_copy,b))-b,2)/np.linalg.norm(b,2)))
        print()
    print("------------------------------------------------------------------------------------------------------------------")






# A = np.array([[1,2,-1,1],[-1,1,2,-1],[2,-1,2,2],[1,1,-1,2]],dtype=float) # n*n
# b = np.array([6,3,14,8],dtype=float) # n*1
# print(solve(A,b))