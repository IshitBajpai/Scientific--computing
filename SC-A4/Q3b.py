import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla



def cubicSpline(n):
    x = np.linspace(-1, 1, n)

     
    h = []
    for i in range(1,n):
        h.append(x[i] - x[i-1])
    # print("h = ",len(h))
    
    # get a_i's
    a = []
    for i in range(n):
        a.append(f(x[i]))

    # print("a = ",len(a))

    A = np.zeros((n,n),dtype=float)
    A[0][0] = 1. 
    A[n-1][n-1] = 1.
    for i in range(1,n-1):
        A[i][i-1] = h[i-1]
        A[i][i] = 2*(h[i-1]+h[i])
        A[i][i+1] = h[i]

    # print(A.shape)

    y = [0.]
    for i in range(1,n-1):
        y_i  = ((3*(a[i+1]-a[i]))/h[i]) - ((3*(a[i]-a[i-1]))/h[i-1])
        y.append(y_i)
    y.append(0.)


    c = np.linalg.solve(A,y)
    # print("c = ",len(c))

    b = []
    d = []

    for i in range(len(c)-1):
        b_i = ((a[i+1]-a[i])/h[i]) - ((h[i]*(2*c[i]+c[i+1]))/3) 
        d_i = (c[i+1]-c[i])/(3*h[i])

        b.append(b_i)
        d.append(d_i)

    # print("b = ",len(b))
    # print("d = ",len(d))

    for i in range(1,n-1):
        plotSpline(x[i-1],a[i-1],b[i-1],c[i-1],d[i-1],x[i-1],x[i])

    # plt.show()
    plotRunge(n)
    
def plotSpline(x,a,b,c,d,low,high):

    x_i = np.linspace(low,high,5)  # plotting for 5 intermediate points
    y_i = []
    for i in range(len(x_i)):
        f_x = a + b*(x_i[i] - x) + c*((x_i[i]-x)**2) + d*((x_i[i]-x)**3)
        y_i.append(f_x)
    
    plt.plot(x_i,y_i)

def f(t):
    return 1/(1 + 25 * t**2)

def plotRunge(n):
    T = np.linspace(-1, 1, 100)
    f_at_T = np.array([f(t_) for t_ in T])
    plt.plot(T, f_at_T,label='Runge',linestyle='--')
    plt.grid(True)
    plt.xlabel('t', fontsize=12)
    plt.ylabel('f(t)', fontsize=12)
    plt.title('Plot of Runge\'s function for n = '+str(n), fontsize=14)
    plt.legend()
    plt.show()

cubicSpline(11)
cubicSpline(21)