import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla



def polynomial_interpolation(n):

    t_i =  np.linspace(-1, 1, n)
    Vn = np.vander(t_i,increasing=True)
    y_i = [f(t_) for t_ in t_i]
    coefficients = np.linalg.solve(Vn,y_i)
    # print(coefficients)


    T = np.linspace(-1, 1, 100)
    p_at_T = []
    for i in range(len(T)):
        x = T[i]
        p_at_T.append(polynomial_function(coefficients,x))

    plt.plot(T, p_at_T,label="Polynomaial interpolated "+" n = "+str(n))
    plt.legend()
    plt.grid(True)
    plt.xlabel('t', fontsize=12)
    plt.ylabel('f(t)', fontsize=12)
    plt.title('Plot of Polynomial Interpolation function', fontsize=14)

    plotRunge(n)
    # plt.show()


def polynomial_function(coeff,x):
    polynomial = 0
    for i in range(len(coeff)):
        polynomial+=coeff[i]*(x**i)
    return polynomial



def f(t):
    return 1/(1 + 25 * t**2)

def plotRunge(n):
    T = np.linspace(-1, 1, 100)
    f_at_T = np.array([f(t_) for t_ in T])
    plt.plot(T, f_at_T,label='Runge')
    plt.grid(True)
    plt.xlabel('t', fontsize=12)
    plt.ylabel('f(t)', fontsize=12)
    plt.title('Plot of Runge\'s function for n = '+str(n), fontsize=14)
    plt.legend()
    plt.show()

polynomial_interpolation(11)
polynomial_interpolation(21)
