import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import numpy.linalg as npla
import scipy.linalg as spla

def f(x):
    f_x = np.exp(-np.sin(x**3)/4.)
    return f_x

def exact_derivative(x):
    d_x = (-3./4)*(x**2)*(np.cos(x**3))*f(x)
    return d_x

def forward_approx(x,h):
    fapprox = (f(x+h)-f(x))/(h)
    return fapprox

def compute_error(h):
    error = []
    x = 1.
    for i in range(len(h)):
        abs_error = abs(exact_derivative(x) - forward_approx(x,h[i]))
        error.append(abs_error)
    return error

def plot(x,y):
    plt.loglog(x,y)
    plt.xlabel('h')
    plt.ylabel('error')
    plt.show()


h = []
for i in range(1,15+1):
    h.append(10**(-i))
h.reverse()
error = compute_error(h)
# print("errors:" , error)
plot(h,error)
for i in range(len(error)):
    print("----------------------------")
    print("h = ",h[i]," error = ",error[i])
    print("----------------------------")

# plt.loglog(h, error)
# plt.show()
