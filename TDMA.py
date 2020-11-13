# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 21:42:57 2020

@author: yujie
"""

'''
python TDMA solver
'''
import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize = (8,5))
def TDMA(N, A, d, phi):
    """
    @N: element number of the solution vector
    @A: the coeff. matrix (2-D vector)
    @d: the vector b in Ax = b
    @phi: solution vector
    """
    Ns = N
    P = np.zeros(Ns)
    Q = np.zeros(Ns)
    a = np.zeros(Ns)
    b = np.zeros(Ns)
    c = np.zeros(Ns)
    
    a[0] = A[0][0]
    b[0] = A[1][0]
    c[0] = 0
    a[Ns-1] = A[Ns-1][Ns-1]
    b[Ns-1] = 0
    c[Ns-1] = A[Ns-2][Ns-1]
    
    for ii in np.arange(1, (Ns-1), 1):
        a[ii] = A[ii][ii]
        b[ii] = A[ii+1][ii]
        c[ii] = A[ii-1][ii]
    
    #compute the value of P1 and Q1:
    P[0] = -(b[0]/a[0])
    Q[0] = d[0]/a[0]
    #forward recursion to compute the value of Pi and Qi:
    for i in np.arange(1, Ns, 1):
        P[i] = -(b[i]/(a[i]+c[i]*P[i-1]))
        Q[i] = (d[i]-c[i]*Q[i-1])/(a[i]+c[i]*P[i-1])
    
    #set P[N] = Q[N]:
    phi[Ns-1] = Q[Ns-1]
    #back recursion to compute the values of phi[i]:
    for i in np.arange((Ns-2), -1, -1):
        phi[i] = P[i]*phi[i+1]+Q[i]
    
    return phi

#construct the 2-D array using numpy:
    ###user define###
startPoint = 0.0
endPoint = 1.0
Np = 11
N = Np - 2
    #################
dx = (endPoint - startPoint)/(Np-1)
A = np.zeros((N, N))
d = np.zeros(N)
phi = np.zeros(N)
"""
The coeff. matrix needs to be defined as:
    [A00, A10, A20, A30, ..., Ai0]
    [A01, A11, A21, A31, ..., Ai1]
    ...          ...          ...
    [A0i, A1i, A2i, A3i, ..., Aii]
which is corresponding to the index of real matrix
There is one example:
"""
"""
A[0][0] = 3
A[0][1] = -2
A[1][0] = -1
A[1][1] = 6
A[1][2] = -2
A[2][1] = -1
A[2][2] = 6
A[2][3] = -2
A[3][2] = -1
A[3][3] = 7

d[0] = 3
d[1] = 4
d[2] = 5
d[3] = -3

x = TDMA(N, A, d, phi)    
print(x)
"""
for ii in np.arange(0, N, 1):
    A[ii][ii] = -(2/(dx*dx)+2)
for i in np.arange(1, N, 1):
    A[i][i-1] = 1/(dx*dx)+3/(2*dx)
for j in np.arange(0, N-1, 1):
    A[j][j+1] = 1/(dx*dx)-3/(2*dx)

d[0] = 2*dx + 3 -2*(1/(dx*dx)-3/(2*dx))
d[-1] = 2*(endPoint - dx) + 3 - (1/(dx*dx)+3/(2*dx))
for jj in np.arange(1, N-1, 1):
    d[jj] = 2*(jj+1)*dx + 3

v = TDMA(N, A, d, phi)
y = np.zeros(Np)
for yy in np.arange(1, Np-1, 1):
    
    y[yy] = v[yy-1]    
y[0] = 2.0
y[-1] = 1

for x in np.arange(0, Np, 1):
    plt.scatter(x/10, y[x], c="r", marker="o")
plt.xlabel("x")
plt.ylabel("y")
  
