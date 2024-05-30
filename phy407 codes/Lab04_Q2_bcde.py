# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 19:33:50 2020

@author: yonatan eyob
"""

# in this code we will find the 10 eigenvalues of a 10x10 Hamiltonian 
# matrix and then find the first 10 eigenvalues of a 100x100 Hamiltonian 
# matrix. Then we will plot the probability density graphs of the ground 
# state and the first two excited states.

import matplotlib.pyplot as plt
import numpy as np

figure,(a)=plt.subplots(nrows=1,ncols=1)
#2b
def H(m,n):     # this is the elements of the hamiltonian matrix. m is the row index and n is the column index
    L=5*100*10**(-12) #meters
    a= 10*1.6022*10**(-19) #joules
    M=9.1094*10**(-31) #kg
    h=1.0545718*10**(-34)
    if m%2==0:
        if m==n:
            return (1/2)*a+np.pi**2*h**2*m**2/(2*M*L**2)
        elif n%2==0:
            return 0
        else:
           return (-8)*a*m*n/(np.pi**2*((m**2)-(n**2))**2) 
    else:
        if m==n:
            return (1/2)*a+np.pi**2*h**2*m**2/(2*M*L**2)
        elif n%2!=0:
            return 0
        else:
            return (-8)*a*m*n/(np.pi**2*((m**2)-(n**2))**2)

#2c
def H_matrix(m,n): #this is the hamiltonian matrix. m is the number of rows and n is the number of columns
    z= []
    for i in range(m):
        row=[]
        for j in range(n):
            row.append(0)
        z.append(row)
    for i in range(m):
        for j in range(n):
            z[i][j]=H(i+1,j+1)
    return z
# note to self: H(a,b)=H_matrix(n,m)[a-1][b-1] if a<=n and b<=m
m=10        #number of rows
n=10        #number of columns
eigan_value,eigan_vector=np.linalg.eigh(H_matrix(m,n))
print('the eiganvalues of the 10x10 hamiltonian matrix expressed in eV: ',
      eigan_value/(1.6022*10**(-19)))

#2d
m=100   #new number of rows
n=100   #new number of columns
eigan_value,eigan_vector=np.linalg.eigh(H_matrix(m,n))
print('\nthe first ten eiganvalues of the 100x100 hamiltonian matrix expressed in eV: ',
      eigan_value[0:10]/(1.6022*10**(-19)))

#2e
def phi(x,energy_level):        #phi is the non-normalized wavefunction. x is the position and energy_level is the energy level you want to consider.
    L=5           #width of the well in angstroms
    if hasattr(x, "__len__")==True:  #checks if x is an array
        z=np.zeros(len(x))
    else:
        z=0
    for i in range(len(eigan_vector[:,energy_level-1])):
        z+= eigan_vector[:,energy_level-1][i]*np.sin((i+1)*np.pi*x/L)
    return z

def simpson(begin,end,N,energy_level,integrand):    #this function takes inputs begin (integral starts) and end (integral ends) to an integral, number of slices (N), and an array or scalar of x for Intensity(x). this outputs the intensity value at the given x input.
    h=(end-begin)/N            #step size
    s=(integrand(begin,energy_level))**2+(integrand(end,energy_level))**2
    for k in range(1,N,2):          #for odd terms
        s+=4*(integrand(begin+k*h,energy_level))**2
    for k in range(2,N,2):          #for even terms
        s+=2*(integrand(begin+k*h,energy_level))**2
    return abs((h/3)*s)

# these are the inputs we will need for the graph:
L=5      #width of well in angstroms
x=np.linspace(0,L,1000)  #position values
N=1000                  # number of steps in the simpson method integral
energy1,energy2,energy3=1,2,3   #energy levels we are graphing the wave functions for
A1=simpson(0,L,N,energy1,phi)   #normalization constant squared for phy(x,1)
A2=simpson(0,L,N,energy2,phi)   #normalization constant squared for phy(x,2)
A3=simpson(0,L,N,energy3,phi)   #normalization constant squared for phy(x,3)
# end of inputs

a.plot(x,abs(phi(x,1))**2/A1,label='$\parallel \psi_{1}(x)\parallel^{2}$')
a.plot(x,abs(phi(x,2))**2/A2,label='$\parallel \psi_{2}(x)\parallel^{2}$')
a.plot(x,abs(phi(x,3))**2/A3,label='$\parallel \psi_{3}(x)\parallel^{2}$')
a.legend(loc='upper right')
a.set_xlabel('x $(angstroms)$')
a.set_ylabel('probability $\parallel \psi_{n}(x)\parallel^{2}$')
a.set_title('probability density functions for Energy Levels n=1,2,3')


