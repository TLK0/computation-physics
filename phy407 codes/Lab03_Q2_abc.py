# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 02:40:58 2020

@author: yonatan eyob
"""

# in this code we will plot the wave functions for energy levels 0,1,2,3, and 30, 
# then we will find the total energy, momentum uncertainty and position uncertainty
# at energy levels 0 to 15 to find a simple rule to calculate total energy, and 
# find a relation between momentum uncertainty and position uncertainty. 

import matplotlib.pyplot as plt
import numpy as np
import gaussxw as gsx #i explain what this does on the 'gaussxw.py' file. make sure this is in the same directory as that file (note to self)
import math as ma

figure,(a)=plt.subplots(nrows=1,ncols=1)
figure,(b)=plt.subplots(nrows=1,ncols=1)

#2a
def H(n: int,x):           #Hermite polynomial function. n indicates the energy level x can be ascalar or an array.
    if hasattr(x, "__len__")==True:  #checks if x is an array
        H0=np.ones(len(x))
    else:
        H0=1
    H1=2*x
    H=2*x*H1-2*1*H0
    if n==0:
        return H0
    if n==1:
        return H1
    else:
        for i in range(1,n):
            H= 2*x*H1-2*(i)*H0
            H0=H1
            H1=H
        return H

def wave_function(n: int,x):       # n indicates the energy level, x is the position. outputs the wave function for the given energy level
    return np.exp(-(1/2)*x**2)*H(n,x)/np.sqrt((float(2**n)*float(ma.factorial(n))*np.sqrt(np.pi)))

n=np.arange(0,4)                #energy levels for question 2a
x1=np.linspace(-4,4,1000)       # array of position values for energy levels 0-3

a.plot(x1,wave_function(n[0],x1), label='n= 0')
a.plot(x1,wave_function(n[1],x1), label='n= 1')
a.plot(x1,wave_function(n[2],x1), label='n= 2')
a.plot(x1,wave_function(n[3],x1), label='n= 3')
a.legend(loc="upper left") 
a.set_ylabel('$\Psi$(x)')
a.set_xlabel('x $(position)$')
a.set_title('Wave Functions for Energy levels n = 0,1,2,3')

#2b
n=30            #energy level= 30
x2=np.linspace(-10,10,1000)     #array of position values for energy level 30
b.plot(x2,wave_function(n,x2), label='n = {0}'.format(n))
b.legend(loc="lower left") 
b.set_ylabel('$\Psi$(x)')
b.set_xlabel('x $(position)$')
b.set_title('Wave Function for Energy level n = 30')

#2c
def position_uncertainty_integrand(n: int,z):    #integrand of the position uncertainty
    return (1+z**2)/((1-z**2)**2)*((z/(1-z**2))**2*abs(wave_function(n,z/(1-z**2)))**2)

def momentum_uncertainty_integrand(n:int,z):    #integrand of the momentum uncertainty
    return (1+z**2)/((1-z**2)**2)*abs(np.exp(-(1/2)*(z/(1-z**2))**2)*(-(z/(1-z**2))*H(n,z/(1-z**2))+2*n*H(n-1,z/(1-z**2)))/np.sqrt((float(2**n)*float(ma.factorial(n))*np.sqrt(np.pi))))**2

def gaussian_quadrature(N,n,f):       # takes the guassian quadrature of function f. N is the number of steps, n is the energy level, f is the integrand we wish to integrate.
    x,w=gsx.gaussxw(N)
    I=0
    for k in range(N):
        I += w[k]*f(n,x[k])
    return I

def energy(N,n):    #energy at energy level n. N is for the number of steps in the gaussian quadrature.
    return (1/2)*(gaussian_quadrature(N,n,position_uncertainty_integrand)+gaussian_quadrature(N,n,momentum_uncertainty_integrand))

N=100       #number of steps in guassian quadrature
                                    

for n in range(16):   #n is the energy levels
    print("For N={a} slices and energy level= {b},\nEnergy= {c}\nposition uncertainty(n={b}) ={d}\n\
momentum uncertainty(n={b}) = {e}\n".format(
a=N,
b=n,
c=energy(N,n),
d=np.sqrt(gaussian_quadrature(N,n,position_uncertainty_integrand)),
e=np.sqrt(gaussian_quadrature(N,n,momentum_uncertainty_integrand))
))