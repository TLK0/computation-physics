# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 17:29:39 2020

@author: Yonatan Eyob
"""
# In this code i will solve for the volume of a hypersphere with radius=1 using 
# Monte Carlo integration.

import numpy as np
import matplotlib.pyplot as plt

def f(q,w,e,r,t,y,u,i,o,p):
    if q**2+w**2+e**2+r**2+t**2+y**2+u**2+i**2+o**2+p**2<=1:
        return 1
    else:
        return 0

N = 1000000 #number of random points 
a,b = -1,1 #the min,max of one of the sides of the 10d hypercube. the sidelength is (b-a)

k = 0 # this will contain the points that land in the hypersphere
k2 = 0 # this will be used for variance
for i in range(N):
    q,w,e,r,t,y,u,i,o,p = (b-a)*np.random.random(10)+a*np.ones(10)
    k += f(q,w,e,r,t,y,u,i,o,p)
    k2 += f(q,w,e,r,t,y,u,i,o,p)**2
    
dim=10 #dimensions
I = k * (b-a)**dim / N # solution to the integral
print('solution=',I)
# below is the error calculation:
var = k2/N - (k/N)**2 # variance <f**2> - <f>**2
sigma = (b-a)**dim*np.sqrt(var/N) # sigma=error
print('error = ', sigma)




