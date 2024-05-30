# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 00:08:55 2020

@author: yonatan eyob
"""

from pylab import *
import numpy as np
#this was used in the class lecture notes and in the textbook. this can also be found in
#http://www-personal.umich.edu/~mejn/cp/programs/gaussxw.py
#i reference this throughout my code for questions 1 and 2
def gaussxw(N):
    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))
    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))
    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    return x, w #outputs the sample points (x) and weights (w) needed to use guassian quadrature
def gaussxwab(N,a,b):
    x, w = gaussxw(N)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w #this allows us to use gaussian quadrature to integrate from a to b

