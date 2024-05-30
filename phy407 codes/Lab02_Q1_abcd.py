# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 20:56:43 2020

@author: yonatan eyob
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as spc

#this code will be for question 1abcd. i will be testing the forward and central differences and compare
#their errors.
figure,(a)=plt.subplots(nrows=1,ncols=1)
figure,(b)=plt.subplots(nrows=1,ncols=1)

#1a
h=10.**(np.arange(-16,1))   #an array of step sizes.
x=0.5                       #x value we are considering.
def e(x,h):                    #being lazy
    return np.exp(-(x+h)**2)
def forward_dervivative(x,h):  #this is the forward derivative function 
    return(e(x,h)-e(x,0))/h

#1b
#to get the analytical value i took the derivative of f(x) (by hand) and plugged in x=0.5    
analytical_value=-2*x*(np.exp(-(x)**2))
forward_error=abs(forward_dervivative(x,h)-analytical_value) #absolute value of the error in forward derivative calculations
print("this is the forward derivative of f(x)=e^{-x^{2}}:",forward_dervivative(x,h))
print("this is the forward derivative error:",forward_error)

#1c
a.loglog(h,forward_error)   #this plot is replicated on question 1d. this is just for a visual for 1c
a.loglog(h,forward_error,".-", label='forward error')
a.legend(loc="lower left")         
a.set_ylabel("absolute value of error")
a.set_xlabel("step size")
a.set_title("the errors of forward derivatives")

#1d
def central_dervivative(x,h):   #this central difference function 
    return(e(x,h)-e(x,-h))/((x+h)-(x-h))
central_error=abs(central_dervivative(x,h)-analytical_value)
b.loglog(h,central_error,".-", label='central error')
b.loglog(h,forward_error,".-", label='forward error')
b.legend(loc="lower left")         
b.set_ylabel("absolute value of error")
b.set_xlabel("step size")
b.set_title("comparing the errors of central and forward derivatives")


