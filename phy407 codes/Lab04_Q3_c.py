# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 03:11:45 2020

@author: yonatan eyob
"""
# in this code we will we will compare the relaxation method, newton method, and bisection
# method with eachother  then we will find the wien displacement constant 
# and find an estimate to the surface temperature of the sun
import matplotlib.pyplot as plt
import numpy as np

#excercise 6.13 parts b,c:
#b)
def F(x):
    return 5*np.exp(-x)+x-5
def derv_F(x):
    return -5*np.exp(-x)+1

x=np.linspace(-1,10,100)
plt.plot(x,F(x),label='F(x)=$5e^{-x}+x-5$')
plt.plot(x,x*0,label='F(x)=0')
plt.title('F(x)=$5e^{-x}+x-5$')
plt.ylabel('F(x)')
plt.xlabel('x')
plt.legend(loc='upper right')
#bi is the bisection method
def bi(i,f,delta):
    x1=i
    x2=f
    m=0
    while abs(x2-x1)>= delta:
        m+=1
        x=(x2+x1)/2
        p= F(x2)*F(x1)
        if p>delta:
            x1=x
        else:
            if p<delta:
                x2=x
    return 'iterations={a}, root: x={b}'.format(a=m,b=x)
#Newt is Newton's method
def Newt(x,delta):
    xn=x
    m=0
    while abs(F(xn))>delta:
        m+=1
        xn=x-F(x)/derv_F(x)
        x=xn
    return 'iterations={a}, root: x={b}'.format(a=m,b=xn)       
#relax is the relaxation method
def relax(x,delta):
    for i in range(10000): 
        if abs(x-(5-5*np.exp(-x)))>delta:   
            x=5-5*np.exp(-x)
        else:
            return 'iterations={a}, root: x={b}'.format(a=i,b=x)
    
#these are the inputs needed to run the functions bi,relax,Newt.      
initial_guess=4 #in bi() this is the lowest value in the initial (x1,x2) range. in Newt() and relax() this is the initial guess.
end_of_range=6 #in bi() this is the largest value in the initial (x1,x2) range.
delta=1e-25     #this value substitutes zero in our functions because python will run into errors if we actually run our function with zeros
#end of inputs. note that we change these input values to get all the values listed in the pdf for 3c.

print('bisection method: {a} \nrelaxation method: {b}\n\
Newtons method: {c}'.format(a=bi(initial_guess,end_of_range,delta),b=relax(initial_guess,delta),c=Newt(initial_guess,delta)))

#now that we know the non-obvious root is about x=4.965114231744276,
# lets solve for the Wien displacement constant using x=4.965114231744276.
x=4.965114              
h=6.62607004*10**(-34)  #Planck's constant (SI units)
kB=1.38064852*10**(-23) #Boltzmann constant (SI units)
c=299792458             #speed of light (SI units)
b=h*c/(kB*x)            #Wien displacement constant formula
print('\nthe Wien displacement constant is {0}meters*Kelvin'.format(b))

#c)
wavelength=502*10**(-9) #peak wavelength in the sun's emitted radiation
T=b/wavelength          #Temperature formula
print('\nby our calculations, the surface temperature of the sun is estimated to be {0} kelvin'.format(T))