# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 18:57:05 2020

@author: Yonatan Eyob
"""
# In this code we will use simulated annealing to find the global minimum of 
# two functions and plot x vs time and y vs time of the simulated annealing process.
import numpy as np
from random import random,randrange,seed
import matplotlib.pyplot as plt
# 1bi


Tmax = 100 # max temperature
Tmin = 1e-4 # minimum temperature
tau = 1e4 # tau time constant
x_initial=np.array([2,2]) # initial (x,y) for f(x,y)

def gaussian():
    sigma=1
    r = np.sqrt(-2*sigma*sigma*np.log(1-random()))
    theta = 2*np.pi*random()
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return np.array([x,y])


def f(x):
    return x[0]**2-np.cos(4*np.pi*x[0])+(x[1]-1)**2

fn = f(x_initial)
T = Tmax
x = x_initial
z=[]
time=[]
t = 0

while T>Tmin:
    t+=1 
    T=Tmax*np.exp(-t/tau) # cooling rate
    
    previous_x=x # previous x
    previous_f=fn # previous f(x,y) (before (x,y) is updated)
    r=gaussian() # (dx,dy)
    x = x+r # updating x
    fn=f(x) # updating f(x,y)
    if random()>np.exp((previous_f-fn)/T):
        x=previous_x
        fn=previous_f
    time.append(t)
    z.append(x)

print('The global minimum of this function is approximately ',x)
z=np.array(z)
x,y=z[:,0],z[:,1]

plt.plot(time,x,'.')
plt.title('x vs time for Simulated annealing on \n f(x,y)=$x^{2}-cos(4\pi x)+(y-1)^{2}$')
plt.xlabel('time')
plt.ylabel('x')
plt.show()
plt.plot(time,y,'.')
plt.title('y vs time for Simulated annealing on \n f(x,y)=$x^{2}-cos(4\pi x)+(y-1)^{2}$')
plt.xlabel('time')
plt.ylabel('y')
plt.show()
       



 
# 1bii
   
def f(x):
    return np.cos(x[0])+np.cos(np.sqrt(2)*x[0])+np.cos(np.sqrt(3)*x[0]) + (x[1]-1)**2
fn = f(x_initial)
T = Tmax
x = x_initial
z=[]
time=[]
t = 0

while T>Tmin:
    t+=1
    T=Tmax*np.exp(-t/tau)
    
    previous_x=x
    previous_f=fn
    r=gaussian()
    x = x+r
    fn=f(x)
    if random()>np.exp((previous_f-fn)/T) or x[0]<0 or x[0]>50 or x[1]<-20 or x[1]>20:
        x=previous_x
        fn=previous_f
    time.append(t)
    z.append(x)

print('The global minimum of this function is approximately ',x)
z=np.array(z)
x,y=z[:,0],z[:,1]

plt.plot(time,x,'.')
plt.title('x vs time for Simulated annealing on \n f(x,y)=$cos(x)+cos(\sqrt{2}x)+cos(\sqrt{3}x)+(y-1)^{2}$')
plt.xlabel('time')
plt.ylabel('x')
plt.show()
plt.plot(time,y,'.')
plt.title('y vs time for Simulated annealing on \n f(x,y)=$cos(x)+cos(\sqrt{2}x)+cos(\sqrt{3}x)+(y-1)^{2}$')
plt.xlabel('time')
plt.ylabel('y')
plt.show()  
        
        
        
        
        
        
        
        
        
        
        
        
        
