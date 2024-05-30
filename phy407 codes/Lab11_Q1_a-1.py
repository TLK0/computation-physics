# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 03:55:45 2020

@author: yonat
"""

from math import sqrt,exp
from numpy import empty
from random import random,randrange,seed
import matplotlib.pyplot as plt

seed(10)
N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e5

# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0]**2+x[1]**2)

# Function to calculate the total length of the tour
def distance():
    s = 0.0
    for i in range(N):
        s += mag(r[i+1]-r[i])
    return s

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()
plt.plot(r[:,0],r[:,1],'.',markersize=10)

g=100
seed(g)  # change this to get new paths

# Main loop
t = 0
T = Tmax
while T>Tmin:

    # Cooling
    t += 1
    T = Tmax*exp(-t/tau)

    # Choose two cities to swap and make sure they are distinct
    i,j = randrange(1,N),randrange(1,N)
    while i==j:
        i,j = randrange(1,N),randrange(1,N)

    # Swap them and calculate the change in distance
    oldD = D
    r[i,0],r[j,0] = r[j,0],r[i,0]
    r[i,1],r[j,1] = r[j,1],r[i,1]
    D = distance()
    deltaD = D - oldD

    # If the move is rejected, swap them back again
    if random()>exp(-deltaD/T):
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = oldD
print('with seed({}) D='.format(g),D)
plt.plot(r[:,0],r[:,1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('simulated annealing with seed({0}) and seed({1})'.format(10,g))




















