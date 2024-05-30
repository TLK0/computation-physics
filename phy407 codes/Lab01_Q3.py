# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 17:58:05 2020

@author: yonatan eyob
"""


import numpy as np
import matplotlib.pyplot as plt
import time

# in this code i will be comparing the matrix multiplication operation speed of 
#the np.dot function and python manuelly matrix multiplying by timing and graphing them.

#i got the "figure,(ax1,ax2)=plt.subplots(nrows=2,ncols=1)" idea from: http://jonathansoma.com/lede/data-studio/classes/small-multiples/long-explanation-of-using-plt-subplots-to-create-small-multiples/
#this figure allows me to plot multiple graphs in 1 run without them stacking.
figure,(ax1,ax2)=plt.subplots(nrows=2,ncols=1)
# the dot() function takes input N  and multiplies 2 NxN matrices using np.dot() and outputs the 
#amount of time this operation took
def dot(N):
    start = time.time()
    A = np.ones([N, N], float)*10
    B = np.ones([N, N], float)*4
    C=np.dot(A,B)
    z=time.time() - start
    return z
# the mult() function takes input N  and multiplies 2 NxN matrices by manuelly matrix multiplying 
#the two matrices (columns dot rows, etc...) and outputs the amount of time this operation took
def mult(N):
    start = time.time()
    C = np.zeros([N,N],float)
    A = np.ones([N, N], float)*10
    B = np.ones([N, N], float)*4
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i,j] += A[i,k]*B[k,j]
    z=time.time() - start
    return z

#the Tdot(), Tmult(), and N() functions takes 3 inputs that works similar to np.linspace. a is the 
#initial N in the list of NxN matrices, b is the final N, c is the number of intervals between a and b (very similar to np.linspace)
#N() outputs N=np.linspace(a,b,c), Tdot() outputs an array such that Tdot()[i]=dot(N[i]),
#Tmult() outputs an array such that Tmult()[i]=mult(N[i]) 
def Tdot(a,b,c):
    i=0
    t=np.zeros(c)
    N=np.linspace(a,b,len(t))
    while i<len(t):
        t[i]=dot(int(N[i]))
        i+=1
    return t
def N(a,b,c):
    N=np.linspace(a,b,c)
    return N
def Tmult(a,b,c):
    i=0
    t=np.zeros(c)
    N=np.linspace(a,b,len(t))
    while i<len(t):
        t[i]=mult(int(N[i]))
        i+=1
    return t

#Ni,Nf,t are the inputs
Ni,Nf,t=0,100,30 #Ni=first N, Nf=last N, t=number of intervals
#end of inputs

#plotting N() vs Tdot()
ax2.plot(N(Ni,Nf,t),Tdot(Ni,Nf,t),label='numpy dot product')
#plotting N() vs Tmult()
ax2.plot(N(Ni,Nf,t),Tmult(Ni,Nf,t),label='python matrix mult.')
#plotting (N())**3 vs Tmult()
ax1.plot(N(Ni,Nf,t)**3,Tmult(Ni,Nf,t),label='python matrix mult. N^3')
#plotting (N())**3 vs Tdot() is unnessesary.
ax2.legend(loc="upper left")
ax2.set_ylabel("time (s)")
ax2.set_xlabel("Number of rows in NxN matrix")

ax1.set_title("matrix multiplication operation vs time")
ax1.legend(loc="upper left")
ax1.set_ylabel("time (s)")
ax1.set_xlabel("Number of rows in NxN matrix cubed")

#bonus (attempting to make each time measurement more accurate)
#i noticed that the times were slightly different everytime i ran the timer function so i made
#these two improved functions (Idot and Imult) that does the same thing as Tdot and Tmult but it takes
# each time measurement v amount of times (v is one of the input values) and takes the average value
# from them for each element element in the array it outputs. so Idot and Imult is the same idea as
#Tdot and Tmult but the difference is that each element in the output array of Idot/Imult is an average 
#of v measurements and Tdot/Tmult is just 1 measurement per element. I was suggested not to hand this in
# (because it is slow) but was told to reference that i attempted to fix the random timing varience.
def Idot(a,b,c,v):
    i=0
    j=0
    add=np.zeros(v)
    t=np.zeros(c)
    N=np.linspace(a,b,len(t))
    while i<len(t):
        while j<v:
            add[j]=dot(int(N[i]))
            j+=1
        avg=np.sum(add)/v    
        t[i]=avg
        i+=1
        j=0
    return t
def Imult(a,b,c,v):
    i=0
    j=0
    add=np.zeros(v)
    t=np.zeros(c)
    N=np.linspace(a,b,len(t))
    while i<len(t):
        while j<v:
            add[j]=mult(int(N[i]))
            j+=1
        avg=np.sum(add)/v    
        t[i]=avg
        i+=1
        j=0
    return t
    
#sample plot
#v=20
#plt.plot(N(Ni,Nf,t),Imult(Ni,Nf,t,v),label='python matrix mult. avg')
#plt.legend(loc="upper left")
#plt.ylabel("time (s)")
#plt.xlabel("Number of rows in NxN matrix")
#plt.title("matrix multiplication operation vs time")




