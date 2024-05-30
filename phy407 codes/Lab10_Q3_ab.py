# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 18:48:13 2020

@author: Yonatan Eyob
"""
# In this code i will be solving an integral using Mean value MC and
# Importance sampling MC 100 times each and plotting the solutions on a
# histograph to see which method gives the most consistent solution. I will do
# this twice; one for 3a and again for 3b with one different integrand per 
# question. There should be four histographs total.
import numpy as np
import matplotlib.pyplot as plt

# 3a)
# Mean value MC:
def f(x):
    return x**(-1/2)/(1+np.exp(x))

N = 10000 # number of random points
a = 0. #beginning of integral
b = 1. # end of integral
k = 0 # will be used to hold the sum of all f(x)
I=[] # this will hold all 100 of the the solutions as a list 

for j in range(100):
    for i in range(N):
        x = np.random.random() #random number between 0 and 1
        k += f(x)
    I.append(k * (b-a) / N)
    k=0 #resetting k
plt.hist(I, 10, range=[0.8, 0.88])
plt.xlabel('x values')
plt.ylabel('number of solutions in bin')
plt.title('Mean Value MC Solutions')
plt.show()




# Importance sampling:
def f(x):
    return x**(-1/2)/(1+np.exp(x))
def w(x):
    return x**(-1/2)


N = 10000 # number of random points
k = 0 # will be used to hold the sum of all f(x)/w(x)

I=[] # this will hold all 100 of the the solutions as a list 
for j in range(100):
    for i in range(N):
        x = (np.random.random())**2 # x(z)=z^2 using the transformation formula
        k += f(x)/w(x)
    int_w=2 # integral of w(x) from x=0 to x=1
    I.append(k * int_w / N)
    k=0 # resetting k 
    
plt.hist(I, 10, range=[0.8, 0.88])
plt.xlabel('x values')
plt.ylabel('number of solutions in bin')
plt.title('Importance sampling MC Solutions')
plt.show()





# 3b)
# Mean value MC:
def f(x):
    return np.exp(-2*abs(x-5))


N = 10000 # number of random points
a = 0. #beginning of integral
b = 10. # end of integral
k = 0 # will be used to hold the sum of all f(x)
I=[] # this will hold all 100 of the the solutions as a list 


for j in range(100):
    for i in range(N):
        x = b*np.random.random() # random number from 0 to 10
        k += f(x)
    I.append(k * (b-a) / N)
    k=0 #resets k

plt.hist(I, 10,range=[0.94,1.06])
plt.xlabel('x values')
plt.ylabel('number of solutions in bin')
plt.title('Mean Value MC Solutions')
plt.show()


# Importance sampling:
def f(x):
    return np.exp(-2*abs(x-5))
def w(x):
    return np.exp(-(x-5)**2/2)/np.sqrt(2*np.pi)



def simpson(begin,end,N,integrand):    #this function takes inputs beginning (integral starts) and end (integral ends) to an integral, number of slices (N), andan integrand function. this outputs the integral of the integrand.
    h=(end-begin)/N            #step size
    s=integrand(begin)+integrand(end)
    for k in range(1,N,2):          #for odd terms
        s+=4*integrand(begin+k*h)
    for k in range(2,N,2):          #for even terms
        s+=2*integrand(begin+k*h)
    return (h/3)*s

N = 10000 # number of random points
k = 0 # will be used to hold the sum of all f(x)/w(x)
a,b=0,10 # beginning,end of integral

I=[] # this will hold all 100 of the the solutions as a list 
for j in range(100):
    for i in range(N):
        x = np.random.normal(5,1) # x input value in f(x) and w(x)
        k += f(x)/w(x)

    int_w=simpson(a,b,300,w) # integral of w(x) from x=0 to x=10
    I.append(k * int_w / N)
    k=0 #resets k
    
plt.hist(I, 10,range=[0.94,1.06])
plt.xlabel('x values')
plt.ylabel('number of solutions in bin')
plt.title('Importance sampling MC Solutions')
plt.show()
