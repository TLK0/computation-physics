# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:38:36 2020

@author: yonatan Eyob
"""
import matplotlib.pyplot as plt
import numpy as np

# In this code we will be finding the Energy and plotting the eigenfunction of 
# R(r) with initial conditions of n=1 and l=0, n=2 and l=0, n=2 and l=1 using RK4
# and the shooting method, and then compare our numerically calculated solutions
# with the explicit solutions given to us on the lab handout. 
# This code is an edited version of the squarewell.py code given by Newman.

# these figures b,c,d are used to plot three graphs in the same code. I prefer this set up.
figure,(b)=plt.subplots(nrows=1,ncols=1)
figure,(c)=plt.subplots(nrows=1,ncols=1)
figure,(d)=plt.subplots(nrows=1,ncols=1)

# questions 3a,3b, and 3c are all done at the same time:

# Constants
m = 9.1094e-31     # Mass of electron
hbar = 1.0546e-34  # Planck's constant over 2*pi
e = 1.6022e-19     # Electron charge
a = 5.2918e-11     # Bohr radius
e0=8.854187817*10**(-12) #vacuum permittivity
L = 20*a    #width of the well
h = 0.0001*a    #step size in the RK4 method
# end of constants


# functions for questions 3a,3b,3c code:
def V(x):   # Potential function
    return (-e**2)/(4*np.pi*e0*x)

def dSR(update,x,E,l):  # this outputs the derivative of R and S as a vector (dR,dS)
    R=update[0]
    S=update[1]
    dR=S/(x**2)
    dS=(l*(l+1)+2*m*x**2/hbar**2*(V(x)-E))*R
    return np.array([dR,dS])
    
def solve1(E,l): # uses RK4 to solve the system of odes in dSR
    R_initial = 0.0 # initial R value
    S_initial = 1.0 # initial S value
    r = np.array([R_initial,S_initial],float) #initial R(r) and S(r) values
    R = np.zeros(int(L//h)+1) # this will be used to store all R(r) values
    i=1
    for x in np.arange(h,L,h): #RK4 loop 
        k1 = h*dSR(r,x,E,l)
        k2 = h*dSR(r+0.5*k1,x+0.5*h,E,l)
        k3 = h*dSR(r+0.5*k2,x+0.5*h,E,l)
        k4 = h*dSR(r+k3,x+h,E,l)
        r += (k1+2*k2+2*k3+k4)/6
        R[i]=r[0] #stores r[0] at each step
        i+=1
    return r[0],R

def simpson(N,E,l):    #this function takes inputs number of slices (N), Energy (E), and l value (l). this outputs the integral of R(r)^2.
    R=(solve1(E,l)[1])**2 # this takes all the R(r) elements and squares them
    s=R[0]+R[len(R)-1] 
    for k in range(1,N,2):          #for odd terms
        s+=4*R[k]
    for k in range(2,N,2):          #for even terms
        s+=2*R[k]
    return abs(((L/a)/(3*N))*s) #stepsize= (L/a-0)/N=(L/a)/N

# End of functions=============================================================


# now we will find the energy using the secant method,
# and then plot the corresponding wavefunctions:
# first we will find the energy for n=1,l=0, and plot the corresponding wavefunction:
l=0 #degree l
n=1 #energy level

# shooting method working its magic
E1=-15*e/n**2   # reasonable guess of Energy
E2=-13*e/n**2   # reasonable guess of Energy
psi2 = solve1(E1,l)[0]
target = e/1000

while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve1(E2,l)[0]
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

x=np.linspace(0,(L/a),len(solve1(E2,l)[1])) # this will be the x-axis of the graphs

integral_R_squared=simpson(len(solve1(E2,l)[1]),E2,l) #this is the solution to the integral of R(r)^2 from r=0 to r=20a
b.plot(x,solve1(E2,l)[1]/integral_R_squared**(1/2),'k',label='numerically  calculated  solution') 
b.set_xlabel("position r in units of $a_{0}$ (bohr radius)")
b.set_ylabel('normalized R(r)')
b.grid(None)
b.set_title("normalized eigenfunction for n={a}, l={b}".format(a=n,b=l))
print("E(n={a},l={b}) ={c}eV".format(c=E2/e,b=l,a=n))



# now we will find the energy for n=2,l=0, and plot the corresponding wavefunction:
l=0 #degree l
n=2 #energy level

# shooting method go again
E1=-15*e/n**2   # reasonable guess of Energy
E2=-13*e/n**2   # reasonable guess of Energy
psi2 = solve1(E1,l)[0]
target = e/1000

while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve1(E2,l)[0]
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

integral_R_squared=simpson(len(solve1(E2,l)[1]),E2,l) #this is the solution to the integral of R(r)^2 from r=0 to r=20a
c.plot(x,solve1(E2,l)[1]/integral_R_squared**(1/2),'k',label='numerically  calculated  solution') 
c.set_xlabel("position r in units of $a_{0}$ (bohr radius)")
c.set_ylabel('normalized R(r)')
c.grid(None)
c.set_title("normalized eigenfunction for n={a}, l={b}".format(a=n,b=l))
print("E(n={a},l={b}) ={c}eV".format(c=E2/e,b=l,a=n))



# now we will find the energy for n=2,l=1, and plot the corresponding wavefunction:
l=1 #degree l
n=2 #energy level

# shooting method go aGAIN
E1=-15*e/n**2   # reasonable guess of Energy
E2=-13*e/n**2   # reasonable guess of Energy
psi2 = solve1(E1,l)[0]
target = e/1000

while abs(E1-E2)>target:
    psi1,psi2 = psi2,solve1(E2,l)[0]
    E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

integral_R_squared=simpson(len(solve1(E2,l)[1]),E2,l) #this is the solution to the integral of R(r)^2 from r=0 to r=20a
d.plot(x,solve1(E2,l)[1]/integral_R_squared**(1/2),'k',label='numerically  calculated  solution') 
d.set_xlabel("position r in units of $a_{0}$ (bohr radius)")
d.set_ylabel('normalized R(r)')
d.grid(None)
d.set_title("normalized eigenfunction for n={a}, l={b}".format(a=n,b=l))
print("E(n={a},l={b}) ={c}eV".format(c=E2/e,b=l,a=n))



# question 3d, plotting the explicit solutions:
def R10(x): # this is the explicit solution for n=1, l=0
    return 2*np.exp(-x/a)/a**(3/2)

def R20(x): # this is the explicit solution for n=2, l=0
    return (2-x/a)*np.exp(-x/(2*a))/(2*np.sqrt(2)*a**(3/2))

def R21(x): # this is the explicit solution for n=2, l=1
    return np.exp(-x/(2*a))*x/(2*np.sqrt(6)*a**(5/2))
    
def simpson(Rn,N):    # this new simpson function will be used to scale down the explicit solutions
    R=Rn**2
    s=R[0]+R[len(R)-1]
    for k in range(1,N,2):          #for odd terms
        s+=4*R[k]
    for k in range(2,N,2):          #for even terms
        s+=2*R[k]
    return abs(((L/a)/(3*N))*s)

x=np.linspace(0,(L/a),100)*a # x values for the R(x) function
I10=simpson(R10(x),len(R10(x))) # normalization constant of R10
I20=simpson(R20(x),len(R20(x))) # normalization constant of R20
I21=simpson(R21(x),len(R21(x))) # normalization constant of R21


b.plot(x/a,R10(x)/I10**(1/2),'r-.',label='explicit solution')
c.plot(x/a,R20(x)/I20**(1/2),'r-.',label='explicit solution')
d.plot(x/a,R21(x)/I21**(1/2),'r-.',label='explicit solution')

b.legend(loc='upper right')
c.legend(loc='upper right')
d.legend(loc='upper right')

