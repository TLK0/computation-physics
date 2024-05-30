# -*- coding: utf-8 -*-
"""


@author: Yonatan Eyob
"""
# in this code i use the Crank-Nicolson scheme to solve the time-dependent 
# Schrödinger equation in the square well, and then ill plot phi*phi  
# for T=(0,T/4,T/2,T). and then Ill plot the energy vs time and position vs time
# graphs.  
import numpy as np
import matplotlib.pyplot as plt
L=10**(-8)
m=9.109*10**(-31)
k=500/L
si=L/25
x0=L/5
x=np.linspace(-L/2,L/2,1024) # x values
dt=10**(-18)
a=abs(x[0]-x[1]) #step size
h=1.0545718*10**(-34)




def phii(x): #this is the non-normalized phi function. i will use this to find the norm constant with integration
    return np.exp(-(x-x0)**2/(4*si**2)+1j*k*x)

    
def simpson(N,x,phii):    #this function takes inputs number of slices (N), this outputs the integral of phii(x)^2.
    phi=abs(phii(x))**2
    hs=L/N
    s=phi[0]+phi[len(phi)-1] 
    for k in range(1,N,2):          #for odd terms
        s+=4*phi[k]
    for k in range(2,N,2):          #for even terms
        s+=2*phi[k]
    return (hs/3)*s

phi0=1/np.sqrt(simpson(len(x),x,phii)) #thisis the norm constant

def phi(x): # thisis the normalized phi function
    return np.exp(-(x-x0)**2/(4*si**2)+1j*k*x)*phi0






def V(x): #potential function for 1a
    return 0

def B(V,p): # B function for the diagonal elements of the Hamiltonian matrix
    A=-h**2/(2*m*a**2) 
    return V(p*a-L/2)-2*A

def H(N,V): # hamiltonian matrix
    A=-h**2/(2*m*a**2) # A value 
    H=np.zeros([N,N])
    for i in range(N):
        H[i][i]=B(V,i+1) # setting diagonal of the hamiltonian
    for i in range(1,N):
        H[i-1][i]=A # setting the subdiagonals
        H[i][i-1]=A
    return H


def I(N): #identity matrix
    return 1*np.eye(N, k=0)# NxN matrix, with ones on 1st super-diagonal

def LL(N,V): # L matrix in the L*phi_n+1=v equation
    return I(N)+1j*dt/(2*h)*H(N,V)

def v(N,V,phi_array):#  v vector in the L*phi_n+1=v equation
    R=I(N)-1j*dt/(2*h)*H(N,V)
    return np.matmul(R,phi_array)

def phi_t(x,T,N,V): #phi that will iterate through time. this function is very SLOW on my laptop...
    phi_array=phi(x)
    pht=phi(x)
    for i in range(int(round(T))):
        pht=np.linalg.solve(LL(N,V),v(N,V,phi_array))
        #pht=np.matmul(np.linalg.inv(LL(N,V)),v(N,V,phi_array))
        phi_array=pht
    return pht

N=1024  # Total position steps
T=3000  # Total timesteps

T=(0,T/4,T/2,T) # I'll plot the wavefunction at these four times. time = T*dt 
for i in T:
    plt.plot(x,np.conj(phi_t(x,i,N,V))*phi_t(x,i,N,V),label='timestep={0}'.format(int(i)))
    plt.title('square well wave function at time step={0}'.format(int(i)))
    plt.xlabel('position (m)')
    plt.ylabel('probability density')
plt.legend(loc='upper left')
plt.show()
    
    
    
    
    
    
#%%   
    
    
def simpson(x,T,N,V,integrand):    #this function takes inputs number of slices (N), this outputs the integral of phii(x)^2.
    phi=integrand(x,T,N,V)
    hs=L/N
    s=phi[0]+phi[len(phi)-1] 
    for k in range(1,N,2):          #for odd terms
        s+=4*phi[k]
    for k in range(2,N,2):          #for even terms
        s+=2*phi[k]
    return (hs/3)*s



def Energy_integrand(x,T,N,V): 
    return np.conj(phi_t(x,T,N,V))*np.matmul(H(N,V),phi_t(x,T,N,V))

def Energy(x,T,N,V): #outputs the energy at each time T that is inputted
    z=[]
    for i in T:
        z.append(simpson(x,i,N,V,Energy_integrand))
        print(i)
    return z

T=np.arange(0,3000,100) #these are the Times we will consider for the energy and position vs time graphs.
# the graph will look more choppy than expected because the time is going in steps of 100.
# this was done to save time because my tablet is too slow to run this code in a reasonable amount of time.
# if you want to plot more smooth graphs, let T=np.arange(0,3000,1).
plt.plot(T,Energy(x,T,N,V))
plt.ylabel('Energy (Joules)')
plt.xlabel('time steps')
plt.title('Square Well Energy vs Time plot')
plt.show()

#%%

def x_integrand(x,T,N,V):
    return np.conj(phi_t(x,T,N,V))*x*phi_t(x,T,N,V)

def X(x,T,N,V): #outputs the position at each time T that is inputted
    z=[]
    for i in T:
        z.append(simpson(x,i,N,V,x_integrand))
        print(i)
    return z


plt.plot(T,X(x,T,N,V))
plt.ylabel('position (m)')
plt.xlabel('time steps')
plt.title('Square Well Position vs Time plot')
plt.show()







    
        
























