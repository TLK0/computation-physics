# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 00:19:50 2020

@author: Yonatan Eyob
"""

# in this code we will plot the kinetic, potential, and total energy at each time step
# of the system in question 3a

import numpy as np
import matplotlib.pyplot as plt

figure,(b)=plt.subplots(nrows=1,ncols=1)


#this is given=================================================================
N = 16
Lx = 4.0
Ly = 4.0
dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)
x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)
xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()
#==============================================================================


def accel_of_particles(r,N): #input is position of particles output is acceleration on each particle 
    accel=np.zeros((N,2))
    rn=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if i==j:
                continue
            ri=np.zeros(2)
            rx=r[i,0]-r[j,0]    #distance particles are apart in the x direction
            ry=r[i,1]-r[j,1]    #distance particles are apart in the y direction
            ri=(rx,ry)
            rn[i][j]=sum(x*x for x in ri)**(1/2) #distance (r) that particles are apart
            accel[i][0]+=12*rx*rn[i][j]**(-8)*(2*rn[i][j]**(-6)-1)  #acceleration in the x direction
            accel[i][1]+=12*ry*rn[i][j]**(-8)*(2*rn[i][j]**(-6)-1)  #acceleration in the y direction
    return accel


def energy(r,N):    # r is the initial positions of the particle, N is the number of particles. this function outputs the total energy, kinetic, and potential energy at each timestep. 
    dt=0.01     #step size
    time=np.arange(0,10,dt)
    v0=1/2*dt*accel_of_particles(r,N)    # initial velocity
    r01=r   # original position of particles
    kinetic=np.ones(len(time))       
    sub_potential=np.zeros(N)   #placeholder
    potential=np.zeros(len(time))
    total_energy=np.zeros(len(time))    
    y=np.zeros(N)   #placeholder
    rdist=np.zeros((N,N)) #placeholder
    i=0
    for h in time:
        r1new=r01+dt*v0 # equation 8 on handout, update position
        k=dt*accel_of_particles(r1new,N)   # equation 9 on handout
        v0=v0+k     # equation 11 on handout, update velocity
        r01=r1new   # updating r(t)
        for j in range(N):
            y[j]=sum(n*n for n in v0[j])*0.5
        kinetic[i]=np.sum(y) # total kinetic energy for i^th timestep

        for a in range(N):
            for b in range(N):
                if a==b:
                    continue
                ri=np.zeros(2)
                rx=r01[a,0]-r01[b,0]    #distance particles are apart in the x direction
                ry=r01[a,1]-r01[b,1]    #distance particles are apart in the y direction
                ri=(rx,ry)
                rdist[a][b]=sum(n*n for n in ri)**(1/2) #distance (r) that particles are apart
                sub_potential[a]+=4*(rdist[a][b]**(-12)-rdist[a][b]**(-6))
        potential[i]= np.sum(sub_potential)/4 # total potential energy for i^th timestep
        sub_potential=sub_potential*0 # resetting placeholder
        
        total_energy[i]=kinetic[i]+potential[i] # total Energy for i^th timestep
        i+=1
    return kinetic,potential,total_energy


# these are the inputs=========================================================
N=16    # number of particles

r=np.zeros((N,2)) # this will be the x,y positions of the particles.
for i in range(N):
    r[i]=x_initial[i],y_initial[i]  # now we have the x,y positions of each particle in a nice single array r

time=np.arange(0,10,0.01) # i want the x axis to span from 0 to 10
# end of inputs================================================================


b.plot(energy(r,N)[0],label='kinetic energy')
b.plot(energy(r,N)[1],label='potential energy')
b.plot(energy(r,N)[2],label='total energy')
b.legend(loc='lower right')
b.set_xlabel('time step')
b.set_ylabel('Energy')
b.set_title('Energy vs Time graph')