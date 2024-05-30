# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:45:20 2020

@author: yonat, kivancaykac
"""

# in this code we will compare the relative and approximate errors of the gaussian quadrature,
# simpson rule, and trapezoidal rule then plot a blowing snow vs average 
# hourly temperature graph using the equation given to us in question 1b.
# %%
import numpy as np
import matplotlib.pyplot as plt
import gaussxw as gsx
from scipy import special

#1a 
#this is similar to the code our group wrote in the previous lab (not the solution code).
def dawson(x, t):  # Dawson function without the integral
    return np.exp(-1.*x**2)*np.exp(t**2)

def trap_dawson(N, a, b, integrand):  # Trapezoidal method
    h = (b-a)/N
    x = b
    s = 0.5*integrand(x, a) + 0.5*integrand(x, b)
    for k in range(1,N):
        s += integrand(x, a+k*h)
    return s*h

def simps_int(N, a, b, integrand):  # Simpson's Rule method
    h = (b-a)/N
    x = b
    s = integrand(x, a) + integrand(x, b)
    for i in range(1, N, 2):  # odds
        s += 4*integrand(x, a+i*h)
    for i in range(2, N, 2):  # evens
        s += 2*integrand(x, a+i*h)
    return s*h/3

def gaussian_quadrature(N,a,b,integrand):   #gaussian quadrature method
    t, w = gsx.gaussxwab(N,a,b)     #outputs the sample points (t) and weights (w) needed to use guassian quadrature
    I = 0
    for k in range(N):
        I += w[k]*integrand(4,t[k])
    return I

#1ai
for N in 2**np.arange(3,12):   # N is the number of slices
    print("For N={a} slices and x=4,\nTrapezoidal Rule gives: {b}\n\
Simpson's Rule gives: {c}\nThe Gaussian quadrature gives:{d}\n".format(a=N,b=trap_dawson(N, 0., 4., dawson), c=simps_int(N, 0., 4.,
dawson),d=gaussian_quadrature(N,0,4,dawson)))

    
#1aii
scipy_value = special.dawsn(4.0)
N = 2**np.arange(3,12)

def rel_error(true_val, exp_val):  # Magnitude of the relative error
    return np.abs(exp_val - true_val)/true_val

relerr_gaussian,apperr_gaussion,relerr_trap,apperr_trap,relerr_simp,apperr_simp = [np.zeros(len(N)) for _ in range(6)]
for i in range(len(N)):
    relerr_gaussian[i] = rel_error(scipy_value,
                                   gaussian_quadrature(N[i],0.,4.,dawson))
    apperr_gaussion[i] = np.abs(gaussian_quadrature(2*N[i],0,4,dawson
                                             )-gaussian_quadrature(N[i],0,4,
                                                                   dawson))
    relerr_trap[i] = rel_error(scipy_value, trap_dawson(N[i], 0, 4, dawson))
    apperr_trap[i] = np.abs(trap_dawson(2*N[i],0,4,dawson
                                        )-trap_dawson(N[i],0,4,dawson))/3
    relerr_simp[i] = rel_error(scipy_value, simps_int(N[i], 0, 4, dawson))
    apperr_simp[i] = np.abs(simps_int(2*N[i], 0, 4, dawson
                                      )-simps_int(N[i], 0, 4, dawson))/15

plt.style.use("seaborn-whitegrid")
plt.plot(N, relerr_gaussian, "k.",
         label="Gaussian Quadrature\n"+"Error, Relative")
plt.plot(N, apperr_gaussion, "b.",
         label="Gaussian Quadrature\n"+r"Error, $\epsilon_N$")
plt.plot(N, relerr_trap, "y.", label="Trapezoidal Method\n"+"Error, Relative")
plt.plot(N, apperr_trap, "g.",
         label="Trapezoidal Method\n"+r"Error, $\epsilon_N$")
plt.plot(N, relerr_simp, "r.", label="Simpson's Method\n"+"Error, Relative")
plt.plot(N, apperr_simp, "c.",
         label="Simpson's Method\n"+r"Error, $\epsilon_N$")
plt.title("Three Methods\n" + r"$|Error|$ and $\epsilon_N$ vs $N$")
plt.ylabel(r"$|Error|$")
plt.xlabel(r"$N$")
plt.yscale('log')  # on logarithmic scale for the y-axis
plt.xscale('log')  # on logarithmic scale for the x-axis
lgd = plt.legend(bbox_to_anchor=(.6, -.15), loc='upper left')
save = True
if(save): plt.savefig("fig1a.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()

# %%
#1b
def P_integrand(u,Ta,th):
    u_mean=11.2+0.365*Ta + 0.00706*Ta**2 + 0.9*np.log(th)
    S=4.3+0.145*Ta+0.00196*Ta**2
    return np.exp(-(u_mean-u)**2/(2*S**2))

def P(N,u10,Ta,th,P_integrand):
    u, w = gsx.gaussxwab(N,0,u10)
    S=4.3+0.145*Ta+0.00196*Ta**2
    I = 0
    for k in range(N):
        I += w[k]*P_integrand(u[k],Ta,th)
    return I*1/(np.sqrt(2*np.pi)*S)
 
Ta=np.arange(-50,30)    # temperature
th1=24                  # snow age
th2=48                  # snow age
th3=72                  # snow age
u10_1=6                 # wind speed
u10_2=8                 # wind speed
u10_3=10                # wind speed
N=100                   #number of slices/sample points in the gaussian quadrature
plt.plot(Ta,P(N,u10_1,Ta,th1,P_integrand),'r',label='u10={0},th={1}'.format(u10_1,th1))
plt.plot(Ta,P(N,u10_1,Ta,th2,P_integrand),'r-.',label='u10={0},th={1}'.format(u10_1,th2))
plt.plot(Ta,P(N,u10_1,Ta,th3,P_integrand),'r:',label='u10={0},th={1}'.format(u10_1,th3))

plt.plot(Ta,P(N,u10_2,Ta,th1,P_integrand),'y',label='u10={0},th={1}'.format(u10_2,th1))
plt.plot(Ta,P(N,u10_2,Ta,th2,P_integrand),'y-.',label='u10={0},th={1}'.format(u10_2,th2))
plt.plot(Ta,P(N,u10_2,Ta,th3,P_integrand),'y:',label='u10={0},th={1}'.format(u10_2,th3))

plt.plot(Ta,P(N,u10_3,Ta,th1,P_integrand),'k',label='u10={0},th={1}'.format(u10_3,th1))
plt.plot(Ta,P(N,u10_3,Ta,th2,P_integrand),'k-.',label='u10={0},th={1}'.format(u10_3,th2))
plt.plot(Ta,P(N,u10_3,Ta,th3,P_integrand),'k:',label='u10={0},th={1}'.format(u10_3,th3))
plt.legend(loc="upper right")
plt.ylabel('blowing snow probability')
plt.xlabel(r'average hourly temperature $^\circ C$')
plt.title('Temperature vs Probability of Blowing Snow')
save = True
if(save): plt.savefig("fig1b.png")
plt.show()
