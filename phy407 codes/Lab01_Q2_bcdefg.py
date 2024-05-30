# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 11:21:19 2020

@author: yonatan Eyob
"""

import numpy as np
import matplotlib.pyplot as plt

#there are two parts to this code. the first half will answer question 2b,2c,2f,2g. the second half will answer question 2d,2e.
#this half will answer question 2b,2c,2f,2g
#the pseudocode for 2a is done on the pdf
#in this half of the code i will be graphing  the xp vs p graph for different growth rates, graphing slightly different
# initial conditions to observe the chaos, and plotting the line of best fit of a semilogy graph to 
#solve for the Lyapunov exponent.

#i got the "figure,(ax1,ax2)=plt.subplots(nrows=2,ncols=1)" idea from: http://jonathansoma.com/lede/data-studio/classes/small-multiples/long-explanation-of-using-plt-subplots-to-create-small-multiples/
#these figures allow me to plot multiple graphs in 1 run without them stacking.
figure,(a)=plt.subplots(nrows=1,ncols=1)
figure2,(b)=plt.subplots(nrows=1,ncols=1)
figure3,(c)=plt.subplots(nrows=1,ncols=1)
figure4,(d)=plt.subplots(nrows=1,ncols=1)


#2b)
def pop(x0,r,pmax):     #input x0,r,pmax (all scalar values). outputs xp values for all p from 0 to pmax (array)
    i=0
    z=np.zeros(pmax)
    z[0]=x0
    while i<pmax-1:
        z[i+1]=r*(1-z[i])*z[i]
        i+=1
    return z



#2c)
x0=0.1          #initial normalized population
xf=x0+x0/10000  #x+e for question (2f)

r0=2            #growth rate = 2
r1=2.5          #growth rate = 2.5
r2=3            #growth rate = 3
r3=3.5          #growth rate = 3.5
r4=3.58         #growth rate = 3.58
r5=3.6          #growth rate = 3.6
r6=3.7          #growth rate = 3.7
rf=3.8          #choatic growth rate for question (2f)

pmax=50         #max years
pff=25          #max years for question (2g) and (2f)

x=np.arange(pmax)   #array with elements 0 to pmax used for the x axis values in question (2d)
xff=np.arange(pff)  #array with elements 0 to pff used for the x axis values in questions (2g) and (2f)

a.plot(x,pop(x0,r0,pmax),".", label='growth={0}'.format(r0))     #xp vs p plot with r=2
a.plot(x,pop(x0,r1,pmax),".", label='growth={0}'.format(r1))     #xp vs p plot with r=2.5  
a.plot(x,pop(x0,r2,pmax),".", label='growth={0}'.format(r2))     #xp vs p plot with r=3
a.plot(x,pop(x0,r3,pmax),".", label='growth={0}'.format(r3))     #xp vs p plot with r=3.5
a.legend(loc="lower right")         
a.set_ylabel("normalized population")
a.set_xlabel("years")
a.set_title("norm. population vs years")

d.plot(x,pop(x0,r4,pmax),".", label='growth={0}'.format(r4))     #xp vs p plot with r=3.58
d.plot(x,pop(x0,r5,pmax),".", label='growth={0}'.format(r5))     #xp vs p plot with r=3.6 
d.plot(x,pop(x0,r6,pmax),".", label='growth={0}'.format(r6))     #xp vs p plot with r=3.7
d.legend(loc="lower right")         
d.set_ylabel("normalized population")
d.set_xlabel("years")
d.set_title("norm. population vs years")

#2f)
b.plot(xff,pop(x0,rf,pff),".", label='x0={0}, growth={1}'.format(x0,rf)) #plots xp vs p using x0,rf,pff,xff
b.plot(xff,pop(xf,rf,pff),".", label='x0={0}, growth={1}'.format(xf,rf)) #plots xp vs p using xf,rf,pff,xff   
b.legend(loc="lower right")
b.set_ylabel("normalized population")
b.set_xlabel("years")
b.set_title("norm. population vs years")

#2g)
s=abs(pop(x0,rf,pff)-pop(xf,rf,pff)) #difference of xp(1)-xp2

def expo(a,n,y):   #exponential function: a*e**(n*y)
    m=a*np.exp(n*y)
    return m
c.plot(xff,expo(s[0],0.45,xff),label='exponential function {0}*e**({1}x)'.format(10**-5,0.45)) #line of best fit
c.semilogy(xff,s,"k.",label='semilogy plot')  #semilogy plot of s vs xff
c.legend(loc="lower right")
c.set_ylabel("normalized population")
c.set_xlabel("years")
c.set_title("norm. population vs years")




#======================================================================================================

#this half of the code answers question 2d,2e
# in this half of the code i will be graphing the bifurcation diagram from 4 different perspectives.

#i got the "figure,(ax1,ax2)=plt.subplots(nrows=2,ncols=1)" idea from: http://jonathansoma.com/lede/data-studio/classes/small-multiples/long-explanation-of-using-plt-subplots-to-create-small-multiples/
#these figures allow me to plot multiple graphs in 1 run without them stacking.
figure,(a2)=plt.subplots(nrows=1,ncols=1)
figure2,(b2)=plt.subplots(nrows=1,ncols=1)
figure3,(c2)=plt.subplots(nrows=1,ncols=1)
figure4,(d2)=plt.subplots(nrows=1,ncols=1)

#this function takes 4 inputs: xp (normalized population at year p), r (growth rate),
# pmax (max years), and cut (the last 100 or 1000 years) and graphs xp vs r
def xpVSr(r, xp,pmax,cut,axes):
    for i in range(pmax):
        xp= r*(1-xp)*xp
        if i>(pmax-cut):
            if axes==1:
                a2.plot(r, xp, 'k.',markersize=0.1)
                a2.set_xlim(min(r), max(r))
            if axes==2:
                b2.plot(r, xp, 'k.',markersize=0.1)
                b2.set_xlim(min(r), max(r))
            if axes==3:
                c2.plot(r, xp, 'k.',markersize=0.1)
                c2.set_xlim(min(r), max(r))
            if axes==4:
                d2.plot(r, xp, 'k.',markersize=0.1)
                d2.set_xlim(min(r), max(r))
    a2.set_title("Bifurcation")
    a2.set_xlabel("growth rate")
    a2.set_ylabel("normalized population")
    b2.set_title("Bifurcation")
    b2.set_xlabel("growth rate")
    b2.set_ylabel("normalized population")
    c2.set_title("Bifurcation")
    c2.set_xlabel("growth rate")
    c2.set_ylabel("normalized population")
    d2.set_title("Bifurcation")
    d2.set_xlabel("growth rate")
    d2.set_ylabel("normalized population")
       

#question 2d, 2e
#these are the inputs required for the function to work    
 
#necessary scalar inputs
pmax=2000      #total number of years considered 
cutl3=100      #considering the last 100 xp for r>3
cutg3=1000     #considering the last 1000 xp for r>3
x0=0.1         #initial normalized population 
dr= 0.005      #r increments are 0.005 for (2d)
dre=10**(-5)   #r increments are 10^-5 for (2e)
#end of scalar inputs

#necessary list inputs
r=np.arange(2,4,dr)             #2 < r < 4, increments of dr
rl3=np.arange(2,3,dr)           # 2 ≤ r < 3, increments of dr
rg3=np.arange(3,4,dr)           # 3 ≤ r < 4, increments of dr    
re=np.arange(3.738,3.745,dre)   #3.738 ≤ r < 3.745, increments of dre

x = x0*np.ones(len(r))          # a list of length r of numbers with the initial value being x0
xl3 = x0*np.ones(len(rl3))      # a list of length rl3 of numbers with the initial value being x0
xg3 = x0*np.ones(len(rg3))      # a list of length rg3 of numbers with the initial value being x0
xe = x0*np.ones(len(re))        # a list of length re of numbers with the initial value being x0
#end of list inputs

xpVSr(re,xe,pmax,cutg3,1)     #this function answers (2e)
xpVSr(rg3,xg3,pmax,cutg3,2)   #this function answers (2d) for r>3, considering last 1000 xp (xp being the normalized population during year p)
xpVSr(rl3,xl3,pmax,cutl3,3)   #this function answers (2d) for r<3, considering last 100 xp
xpVSr(r,x,pmax,cutg3,4)       #this is putting r>3 and r<3 into one graph while considering the last 1000 xp and 2<r<4 

