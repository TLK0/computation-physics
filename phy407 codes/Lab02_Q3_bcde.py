# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 18:06:56 2020

@author: yonatan eyob
"""
#this code will cover question 3. in this code, we use the intensity function given to us on page 206 in the
#textbook to make density plots for 3 dufferent transmission functions.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib as mpl

#every measurement that is measuring distance will be in meters in this code.
#these figure arguments are used to make 7 different graphs. figure 5, 6, and 7 are for color bars
figure,(a)=plt.subplots(nrows=1,ncols=1)
figure2,(b)=plt.subplots(nrows=1,ncols=1)
figure3,(c)=plt.subplots(nrows=1,ncols=1)
figure4,(d)=plt.subplots(nrows=1,ncols=1)
figure5,(e)= plt.subplots(figsize=(10, 1))
figure5.subplots_adjust(bottom=0.6)
figure6,(f)= plt.subplots(figsize=(10, 1))
figure6.subplots_adjust(bottom=0.6)
figure7,(g)= plt.subplots(figsize=(10, 1))
figure7.subplots_adjust(bottom=0.6)


#3b==================================================================================
def q(u):                          #transmission function 
    seperation=20*10**(-6)         #given seperation is 20 micrometers 
    a=np.pi/seperation             #seperation=pi/a
    return np.sin(a*u)**2

#3c==================================================================================
#we want to integrate with respect to du
def integrand(u,x):          #this is the function we want to integrate with respect to du. x is the distance from the center of the screen
    f=1                      #focal length
    wave=500*10**(-9)        #wavelength
    return np.sqrt(q(u))*np.exp((2j*np.pi*x*u)/(wave*f))

def simpson(begin,end,N,x):    #this function takes inputs beginning (integral starts) and end (integral ends) to an integral, number of slices (N), and an array or scalar of x for Intensity(x). this outputs the intensity value at the given x input.
    h=(end-begin)/N            #step size
    s=integrand(begin,x)+integrand(end,x)
    for k in range(1,N,2):          #for odd terms
        s+=4*integrand(begin+k*h,x)
    for k in range(2,N,2):          #for even terms
        s+=2*integrand(begin+k*h,x)
    return abs((h/3)*s)

def array(x,N,w):              #this function puts simpson values in an array with length equal to len(x).
                               #this was made so i could use the output array on imshow(). 
    begin=-w/2
    end=w/2
    z=np.zeros(len(x))
    i=0
    while i < len(x):
        z[i]=simpson(begin,end,N,x[i])
        i+=1
    return z


#these are the inputs needed for 3d,3ei###
screen_width=0.1               #this is the width of the screen
Num_of_x_intervals=250         #this is how many x values we will find the intensities for
x1=np.linspace(-screen_width/2,screen_width/2,Num_of_x_intervals)       #choices of x in intensity function I(x) 
x=np.concatenate((x1,np.zeros(1))) #putting 0 in the array because i know a max intensity occurs at x=0
x.sort()
seperation=20*10**(-6)                      #slit seperation length
w=10*seperation                             #total width of grating
N=500                                      #number of slices in the simpson function  
#end of inputs###


# Make an Intensity vs position plot for a visual:
a.plot(x,simpson(-w/2,w/2,N,x))
a.set_xlabel("x(meters)")
a.set_ylabel("intensity")
a.set_title("difraction for transmission function $sin(au)^{2}$")
#3d===================================================================================
# Making a density plot
z=np.array([array(x,N,w)])
b.imshow(z,extent=[-screen_width/2,screen_width/2,0,screen_width/4],cmap=cm.hot)
b.set_xlabel("x(meters)")
b.grid(None)
b.set_title("difraction for transmission function $sin(au)^{2}$")


#3ei
def q2(u):     #new transmission function for question 3ei
    return q(u)*q(u/2)

def integrand2(u,x):          #this is the function we want to integrate over
    f=1                      #focal length
    wave=500*10**(-9)        #wavelength
    return np.sqrt(q2(u))*np.exp((2j*np.pi*x*u)/(wave*f))

def simpson2(begin,end,N,x):    #this function takes inputs beginning and end to an integral, number of slices, and an array of x for I(x). this outputs the intensity value at the given x
    h=(end-begin)/N            #step size
    s=integrand2(begin,x)+integrand2(end,x)
    for k in range(1,N,2):          #for the odd terms
        s+=4*integrand2(begin+k*h,x)
    for k in range(2,N,2):          #for the even terms
        s+=2*integrand2(begin+k*h,x)
    return abs((h/3)*s)

def array2(x,N,w):              #this function puts simpson values in an array with length equal to len(x) 
    begin=-w/2
    end=w/2
    z=np.zeros(len(x))
    i=0
    while i < len(x):
        z[i]=simpson2(begin,end,N,x[i])
        i+=1
    return z

#the only input we have to update from the previous density plots is z (for the new array2 function)
z2=np.array([array2(x,N,w)])
c.imshow(z2,extent=[-screen_width/2,screen_width/2,0,screen_width/4],cmap=cm.hot)
c.set_xlabel("x(meters)")
c.grid(None)
c.set_title("difraction for transmission function $sin(au)^{2}*sin(au/2)^{2}$")

#=================================================================================================
#3eii


def q3(u):   #new transmission function for 3eii. this function is piecewise. it is 1 if u is within the "slit range" and 0 everywhere else
    seperation=60*10**(-6)                       #seperation of the slits in meters
    width_of_slit1=10*10**(-6)                   #width of the 10 micrometer slit
    width_of_slit2=20*10**(-6)                   #width of the 20 micrometer slit
    if -seperation/2 - width_of_slit1 < u < -seperation/2:
        return 1
    elif seperation/2 < u < seperation/2 + width_of_slit2:
        return 1
    else:
        return 0
def integrand3(u,x):          #this is the function we want to integrate over
    f=1                      #focal length
    wave=500*10**(-9)        #wavelength

    return np.sqrt(q3(u))*np.exp((2j*np.pi*x*u)/(wave*f))

def simpson3(begin,end,N,x):    #this function takes inputs beginning and end to an integral, number of slices, and an array of x for I(x). this outputs the intensity value at the given x
    h=(end-begin)/N            #step size
    s=integrand3(begin,x)+integrand3(end,x)
    for k in range(1,N,2):          #for the odd terms
        s+=4*integrand3(begin+k*h,x)
    for k in range(2,N,2):          #for the even terms
        s+=2*integrand3(begin+k*h,x)
    return abs((h/3)*s)

def array3(x,N,w):              #this function puts simpson values in an array with length equal to len(x) 
    begin=-w/2
    end=w/2
    z=np.zeros(len(x))
    i=0
    while i < len(x):
        z[i]=simpson3(begin,end,N,x[i])
        i+=1
    return z

#the only input we have to update from the previous density plots are w (for the new total width of the screen)
#and z (for the new array3 function)
    
#when the seperation of the slits is 60 micrometers: 
seperation=60*10**(-6)                       #seperation of the slits in meters
width_of_slit1=10*10**(-6)                   #width of the 10 micrometer slit
width_of_slit2=20*10**(-6)                   #width of the 20 micrometer slit
w=seperation + width_of_slit1 + width_of_slit2     #total width of grating
z3=np.array([array3(x,N,w)])
d.imshow(z3,extent=[-screen_width/2,screen_width/2,0,screen_width/4],cmap=cm.hot)
d.set_xlabel("x(meters)")
d.grid(None)
d.set_title("difraction from 2 non identical slits 60 micrometers apart")


#this is the colour bar, used to indicate which color on the density graph refers to what Intensity.
#i made these bar graphs seperately from the density plot so i can adjust the size and orientation on latex 
#to what I prefer.
# i found info on how to set these color bars up from this website:https://matplotlib.org/3.1.0/tutorials/colors/colorbar_only.html
#the width of the grating is different for question 3d/3ei and question 3eii so I'll update the width.
#this color bar is for question 3d, with transmission function q(u)
seperation=20*10**(-6)                      #slit seperation length for quuestion 3d
w=10*seperation                             #total width of grating for quesion 3d
cmap = mpl.cm.hot
norm = mpl.colors.Normalize(vmin=0, vmax=max(simpson(-w/2,w/2,N,x)))
cb1 = mpl.colorbar.ColorbarBase(e, cmap=cmap,norm=norm, orientation='horizontal')
cb1.set_label('Intensity($W/m^{2}$) (question 3d)')
print("this is the max intensity in the density graph for quesion 3d: ",max(simpson(-w/2,w/2,N,x)))

#this color bar is for question 3ei, with transmission function q2(u)
seperation=20*10**(-6)                      #slit seperation length for quuestion 3ei
w=10*seperation                             #total width of grating for quesion 3ei
cmap = mpl.cm.hot
norm = mpl.colors.Normalize(vmin=0, vmax=max(simpson2(-w/2,w/2,N,x)))
cb1 = mpl.colorbar.ColorbarBase(f, cmap=cmap,norm=norm, orientation='horizontal')
cb1.set_label('Intensity($W/m^{2}$) (question 3ei)')
print("this is the max intensity in the density graph for quesion 3ei: ",max(simpson2(-w/2,w/2,N,x)))

#this color bar is for question 3eii, with transmission function q3(u)
seperation=60*10**(-6)                       #seperation of the slits in meters for question 3eii
width_of_slit1=10*10**(-6)                   #width of the 10 micrometer slit for question 3eii
width_of_slit2=20*10**(-6)                   #width of the 20 micrometer slit for question 3eii
w=seperation + width_of_slit1 + width_of_slit2     #total width of grating for question 3eii
cmap = mpl.cm.hot
norm = mpl.colors.Normalize(vmin=0, vmax=max(simpson3(-w/2,w/2,N,x)))
cb1 = mpl.colorbar.ColorbarBase(g, cmap=cmap,norm=norm, orientation='horizontal')
cb1.set_label('Intensity($W/m^{2}$) (question 3eii)')
print("this is the max intensity in the density graph for quesion 3eii: ",max(simpson3(-w/2,w/2,N,x)))

