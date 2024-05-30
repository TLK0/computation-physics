# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:35:30 2020

@author: Yonatan Eyob
"""
# in this code we will load and plot data from dow.txt and then take the FFT of
# the data, set the last 90% of the FFT elements to zero, the plot the Inverse
# FFT, and then do the same thing but instead we will set the last 98% of the 
# FFT elements to zero. We will plot the 3 plots on the same graph to observe
# the differences in the plots.

import numpy as np
import matplotlib.pyplot as plt

figure,(c)=plt.subplots(nrows=1,ncols=1)

# 1b)
# Exercise 7.4
# 7.4,a)

dow=np.loadtxt('dow.txt') # load in the dow file
x=np.linspace(0,48,len(dow))     # x axis to plot the dow file elements (in months)
c.plot(x,dow,label='original dow data')

#7.4,b)
Fourier_dow=np.fft.rfft(dow)   # taking the FFT of the dow file

z=np.zeros(len(Fourier_dow)*9//10+1)   #making an array of zeros equal to 90% the length of Fy 
Fourier_dow=np.concatenate([Fourier_dow[:len(Fourier_dow)//10],z]) #combining the first 10% of Fourier_dow with z
# now  Fourier_dow is the same length as before but only the first 10% of the elements are non-zero

#7.4,c)
Iy=np.fft.irfft(Fourier_dow)    #inverse fourier transform
c.plot(x,Iy,'r-.',label='10% of the FFT')


#7.4,e)
Fourier_dow=np.fft.rfft(dow)    # taking the FFT of the dow file
z=np.zeros(len(Fourier_dow)*98//100+1)    #making an array of zeros equal to 98% the length of Fourier_dow
Fourier_dow=np.concatenate([Fourier_dow[:len(Fourier_dow)*2//100],z])   #combining the first 2% of Fourier_dow with z
# now  Fourier_dow is the same length as before but only the first 2% of the elements are non-zero

Iy=np.fft.irfft(Fourier_dow)    #inverse fourier transform
c.plot(x,Iy,'k-.',markersize=0.8,label='2% of the FFT')


c.legend(loc="lower left")
c.set_xlabel('months')
c.set_ylabel('value of the Dow')
c.set_title('The Dow from late 2006 to the end of 2010')




 