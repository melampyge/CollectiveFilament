# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:33:09 2015

@author: jaeger
"""

import numpy as np
import matplotlib.pyplot as plt

#Extract angular orientation of rod from coordinates of atoms and plot it 
#Then the general swimming direction is substracted and the signal fourier transformed

#angle zero for initial orientation along negative x axis, then increases clockwise and also increases further for 2Pi rotations 
    #xs= x-coordinates of all atoms over the whole simulation: xs[time][atom]
    #ys= y-coordinates of all atoms over the whole simulation: ys[time][atom]
    #sets= number of timesteps in simulation
    #rod= length of the rod in number of atoms
    #loc= path to where the graphic should be saved

def angle_of_rod(xs,ys,sets,rod,loc):
    """ compute the rod orientation and make a Fourier transform"""
    # compute rod angle, account for complete rotations
    r_x = xs[:,0] - xs[:,rod-1]
    r_y = ys[:,0] - ys[:,rod-1]
    theta = np.arctan2(r_x, r_y)
    dtheta = theta[1:] - theta[:-1]
    for i in range(sets-1):
        if dtheta[i] > 0.5*np.pi:
            dtheta[i] -= np.pi
        if dtheta[i] < -0.5*np.pi:
            dtheta[i] += np.pi
        theta[i+1] = theta[i] + dtheta[i]

    plt.figure()
    plt.plot(theta)
    #data=[theta[i]-smoothList(theta)[i] for i in range(len(theta))]
    #plt.plot(data)
    plt.title('Angle of rod from negative x-axis(increasing clockwise)')
    plt.savefig(loc+'angle.png')
    # plt.show()

    FTangle=np.fft.fft(theta)
    k=float(np.argmax(np.abs(FTangle)[:sets/2.0]))
    f=k/float(sets)
    #plt.figure()
    #plt.plot(np.abs(FTangle))   
    #plt.xlim(0,500)
    #plt.title('FT of straightened angle data')
    #plt.show()
    return theta, f
    
    
def smoothList(list,degree=200):   #use 3000 for noice case

     smoothed=np.zeros((len(list)))
     for i in range(len(smoothed)):  
         if i<degree/2:
             smoothed[i]=np.sum(list[i:i+degree])/float(degree)          
         elif i>=len(list)-degree/2:
             smoothed[i]=np.sum(list[i-degree:i])/float(degree)  
         else:
             smoothed[i]=np.sum(list[i-degree/2:i+degree/2])/float(degree)  

     return smoothed      
    
