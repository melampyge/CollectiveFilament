# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:47:22 2015

@author: jaeger

Calculate the mean square rotation and plot it on a log-log graph
"""
import numpy as np
import matplotlib.pyplot as plt

      
#Plot a graph of the MSD averaged over the whole simulation    
    #xs= x-coordinates of all atoms over the whole simulation: xs[time][atom]
    #ys= y-coordinates of all atoms over the whole simulation: ys[time][atom]
    #sets= number of time steps in simulation
    #loc= path to where the graphic should be saved
    
def plot_msr(xs,ys,sets,rod,loc):
    
  #Calculate angle of rod for whole simulation
    x_f=[]                        
    y_f=[]
    x_b=[]
    y_b=[]
    r_x=[]
    r_y=[]
    theta=[0]
    angle=0
    turn=0
    for i in range(sets):
        x_f.append(xs[i][0])
        y_f.append(ys[i][0])
        x_b.append(xs[i][rod-1])
        y_b.append(ys[i][rod-1])
        r_x.append(xs[i][0]-xs[i][rod-1])
        r_y.append(ys[i][0]-ys[i][rod-1])
        l=np.sqrt(r_y[-1]**2+r_x[-1]**2)
        if r_y[-1]>0:
            if r_x[-1]<0:
                angle=np.arcsin(r_y[-1]/l)
            elif r_x[-1]==0:
                angle=np.pi/2
            elif r_x[-1]>0:
                angle=np.pi/2+np.arcsin(r_x[-1]/l)
        elif r_y[-1]<0:
            if r_x[-1]>0:
                angle=np.pi+np.arccos(r_x[-1]/l)
            elif r_x[-1]==0:
                angle=3*np.pi/2
            elif r_x[-1]<0:
                angle=2*np.pi-np.arccos(-r_x[-1]/l)
        elif r_y[-1]==0:
            if r_x[-1]>0:
                angle=np.pi
            elif r_x[-1]<0:
                angle=0
        if np.mod(theta[-1],2*np.pi)>7*np.pi/4 and angle<np.pi/4:
            turn+=1
        if np.mod(theta[-1],2*np.pi)<np.pi/4 and angle>7*np.pi/4:
            turn-=1
        theta.append(turn*2*np.pi+angle)
        
        
    MSR=[]
    error=[]
  #list of time intervals to calcuate the msd for 
    ts=[1,int(sets/50000.0),int(sets/30000.0),int(sets/10000.0),int(sets/6000.0),int(sets/3000.0),int(sets/1000.0),int(sets/500.0),int(sets/200.0),int(sets/100.0),int(sets/50.0),int(sets/20.0),int(sets*2/10.0),int(sets*3/10.0),int(sets*6/10.0),int(sets*8/10.0)]#2500,3000,3500,4000,4500,5000,5500,6000]
    for t in ts:
        av,std, data=msr(theta,sets,t)
        MSR.append(av)
        error.append(std/np.sqrt(len(data)))
        #if t==int(sets*2/10.0):
        #    print data, MSD[-1], error[-1]
        
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.errorbar(ts,MSR,yerr=error)
    ax.set_xscale('log')  
    ax.set_yscale('log')  
    ax.set_xlabel('time', fontsize=13)
    ax.set_ylabel('MSR', fontsize=13)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    fig.suptitle('MSR (log-log plot)',fontsize=16)
    plt.savefig(loc+'msr.png')
    plt.show()
  
    return MSR , error
    
    
    
    
#Calculate the msr for a given dt
#For the average over the simulation not every time is used for the interval but we jump by a number of timesteps that is given by step
    #xs= x-coordinates of all atoms over the whole simulation: xs[time][atom]
    #ys= y-coordinates of all atoms over the whole simulation: ys[time][atom]
    #sets= number of time steps in simulation
    #step
    
def msr(theta,sets,dt,step=1): #x is a list of x coordinates of one atom
    t=0
    msr=[]
    while (t+dt) < sets:
        msr.append(np.square(theta[t+dt]-theta[t]))
        t+=step 
    return np.mean(msr), np.std(msr), msr
    
    
    
    




















