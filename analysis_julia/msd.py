# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:47:22 2015

@author: jaeger

Calculate the mean square displacement and plot it on a log-log graph
"""
import numpy as np
import matplotlib.pyplot as plt

      
#Plot a graph of the MSD averaged over the whole simulation    
    #xs= x-coordinates of all atoms over the whole simulation: xs[time][atom]
    #ys= y-coordinates of all atoms over the whole simulation: ys[time][atom]
    #sets= number of time steps in simulation
    #loc= path to where the graphic should be saved
    
def plot_msd(xs,ys,sets,loc):
    MSD=[]
    error=[]
    x=[]
    y=[]
    for i in range(sets):
        x.append(np.average(xs[i]))
        y.append(np.average(ys[i]))
  #list of time intervals to calcuate the msd for 
    ts=[1,int(sets/50000.0),int(sets/30000.0),int(sets/10000.0),int(sets/6000.0),int(sets/3000.0),int(sets/1000.0),int(sets/500.0),int(sets/200.0),int(sets/100.0),int(sets/50.0),int(sets/20.0),int(sets*2/10.0),int(sets*3/10.0),int(sets*6/10.0),int(sets*8/10.0)]#2500,3000,3500,4000,4500,5000,5500,6000]
    for t in ts:
        av,std, data=msd(x,y,sets,t)
        MSD.append(av)
        error.append(std/np.sqrt(len(data)))
        #if t==int(sets*2/10.0):
        #    print data, MSD[-1], error[-1]
        
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.errorbar(ts,MSD,yerr=error)
    ax.set_xscale('log')  
    ax.set_yscale('log')  
    ax.set_xlabel('time', fontsize=13)
    ax.set_ylabel('MSD', fontsize=13)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    fig.suptitle('MSD (log-log plot)',fontsize=16)
    plt.savefig(loc+'msd.png')
    plt.show()
  
    return MSD , error
    
    
    
    
#Calculate the msd for a given dt
#For the average over the simulation not every time is used for the interval but we jump by a number of timesteps that is given by step
    #xs= x-coordinates of all atoms over the whole simulation: xs[time][atom]
    #ys= y-coordinates of all atoms over the whole simulation: ys[time][atom]
    #sets= number of time steps in simulation
    #step
    
def msd(x,y,sets,dt,step=1): #x is a list of x coordinates of one atom
    t=0
    msd=[]
    while (t+dt) < sets:
        msd.append(np.square(x[t+dt]-x[t])+np.square(y[t+dt]-y[t]))
        t+=step 
    return np.mean(msd), np.std(msd), msd
    
    
    
    
#Calculate speed of first atom and smooth out the oscillations    
    
def plot_velocity(xs,ys,sets,timeunit,loc):    
    x=[]
    y=[]
    vx=[]
    vy=[]
    v=[]
    x0=xs[0][0]
    y0=ys[0][0]
    vav=0
    for i in range(sets):
        x.append(xs[i][0])
        y.append(ys[i][0]) 
        vx.append((x[-1]-x0)/float(timeunit))
        x0=x[-1]
        vy.append((y[-1]-y0)/float(timeunit))
        y0=y[-1]
        vel=np.sqrt(vx[-1]**2+vy[-1]**2)
        v.append(vel)
        vav+=vel        
    av=vav/float(sets)
  #Angle of swimming direction
    a=[]
    for i in range(len(vx)):
        alpha=np.arctan(np.absolute(vy[i]/float(vx[i])))
        if vx[i]<0 and vy[i]<0:
            a.append(-alpha)
        elif vx[i]<0 and vy[i]>0:
            a.append(alpha)
        elif vx[i]>0 and vy[i]>0:
            a.append(np.pi-alpha)
        elif vx[i]>0 and vy[i]<0:
            a.append(-np.pi+alpha)
    plt.plot(vx[1:],label='x-component')
    plt.title('Velocity')
    plt.plot(vy[1:],label='y-component')
    vx_s=smoothList(vx)
    vy_s=smoothList(vy)
    plt.plot(vy_s,label='y-smooth')
    plt.plot(vx_s,label='x-smooth')
    plt.plot(v[1:],label='total speed')
    plt.legend()
    plt.savefig(loc+'velocity.png')
    plt.show()

    b=smoothList(a,degree=100)
    data=[a[i]-b[i] for i in range(len(a))]
    plt.plot(a)
    plt.plot(b)
    plt.show()

    FTangle=np.fft.fft(a)
    plt.plot(FTangle.real, label='Re')
    plt.plot(FTangle.imag, label='Im')    
    plt.xlim(0,500)
    plt.title('FT of angle of swimming direction')
    plt.show()      
    
    FTangle=np.fft.fft(data)
    plt.plot(FTangle.real, label='Re')
    plt.plot(FTangle.imag, label='Im')    
    #plt.plot(FTangle.real+FTangle.imag)
    plt.xlim(0,500)
    plt.title('FT of smoothed angle of swimming direction')
    plt.savefig(loc+'FT_smoothed_velocity.png')
    plt.show()
    
    return av,vx_s,vy_s
    
    
    
    
    
#Smoothes out the data given in list by averaging over the number of adjacent entries given in degree    

def smoothList(list,degree=300):   #use 3000 for noice case

     smoothed=[0]*len(list)  

     for i in range(len(smoothed)):  
         if i<degree/2:
             smoothed[i]=sum(list[i:i+degree])/float(degree)          
         elif i>=len(list)-degree/2:
             smoothed[i]=sum(list[i-degree:i])/float(degree)  
         else:
             smoothed[i]=sum(list[i-degree/2:i+degree/2])/float(degree)  

     return smoothed  




#Plot the x and y displacement of the filament(first atom) over the whole simulation

def plot_displacements(xs,ys,sets):    
    x=[]
    y=[]
    for i in range(sets):
        x.append(xs[i][0])
        y.append(ys[i][0]) 
    plt.plot(x)
    plt.title('Position x')
    plt.show()
    plt.plot(y)
    plt.title('Position y')
    plt.show()



















