# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:47:22 2015

@author: jaeger
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from data_from_xyz import *
      
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
	if vx[i]==0 and vy[i]>0:
		alpha=np.pi/2.0
	elif vx[i]==0 and vy[i]<0:
		alpha=-np.pi/2.0
	elif vy[i]==0 and vx[i]>0:
		alpha=np.pi
	elif vy[i]==0 and vx[i]<0:
		alpha=0
	elif vx[i]==0 and vy[i]==0:
		alpha=0
	else:
        	alpha=np.arctan(np.absolute(vy[i]/float(vx[i])))

        if vx[i]<0 and vy[i]<0:
            a.append(-alpha)
        elif vx[i]<0 and vy[i]>0:
            a.append(alpha)
        elif vx[i]>0 and vy[i]>0:
            a.append(np.pi-alpha)
        elif vx[i]>0 and vy[i]<0:
            a.append(-np.pi+alpha)
    
    GF=lambda x:ndimage.filters.gaussian_filter(x,2)    
    
    vx_s=smoothList(vx)
    vy_s=smoothList(vy)
    b=ndimage.filters.gaussian_filter(a,50)
    data=[a[i]-b[i] for i in range(len(a))]
    FTangle=np.fft.fft(a)
    FTangleb=np.fft.fft(data)
    
    sig=50
    f, axarr = plt.subplots(3, 2)

    axarr[0, 0].plot(vy[1:],label='y')
    axarr[0, 0].plot(vx[1:],label='x')
    #axarr[0, 0].plot(v[1:],label='total speed-smooth')
    lx,lX,ly,lY=axarr[0, 0].axis()
    YY=max((abs(ly),abs(lY)))
    axarr[0, 0].set_ylim([-YY,1+YY])
    axarr[0, 0].legend()
    axarr[0, 0].set_title('Velocity')

    axarr[0, 1].plot(ndimage.filters.gaussian_filter(vy,sig),label='y-smooth')
    axarr[0, 1].plot(ndimage.filters.gaussian_filter(vx,sig),label='x-smooth')
    axarr[0, 1].plot(ndimage.filters.gaussian_filter(v[1:],sig),label='total speed-smooth')
    lx,lX,ly,lY=axarr[0, 1].axis()
    YY=max((abs(ly),abs(lY)))
    axarr[0, 1].set_ylim([-YY,2+YY])
    axarr[0, 1].legend()
    axarr[0, 1].set_title('Velocity')

    axarr[1, 0].plot(a)
    lx,lX,ly,lY=axarr[1, 0].axis()
    YY=max((abs(ly),abs(lY)))
    axarr[1, 0].set_ylim([-YY,YY])
    axarr[1, 0].set_title('Swimming direction-Angle')

    axarr[1, 1].plot(data)
    lx,lX,ly,lY=axarr[1, 1].axis()
    YY=max((abs(ly),abs(lY)))
    axarr[1, 1].set_ylim([-YY,YY])
    axarr[1, 1].set_title('Oscillations in swimming')

    axarr[2, 0].plot(GF(np.abs(FTangle)))
    axarr[2, 0].set_xlim([0,500])
    axarr[2, 0].set_title('FT of swimming direction')
    axarr[2, 1].plot(GF(np.abs(FTangleb)))
    axarr[2, 1].set_xlim([0,500])
    axarr[2, 1].set_title('FT of oscillations')   

    plt.setp([j.get_xticklabels() for j in axarr[0, :]], visible=False)
    plt.setp([j.get_xticklabels() for j in axarr[1, :]], visible=False)
    plt.suptitle('Analysis of swimming', fontsize=16)
    plt.savefig(loc+'velocity_analysis.png')
    #plt.show()    
    
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






'''location='/work/iff_th2/jaeger/simulations/noice_1/tail_100/head_30/kappa_200/'
rod=30
run='run_p4.xyz'
folder='p4_new/'
name='Force=0.4'
data=read_coordinates(location+run, 0)
xs=data[0]
ys=data[1]
sets=data[2]
n=data[3]

plot_velocity(xs,ys,sets,5,location+folder,name)'''









