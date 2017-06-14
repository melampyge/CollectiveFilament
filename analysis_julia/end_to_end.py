# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 10:24:14 2015

@author: jaeger

Here we calculate the angle between the rod and the end-to-end vector of the tail and analyse it for periodic signals
"""
import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from data_from_xyz import *

def end_to_end(xs,ys,sets,n,rod,loc, PLOT=True):
    sigma=40
    r_rod_x=[]
    r_rod_y=[]
    r_rod=[]
    r_tail_x=[]
    r_tail_y=[]
    r_tail=[]
    n_x=[]
    n_y=[]
    theta=[]
    h=[]
    for i in range(sets):
        r_rod_x.append(xs[i][0]-xs[i][rod-1])
        r_rod_y.append(ys[i][0]-ys[i][rod-1])
        r_rod.append(np.sqrt(r_rod_y[-1]**2+r_rod_x[-1]**2))
        n_x.append(-r_rod_y[-1]/r_rod[-1])
        n_y.append(r_rod_x[-1]/r_rod[-1])
        r_tail_x.append(xs[i][n-1]-xs[i][rod-1])
        r_tail_y.append(ys[i][n-1]-ys[i][rod-1])
        r_tail.append(np.sqrt(r_tail_y[-1]**2+r_tail_x[-1]**2))
        if (n_x[-1]*r_tail_x[-1]+n_y[-1]*r_tail_y[-1])<0:
            sign=1
        else:
            sign=-1
        #Calculation of angle and the height of the last atom
        theta.append(sign*(np.pi-np.arccos((r_rod_x[-1]*r_tail_x[-1]+r_rod_y[-1]*r_tail_y[-1])/float(r_rod[-1]*r_tail[-1]))))
        h.append(r_tail[-1]*np.sin(np.pi-theta[-1]))
    FTangle=np.fft.fft(theta-np.mean(theta))

    GF=lambda x:ndimage.filters.gaussian_filter(x,2)

    gauss_theta=ndimage.filters.gaussian_filter(theta,sigma)
    v=np.diff(gauss_theta)
    a=v>0
    changes=np.diff(a)
    number=np.sum(changes)

    if PLOT==True:
        f, axarr = plt.subplots(2, 2)
        axarr[0, 0].plot([0,6000],[0,0])
        axarr[0, 0].plot(theta)
        axarr[0, 0].set_ylim([-2.5,2.5])
        axarr[0, 0].set_title('Angle end-to-end')
        axarr[0, 1].plot([0,6000],[0,0])
        axarr[0, 1].plot(h)
        axarr[0, 1].set_ylim([-100,100])
        axarr[0, 1].set_title('Height of last atom')
        axarr[1, 0].plot(GF(np.abs(FTangle)))
        axarr[1, 0].set_xlim([0,500])
        axarr[1, 0].set_title('FT of angle')
        axarr[1, 1].hist(theta, bins=50)
        axarr[1, 1].set_xlim([-2.5,2.5])
        axarr[1, 1].set_title('Distribution of angle')   
        '''
        axarr[2, 0].plot(v) 
        lx,lX,ly,lY=axarr[2, 0].axis()
        YY=max((abs(ly),abs(lY)))
        axarr[2, 0].set_ylim([-YY,YY])
        axarr[2, 1].plot(changes,'.') 
        axarr[2, 1].set_ylim([-1,2])
        '''
        plt.setp([aa.get_xticklabels() for aa in axarr[0, :]], visible=False)
        plt.suptitle('End-to-end vector of tail', fontsize=24)
        f.tight_layout(renderer=None, pad=4, h_pad=None, w_pad=None, rect=None)
        plt.savefig(loc+'end_to_end.png',dpi=200)
        plt.show()

    return number
    
        
'''
location='/work/iff_th2/jaeger/simulations/noise_p001/tail_100/head_30/kappa_200/'
rod=30

run='run_p3.xyz'
folder='p3_new/'
name='Force=0.3'
data=read_coordinates(location+run, 0)
xs=data[0]
ys=data[1]
sets=data[2]
n=data[3]

end_to_end(xs,ys,sets,n,rod,location+folder,PLOT=True)
'''
