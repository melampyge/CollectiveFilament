# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 08:59:28 2015

@author: jaeger

This script reads in the eigenvalues from npy-files  
"""
import matplotlib.pyplot as plt
import numpy as np

force=['p05','p1','p2','p3','p4','p6', 'p8','1']

ev_head_5=np.zeros((8,103))
ev_head_5_noise=np.zeros((8,103))
ev_head_10=np.zeros((8,108))
ev_head_10_noise=np.zeros((8,108))
ev_head_15=np.zeros((8,113))
ev_head_15_noise=np.zeros((8,113))
ev_head_20=np.zeros((8,118))
ev_head_20_noise=np.zeros((8,118))
ev_head_30=np.zeros((8,128))
ev_head_30_noise=np.zeros((8,128))
ev_head_40=np.zeros((8,138))
ev_head_40_noise=np.zeros((8,138))
ev_head_60=np.zeros((8,158))
ev_head_60_noise=np.zeros((8,158))


for i in range(8):
    ev_head_5[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_5/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_5_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_5/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    ev_head_10[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_10/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_10_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_10/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    ev_head_15[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_15/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_15_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_15/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    ev_head_20[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_20/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_20_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_20/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    ev_head_30[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_30/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_30_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_30/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    ev_head_40[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_40/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_40_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_40/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    ev_head_60[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_60/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    ev_head_60_noise[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_60/kappa_200/'+force[i]+'_new/eigenvalues_noise.npy')
    
colors=plt.cm.gist_heat(np.linspace(0,0.8,8))    
    
for i in range(8):
    plt.plot(ev_head_5[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=5')
plt.show()

for i in range(8):
    plt.plot(ev_head_10[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=10')
plt.show()

for i in range(8):
    plt.plot(ev_head_15[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=15')
plt.show()

for i in range(8):
    plt.plot(ev_head_20[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=20')
plt.show()

for i in range(8):
    plt.plot(ev_head_30[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=30')
plt.show()

for i in range(8):
    plt.plot(ev_head_40[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=40')
plt.show()

for i in range(8):
    plt.plot(ev_head_60[i],'--o',markersize=3, label=force[i], color=colors[7-i])
plt.legend()
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.title('Rod=60')
plt.show()