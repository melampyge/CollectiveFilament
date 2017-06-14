# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 08:59:28 2015

@author: jaeger

Compare the cumulative sum of the eigenvalues for differnt noise level
"""
import matplotlib.pyplot as plt
import numpy as np

force=['p05','p1','p2','p3','p4','p6', 'p8','1']
f=[0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0]
flexure=[i*np.power(100,3)/200.0 for i in f]

colors=plt.cm.gist_heat(np.linspace(0,0.8,8))  


ev_head_30_1=np.zeros((8,128))
ev_head_30_p1=np.zeros((8,128))
ev_head_30_p01=np.zeros((8,128))

cu_ev_head_30_1=np.zeros((8,128))
cu_ev_head_30_p1=np.zeros((8,128))
cu_ev_head_30_p01=np.zeros((8,128))

noise_p1=[]
noise_p01=[]
noise_1=[]
noise_p001=[0.0994437888534,0.951114624964,0.981635679899,0.983459197284,0.984873885266,0.98217970712,0.984055956552,0.985333983288]

for i in range(8):
    
    ev_head_30_1[i]=np.load('/work/iff_th2/jaeger/simulations/noise_1/tail_100/head_30/kappa_200/'+force[i]+'_new/eigenvalues.npy')
    cu_ev_head_30_1[i]=ev_head_30_1[i].cumsum()/ev_head_30_1[i].sum()

    ev_head_30_p1[i]=np.load('/work/iff_th2/jaeger/simulations/noise_p1/tail_100/head_30/kappa_200/'+force[i]+'/eigenvalues.npy')
    cu_ev_head_30_p1[i]=ev_head_30_p1[i].cumsum()/ev_head_30_p1[i].sum()
    
    ev_head_30_p01[i]=np.load('/work/iff_th2/jaeger/simulations/noise_p01/tail_100/head_30/kappa_200/'+force[i]+'/eigenvalues.npy')
    cu_ev_head_30_p01[i]=ev_head_30_p01[i].cumsum()/ev_head_30_p01[i].sum()
    
    noise_p01.append(cu_ev_head_30_p01[i][1])
    noise_p1.append(cu_ev_head_30_p1[i][1])
    noise_1.append(cu_ev_head_30_1[i][1])
    
 
  
fig= plt.figure(figsize=(12,8))
ax= fig.add_subplot(1,1,1)
ax.plot(flexure,noise_p001,'--o', label='kT=0.001',markersize=9, color=colors[0])
ax.plot(flexure,noise_p01,'--o', label='kT=0.01',markersize=9,color=colors[2])
ax.plot(flexure,noise_p1,'--o', label='kT=0.1',markersize=9,color=colors[4])
ax.plot(flexure,noise_1,'--o', label='kT=1',markersize=9,color=colors[6])

ax.set_xlabel(r'$f_p l^3/ \kappa$', fontsize=26)
ax.set_ylabel(r'$\frac{ev(1)+ ev(2)}{total \ sum}$', fontsize=28)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)
fig.suptitle('Sum of 2 largest eigenvalues for differnt noise levels',fontsize=24)
fig.tight_layout(renderer=None, pad=3, h_pad=None, w_pad=None, rect=None)
plt.legend(loc=4)
plt.show()

