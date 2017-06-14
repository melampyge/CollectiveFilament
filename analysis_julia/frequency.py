# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 11:03:29 2015

@author: jaeger
"""
import numpy as np
import matplotlib.pyplot as plt

colors=plt.cm.gist_heat(np.linspace(0,0.8,6))
'''
plt.plot(force[7:],f_10[7:],'--',label='rod=10',color=colors[0])
plt.plot(force[4:],f_15[4:],'--',label='rod=15',color=colors[1])
plt.plot(force[4:],f_20[4:],'--',label='rod=20',color=colors[2])
plt.plot(force[1:],f_30[1:],'--',label='rod=30',color=colors[3])
plt.plot(force[1:],f_40[1:],'--',label='rod=40',color=colors[4])
plt.plot(force[1:],f_60[1:],'--',label='rod=60',color=colors[5])
plt.plot(force[1:],f_60_05[1:],color=colors[5])
plt.plot(force[1:],f_40_05[1:],color=colors[4])
plt.plot(force[1:],f_30_05[1:],color=colors[3])
plt.plot(force[4:],f_20_05[4:],color=colors[2])
plt.plot(force[4:],f_15_05[4:],color=colors[1])
plt.plot(force[7:],f_10_05[7:],color=colors[0])

plt.legend(loc=1)
plt.title('Frequency of beating')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$f$',fontsize=18)
plt.savefig(loc+'vary_rod_frequency1.png')
plt.show()
'''
fig= plt.figure(figsize=(12,8))
ax= fig.add_subplot(1,1,1)

ax.plot(force[7:],f_10_scaled[7:],'--',label='rod=10',linewidth=2,color=colors[0])
ax.plot(force[4:],f_15_scaled[4:],'--',label='rod=15',linewidth=2,color=colors[1])
ax.plot(force[4:],f_20_scaled[4:],'--',label='rod=20',linewidth=2,color=colors[2])
ax.plot(force[1:],f_30_scaled[1:],'--',label='rod=30',linewidth=2,color=colors[3])
ax.plot(force[1:],f_40_scaled[1:],'--',label='rod=40',linewidth=2,color=colors[4])
ax.plot(force[1:],f_60_scaled[1:],'--',label='rod=60',linewidth=2,color=colors[5])
ax.plot(force[1:],f_60_05_scaled[1:],color=colors[5])
ax.plot(force[1:],f_40_05_scaled[1:],color=colors[4])
ax.plot(force[1:],f_30_05_scaled[1:],color=colors[3])
ax.plot(force[4:],f_20_05_scaled[4:],color=colors[2])
ax.plot(force[4:],f_15_05_scaled[4:],color=colors[1])
ax.plot(force[7:],f_10_05_scaled[7:],color=colors[0])


ax.set_xlim([400,5100])
ax.set_xlabel(r'$f_p l^3/ \kappa$',fontsize=24)
ax.set_ylabel(r'$f \ \tau$',fontsize=24)
ax.legend(loc=4)
fig.suptitle(r'Beating frequency (dashed=$f_p$ is changed)',fontsize=32)
fig.tight_layout(renderer=None, pad=4, h_pad=None, w_pad=None, rect=None)
plt.savefig(loc+'vary_rod_frequency_scaled1.png')
plt.show()



'''
plt.plot(force2[7:],f_10_scaled[7:],'--',label='rod=10',color=colors[0])
plt.plot(force2_fp[4:],f_15_scaled[4:],'--',label='rod=15',color=colors[1])
plt.plot(force2_fp[4:],f_20_scaled[4:],'--',label='rod=20',color=colors[2])
plt.plot(force2_fp[1:],f_30_scaled[1:],'--',label='rod=30',color=colors[3])
plt.plot(force2_fp[1:],f_40_scaled[1:],'--',label='rod=40',color=colors[4])
plt.plot(force2_fp[1:],f_60_scaled[1:],'--',label='rod=60',color=colors[5])


factor=0.007
x=np.array([19531250.0,
 312500000.0,
 1582031250.0,
 5000000000.0,
 25312500000.0,
 80000000000.0,
 405000000000.0,
 1280000000000.0,
 3125000000000.0])
#x=np.array([0,0.1*10**12,0.2*10**12,0.25*10**12,0.5*10**12,0.75*10**12,1.0*10**12,1.25*10**12,1.5*10**12,1.75*10**12,2.0*10**12,2.25*10**12,2.5*10**12,2.75*10**12,3.0*10**12])
y=np.zeros(x.shape)
for i in range(len(x)):
    y[i]=x[i]**(1/3.0)/2.0*factor
plt.plot(x,y,'k.')
'''
