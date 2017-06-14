# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 11:03:29 2015

@author: jaeger
"""
import numpy as np
import matplotlib.pyplot as plt

colors=plt.cm.coolwarm(np.linspace(0,1,6))

plt.plot(force[7:],v_10[7:],'--',label='rod=10',color=colors[0])
plt.plot(force[4:],v_15[4:],'--',label='rod=15',color=colors[1])
plt.plot(force[4:],v_20[4:],'--',label='rod=20',color=colors[2])
plt.plot(force[1:],v_30[1:],'--',label='rod=30',color=colors[3])
plt.plot(force[1:],v_40[1:],'--',label='rod=40',color=colors[4])
plt.plot(force[1:],v_60[1:],'--',label='rod=60',color=colors[5])
plt.plot(force[1:],v_60_05[1:],color=colors[5])
plt.plot(force[1:],v_40_05[1:],color=colors[4])
plt.plot(force[1:],v_30_05[1:],color=colors[3])
plt.plot(force[4:],v_20_05[4:],color=colors[2])
plt.plot(force[4:],v_15_05[4:],color=colors[1])
plt.plot(force[7:],v_10_05[7:],color=colors[0])

plt.legend(loc=1)
plt.xlim(0,7500)
plt.title('Swimming speed')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$v$')
plt.savefig(loc+'vary_rod_speed.png')
plt.show()




plt.plot(force[7:],v_10_scaled[7:],'--',label='rod=10',color=colors[0])
plt.plot(force[4:],v_15_scaled[4:],'--',label='rod=15',color=colors[1])
plt.plot(force[4:],v_20_scaled[4:],'--',label='rod=20',color=colors[2])
plt.plot(force[1:],v_30_scaled[1:],'--',label='rod=30',color=colors[3])
plt.plot(force[1:],v_40_scaled[1:],'--',label='rod=40',color=colors[4])
plt.plot(force[1:],v_60_scaled[1:],'--',label='rod=60',color=colors[5])
plt.plot(force[1:],v_60_05_scaled[1:],color=colors[5])
plt.plot(force[1:],v_40_05_scaled[1:],color=colors[4])
plt.plot(force[1:],v_30_05_scaled[1:],color=colors[3])
plt.plot(force[4:],v_20_05_scaled[4:],color=colors[2])
plt.plot(force[4:],v_15_05_scaled[4:],color=colors[1])
plt.plot(force[7:],v_10_05_scaled[7:],color=colors[0])

plt.legend(loc=1)
plt.xlim(0,7500)
plt.title('Swimming speed')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$v/v_{ref}$',fontsize=18)
plt.savefig(loc+'vary_rod_speed_scaled_beating.png')
plt.show()

fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.plot(force,v_10_scaled,'--',label='s/l=0.1',color=colors[0])
ax.plot(force,v_15_scaled,'--',label='s/l=0.15',color=colors[1])
ax.plot(force,v_20_scaled,'--',label='s/l=0.20',color=colors[2])
ax.plot(force,v_30_scaled,'--',label='s/l=0.30',color=colors[3])
ax.plot(force,v_40_scaled,'--',label='s/l=0.40',color=colors[4])
ax.plot(force,v_60_scaled,'--',label='s/l=0.60',color=colors[5])
ax.plot(force,v_60_05_scaled,color=colors[5])
ax.plot(force,v_40_05_scaled,color=colors[4])
ax.plot(force,v_30_05_scaled,color=colors[3])
ax.plot(force,v_20_05_scaled,color=colors[2])
ax.plot(force,v_15_05_scaled,color=colors[1])
ax.plot(force,v_10_05_scaled,color=colors[0])


ax.set_xlim([0,7500])
fig.suptitle('Swimming speed (dashed line=$f_p$ is changed)',fontsize=16)
ax.set_xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
ax.set_ylabel(r'$v/v_{ref}$',fontsize=18)
plt.legend(loc=1)
fig.tight_layout(renderer=None, pad=2.2, h_pad=None, w_pad=None, rect=None)
plt.savefig(loc+'vary_rod_speed_scaled.png')
plt.show()
