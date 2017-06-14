# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 14:27:13 2015

@author: jaeger
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def read_analysis(name,loc):
    msd_ratio=[]
    speed=[]
    beating=[]
    frequency=[]

    f= open(loc + name,'r')
    t=f.readlines()
    for i in t:
        text=i.split()
        if text==[]:
            text=[0]
        #if text[0]=='name':
        if text[0]=='MSD-ratio':
            msd_ratio.append(float(text[1]))
        if text[0]=='speed':
            speed.append(float(text[1]))
        if text[0]=='PMA-2-Eigenvalues' or text[0]=='beating':
            beating.append(float(text[1]))
        if text[0]=='frequency':
            frequency.append(float(text[1]))
    f.close()
    return msd_ratio, speed, beating, frequency

l=100    
loc='/work/iff_th2/jaeger/simulations/noise_1/tail_100/'
name='analysis_new.txt'

force_=[0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0]
force=[i*l**3/200.0 for i in force_]
force2_fp=[i**4/200.0 for i in force_]
    
location=loc+'head_5/kappa_200/'
msd_5, v_5, beating_5, f_5 = read_analysis(name, location)
location=loc+'head_10/kappa_200/'
msd_10, v_10, beating_10, f_10 = read_analysis(name, location)
location=loc+'head_15/kappa_200/'
msd_15, v_15, beating_15, f_15 = read_analysis(name, location)
location=loc+'head_20/kappa_200/'
msd_20, v_20, beating_20, f_20 = read_analysis(name, location)
location=loc+'head_30/kappa_200/'
msd_30, v_30, beating_30, f_30 = read_analysis(name, location)
location=loc+'head_40/kappa_200/'
msd_40, v_40, beating_40, f_40 = read_analysis(name, location)
location=loc+'head_60/kappa_200/'
msd_60, v_60, beating_60, f_60 = read_analysis(name, location)



colormap=matplotlib.cm.brg
f=[]
r=[]
c=[]
for i in range(len(force)):
    ratio=5/100.0
    b=beating_5[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)

        
    ratio=10/100.0
    b=beating_10[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=15/100.0
    b=beating_15[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=20/100.0
    b=beating_20[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=30/100.0
    b=beating_30[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=40/100.0
    b=beating_40[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=60/100.0
    b=beating_60[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
   
   
f.append(10000)
r.append(1)
c.append(1)
f.append(11000)
r.append(1)
c.append(0)

fig= plt.figure()
ax= fig.add_subplot(1,1,1)
p=ax.scatter(f,r,s=40, c=c, cmap=colormap, edgecolors='None')
        
x=np.delete(range(6000),np.s_[:200], 0)
y=np.zeros(x.shape)
for i in range(len(x)):
    y[i]=1/float((float(x[i])/100.0)-1)
    

ax.plot(x,y,'k')
ax.set_xlabel(r'$f_p l^3/ \kappa$', fontsize=18)
ax.set_ylabel(r'$s/l$', fontsize=18)
ax.set_xlim([0,6000])
ax.set_ylim([0,0.7])
p=fig.colorbar(p, orientation = 'vertical')
p.set_label(r'$\frac{ev(1)+ ev(2)}{total \ sum}$',fontsize=18)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('Phase diagram (vary $f_p$)',fontsize=18)
plt.savefig(loc+'vary_rod_fp_beating.png')
plt.show()
    
#plt.plot(f,c,'-o')
#plt.show()