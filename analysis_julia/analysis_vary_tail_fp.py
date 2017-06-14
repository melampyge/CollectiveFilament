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
loc='/work/iff_th2/jaeger/simulations/noise_p001/head_10/'
name='analysis.txt'

force=[0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0]
force=[i*l**3/200.0 for i in force]
force2=[i**4/200.0 for i in force]
    
location=loc+'tail_200/kappa_200/'
msd_5, v_5, beating_5, f_5 = read_analysis(name, location)
location=loc+'tail_100/kappa_200/'
msd_10, v_10, beating_10, f_10 = read_analysis(name, location)
location=loc+'tail_67/kappa_200/'
msd_15, v_15, beating_15, f_15 = read_analysis(name, location)
location=loc+'tail_50/kappa_200/'
msd_20, v_20, beating_20, f_20 = read_analysis(name, location)
location=loc+'tail_33/kappa_200/'
msd_30, v_30, beating_30, f_30 = read_analysis(name, location)
location=loc+'tail_25/kappa_200/'
msd_40, v_40, beating_40, f_40 = read_analysis(name, location)
location=loc+'tail_17/kappa_200/'
msd_60, v_60, beating_60, f_60 = read_analysis(name, location)

factor=5000
v_5_scaled=[v_5[i]/float(force[i]*200/2.0/210)*factor for i in range(len(force))]
v_10_scaled=[v_10[i]/float(force[i]*100/2.0/110)*factor for i in range(len(force))]
v_15_scaled=[v_15[i]/float(force[i]*67/2.0/77)*factor for i in range(len(force))]
v_20_scaled=[v_20[i]/float(force[i]*50/2.0/60)*factor for i in range(len(force))]
v_30_scaled=[v_30[i]/float(force[i]*33/2.0/43)*factor for i in range(len(force))]
v_40_scaled=[v_40[i]/float(force[i]*25/2.0/35)*factor for i in range(len(force))]
v_60_scaled=[v_60[i]/float(force[i]*17/2.0/27)*factor for i in range(len(force))]

plt.plot(force,v_5_scaled,label='tail=200')
plt.plot(force,v_10_scaled,label='tail=100')
plt.plot(force,v_15_scaled,label='tail=67')
plt.plot(force,v_20_scaled,label='tail=50')
plt.plot(force,v_30_scaled,label='tail=33')
plt.plot(force,v_40_scaled,label='tail=25')
plt.plot(force,v_60_scaled,label='tail=17')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.title('Comparison of swimming speed for fixed kappa')
plt.xlabel(r'$\frac{f_p l^3}{\kappa}$',fontsize=18)
plt.ylabel(r'$v/v_{ref}$')
plt.savefig(loc+'vary_rod_fp_speed_scaled.png')
plt.show()

plt.plot(force,v_5,label='tail=200')
plt.plot(force,v_10,label='tail=100')
plt.plot(force,v_15,label='tail=67')
plt.plot(force,v_20,label='tail=50')
plt.plot(force,v_30,label='tail=33')
plt.plot(force,v_40,label='tail=25')
plt.plot(force,v_60,label='tail=17')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.title('Comparison of swimming speed for fixed kappa')
plt.xlabel(r'$\frac{f_p l^3}{\kappa}$',fontsize=18)
plt.ylabel('$v$')
plt.savefig(loc+'vary_rod_fp_speed.png')
plt.show()

tau=2.0*100**3/200.0
f_5_scaled=[f_5[i]*tau for i in range(len(force))]
f_10_scaled=[f_10[i]*tau for i in range(len(force))]
f_15_scaled=[f_15[i]*tau for i in range(len(force))]
f_20_scaled=[f_20[i]*tau for i in range(len(force))]
f_30_scaled=[f_30[i]*tau for i in range(len(force))]
f_40_scaled=[f_40[i]*tau for i in range(len(force))]
f_60_scaled=[f_60[i]*tau for i in range(len(force))]

plt.plot(force,f_5_scaled,label='tail=200')
plt.plot(force,f_10_scaled,label='tail=100')
plt.plot(force,f_15_scaled,label='tail=67')
plt.plot(force,f_20_scaled,label='tail=50')
plt.plot(force,f_30_scaled,label='tail=33')
plt.plot(force,f_40_scaled,label='tail=25')
plt.plot(force,f_60_scaled,label='tail=17')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.title('Frequency of beating for fixed kappa')
plt.xlabel(r'$\frac{f_p l^3}{\kappa}$',fontsize=18)
plt.ylim(0,90)
plt.ylabel(r'$f \ \tau$')
plt.savefig(loc+'vary_rod_fp_frequency_scaled.png')
plt.show()

plt.plot(force,f_5,label='tail=200')
plt.plot(force,f_10,label='tail=100')
plt.plot(force,f_15,label='tail=67')
plt.plot(force,f_20,label='tail=50')
plt.plot(force,f_30,label='tail=33')
plt.plot(force,f_40,label='tail=25')
plt.plot(force,f_60,label='tail=17')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.title('Frequency of beating for fixed kappa')
plt.xlabel(r'$\frac{f_p l^3}{\kappa}$',fontsize=18)
plt.ylim(0,0.01)
plt.ylabel('$f$')
plt.savefig(loc+'vary_rod_fp_frequency.png')
plt.show()

'''plt.plot(force,msd_5,label='rod=5')
plt.plot(force,msd_10,label='rod=10')
plt.plot(force,msd_15,label='rod=15')
plt.plot(force,msd_20,label='rod=20')
plt.plot(force,msd_30,label='rod=30')
plt.plot(force,msd_40,label='rod=40')
plt.plot(force,msd_60,label='rod=60')
plt.legend(loc=2)
plt.title('MSDy/MSDx')
plt.xlabel('f_p*l^3/kappa')
plt.ylabel('ratio')
plt.savefig('noice_vary_rod_fp_msd.png')
plt.show()'''


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
    
plt.scatter(f,r,s=40, c=c, cmap=colormap, edgecolors='None')
        
x=np.delete(range(6000),np.s_[:200], 0)
y=np.zeros(x.shape)
for i in range(len(x)):
    y[i]=1/float((float(x[i])/100.0)-1)
plt.plot(x,y,'k')
plt.xlim(0,6000)
plt.ylim(0,0.7)
plt.colorbar()
plt.xlabel(r'$\frac{f_p l^3}{\kappa}$',fontsize=18)
plt.ylabel(r'$s/l$',fontsize=18)
plt.title('Beating determined by eigenvalues',fontsize=14)
plt.savefig(loc+'vary_rod_ka_beating.png')
plt.show()
    

