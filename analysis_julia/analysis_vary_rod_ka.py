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
        if text[0]=='PMA-2-Eigenvalues':
            beating.append(float(text[1]))
        if text[0]=='frequency':
            frequency.append(float(text[1]))
    f.close()
    return msd_ratio, speed, beating, frequency
    
r0=1
l=100*r0
kappa=[2000,1000,667,500,333,250,167,125,100]
force=[0.5*l**3/float(i) for i in kappa]
force2_ka=[0.5**4/float(i) for i in kappa]

loc='/work/iff_th2/jaeger/simulations/noise_p001/tail_100/'
name='analysis.txt'    

   
location=loc+'head_5/fp_p5/'
msd_5_05, v_5_05, beating_5_05, f_5_05 = read_analysis(name, location)
location=loc+'head_10/fp_p5/'
msd_10_05, v_10_05, beating_10_05, f_10_05 = read_analysis(name, location)
location=loc+'head_15/fp_p5/'
msd_15_05, v_15_05, beating_15_05, f_15_05 = read_analysis(name, location)
location=loc+'head_20/fp_p5/'
msd_20_05, v_20_05, beating_20_05, f_20_05 = read_analysis(name, location)
location=loc+'head_30/fp_p5/'
msd_30_05, v_30_05, beating_30_05, f_30_05 = read_analysis(name, location)
location=loc+'head_40/fp_p5/'
msd_40_05, v_40_05, beating_40_05, f_40_05 = read_analysis(name, location)
location=loc+'head_60/fp_p5/'
msd_60_05, v_60_05, beating_60_05, f_60_05 = read_analysis(name, location)

v_5_05_scaled=[i/(0.5*100/2.0/105) for i in v_5_05]
v_10_05_scaled=[i/(0.5*100/2.0/110) for i in v_10_05]
v_15_05_scaled=[i/(0.5*100/2.0/115) for i in v_15_05]
v_20_05_scaled=[i/(0.5*100/2.0/120) for i in v_20_05]
v_30_05_scaled=[i/(0.5*100/2.0/130) for i in v_30_05]
v_40_05_scaled=[i/(0.5*100/2.0/140) for i in v_40_05]
v_60_05_scaled=[i/(0.5*100/2.0/160) for i in v_60_05]

plt.plot(force,v_5_05_scaled,label='rod=5')
plt.plot(force,v_10_05_scaled,label='rod=10')
plt.plot(force,v_15_05_scaled,label='rod=15')
plt.plot(force,v_20_05_scaled,label='rod=20')
plt.plot(force,v_30_05_scaled,label='rod=30')
plt.plot(force,v_40_05_scaled,label='rod=40')
plt.plot(force,v_60_05_scaled,label='rod=60')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.title('Comparison of swimming speed for fixed force')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$v/v_{ref}$',fontsize=18)
plt.savefig(loc+'vary_rod_ka_speed_scaled.png')
plt.show()

plt.plot(force,v_5_05,label='rod=5')
plt.plot(force,v_10_05,label='rod=10')
plt.plot(force,v_15_05,label='rod=15')
plt.plot(force,v_20_05,label='rod=20')
plt.plot(force,v_30_05,label='rod=30')
plt.plot(force,v_40_05,label='rod=40')
plt.plot(force,v_60_05,label='rod=60')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.title('Comparison of swimming speed for fixed force')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$v$',fontsize=18)
plt.savefig(loc+'vary_rod_ka_speed.png')
plt.show()

time=2.0/0.5
f_5_05_scaled=[f_5_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]
f_10_05_scaled=[f_10_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]
f_15_05_scaled=[f_15_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]
f_20_05_scaled=[f_20_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]
f_30_05_scaled=[f_30_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]
f_40_05_scaled=[f_40_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]
f_60_05_scaled=[f_60_05[i]*2.0*100**3/float(kappa[i]) for i in range(len(force))]

plt.plot(force,f_5_05_scaled,label='rod=5')
plt.plot(force,f_10_05_scaled,label='rod=10')
plt.plot(force,f_15_05_scaled,label='rod=15')
plt.plot(force,f_20_05_scaled,label='rod=20')
plt.plot(force,f_30_05_scaled,label='rod=30')
plt.plot(force,f_40_05_scaled,label='rod=40')
plt.plot(force,f_60_05_scaled,label='rod=60')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.ylim(0,60)
plt.title('Frequency of beating for fixed force')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$f \ \tau$',fontsize=18)
plt.savefig(loc+'vary_rod_ka_frequency_scaled.png')
plt.show()


plt.plot(force,f_5_05,label='rod=5')
plt.plot(force,f_10_05,label='rod=10')
plt.plot(force,f_15_05,label='rod=15')
plt.plot(force,f_20_05,label='rod=20')
plt.plot(force,f_30_05,label='rod=30')
plt.plot(force,f_40_05,label='rod=40')
plt.plot(force,f_60_05,label='rod=60')
plt.legend(loc=1)
plt.xlim(0,7000)
plt.ylim(0,0.005)
plt.title('Frequency of beating for fixed force')
plt.xlabel(r'$f_p l^3/ \kappa$',fontsize=18)
plt.ylabel(r'$f$',fontsize=18)
plt.savefig(loc+'vary_rod_ka_frequency.png')
plt.show()

'''plt.plot(force,msd_5,label='rod=5')
plt.plot(force,msd_10,label='rod=10')
plt.plot(force,msd_15,label='rod=15')
plt.plot(force,msd_20,label='rod=20')
plt.plot(force,msd_30,label='rod=30')
plt.plot(force,msd_40,label='rod=40')
plt.plot(force,msd_60,label='rod=60')
plt.legend(loc=2)
plt.xlim(2000,100)
plt.title('MSDy/MSDx')
plt.xlabel('kappa')
plt.ylabel('ratio')
plt.savefig('vary_rod_ka_msd.png')
plt.show()'''


colormap=matplotlib.cm.brg
f=[]
r=[]
c=[]
for i in range(len(force)):
    ratio=5/100.0
    b=beating_5_05[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)

        
    ratio=10/100.0
    b=beating_10_05[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=15/100.0
    b=beating_15_05[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=20/100.0
    b=beating_20_05[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=30/100.0
    b=beating_30_05[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=40/100.0
    b=beating_40_05[i]
    f.append(force[i])
    r.append(ratio)
    c.append(b)
        
    ratio=60/100.0
    b=beating_60_05[i]
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
fig.suptitle('Phase diagram for beating (no noise)',fontsize=18)
plt.savefig(loc+'vary_rod_ka_beating.png')
plt.show()