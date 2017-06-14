# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:18:47 2015

@author: jaeger

Script to analyse the the mean sqared displacement of the filaments in the simulations
The data is saved in msd.npy and can therefore be read in and summarised in plots which is done below it is important that the files are 
at the specified locations.
"""
from msd import *
from data_from_xyz import *
import numpy as np

def analysis_msd(location, run, folder, t_s, t_e, a_s, a_e, rod):
    
    #Initialyse data
    data=read_coordinates(location+run,t_e)
    xs=data[0]
    ys=data[1]
    sets=data[2]  
    #n=data[3]   
    
    msd=plot_msd(xs,ys,sets,location+folder)
    
    np.save(location+folder+'msd.npy',msd)
    
#################################################################################################################################
#The analysis can be used in the following way
'''
location='/work/iff_th2/jaeger/simulations/msd/head_5/kappa_200/' 

rod=5
deltaT=5
a_s=0
a_e=159
t_s=0
t_e=0

  
analysis_msd(location, 'run_p05.xyz', 'p05_msd/', t_s, t_e, a_s, a_e, rod)
print '1 done'
analysis_msd(location, 'run_p1.xyz', 'p1_msd/', t_s, t_e, a_s, a_e, rod)
print '2 done'
analysis_msd(location, 'run_p2.xyz', 'p2_msd/', t_s, t_e, a_s, a_e, rod)
print '4 done'
analysis_msd(location, 'run_p3.xyz', 'p3_msd/', t_s, t_e, a_s, a_e, rod)
print '5 done'
analysis_msd(location, 'run_p4.xyz', 'p4_msd/', t_s, t_e, a_s, a_e, rod)
print '6 done'
analysis_msd(location, 'run_p6.xyz', 'p6_msd/', t_s, t_e, a_s, a_e, rod)
print '7 done'
analysis_msd(location, 'run_p8.xyz', 'p8_msd/', t_s, t_e, a_s, a_e, rod)
print '8 done'
analysis_msd(location, 'run_1.xyz', '1_msd/', t_s, t_e, a_s, a_e, rod)
print '9 done' 

'''
################################################################################################################################
#Once the msd for all simulations has been found and saved by the above function, the code below summarises the data in in plots

tau=float(2*np.power(100,3)/200.0)
sets=60001
colors=plt.cm.gist_heat(np.linspace(0,0.8,8))
t=[1,int(sets/50000.0),int(sets/30000.0),int(sets/10000.0),int(sets/6000.0),int(sets/3000.0),int(sets/1000.0),int(sets/500.0),int(sets/200.0),int(sets/100.0),int(sets/50.0),int(sets/20.0),int(sets*2/10.0),int(sets*3/10.0),int(sets*6/10.0),int(sets*8/10.0)]
ts=[i/tau for i in t]

########################################
#Rod=5

location='/work/iff_th2/jaeger/simulations/msd/head_5/kappa_200/' 

msd_head_5_p05=np.load(location+'p05_msd/'+'msd.npy')
msd_head_5_p1=np.load(location+'p1_msd/'+'msd.npy')
msd_head_5_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_5_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_5_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_5_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_5_p8=np.load(location+'p8_msd/'+'msd.npy')
msd_head_5_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p05[0],yerr=msd_head_5_p05[1],label='F=0.05', color=colors[7])
ax.errorbar(ts,msd_head_5_p1[0],yerr=msd_head_5_p1[1],label='F=0.1', color=colors[6])
ax.errorbar(ts,msd_head_5_p2[0],yerr=msd_head_5_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_5_p3[0],yerr=msd_head_5_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_5_p4[0],yerr=msd_head_5_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_5_p6[0],yerr=msd_head_5_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_5_p8[0],yerr=msd_head_5_p8[1],label='F=0.8', color=colors[1])
ax.errorbar(ts,msd_head_5_1[0],yerr=msd_head_5_1[1],label='Force=1', color=colors[0])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=5',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
#plt.show()
'''


########################################
#Rod=10

location='/work/iff_th2/jaeger/simulations/msd/head_10/kappa_200/' 

msd_head_10_p05=np.load(location+'p05_msd/'+'msd.npy')
msd_head_10_p1=np.load(location+'p1_msd/'+'msd.npy')
msd_head_10_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_10_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_10_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_10_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_10_p8=np.load(location+'p8_msd/'+'msd.npy')
msd_head_10_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_10_p05[0],yerr=msd_head_10_p05[1],label='F=0.05', color=colors[7])
ax.errorbar(ts,msd_head_10_p1[0],yerr=msd_head_10_p1[1],label='F=0.1', color=colors[6])
ax.errorbar(ts,msd_head_10_p2[0],yerr=msd_head_10_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_10_p3[0],yerr=msd_head_10_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_10_p4[0],yerr=msd_head_10_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_10_p6[0],yerr=msd_head_10_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_10_p8[0],yerr=msd_head_10_p8[1],label='F=0.8', color=colors[1])
ax.errorbar(ts,msd_head_10_1[0],yerr=msd_head_10_1[1],label='Force=1', color=colors[0])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=18)
ax.set_ylabel('MSD', fontsize=18)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=10',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
#plt.show()
'''

########################################
#Rod=15

location='/work/iff_th2/jaeger/simulations/msd/head_15/kappa_200/' 


msd_head_15_p05=np.load(location+'p05_msd/'+'msd.npy')
msd_head_15_p1=np.load(location+'p1_msd/'+'msd.npy')
#msd_head_15_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_15_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_15_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_15_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_15_p8=np.load(location+'p8_msd/'+'msd.npy')
msd_head_15_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_15_p05[0],yerr=msd_head_15_p05[1],label='F=0.05', color=colors[7])
ax.errorbar(ts,msd_head_15_p1[0],yerr=msd_head_15_p1[1],label='F=0.1', color=colors[6])
#ax.errorbar(ts,msd_head_15_p2[0],yerr=msd_head_15_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_15_p3[0],yerr=msd_head_15_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_15_p4[0],yerr=msd_head_15_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_15_p6[0],yerr=msd_head_15_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_15_p8[0],yerr=msd_head_15_p8[1],label='F=0.8', color=colors[1])
ax.errorbar(ts,msd_head_15_1[0],yerr=msd_head_15_1[1],label='Force=1', color=colors[0])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=15',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
#plt.show()
'''

########################################
#Rod=20

location='/work/iff_th2/jaeger/simulations/msd/head_20/kappa_200/' 


msd_head_20_p05=np.load(location+'p05_msd/'+'msd.npy')
msd_head_20_p1=np.load(location+'p1_msd/'+'msd.npy')
#msd_head_15_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_20_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_20_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_20_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_20_p8=np.load(location+'p8_msd/'+'msd.npy')
msd_head_20_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_20_p05[0],yerr=msd_head_20_p05[1],label='F=0.05', color=colors[7])
ax.errorbar(ts,msd_head_20_p1[0],yerr=msd_head_20_p1[1],label='F=0.1', color=colors[6])
#ax.errorbar(ts,msd_head_20_p2[0],yerr=msd_head_20_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_20_p3[0],yerr=msd_head_20_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_20_p4[0],yerr=msd_head_20_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_20_p6[0],yerr=msd_head_20_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_20_p8[0],yerr=msd_head_20_p8[1],label='F=0.8', color=colors[1])
ax.errorbar(ts,msd_head_20_1[0],yerr=msd_head_20_1[1],label='Force=1', color=colors[0])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=20',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
#plt.show()
'''


########################################
#Rod=30

location='/work/iff_th2/jaeger/simulations/msd/head_30/kappa_200/' 


msd_head_30_p05=np.load(location+'p05_msd/'+'msd.npy')
msd_head_30_p1=np.load(location+'p1_msd/'+'msd.npy')
#msd_head_30_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_30_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_30_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_30_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_30_p8=np.load(location+'p8_msd/'+'msd.npy')
msd_head_30_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_30_p05[0],yerr=msd_head_30_p05[1],label='F=0.05', color=colors[7])
ax.errorbar(ts,msd_head_30_p1[0],yerr=msd_head_30_p1[1],label='F=0.1', color=colors[6])
#ax.errorbar(ts,msd_head_15_p2[0],yerr=msd_head_30_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_30_p3[0],yerr=msd_head_30_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_30_p4[0],yerr=msd_head_30_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_30_p6[0],yerr=msd_head_30_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_30_p8[0],yerr=msd_head_30_p8[1],label='F=0.8', color=colors[1])
ax.errorbar(ts,msd_head_30_1[0],yerr=msd_head_30_1[1],label='Force=1', color=colors[0])

ax.plot([2,20],[1000,100000], label='slope=2')
ax.plot([np.power(10,2),np.power(10,3)],[2*np.power(10,6),2*np.power(10,7)], label='slope=1')

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=30',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
#plt.show()
'''

########################################
#Rod=40

location='/work/iff_th2/jaeger/simulations/msd/head_40/kappa_200/' 


msd_head_40_p05=np.load(location+'p05_msd/'+'msd.npy')
#msd_head_40_p1=np.load(location+'p1_msd/'+'msd.npy')
msd_head_40_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_40_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_40_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_40_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_40_p8=np.load(location+'p8_msd/'+'msd.npy')
msd_head_40_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_40_p05[0],yerr=msd_head_40_p05[1],label='F=0.05', color=colors[7])
#ax.errorbar(ts,msd_head_40_p1[0],yerr=msd_head_40_p1[1],label='F=0.1', color=colors[6])
ax.errorbar(ts,msd_head_40_p2[0],yerr=msd_head_40_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_40_p3[0],yerr=msd_head_40_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_40_p4[0],yerr=msd_head_40_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_40_p6[0],yerr=msd_head_40_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_40_p8[0],yerr=msd_head_40_p8[1],label='F=0.8', color=colors[1])
ax.errorbar(ts,msd_head_40_1[0],yerr=msd_head_40_1[1],label='Force=1', color=colors[0])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=40',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
#plt.show()
'''

########################################
#Rod=60

location='/work/iff_th2/jaeger/simulations/msd/head_60/kappa_200/' 


msd_head_60_p05=np.load(location+'p05_msd/'+'msd.npy')
msd_head_60_p1=np.load(location+'p1_msd/'+'msd.npy')
msd_head_60_p2=np.load(location+'p2_msd/'+'msd.npy')
msd_head_60_p3=np.load(location+'p3_msd/'+'msd.npy')
msd_head_60_p4=np.load(location+'p4_msd/'+'msd.npy')
msd_head_60_p6=np.load(location+'p6_msd/'+'msd.npy')
msd_head_60_p8=np.load(location+'p8_msd/'+'msd.npy')
#msd_head_60_1=np.load(location+'1_msd/'+'msd.npy')
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_60_p05[0],yerr=msd_head_60_p05[1],label='F=0.05', color=colors[7])
ax.errorbar(ts,msd_head_60_p1[0],yerr=msd_head_60_p1[1],label='F=0.1', color=colors[6])
ax.errorbar(ts,msd_head_60_p2[0],yerr=msd_head_60_p2[1],label='F=0.2', color=colors[5])
ax.errorbar(ts,msd_head_60_p3[0],yerr=msd_head_60_p3[1],label='F=0.3', color=colors[4])
ax.errorbar(ts,msd_head_60_p4[0],yerr=msd_head_60_p4[1],label='F=0.4', color=colors[3])
ax.errorbar(ts,msd_head_60_p6[0],yerr=msd_head_60_p6[1],label='F=0.6', color=colors[2])
ax.errorbar(ts,msd_head_60_p8[0],yerr=msd_head_60_p8[1],label='F=0.8', color=colors[1])
#ax.errorbar(ts,msd_head_60_1[0],yerr=msd_head_60_1[1],label='Force=1', color=colors[0])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Rod=60',fontsize=16)
plt.legend()
plt.savefig(location+'msd.png')
##plt.show()
'''

############################################################################################
#Compare different rod length at fixed force

########################################
#Force=0.05

location='/work/iff_th2/jaeger/simulations/msd/' 
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p05[0],yerr=msd_head_5_p05[1],label='Rod=5', color=colors[7])
ax.errorbar(ts,msd_head_10_p05[0],yerr=msd_head_10_p05[1],label='Rod=10', color=colors[6])
ax.errorbar(ts,msd_head_15_p05[0],yerr=msd_head_15_p05[1],label='Rod=15', color=colors[5])
ax.errorbar(ts,msd_head_20_p05[0],yerr=msd_head_20_p05[1],label='Rod=20', color=colors[4])
ax.errorbar(ts,msd_head_30_p05[0],yerr=msd_head_30_p05[1],label='Rod=30', color=colors[3])
ax.errorbar(ts,msd_head_40_p05[0],yerr=msd_head_40_p05[1],label='Rod=40', color=colors[2])
ax.errorbar(ts,msd_head_60_p05[0],yerr=msd_head_60_p05[1],label='Rod=60', color=colors[1])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Force=0.05',fontsize=16)
plt.legend()
plt.savefig(location+'msd_p05.png')
##plt.show()


########################################
#Force=0.1

fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p1[0],yerr=msd_head_5_p1[1],label='Rod=5', color=colors[7])
ax.errorbar(ts,msd_head_10_p1[0],yerr=msd_head_10_p1[1],label='Rod=10', color=colors[6])
ax.errorbar(ts,msd_head_15_p1[0],yerr=msd_head_15_p1[1],label='Rod=15', color=colors[5])
ax.errorbar(ts,msd_head_20_p1[0],yerr=msd_head_20_p1[1],label='Rod=20', color=colors[4])
ax.errorbar(ts,msd_head_30_p1[0],yerr=msd_head_30_p1[1],label='Rod=30', color=colors[3])
#ax.errorbar(ts,msd_head_40_p1[0],yerr=msd_head_40_p1[1],label='Rod=40', color=colors[2])
ax.errorbar(ts,msd_head_60_p1[0],yerr=msd_head_60_p1[1],label='Rod=60', color=colors[1])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Force=0.1',fontsize=16)
plt.legend()
plt.savefig(location+'msd_p1.png')
##plt.show()


########################################
#Force=0.3

fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p3[0],yerr=msd_head_5_p3[1],label='Rod=5', color=colors[7])
ax.errorbar(ts,msd_head_10_p3[0],yerr=msd_head_10_p3[1],label='Rod=10', color=colors[6])
ax.errorbar(ts,msd_head_15_p3[0],yerr=msd_head_15_p3[1],label='Rod=15', color=colors[5])
ax.errorbar(ts,msd_head_20_p3[0],yerr=msd_head_20_p3[1],label='Rod=20', color=colors[4])
ax.errorbar(ts,msd_head_30_p3[0],yerr=msd_head_30_p3[1],label='Rod=30', color=colors[3])
ax.errorbar(ts,msd_head_40_p3[0],yerr=msd_head_40_p3[1],label='Rod=40', color=colors[2])
ax.errorbar(ts,msd_head_60_p3[0],yerr=msd_head_60_p3[1],label='Rod=60', color=colors[1])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Force=0.3',fontsize=16)
plt.legend()
plt.savefig(location+'msd_p3.png')
##plt.show()


########################################
#Force=0.4

fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p4[0],yerr=msd_head_5_p4[1],label='Rod=5', color=colors[7])
ax.errorbar(ts,msd_head_10_p4[0],yerr=msd_head_10_p4[1],label='Rod=10', color=colors[6])
ax.errorbar(ts,msd_head_15_p4[0],yerr=msd_head_15_p4[1],label='Rod=15', color=colors[5])
ax.errorbar(ts,msd_head_20_p4[0],yerr=msd_head_20_p4[1],label='Rod=20', color=colors[4])
ax.errorbar(ts,msd_head_30_p4[0],yerr=msd_head_30_p4[1],label='Rod=30', color=colors[3])
ax.errorbar(ts,msd_head_40_p4[0],yerr=msd_head_40_p4[1],label='Rod=40', color=colors[2])
ax.errorbar(ts,msd_head_60_p4[0],yerr=msd_head_60_p4[1],label='Rod=60', color=colors[1])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Force=0.4',fontsize=16)
plt.legend()
plt.savefig(location+'msd_p4.png')
##plt.show()
'''

########################################
#Force=0.6

fig= plt.figure(figsize=(12,8))
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p6[0],yerr=msd_head_5_p6[1],label='s/l=0.05',linewidth=2, color=colors[7])
ax.errorbar(ts,msd_head_10_p6[0],yerr=msd_head_10_p6[1],label='s/l=0.1',linewidth=2, color=colors[6])
ax.errorbar(ts,msd_head_15_p6[0],yerr=msd_head_15_p6[1],label='s/l=0.15',linewidth=2, color=colors[5])
ax.errorbar(ts,msd_head_20_p6[0],yerr=msd_head_20_p6[1],label='s/l=0.2',linewidth=2, color=colors[4])
ax.errorbar(ts,msd_head_30_p6[0],yerr=msd_head_30_p6[1],label='s/l=0.3',linewidth=2, color=colors[3])
ax.errorbar(ts,msd_head_40_p6[0],yerr=msd_head_40_p6[1],label='s/l=0.4',linewidth=2, color=colors[2])
ax.errorbar(ts,msd_head_60_p6[0],yerr=msd_head_60_p6[1],label='s/l=0.6',linewidth=2, color=colors[1])

ax.plot([1/tau,10/tau],[1000,100000],linewidth=2, color='k')
ax.plot([np.power(10,3)/tau,np.power(10,4)/tau],[2*np.power(10,5),2*np.power(10,6)],linewidth=2, color='k')
ax.text(2.5/tau, 1.5*np.power(10,4), '2', fontsize=14)
ax.text(3*np.power(10,3)/tau, 2*np.power(10,5), '1', fontsize=14)

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel(r'$t/\tau$', fontsize=24)
ax.set_ylabel('MSD', fontsize=18)
ax.set_xlim([1/tau,np.power(10,5)/tau])
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
fig.suptitle(r'MSD (log-log plot): $f_p l^3/\kappa$=3000',fontsize=28)
plt.legend(loc=4)
fig.tight_layout(renderer=None, pad=4, h_pad=None, w_pad=None, rect=None)
plt.savefig(location+'msd_p6.png')
#plt.show()


########################################
#Force=0.8
'''
fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_p8[0],yerr=msd_head_5_p8[1],label='Rod=5', color=colors[7])
ax.errorbar(ts,msd_head_10_p8[0],yerr=msd_head_10_p8[1],label='Rod=10', color=colors[6])
ax.errorbar(ts,msd_head_15_p8[0],yerr=msd_head_15_p8[1],label='Rod=15', color=colors[5])
ax.errorbar(ts,msd_head_20_p8[0],yerr=msd_head_20_p8[1],label='Rod=20', color=colors[4])
ax.errorbar(ts,msd_head_30_p8[0],yerr=msd_head_30_p8[1],label='Rod=30', color=colors[3])
ax.errorbar(ts,msd_head_40_p8[0],yerr=msd_head_40_p8[1],label='Rod=40', color=colors[2])
ax.errorbar(ts,msd_head_60_p8[0],yerr=msd_head_60_p8[1],label='Rod=60', color=colors[1])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Force=0.8',fontsize=16)
plt.legend()
plt.savefig(location+'msd_p8.png')
#plt.show()


########################################
#Force=1

fig= plt.figure()
ax= fig.add_subplot(1,1,1)

ax.errorbar(ts,msd_head_5_1[0],yerr=msd_head_5_1[1],label='Rod=5', color=colors[7])
ax.errorbar(ts,msd_head_10_1[0],yerr=msd_head_10_1[1],label='Rod=10', color=colors[6])
ax.errorbar(ts,msd_head_15_1[0],yerr=msd_head_15_1[1],label='Rod=15', color=colors[5])
ax.errorbar(ts,msd_head_20_1[0],yerr=msd_head_20_1[1],label='Rod=20', color=colors[4])
ax.errorbar(ts,msd_head_30_1[0],yerr=msd_head_30_1[1],label='Rod=30', color=colors[3])
ax.errorbar(ts,msd_head_40_1[0],yerr=msd_head_40_1[1],label='Rod=40', color=colors[2])
#ax.errorbar(ts,msd_head_60_1[0],yerr=msd_head_60_1[1],label='Rod=60', color=colors[1])

ax.set_xscale('log')  
ax.set_yscale('log')  
ax.set_xlabel('time', fontsize=13)
ax.set_ylabel('MSD', fontsize=13)
ax.set_xlim([1,np.power(10,7)])
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.suptitle('MSD (log-log plot): Force=1',fontsize=16)
plt.legend()
plt.savefig(location+'msd_1.png')
##plt.show()
'''