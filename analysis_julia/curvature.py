# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:16:56 2015

@author: jaeger

Guglielmo's code is used to do the principle mode analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
#from data_from_xyz import *
from pma import *
import scipy.ndimage as ndimage
import matplotlib.ticker
from scipy import linalg


#Curvature from 3 atoms
    
def local_curvature(r_0,r_1,r_2):     
    #the curvature is computed about r_1, z-components are expected to be zero
    #curvature is calculated from the area of the signed triangle (given by the three atoms)
    #curvature is positive if filament bends to the left
    ax=r_0[0]-r_1[0]
    ay=r_0[1]-r_1[1]
    bx=r_2[0]-r_1[0]
    by=r_2[1]-r_1[1]
    cx=r_2[0]-r_0[0]
    cy=r_2[1]-r_0[1]
    #a, b and c are the 3 sides of the triangle given by the atoms
    a=np.sqrt(ax**2+ay**2)
    b=np.sqrt(bx**2+by**2)
    c=np.sqrt(cx**2+cy**2)
    A=0.5*(ax*by-ay*bx)
    return float(4.0*A/float(a*b*c))
    
    
    
#Find curvature for whole simulation: curv[time][atom]

def curvature_coor(timesteps,atoms,xs,ys):
    curv=np.zeros((timesteps,atoms-2))
    for j in range(timesteps):
        for i in range(atoms-2):
            curv[j][i]=local_curvature([xs[j][i],ys[j][i],0],[xs[j][i+1],ys[j][i+1],0],[xs[j][i+2],ys[j][i+2],0])
    return curv
    
    
    
#Retruns the section of the data between time_start and time_end and only for atoms between atom_start and atom_end
    
def cut(c,time_start,time_end,atom_start,atom_end):
    return np.delete(np.delete(np.delete(np.delete(c, np.s_[time_end+1:], 0), np.s_[:time_start], 0), np.s_[atom_end+1:], 1), np.s_[:atom_start], 1)
    
    
    
#Display curvaure by a colorbar on a graph that shows position against time    
    
def plot_curvature(c,loc):
    #Fix extreme values so that colorbar is the same for all simulations
    c[0][0]=-0.08
    c[0][1]=0.08

    nx,ny=np.shape(c)
    fig= plt.figure(figsize=(6,3))
    ax= fig.add_subplot(1,1,1)
    cs=ax.pcolor(np.transpose(c))
    cb=fig.colorbar(cs, orientation = 'horizontal')
    
    #tick_locs   = [-0.08,-0.04,0,0.04,0.08]
    #cb.locator     = matplotlib.ticker.FixedLocator(tick_locs)
    #cb.update_ticks()
    
    cb.set_label('curvature',fontsize=16)
    #ax.xaxis.tick_top()
    ax.set_xlim([0,nx])
    ax.set_ylim([0,ny])
    ax.set_xlabel('time', fontsize=16)
    ax.set_ylabel('Arclength', fontsize=18)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    fig.suptitle('Curvature',fontsize=16)
    fig.tight_layout(renderer=None, pad=2.2, h_pad=None, w_pad=None, rect=None)
    plt.savefig(loc+'curvature.png',dpi=300) 
    #plt.show()
    
    
    
#The principal mode analysis is done and various quantities extracted and plotted    
#coordinates= 2D array of curvature coordinates with time on the first axis
#angle= angle of the rod to the neative x axis at each timestep
    
def plot_pma(coordinates,angle,loc,rod,CM=0):
  #Solve for the eigensystem
    nt,ns=np.shape(coordinates)
    eigensystem=pma(coordinates)
    eigenvalues=eigensystem[0]
    
  #Generate noise
    Cvar=coordinates.var()**.5
    CNH=np.random.randn(coordinates.shape[0],coordinates.shape[1])*Cvar
    eigensystem_noise=pma(CNH)
    eigenvalues_noise=eigensystem_noise[0]
    
  #Take noise away from data
    norm_eig=eigenvalues-eigenvalues_noise
    cu_norm_eig=norm_eig.cumsum()/eigenvalues.sum()
    
  #Plot cumulative sum of eigenvalues
    cu_eigenvalues=eigenvalues.cumsum()/eigenvalues.sum()
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.plot(cu_eigenvalues,'-o', markersize=5, label='Pure data')
    ax.plot(cu_norm_eig,'-o', markersize=5, label='Noise removed')
    ax.set_xlim([-5,60])
    ax.set_xlabel('Eigenvalues in decreasing order', fontsize=14)
    ax.set_ylabel('Normalized cumulative sum', fontsize=14)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.legend()
    fig.suptitle('Cumulative sum of eigenvalues',fontsize=16)
    plt.savefig(loc+'cusum_eigenvalues.png')
    #plt.show()
    
  #Plot eigenvalues
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.plot(eigenvalues,'-o', markersize=3, label='Pure data')
    ax.plot(norm_eig,'-o', markersize=5, label='Noise removed')
    ax.plot(eigenvalues_noise,'-o', markersize=3, label='Noise')
    ax.set_xlim([-5,60])
    ax.set_xlabel('Eigenvalues in decreasing order', fontsize=14)
    ax.set_ylabel('Values', fontsize=14)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    fig.suptitle('Eigenvalues',fontsize=16)
    plt.legend()
    plt.savefig(loc+'eigenvalues.png')
    #plt.show()
    
  #Plot the two dominant amplitudes against each other
    eigenvectors=eigensystem[1]
    w=eigenvectors[:,:2]
    amplitudes=get_XY_fast(coordinates,w)
    """
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.plot(amplitudes[0],amplitudes[1])
    ax.set_xlabel('Amplitude 1', fontsize=14)
    ax.set_ylabel('Amplitude 2', fontsize=14)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    fig.suptitle('Relation of 1st and 2nd amplitude',fontsize=16)
    ax.set_aspect("equal")
    fig.tight_layout(renderer=None, pad=2.2, h_pad=None, w_pad=None, rect=None)
    plt.savefig(loc+'amplitudes.png')
    """
    #plt.show()
    
  #Plot the amplitudes of the first 2 modes over time
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.plot(amplitudes[0], label='Mode 1',linewidth=2)
    ax.plot(amplitudes[1], label='Mode 2',linewidth=2)    
    ax.set_xlabel('time', fontsize=14)
    ax.set_ylabel('amplitude', fontsize=14)
    ax.set_ylim([-0.8,0.8])
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.legend()
    fig.suptitle('Development of amplitudes over time',fontsize=16)
    fig.tight_layout()
    plt.savefig(loc+'amp.png')
    #plt.show()
    
  #Find frequency from fourier transform of the 2 dominant amplitudes
    ft0=np.fft.fft(amplitudes[0])
    ft1=np.fft.fft(amplitudes[1])
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.plot(np.abs(ft0), label='1')  
    ax.plot(np.abs(ft1), label='2')      
    ax.set_xlabel('wavenumber', fontsize=14)
    ax.set_ylabel('fourier transform', fontsize=14)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xlim([0,200])
    plt.legend()
    fig.suptitle('Fourier Transforms of amplitudes -> frequency',fontsize=16)
    plt.savefig(loc+'ft_amplitude.png')
    #plt.show()
    k=float(np.argmax(np.abs(ft0)[:nt/2.0]))
    f=k/float(nt)
    if f==0:
        T=nt
        print 'Period not detected'
    else:
        T=1/float(f)
        print T, f


  #Reconstruct shape of principal components of swimmer 
    xs1, ys1=reconstruct_pricipal_modes(eigenvectors[:,0],amplitudes[0])
    xs2, ys2=reconstruct_pricipal_modes(eigenvectors[:,1],amplitudes[1])
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    ax.plot(xs1,ys1,linewidth=2)
    ax.plot(xs2,ys2,linewidth=2)         
    ax.set_xlabel(r'$x$', fontsize=18)
    ax.set_ylabel(r'$y$', fontsize=18)
    #ax.set_xlim([-10,140])
    #ax.set_ylim([-70,70])
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    fig.suptitle('Shape of first two eigenmodes',fontsize=16)
    ax.set_aspect('equal')
    plt.savefig(loc+'eigenmodes.png',dpi=200)
    #plt.show()
    
  #Show development of motion over one cycle
    fig= plt.figure()
    ax= fig.add_subplot(1,1,1)
    instances=22
    dt=T/float(instances)
    colors=plt.cm.gist_heat(np.linspace(0,0.8,instances))
    start_time=find_start(amplitudes[0],T,20)
    for i in range(instances):
        xs, ys=reconstruct_shape(eigenvectors[:,0],eigenvectors[:,1],amplitudes[0],amplitudes[1],int(start_time+i*dt),angle,rod)
        ax.plot(xs,ys,color=colors[i],linewidth=2, label=str(i+1))
    ax.set_xlabel(r'$x$', fontsize=18)
    ax.set_ylabel(r'$y$', fontsize=18)
    #ax.set_xlim([-10,140])
    #ax.set_ylim([-70,70])
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    fig.suptitle('Development of shape over one period',fontsize=16)
    ax.set_aspect('equal')
    plt.savefig(loc+'shape.png')
    #plt.show()        
    
  #Plot curvature in new coordinates
    if CM==1:
        cm=eigensystem[2]
        cs=plt.pcolor(np.transpose(cm))
        nx,ny=np.shape(cm)
        cb= plt.colorbar(cs, orientation = 'horizontal')
        cb.set_label('correlation')
        plt.xlim(0,nx)
        plt.ylim(0,ny)
        plt.xlabel('time')
        plt.ylabel('position')
        plt.show()
        
    GF=lambda x:ndimage.filters.gaussian_filter(x,50)        
        
  #Find the number of extrema along the tail (need to adjust the Gaussian filter so that the noise is not counted!!!!)
    start=0
    end=nt
    ks=[]
    for j in range(end-start):
        a=GF(coordinates[start+j])[rod:ns-10]
        maxima=np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True] 
        minima=np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
        extrema= maxima + minima
        k=-2
        for i in extrema:
            if i==True:
                k+=1
        ks.append(k)
    k=np.average(ks)    
    print 'Extrema: ', k    
        
        
  #The normalized sum of the first two eienvalues and the frequency are returned
    return cu_eigenvalues[1], f, k, eigenvalues, eigenvalues_noise

#Convert the principal modes which are given in curvature coordinates into the acutal shape in x and y

def reconstruct_pricipal_modes(ev,amp):
    angle=[]
    N=ev.shape[0]
    xs=np.zeros(N)
    ys=np.zeros(N)
    angle=amp.max()*np.cumsum(ev)
    C=np.cos(angle)
    S=np.sin(angle)
    xs=np.cumsum(C)
    ys=np.cumsum(S)
    return xs, ys


#Reconstruct the shape of the filament from principal modes at a given time. 
#Curvature coordinates are converted to angles and then positions
#the angle is used to reintroduce the orientation of the rod

def reconstruct_shape(ev1,ev2,amp1,amp2,time,angle,rod):
    phi=[]
    N=ev1.shape[0]
    xs=np.zeros(N)
    ys=np.zeros(N)
    phi=amp1[time]*np.cumsum(ev1)+amp2[time]*np.cumsum(ev2)
    C=np.cos(phi)
    S=np.sin(phi)
    xs=np.cumsum(C)
    ys=np.cumsum(S)
    # shift by the center of mass
    comx = np.average(xs)
    comy = np.average(ys)
    xs -= comx
    ys -= comy
    # compute the prinicipal moments
    I = np.zeros((2,2))
    I[0,0] = np.sum(ys**2)
    I[1,1] = np.sum(xs**2)
    I[0,1] = np.sum(xs*ys)
    I[1,0] = I[0,1]
    # find the largest eigenvalue and eigenvector of the moment of intertia tensor
    w, v = linalg.eig(I)
    idmax = np.argmax(w)
    wmax = w[idmax]
    vmax = v[idmax]
    # rotate xs and ys such that vmax points in the x direction
    phi = np.arctan2(vmax[0],vmax[1])
    cos = np.cos(phi)
    sin = np.sin(phi)
    xnew = xs*cos - ys*sin
    ynew = xs*sin + ys*cos
    # mirror image, if required
    if xnew[0] > 0:
        xnew *= -1
        ynew *= -1
    
    '''
  #Reintroduce orientation of rod
    th=-angle[time] #Turn clockwise for positive angle
    rotxs=np.cos(th)*xs-np.sin(th)*ys
    rotys=np.sin(th)*xs+np.cos(th)*ys    
  #Fix at center of mass of rod
    x=np.sum(rotxs[0:rod])/float(rod)
    y=np.sum(rotys[0:rod])/float(rod)
    new_xs=[rotxs[j]-x for j in range(len(rotxs))]
    new_ys=[rotys[j]-y for j in range(len(rotys))] 
    ''' 
    return xnew, ynew
    
    
    
    
#For plotting the motion over one period, we don't want to the start or the end of the simulation but a period in the middle    
    
def find_start(a,T,beats):
    maxima=np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
    i=0
    for j in range(beats):
        if i>len(a)-1-T:
            i=len(a)-2*T
            break
        while maxima[i]==False:
            i+=1
        i+=1
    return i



#Smooth an array

def smooth1DArray(list,degree=10):  

     smoothed=np.zeros(np.shape(list))  

     for i in range(len(smoothed)):  
         if i<degree/2:
             smoothed[i]=sum(list[i:i+degree])/float(degree)          
         elif i>=len(list)-degree/2:
             smoothed[i]=sum(list[i-degree:i])/float(degree)  
         else:
             smoothed[i]=sum(list[i-degree/2:i+degree/2])/float(degree)  

     return smoothed 







        
