#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:18:47 2015

@author: jaeger

Modified on Thu Sep 10

@author: isele-holder

Script to do the analysis of several simulations (rod length fixed but different forces analysed)
Analysed quantities: Angle of rod, velocity of filament, relation of rod to end-to-end vector of tail, curvature (PCA)
"""

import matplotlib as mpl
mpl.use('Agg')
from curvature import *
from velocity import *
from pma import *
from data_from_xyz import *
from angle import *
from end_to_end import *
import sys, os

try:
    xyzfile = sys.argv[1]                    # file with coordinates
    folder = sys.argv[2]                     # folder to store the results
    logfilename = sys.argv[3]                # log file
    head = int(sys.argv[4])                  # number of head atoms
    from_atom = int(sys.argv[5])             # starting atom for analysis
    to_atom = int(sys.argv[6])               # last atom for analysis
    dT = float(sys.argv[7])                  # time distance between snapshots in simulation time units
    
except:
    print 'Usage: ' + sys.argv[0] + '      xyz-file     output folder       log file      #(head atoms)       lower atom ID         uppper atom ID       dT'
    exit()

###############################################################################

def analysis(xyzfile,folder,logfile,rod,deltaT,a_s,a_e,t_s=0,t_e=0): #set t_e=0 for all data points to be read
    
    #Initialyse data
    print 'reading the coordinates'
    data=read_coordinates(xyzfile,t_e)
    xs=data[0]
    ys=data[1]
    sets=data[2]  
    n=data[3]   
    if t_e==0:
        t_e=sets
    elif t_e>sets:
        t_e=sets
        
    #Possibly cut the data set for analysis
    print 'cutting the data'
    newxs=cut(xs,t_s,t_e,a_s,a_e)
    newys=cut(ys,t_s,t_e,a_s,a_e)
    sets=t_e-t_s
    n=a_e-a_s+1
    c=curvature_coor(sets,n,newxs,newys)   
    
    #Plot the angle the rod makes with postitive x axis
    print 'plot rod angle with x axis'
    angle,f_angle=angle_of_rod(newxs,newys,sets,10,folder)
    
    #Plot of velocity in x and y-direction and total speed
    print 'plot velocity'
    v_av,vx_sm,vy_sm=plot_velocity(newxs,newys,sets,deltaT,folder)
    logfile.write('speed ' + str(v_av)+ '\n')
    
    #Analyse the angle between the rod and the end-to-end vector of the tail
    print 'analyse end-2-end results'
    end_to_end(xs,ys,sets,n,rod,folder)
   
    #Plot the curvature: arclength against time
    #print 'plot the curvature'
    #plot_curvature(c,folder)   
   
    #Principle component analysis
    print 'perform the principal component analysis' 
    cu_eig,freq, k, eigenvalues, eigenvalues_noise=plot_pma(c,angle,folder,rod)
    np.save(folder+'eigenvalues.npy',eigenvalues)
    np.save(folder+'eigenvalues_noise.npy',eigenvalues_noise)
    logfile.write('PMA-2-Eigenvalues ' + str(cu_eig) + '\n')
    logfile.write('curvature extrema ' +str(k)+ '\n' + '\n')
    return

###############################################################################
 
def main():
    """ analyse principal components"""
    logfile=open(logfilename,'w') #The data for all different forces is written in this one file
    os.system('mkdir ' + folder)
    analysis(xyzfile, folder, logfile, head, dT, from_atom, to_atom)
    logfile.close()
    return

###############################################################################

if __name__ == '__main__':
    main()
