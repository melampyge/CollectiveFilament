#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

##################################################################

# OrientVel class: compare orientation and velocity

##################################################################

class OrientVel:

    ##############################################################
      
    def __init__(self, nsteps, natoms, line):
        """ initialize: allocate density array"""
        ### set histogram parameters
        line = line.split()
        self.nevery = int(line[1])
        ### allocate arrays to store data
        self.x_store = np.zeros((nsteps/self.nevery, natoms))
        self.y_store = np.zeros((nsteps/self.nevery, natoms))
        self.scprod = np.zeros((nsteps/self.nevery, natoms))
        self.counter = 0
        return

    ##############################################################

    def compute(self,step,x,y,vx,vy,phi,natoms,plot = 'False'):
        """ compute a density distribution and a histogram
            of the density distribution """
        ### check whether quantity needs to be computed
        if step % self.nevery != 0:
            return
        ### compute scalar product of velocity and orientation
        cos = np.cos(phi)
        sin = np.sin(phi)
        ### compute product and increase the counter
        self.scprod[self.counter] = cos*vx + sin*vy
        self.counter = self.counter + 1
        ### generate a plot for testing purposes
        if plot == 'True':
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.scatter(x,y,c = self.scprod[self.counter-1], s = 1.2, linewidth = 0)
            ax.set_aspect('equal')
            plt.show()
            plt.close()
        return

