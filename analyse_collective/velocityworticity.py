#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import density

class VelocityWorticity:

    def __init__(self, nsteps, line):
        """ set some variables and allocate output arrays"""
        ### allocate output arrays and set values
        line = line.split()
        self.nbins = int(line[1])
        self.rho = np.zeros((nsteps, self.nbins, self.nbins))
        self.vx = np.zeros((nsteps, self.nbins, self.nbins))
        self.vy = np.zeros((nsteps, self.nbins, self.nbins))
        self.worticity = np.zeros((nsteps, self.nbins, self.nbins))
        return

    ##############################################################

    def compute(self,step,x,y,vx,vy,natoms,lx,ly,plot = 'False'):
        """ bin the velocities and particles + compute the worticity"""
        ### fill all atoms and velocities into bins
        for i in range(natoms):
            # get coordinates
            xi = x[i]
            yi = y[i]
            # get current bin
            segx = int(xi/lx*self.nbins)
            segy = int(yi/ly*self.nbins)
            # add data to bin
            self.rho[step,segx,segy] += 1
            self.vx[step,segx,segy] += vx[i]
            self.vy[step,segx,segy] += vy[i]
        # normalize velocities
        for i in range(self.nbins):
            for j in range(self.nbins):
                if self.rho[step,i,j] > 1:
                    self.vx[step,i,j] /= self.rho[step,i,j]
                    self.vy[step,i,j] /= self.rho[step,i,j]
        # transform number counts to densities
        wx = lx/self.nbins
        wy = ly/self.nbins
        self.rho[step] /= wx*wy
        ### compute the worticity
        for i in range(self.nbins):
            for j in range(self.nbins):
                # compute velocity gradients using finite differences
                duy_dx = (self.vy[step,(i+1)%self.nbins,j] - self.vy[step,i-1,j])/(2*wx)
                dux_dy = (self.vx[step,i,(j+1)%self.nbins] - self.vx[step,i,j-1])/(2*wy)
                self.worticity[step,i,j] = duy_dx - dux_dy
        ### generate plots for testing purposes
        if plot == 'True':
            fig = plt.figure()
            ax1 = plt.subplot(221)
            ax2 = plt.subplot(222)
            ax3 = plt.subplot(223)
            ax4 = plt.subplot(224)
            ax1.imshow(self.rho[step].transpose(), origin = 'lower')
            ax2.plot(x,y,ls = '', marker = 'o', markersize = 1)
            ax2.set_aspect('equal')
            ax3.imshow(self.vx[step].transpose(), origin = 'lower')
            ax4.imshow(self.vy[step].transpose(), origin = 'lower')
            plt.show()
            plt.close()
        return

