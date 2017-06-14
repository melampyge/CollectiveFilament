#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

##################################################################

# Density class: computes and stores the density as histograms

##################################################################

class Density:


    ##############################################################
      
    def __init__(self, nsteps, line):
        """ initialize: allocate density array"""
        ### set histogram parameters
        line = line.split()
        self.nbins2d = int(line[1])
        self.nbins1d = int(line[2])
        self.rho_min = float(line[3])
        self.rho_max = float(line[4])
        ### allocate arrays to store data
        self.rho2d = np.zeros((nsteps,self.nbins2d,self.nbins2d))
        self.xedges = 0
        self.yedges = 0
        self.hist_rho = np.zeros((nsteps,self.nbins1d))
        self.edges1d = 0
        return

    ##############################################################

    def compute(self,step,x,y,lx,ly,natoms, plot = 'False'):
        """ compute a density distribution and a histogram
            of the density distribution """
        ### compute 2D distribution of the density
        hist2d,self.xedges,self.yedges = np.histogram2d(x, y, bins = self.nbins2d, weights = np.ones(natoms))
        # transfrom number counts to densities
        hist2d *= self.nbins2d**2 / (lx*ly)
        ### compute histogram of the density
        hist1d, self.edges1d = np.histogram(hist2d, bins=self.nbins1d, range = (self.rho_min, self.rho_max))
        ### add data to results arrays
        self.rho2d[step] = np.copy(hist2d)
        self.hist_rho[step] = np.copy(hist1d)

        if plot == 'True':
            fig = plt.figure()
            ax1 = plt.subplot(121)
            ax2 = plt.subplot(122)
            ax1.imshow(self.rho2d[step].transpose(), origin = 'lower')
            ax2.plot(x,y,ls='',marker = 'o', markersize = 1)
            ax2.set_aspect('equal')
            plt.show()
            plt.close()
        return

